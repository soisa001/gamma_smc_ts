"""Existing phylogenetic programs benchmarked on synthetic SV genealogies.

Phases: prepare, iqtree, beast, relate, analyze. All truth is confined to panel
ascertainment, explicitly labelled oracle windows, and evaluation. The SV is
never supplied to sequence-tree inference. Outputs are synthetic only.
"""

import argparse
import hashlib
import io
import json
import os
import re
import subprocess
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import reduce
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import msprime
import numpy as np
import pandas as pd
import tskit
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def seed(cfg, name):
    return (
        int.from_bytes(
            hashlib.sha256(f"{cfg['seed']}:{name}".encode()).digest()[:4], "little"
        )
        % (2**31 - 1)
        + 1
    )


def valid_receipt(directory, name, spec):
    p = directory / (name + ".receipt.json")
    if not p.exists():
        return False
    old = json.loads(p.read_text())
    return old["spec"] == spec and all(
        (directory / f).exists() and sha(directory / f) == h
        for f, h in old["files"].items()
    )


def receipt(directory, name, spec, files, **extra):
    save(
        directory / (name + ".receipt.json"),
        {"spec": spec, "files": {str(f): sha(directory / f) for f in files}, **extra},
    )


def add_chimp(humans, cfg, name):
    """Continue fully coalesced human roots from the assumed species split.

    Human genealogies below the split are held fixed. No human root may exceed
    it. The chimp lineage and ancestral recombination/coalescence are simulated
    with msprime, conditional on those human histories. Re-overlay all SNPs.
    """
    split = cfg["chimp_split_years"] / cfg["generation_years"]
    assert max(humans.tables.nodes.time) < split
    tables = humans.dump_tables()
    tables.migrations.clear()
    tables.sites.clear()
    tables.mutations.clear()
    tables.populations.clear()
    tables.populations.metadata_schema = tskit.MetadataSchema(None)
    tables.populations.add_row()
    tables.nodes.population = np.zeros(len(tables.nodes), dtype=np.int32)
    chimp = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
    ts = msprime.sim_ancestry(
        initial_state=tables.tree_sequence(),
        start_time=split,
        population_size=cfg["human_chimp_ancestor_diploid_ne"],
        recombination_rate=cfg["recombination_rate"],
        random_seed=seed(cfg, name + ":chimp"),
    )
    ts = ts.simplify(samples=list(humans.samples()) + [chimp])
    ts = msprime.sim_mutations(
        ts,
        rate=cfg["mutation_rate"],
        model=msprime.JC69(),
        random_seed=seed(cfg, name + ":snps"),
    )
    assert ts.num_samples == humans.num_samples + 1
    return ts


def prepare_one(item, cfg, source, out):
    rep = source / f"rep{item['replicate']:04d}"
    directory = out / item["id"]
    directory.mkdir(parents=True, exist_ok=True)
    spec = {
        "phase": "prepare/v2",
        "cfg": cfg,
        "ancestry": sha(rep / "ancestry.trees"),
        "observations": sha(rep / "observations.json"),
    }
    if valid_receipt(directory, "prepare", spec):
        return {"id": item["id"], "cached": True}
    obs = next(
        x
        for x in json.loads((rep / "observations.json").read_text())
        if x["id"] == item["id"]
    )
    states = np.array(obs["carriers"])
    rng = np.random.default_rng(seed(cfg, item["id"] + ":sample"))
    n = cfg["modern_haplotypes"]
    nc = min(n // 2, int(states.sum()))
    chosen_c = rng.choice(np.flatnonzero(states), nc, replace=False)
    chosen_n = rng.choice(np.flatnonzero(~states), n - nc, replace=False)
    chosen = np.sort(np.r_[chosen_c, chosen_n]).astype(int)
    ancestry = tskit.load(rep / "ancestry.trees")
    tables = ancestry.dump_tables()
    tables.migrations.clear()
    humans = tables.tree_sequence().simplify(samples=ancestry.samples()[chosen])
    focal = humans.at(obs["position"])
    carriers = [f"h{i:03d}" for i in range(n) if states[chosen[i]]]
    clade = reduce(focal.mrca, np.flatnonzero(states[chosen]))
    lower = focal.time(clade) * cfg["generation_years"]
    upper = focal.time(focal.parent(clade)) * cfg["generation_years"]
    assert lower <= item["age_years"] <= upper
    ts = add_chimp(humans, cfg, item["id"])
    ts.dump(directory / "chimp_snps.trees")
    labels = [f"h{i:03d}" for i in range(n)] + ["chimp"]
    length = int(ts.sequence_length)
    # At invariant/unmutated sites the ancestral state is unobserved. Under JC69
    # a deterministic random uniform base is valid; site-specific simulation
    # ancestral states replace it wherever msprime recorded mutations.
    bases = np.array(list("ACGT"), dtype="U1")
    anc = rng.choice(bases, size=length)
    alignment = np.tile(anc, (n + 1, 1))
    haps = []
    flipped = dropped = 0
    for var in ts.variants():
        pos = int(var.site.position)
        alleles = np.array(var.alleles)
        calls = alleles[var.genotypes]
        alignment[:, pos] = calls
        modern = sorted(set(calls[:-1]))
        if len(modern) != 2:
            continue
        if calls[-1] not in modern:
            dropped += 1
            continue
        ancestral = calls[-1]
        derived = next(a for a in modern if a != ancestral)
        flipped += ancestral != var.site.ancestral_state
        haps.append(
            f"1 s{pos} {pos + 1} {ancestral} {derived} "
            + " ".join(str(int(x != ancestral)) for x in calls[:-1])
        )
    # Mask focal coordinate: the SV is not an observed nucleotide character.
    alignment[:, obs["position"]] = "N"
    haps = [h for h in haps if int(h.split()[2]) != obs["position"] + 1]
    (directory / "human.haps").write_text("\n".join(haps) + "\n")
    (directory / "human.sample").write_text(
        "ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{x} NA 0\n" for x in labels[:-1])
    )
    (directory / "genetic.map").write_text(
        f"pos rate map\n1 1 0\n{length} 1 {(length - 1) * 1e-6}\n"
    )
    (directory / "rate.txt").write_text(str(cfg["mutation_rate"]) + "\n")
    windows = {
        "oracle": [
            int(np.ceil(focal.interval.left)),
            int(np.ceil(focal.interval.right)),
        ]
    }
    for width in cfg["windows_bp"]:
        windows[f"w{width}"] = [
            max(0, obs["position"] - width // 2),
            min(length, obs["position"] + width // 2),
        ]
    winstats = {}
    files = [
        "chimp_snps.trees",
        "human.haps",
        "human.sample",
        "genetic.map",
        "rate.txt",
    ]
    for name, (left, right) in windows.items():
        sub = alignment[:, left:right]
        (directory / f"{name}.fasta").write_text(
            "".join(f">{label}\n{''.join(seq)}\n" for label, seq in zip(labels, sub))
        )
        polymorphic = np.any(sub != sub[0], axis=0) & np.all(sub != "N", axis=0)
        human_poly = np.any(sub[:-1] != sub[0], axis=0) & np.all(sub != "N", axis=0)
        trees = sum(
            t.interval.left < right and t.interval.right > left for t in humans.trees()
        )
        winstats[name] = {
            "left": left,
            "right": right,
            "length": right - left,
            "variable_sites": int(polymorphic.sum()),
            "human_variable_sites": int(human_poly.sum()),
            "true_human_trees": trees,
        }
        files.append(f"{name}.fasta")
    save(
        directory / "input.json",
        {
            "id": item["id"],
            "replicate": item["replicate"],
            "position": obs["position"],
            "selected_modern_haplotypes": chosen.tolist(),
            "carriers": carriers,
            "labels": labels,
            "windows": winstats,
            "relate_snps": len(haps),
            "chimp_polarization_flips": int(flipped),
            "chimp_polarization_dropped": dropped,
        },
    )
    save(
        directory / "evaluation.json",
        {
            **item,
            "subset_lower_years": lower,
            "subset_upper_years": upper,
            "source_truth_sha256": sha(rep / "truth.json"),
        },
    )
    files += ["input.json", "evaluation.json"]
    receipt(directory, "prepare", spec, files)
    return {"id": item["id"], "windows": winstats, "relate_snps": len(haps)}


def run_command(command, directory, log, timeout):
    start = time.monotonic()
    save(
        directory / (log + ".command.json"),
        {
            "command": command,
            "cwd": str(directory),
            "timeout_seconds": timeout,
            "start_unix": time.time(),
        },
    )
    with (directory / (log + ".stdout.txt")).open("w") as f:
        p = subprocess.run(
            command, cwd=directory, stdout=f, stderr=subprocess.STDOUT, timeout=timeout
        )
    if p.returncode:
        raise RuntimeError(
            f"{directory.name}/{log}: exit {p.returncode}: "
            + (directory / (log + ".stdout.txt")).read_text(errors="replace")[-2500:]
        )
    return time.monotonic() - start


def iqtree_one(directory, window, cfg, binary):
    name = "iq_" + window
    spec = {
        "alignment": sha(directory / f"{window}.fasta"),
        "binary": sha(binary),
        "model": "JC",
        "bootstraps": cfg["iqtree_bootstraps"],
        "seed": seed(cfg, directory.name + name),
        "rate": cfg["mutation_rate"],
        "phase": "iqtree/v2",
    }
    if valid_receipt(directory, name, spec):
        return {"task": directory.name + "/" + name, "cached": True}
    labels = json.loads((directory / "input.json").read_text())["labels"]
    (directory / "dates.txt").write_text("".join(f"{label} 0\n" for label in labels))
    common = [
        str(binary),
        "-s",
        f"{window}.fasta",
        "-m",
        "JC",
        "-o",
        "chimp",
        "-nt",
        "1",
        "-seed",
        str(spec["seed"]),
        "-redo",
        "-keep-ident",
    ]
    cmd = common + ["-pre", name, "-b", str(cfg["iqtree_bootstraps"])]
    if not valid_receipt(directory, name + "_ml", spec):
        runtime = run_command(cmd, directory, name, cfg["timeout_seconds"])
        assert (
            len(list(Phylo.parse(directory / (name + ".boottrees"), "newick")))
            == cfg["iqtree_bootstraps"]
        )
        receipt(
            directory,
            name + "_ml",
            spec,
            [name + ".treefile", name + ".boottrees", name + ".iqtree"],
            runtime_seconds=runtime,
        )
    else:
        runtime = json.loads((directory / (name + "_ml.receipt.json")).read_text())[
            "runtime_seconds"
        ]
    # IQ-TREE's standard bootstrap mode does not run the dating stage itself.
    cmd = common + [
        "-te",
        name + ".treefile",
        "-pre",
        name + "_dating",
        "--date",
        "dates.txt",
        "--date-options",
        "-w rate.txt -l -1 -u 0 -U 0",
    ]
    runtime += run_command(cmd, directory, name + "_dating", cfg["timeout_seconds"])
    for ext in ("nwk", "lsd", "nex"):
        (directory / (name + ".timetree." + ext)).write_bytes(
            (directory / (name + "_dating.timetree." + ext)).read_bytes()
        )
    # LSD2's Newick output retains substitution units; its NEXUS tree carries
    # calendar/time units. With rate per generation, NEXUS edges are generations.
    treefile = name + ".timetree.nex"
    tree = read_lsd_time_tree(directory / treefile)
    assert len(tree.get_terminals()) == cfg["modern_haplotypes"] + 1
    receipt(
        directory,
        name,
        spec,
        [
            treefile,
            name + ".treefile",
            name + ".iqtree",
            name + ".boottrees",
            name + ".timetree.lsd",
        ],
        runtime_seconds=runtime,
    )
    return {"task": directory.name + "/" + name, "seconds": runtime}


def beast_xml(directory, window, chain, cfg):
    from Bio import SeqIO

    data = list(SeqIO.parse(directory / f"{window}.fasta", "fasta"))
    # Starting topology is estimated from sequence only by IQ-TREE; date in
    # generations. It supplies no truth and is free to change during MCMC.
    init = read_lsd_time_tree(directory / f"iq_{window}.timetree.nex")
    for node in init.find_clades():
        node.confidence = None
        node.comment = None
        node.name = node.name if node.is_terminal() else None
        if node.branch_length is not None:
            node.branch_length = max(node.branch_length, 0.1)
    # Tiny positive edges above may break ultrametricity: explicitly refit node
    # heights bottom-up before emitting an initial state.
    heights = {}
    for node in init.find_clades(order="postorder"):
        heights[node] = (
            0.0
            if node.is_terminal()
            else max(heights[c] + c.branch_length for c in node.clades)
        )
        for child in node.clades:
            child.branch_length = heights[node] - heights[child]
    init.root.branch_length = 0
    buf = io.StringIO()
    Phylo.write(init, buf, "newick", format_branch_length="%.12g")
    newick = buf.getvalue().strip()
    prefix = f"beast_{window}_{chain}"
    sequences = "\n".join(
        f'<sequence spec="Sequence" taxon="{s.id}" value="{s.seq}"/>' for s in data
    )
    taxa = "".join(
        f'<taxon id="{s.id}" spec="Taxon"/>' for s in data if s.id != "chimp"
    )
    xml = f'''<beast version="2.7" namespace="beast.base.core:beast.base.inference:beast.base.inference.parameter:beast.base.inference.distribution:beast.base.evolution.alignment:beast.base.evolution.tree:beast.base.evolution.tree.coalescent:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood:beast.base.evolution.branchratemodel">
<data id="alignment" spec="Alignment" dataType="nucleotide">{sequences}</data>
<run id="mcmc" spec="MCMC" chainLength="{cfg["beast_chain_length"]}">
<state storeEvery="100000">
<stateNode id="tree" spec="TreeParser" taxa="@alignment" newick="{newick}" IsLabelledNewick="true" adjustTipHeights="false"/>
<stateNode id="popSize" spec="RealParameter" lower="1" upper="10000000">{cfg["beast_haploid_ne_initial"]}</stateNode>
</state>
<distribution id="posterior" spec="CompoundDistribution">
<distribution id="prior" spec="CompoundDistribution">
<distribution id="coalescent" spec="Coalescent">
<treeIntervals spec="TreeIntervals" tree="@tree"/>
<populationModel spec="ConstantPopulation" popSize="@popSize"/>
</distribution>
<distribution spec="Prior" x="@popSize"><distr spec="LogNormalDistributionModel" M="10.308952660644293" S="2.0" meanInRealSpace="false"/></distribution>
<distribution id="humanMonophyly" spec="MRCAPrior" tree="@tree" monophyletic="true"><taxonset spec="TaxonSet">{taxa}</taxonset></distribution>
</distribution>
<distribution id="likelihood" spec="TreeLikelihood" data="@alignment" tree="@tree">
<siteModel spec="SiteModel"><substModel spec="JukesCantor"/></siteModel>
<branchRateModel spec="StrictClockModel" clock.rate="{cfg["mutation_rate"]}"/>
</distribution>
</distribution>
<operator spec="ScaleOperator" parameter="@popSize" scaleFactor="0.75" weight="3"/>
<operator spec="ScaleOperator" tree="@tree" scaleFactor="0.75" weight="3"/>
<operator spec="ScaleOperator" tree="@tree" rootOnly="true" scaleFactor="0.75" weight="3"/>
<operator spec="beast.base.evolution.operator.Uniform" tree="@tree" weight="30"/>
<operator spec="SubtreeSlide" tree="@tree" size="1000" gaussian="true" weight="15"/>
<operator spec="Exchange" tree="@tree" isNarrow="true" weight="15"/>
<operator spec="Exchange" tree="@tree" isNarrow="false" weight="3"/>
<operator spec="WilsonBalding" tree="@tree" weight="3"/>
<logger fileName="{prefix}.log" logEvery="{cfg["beast_log_every"]}"><log idref="posterior"/><log idref="likelihood"/><log idref="prior"/><log idref="popSize"/><log spec="TreeHeightLogger" tree="@tree"/></logger>
<logger fileName="{prefix}.trees" logEvery="{cfg["beast_log_every"]}" mode="tree"><log idref="tree"/></logger>
<logger logEvery="100000"><log idref="posterior"/></logger>
</run></beast>'''
    return prefix, xml


def at_least_requested_chain(existing, requested):
    pattern = r'chainLength="(\d+)"'
    old, new = re.search(pattern, existing), re.search(pattern, requested)
    return bool(
        old
        and new
        and int(old[1]) >= int(new[1])
        and re.sub(pattern, 'chainLength="ANY"', existing)
        == re.sub(pattern, 'chainLength="ANY"', requested)
    )


def beast_one(directory, window, chain, cfg, tools):
    deadline = time.monotonic() + cfg["timeout_seconds"]
    while not (directory / f"iq_{window}.receipt.json").exists():
        if time.monotonic() > deadline:
            raise TimeoutError(
                f"Waiting for IQ-TREE initialization: {directory}/{window}"
            )
        time.sleep(1)
    name, xml = beast_xml(directory, window, chain, cfg)
    jars = [
        tools / "launcher.jar",
        tools / "beast-base/lib/BEAST.base.jar",
        tools / "beast-app/lib/BEAST.app.jar",
    ]
    common_spec = {
        "jars": [sha(p) for p in jars],
        "seed": seed(cfg, directory.name + name),
    }
    prior_receipt = directory / (name + ".receipt.json")
    xml_path = directory / (name + ".xml")
    if prior_receipt.exists() and xml_path.exists():
        old_spec = json.loads(prior_receipt.read_text())["spec"]
        if (
            all(old_spec[k] == v for k, v in common_spec.items())
            and at_least_requested_chain(xml_path.read_text(), xml)
            and valid_receipt(directory, name, old_spec)
        ):
            return {"task": directory.name + "/" + name, "cached": True}
    xml_path.write_text(xml)
    spec = {**common_spec, "xml": sha(xml_path)}
    cmd = [
        "java",
        "-Xmx5g",
        "-Djava.only=true",
        "-cp",
        os.pathsep.join(str(p) for p in jars),
        "beastfx.app.beast.BeastMCMC",
        "-overwrite",
        "-seed",
        str(spec["seed"]),
        "-threads",
        "1",
        name + ".xml",
    ]
    runtime = run_command(cmd, directory, name, cfg["timeout_seconds"])
    stdout = (directory / (name + ".stdout.txt")).read_text(errors="replace")
    # BeastMCMC can return exit code zero after XML parsing errors. Never accept
    # stale output traces left by an earlier run in that case.
    assert "Total calculation time:" in stdout and "End likelihood:" in stdout, (
        f"BEAST did not finish: {directory}/{name}"
    )
    trace = pd.read_csv(directory / (name + ".log"), sep="\t", comment="#")
    assert trace.iloc[-1, 0] == cfg["beast_chain_length"]
    receipt(
        directory,
        name,
        spec,
        [name + ".log", name + ".trees", name + ".xml"],
        runtime_seconds=runtime,
    )
    return {"task": directory.name + "/" + name, "seconds": runtime}


def relate_one(directory, cfg, binary):
    spec = {
        "input": {
            f: sha(directory / f) for f in ("human.haps", "human.sample", "genetic.map")
        },
        "binary": sha(binary),
        "mu": cfg["mutation_rate"],
        "haploid_ne": cfg["relate_haploid_ne"],
        "seed": seed(cfg, directory.name + ":relate"),
    }
    if valid_receipt(directory, "relate", spec):
        return {"task": directory.name + "/relate", "cached": True}
    scratch = directory / "relate"
    if scratch.exists():
        assert scratch.resolve().parent == directory.resolve()
        scratch.rename(directory / f"relate_scratch_{time.time_ns()}")
    cmd = [
        str(binary),
        "--mode",
        "All",
        "--haps",
        "human.haps",
        "--sample",
        "human.sample",
        "--map",
        "genetic.map",
        "-m",
        str(cfg["mutation_rate"]),
        "-N",
        str(cfg["relate_haploid_ne"]),
        "--seed",
        str(spec["seed"]),
        "--memory",
        "1",
        "-o",
        "relate",
    ]
    runtime = run_command(cmd, directory, "relate", cfg["timeout_seconds"])
    assert (directory / "relate.anc").stat().st_size > 100
    receipt(
        directory, "relate", spec, ["relate.anc", "relate.mut"], runtime_seconds=runtime
    )
    return {"task": directory.name + "/relate", "seconds": runtime}


def read_lsd_time_tree(path):
    tree = Phylo.read(path, "nexus")
    date = re.search(r'date="?(-?[\d.eE+]+)', tree.root.comment or "")
    assert date, "Missing LSD2 dated root annotation"
    depth = tree.depths()
    observed = np.mean([depth[t] for t in tree.get_terminals()])
    assert np.isclose(observed, -float(date[1]), rtol=1e-5, atol=0.5), (
        "LSD2 edge/date units disagree"
    )
    return tree


def branch_bounds(tree, carriers, scale=25):
    """Return a branch only for an exact carrier clade; never force placement."""
    carriers = set(carriers)
    tips = {t.name for t in tree.get_terminals()}
    assert carriers < tips
    crown = tree.common_ancestor(list(carriers))
    if {t.name for t in crown.get_terminals()} != carriers:
        return None
    parents = {c: p for p in tree.find_clades() for c in p.clades}
    if crown not in parents:
        return None
    depth = tree.depths()
    terminal_depth = [depth[t] for t in tree.get_terminals()]
    height = float(np.mean(terminal_depth))
    assert np.ptp(terminal_depth) < max(0.5, height * 1e-5), "Not a time tree"
    # Printed LSD2 edges have finite precision. Local bottom-up heights avoid
    # subtracting two large root distances for recent nodes with zero edges.
    heights = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            heights[node] = 0.0
        else:
            heights[node] = max(
                float(
                    np.mean([heights[c] + (c.branch_length or 0) for c in node.clades])
                ),
                max(heights[c] for c in node.clades),
            )
    lo, hi = heights[crown] * scale, heights[parents[crown]] * scale
    assert lo >= 0 and hi >= lo
    return lo, hi


def read_relate(directory, position):
    from Bio.Phylo.BaseTree import Clade, Tree

    mut = pd.read_csv(directory / "relate.mut", sep=";")
    anchor = int(
        np.clip(
            np.searchsorted(mut.pos_of_snp.to_numpy(), position + 1, side="right") - 1,
            0,
            len(mut) - 1,
        )
    )
    index = int(mut.iloc[anchor].tree_index)
    lines = (directory / "relate.anc").read_text().splitlines()
    n = int(lines[0].split()[1])
    assert len(lines) - 2 == int(lines[1].split()[1])
    entries = re.findall(r"(-?\d+):\(([^)]+)\)", lines[index + 2])
    assert len(entries) == 2 * n - 1
    nodes = [Clade(name=f"h{i:03d}" if i < n else None) for i in range(len(entries))]
    roots = []
    for i, (parent, values) in enumerate(entries):
        nodes[i].branch_length = float(values.split()[0])
        if int(parent) == -1:
            roots.append(nodes[i])
        else:
            assert i < int(parent) < len(entries)
            nodes[int(parent)].clades.append(nodes[i])
    assert len(roots) == 1
    return Tree(root=roots[0], rooted=True), int(mut.iloc[anchor].pos_of_snp) - 1


def ess(values):
    """Initial-positive-pair autocorrelation estimate; constant traces flagged."""
    x = np.asarray(values, dtype=float)
    x = x - x.mean()
    if np.var(x) == 0:
        return 0.0
    f = np.fft.rfft(x, n=2 * len(x))
    ac = np.fft.irfft(f * f.conj())[: len(x)]
    ac /= ac[0]
    total = 0.0
    previous = np.inf
    for i in range(1, len(ac) - 1, 2):
        pair = min(previous, ac[i] + ac[i + 1])
        if pair <= 0:
            break
        total += pair
        previous = pair
    return float(min(len(x), len(x) / (1 + 2 * total)))


def split_rhat(chains):
    length = min(map(len, chains)) // 2
    a = np.array(
        [
            part
            for c in chains
            for part in (np.asarray(c)[:length], np.asarray(c)[-length:])
        ]
    )
    within = np.mean(np.var(a, axis=1, ddof=1))
    if within == 0:
        return float("inf")
    between = length * np.var(a.mean(axis=1), ddof=1)
    return float(np.sqrt(((length - 1) / length * within + between / length) / within))


def summarize_bounds(bounds, truth, posterior=False):
    result = {
        "monophyly_probability": len(bounds[0]) / bounds[1],
        "lower_years": np.nan,
        "upper_years": np.nan,
        "estimate_years": np.nan,
        "age_contained": False,
        "fold_error": np.nan,
        "timing_call": "unresolved",
    }
    b = np.asarray(bounds[0])
    if not len(b):
        return result
    result["branch_lower_median_years"] = float(np.median(b[:, 0]))
    result["branch_upper_median_years"] = float(np.median(b[:, 1]))
    if posterior:
        # Deterministic quadrature for uniform placement conditional on each
        # compatible sampled branch, then mix across posterior trees.
        q = (np.arange(100) + 0.5) / 100
        ages = (b[:, 0, None] + (b[:, 1] - b[:, 0])[:, None] * q).ravel()
        lo, point, hi = np.quantile(ages, [0.05, 0.5, 0.95])
    else:
        lo, hi = b[0]
        point = (lo + hi) / 2
    result.update(
        lower_years=lo,
        upper_years=hi,
        estimate_years=point,
        age_contained=bool(lo <= truth["age_years"] <= hi),
        fold_error=max(point / truth["age_years"], truth["age_years"] / point)
        if point > 0
        else np.inf,
    )
    if result["monophyly_probability"] >= 0.95:
        result["timing_call"] = (
            "post" if hi < 50000 else ("pre" if lo > 50000 else "unresolved")
        )
    return result


def analyze(dirs, cfg, out):
    rows, diagnostics, windows = [], [], []
    for d in dirs:
        inp = json.loads((d / "input.json").read_text())
        truth = json.loads((d / "evaluation.json").read_text())
        base = {
            "id": d.name,
            "replicate": truth["replicate"],
            "history": truth["history"],
            "true_age_years": truth["age_years"],
            "true_lower_years": truth["subset_lower_years"],
            "true_upper_years": truth["subset_upper_years"],
        }
        rows.append(
            {
                **base,
                "method": "oracle_tree",
                "window": "focal",
                "runtime_seconds": 0,
                "usable": True,
                **summarize_bounds(
                    ([(truth["subset_lower_years"], truth["subset_upper_years"])], 1),
                    truth,
                ),
            }
        )
        for window, stats in inp["windows"].items():
            windows.append(
                {
                    **base,
                    "window": window,
                    **stats,
                    "relate_snps": inp["relate_snps"],
                    "chimp_polarization_flips": inp["chimp_polarization_flips"],
                }
            )
            name = "iq_" + window
            if (d / (name + ".receipt.json")).exists():
                tree = read_lsd_time_tree(d / (name + ".timetree.nex"))
                b = branch_bounds(tree, inp["carriers"], cfg["generation_years"])
                support = []
                for bt in Phylo.parse(d / (name + ".boottrees"), "newick"):
                    bt.root_with_outgroup("chimp")
                    crown = bt.common_ancestor(inp["carriers"])
                    support.append(
                        {x.name for x in crown.get_terminals()} == set(inp["carriers"])
                    )
                rec = json.loads((d / (name + ".receipt.json")).read_text())
                row = {
                    **base,
                    "method": "IQTREE_LSD2",
                    "window": window,
                    "runtime_seconds": rec["runtime_seconds"],
                    "bootstrap_carrier_clade": np.mean(support),
                    "usable": b is not None,
                    **summarize_bounds(([b] if b else [], 1), truth),
                }
                rows.append(row)
            if window in cfg["beast_windows"]:
                chain_bounds, logs, indicator_chains, crown_chains, parent_chains = (
                    [],
                    [],
                    [],
                    [],
                    [],
                )
                runtime = 0
                for chain in range(cfg["beast_chains"]):
                    name = f"beast_{window}_{chain}"
                    if not (d / (name + ".receipt.json")).exists():
                        break
                    trace = pd.read_csv(d / (name + ".log"), sep="\t", comment="#")
                    burn = int(len(trace) * cfg["beast_burnin_fraction"])
                    logs.append(trace.iloc[burn:])
                    trees = list(Phylo.parse(d / (name + ".trees"), "nexus"))[burn:]
                    one, indicators, crowns, parents = [], [], [], []
                    for tree in trees:
                        b = branch_bounds(
                            tree, inp["carriers"], cfg["generation_years"]
                        )
                        indicators.append(int(b is not None))
                        if b is not None:
                            one.append(b)
                        # Track MRCA age even if carriers fail monophyly, so
                        # convergence diagnostic does not discard bad states.
                        dep = tree.depths()
                        height = np.mean([dep[t] for t in tree.get_terminals()])
                        crown = tree.common_ancestor(inp["carriers"])
                        parentmap = {c: p for p in tree.find_clades() for c in p.clades}
                        crowns.append(height - dep[crown])
                        parents.append(
                            height - dep[parentmap[crown]]
                            if crown in parentmap
                            else height
                        )
                    chain_bounds.extend(one)
                    indicator_chains.append(indicators)
                    crown_chains.append(crowns)
                    parent_chains.append(parents)
                    runtime += json.loads((d / (name + ".receipt.json")).read_text())[
                        "runtime_seconds"
                    ]
                if len(logs) == cfg["beast_chains"]:
                    checks = []
                    all_series = {
                        c: [l[c].to_numpy() for l in logs] for c in logs[0].columns[1:]
                    }
                    all_series.update(
                        carrier_mrca=crown_chains, carrier_parent=parent_chains
                    )
                    indicator_means = [np.mean(x) for x in indicator_chains]
                    checks.append(np.ptp(indicator_means) <= 0.05)
                    for column, values in all_series.items():
                        effective = min(ess(x) for x in values)
                        rhat = split_rhat(values)
                        checks.append(
                            effective >= cfg["beast_min_ess"]
                            and rhat <= cfg["beast_max_rhat"]
                        )
                        diagnostics.append(
                            {
                                "id": d.name,
                                "window": window,
                                "parameter": column,
                                "min_chain_ess": effective,
                                "split_rhat": rhat,
                            }
                        )
                    converged = all(checks)
                    summary = summarize_bounds(
                        (chain_bounds, sum(map(len, indicator_chains))),
                        truth,
                        posterior=True,
                    )
                    row = {
                        **base,
                        "method": "BEAST2",
                        "window": window,
                        "runtime_seconds": runtime,
                        "converged": converged,
                        "usable": converged
                        and summary["monophyly_probability"] >= 0.95,
                        **summary,
                    }
                    if not row["usable"]:
                        row["timing_call"] = "unresolved"
                    rows.append(row)
        if (d / "relate.receipt.json").exists():
            tree, anchor = read_relate(d, inp["position"])
            b = branch_bounds(tree, inp["carriers"], cfg["generation_years"])
            rows.append(
                {
                    **base,
                    "method": "Relate",
                    "window": "500kb",
                    "runtime_seconds": json.loads(
                        (d / "relate.receipt.json").read_text()
                    )["runtime_seconds"],
                    "anchor_distance_bp": abs(anchor - inp["position"]),
                    "usable": b is not None,
                    **summarize_bounds(([b] if b else [], 1), truth),
                }
            )
        print(f"Analyzed {d.name}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "estimates.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(
        out / "beast_diagnostics.tsv", sep="\t", index=False
    )
    pd.DataFrame(windows).to_csv(out / "windows.tsv", sep="\t", index=False)
    summaries = []
    for (method, window), group in frame.groupby(["method", "window"]):
        usable = group[group.usable]
        calls = group[group.timing_call != "unresolved"]
        wrong = sum((calls.timing_call == "pre") != (calls.true_age_years > 50000))
        summaries.append(
            {
                "method": method,
                "window": window,
                "events": len(group),
                "usable_events": len(usable),
                "mean_monophyly_probability": group.monophyly_probability.mean(),
                "median_fold_error_usable": usable.fold_error.median(),
                "age_containment_usable": usable.age_contained.mean(),
                "timing_calls": len(calls),
                "wrong_calls": wrong,
                "median_runtime_seconds": group.runtime_seconds.median(),
            }
        )
    pd.DataFrame(summaries).to_csv(out / "summary.tsv", sep="\t", index=False)
    detailed = []
    for (method, window, history), group in frame.groupby(
        ["method", "window", "history"]
    ):
        usable = group[group.usable]
        detailed.append(
            {
                "method": method,
                "window": window,
                "history": history,
                "n": len(group),
                "usable": len(usable),
                "mean_monophyly_probability": group.monophyly_probability.mean(),
                "median_fold_error_usable": usable.fold_error.median(),
                "age_containment_usable": usable.age_contained.mean(),
                "timing_calls": int((group.timing_call != "unresolved").sum()),
                "median_bootstrap_carrier_clade": group.bootstrap_carrier_clade.median(),
            }
        )
    pd.DataFrame(detailed).to_csv(out / "by_history.tsv", sep="\t", index=False)
    print(pd.DataFrame(summaries).to_string(index=False))


def audit(dirs, cfg, out):
    checks = []
    for d in dirs:
        inp = json.loads((d / "input.json").read_text())
        for p in d.glob("*.receipt.json"):
            rec = json.loads(p.read_text())
            assert valid_receipt(
                d, p.name.removesuffix(".receipt.json"), rec["spec"]
            ), f"Corrupted receipt: {p}"
            checks.append(
                {"event": d.name, "receipt": p.name, "checked_files": len(rec["files"])}
            )
        for w in inp["windows"]:
            tree = read_lsd_time_tree(d / f"iq_{w}.timetree.nex")
            assert {t.name for t in tree.get_terminals()} == set(inp["labels"])
        read_relate(d, inp["position"])
        for w in cfg["beast_windows"]:
            for chain in range(cfg["beast_chains"]):
                name = f"beast_{w}_{chain}"
                trace = pd.read_csv(d / (name + ".log"), sep="\t", comment="#")
                assert len(trace) >= 1001 and np.isfinite(trace.to_numpy()).all()
                assert trace.iloc[-1, 0] >= cfg["beast_chain_length"]
                states = [
                    int(x)
                    for x in re.findall(
                        r"tree STATE_(\d+)", (d / (name + ".trees")).read_text()
                    )
                ]
                assert np.array_equal(states, trace.iloc[:, 0].to_numpy()), name
    save(
        out / "audit.json",
        {
            "passed": True,
            "events": len(dirs),
            "checks": checks,
            "source_sha256": sha(__file__),
        },
    )
    print(f"Audited {len(dirs)} events and {len(checks)} hash receipts")


def restore_bundle(out, bundle):
    manifest_path = bundle.with_name("inputs_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    assert sha(bundle) == manifest["zip_sha256"]
    with zipfile.ZipFile(bundle) as z:
        assert z.testzip() is None
        for item in z.infolist():
            path = out / item.filename
            assert path.resolve().is_relative_to(out.resolve())
            data = z.read(item)
            assert hashlib.sha256(data).hexdigest() == manifest["files"][item.filename]
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists() and path.read_bytes() == data:
                continue
            if path.exists():
                raise RuntimeError(
                    f"Refusing to overwrite different existing input: {path}"
                )
            path.write_bytes(data)
    print(f"Restored and verified {len(manifest['files'])} synthetic input files")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=["prepare", "restore", "iqtree", "beast", "relate", "analyze", "audit"],
        required=True,
    )
    parser.add_argument(
        "--source", type=Path, default=ROOT.parent / "sim/sv_mutation_age_20261024"
    )
    parser.add_argument(
        "--out", type=Path, default=ROOT.parent / "sim/sv_phylogeny_20261028"
    )
    parser.add_argument(
        "--config", type=Path, default=ROOT / "configs/sv_phylogeny.json"
    )
    parser.add_argument("--tools", type=Path, default=ROOT / ".tools/phylo")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--event", help="Run only this existing panel event ID")
    parser.add_argument("--window")
    parser.add_argument("--chain-length", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument(
        "--bundle",
        type=Path,
        default=ROOT / "docs/results/sv_phylogeny_20261028/benchmark_inputs.zip",
    )
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if args.chain_length:
        cfg["beast_chain_length"] = args.chain_length
    if args.workers:
        cfg["workers"] = args.workers
    args.out.mkdir(parents=True, exist_ok=True)
    invocation = {
        "args": {
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        "config": cfg,
        "source_sha256": sha(__file__),
        "time_unix": time.time(),
    }
    (args.out / "invocations").mkdir(exist_ok=True)
    save(args.out / "invocations" / f"{time.time_ns()}_{args.phase}.json", invocation)
    save(args.out / f"invocation_{args.phase}.json", invocation)
    if args.phase == "restore":
        restore_bundle(args.out, args.bundle)
        return
    jobs = []
    if args.phase == "prepare":
        rows = []
        for rep in sorted(args.source.glob("rep*")):
            if int(rep.name[3:]) % 4 == 3:
                rows.extend(json.loads((rep / "truth.json").read_text()))
        selected = []
        for history in sorted({r["history"] for r in rows}):
            selected.extend(
                [r for r in rows if r["history"] == history][
                    : cfg["events_per_history"]
                ]
            )
        selected.sort(key=lambda r: r["id"])
        save(args.out / "panel.json", selected)
        for r in selected[: args.limit]:
            jobs.append((prepare_one, (r, cfg, args.source, args.out)))
    else:
        dirs = [
            args.out / r["id"]
            for r in json.loads((args.out / "panel.json").read_text())
            if not args.event or r["id"] == args.event
        ][: args.limit]
        for d in dirs:
            if args.phase == "iqtree":
                binary = (
                    next(args.tools.glob("iqtree/**/iqtree2.exe"))
                    if os.name == "nt"
                    else Path(os.environ["IQTREE_BIN"])
                )
                for w in json.loads((d / "input.json").read_text())["windows"]:
                    if not args.window or w == args.window:
                        jobs.append((iqtree_one, (d, w, cfg, binary)))
            elif args.phase == "beast":
                for w in cfg["beast_windows"]:
                    if not args.window or w == args.window:
                        for chain in range(cfg["beast_chains"]):
                            jobs.append((beast_one, (d, w, chain, cfg, args.tools)))
            elif args.phase == "relate":
                binary = (
                    args.tools / "relate-build/Relate.exe"
                    if os.name == "nt"
                    else Path(os.environ["RELATE_BIN"])
                )
                jobs.append((relate_one, (d, cfg, binary)))
        if args.phase == "analyze":
            analyze(dirs, cfg, args.out)
            return
        if args.phase == "audit":
            audit(dirs, cfg, args.out)
            return
    failures = []
    with ThreadPoolExecutor(max_workers=cfg["workers"]) as pool:
        futures = [pool.submit(fn, *a) for fn, a in jobs]
        for future in as_completed(futures):
            try:
                print(json.dumps(future.result()), flush=True)
            except Exception as exc:
                failures.append(str(exc))
                print(json.dumps({"error": str(exc)}), flush=True)
    save(args.out / f"{args.phase}_failures.json", failures)
    if failures:
        raise SystemExit(
            f"{len(failures)} failed tasks; see {args.phase}_failures.json"
        )


if __name__ == "__main__":
    main()

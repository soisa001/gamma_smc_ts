"""Build unmodified Relate algorithms with a Windows resource-reporting shim.

Requires a sparse Relate checkout, w64devkit and zlib-1.3.1 under .tools/phylo.
Only getrusage (diagnostic reporting) is supplied locally; algorithm source is
unmodified. The benchmark measures elapsed time outside the executable.
"""

import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / ".tools/phylo"
src = base / "relate-sparse/include"
out = base / "relate-build"
out.mkdir(exist_ok=True)
(out / "sys").mkdir(exist_ok=True)
(out / "sys/resource.h").write_text("""#pragma once
#include <sys/time.h>
#include <windows.h>
#include <psapi.h>
#define RUSAGE_SELF 0
struct rusage { timeval ru_utime; long long ru_maxrss; };
inline int getrusage(int, rusage *r) {
 FILETIME c,e,k,u; GetProcessTimes(GetCurrentProcess(),&c,&e,&k,&u);
 ULARGE_INTEGER v; v.LowPart=u.dwLowDateTime;v.HighPart=u.dwHighDateTime;
 r->ru_utime.tv_sec=v.QuadPart/10000000;
 r->ru_utime.tv_usec=(v.QuadPart%10000000)/10;
 PROCESS_MEMORY_COUNTERS p; GetProcessMemoryInfo(GetCurrentProcess(),&p,sizeof(p));
 r->ru_maxrss=p.PeakWorkingSetSize/1024; return 0;
}
""")
zlib = base / "zlib-1.3.1"
gcc = base / "w64devkit/bin/gcc.exe"
cpp = base / "w64devkit/bin/g++.exe"
cfiles = [
    zlib / (n + ".c")
    for n in "adler32 compress crc32 deflate gzclose gzlib gzread gzwrite infback inffast inflate inftrees trees uncompr zutil".split()
]
cppfiles = [
    src / "src" / (n + ".cpp")
    for n in "cxxopts filesystem plot fast_log collapsed_matrix data sample fast_painting anc mutations tree_builder branch_length_estimator anc_builder tree_comparer".split()
]
cppfiles += [src / "src/gzstream/gzstream.cpp", src / "pipeline/Relate.cpp"]
cppfiles = [p for p in cppfiles if p.exists()]
# Windows forbids unlinking an open file. Upstream leaves a paint stream open
# until process termination. Preserve those scratch files and warn on directory
# cleanup failure; inference and output-writing code are unchanged.
filesystem = (src / "src/filesystem.cpp").read_text()
prefix, suffix = filesystem.split("filesys::RmDir", 1)
suffix = suffix.replace("exit(1);", "return;")
(out / "filesystem.cpp").write_text(prefix + "filesys::RmDir" + suffix)
cppfiles = [
    out / "filesystem.cpp" if p.name == "filesystem.cpp" else p for p in cppfiles
]


def compile_one(p):
    obj = out / (p.stem + ".o")
    cmd = [
        str(cpp if p.suffix == ".cpp" else gcc),
        "-O2",
        "-c",
        str(p),
        "-o",
        str(obj),
        "-I" + str(out),
        "-I" + str(src / "src"),
        "-I" + str(src / "src/gzstream"),
        "-I" + str(zlib),
    ]
    if p.suffix == ".cpp":
        cmd += ["-std=c++14", "-include", "direct.h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    (out / (p.stem + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(p.name + "\n" + result.stderr[-8000:])
    print(p.name, flush=True)
    return str(obj)


with ThreadPoolExecutor(max_workers=4) as pool:
    objects = list(pool.map(compile_one, cfiles + cppfiles))
subprocess.run(
    [str(cpp), *objects, "-static", "-lpsapi", "-o", str(out / "Relate.exe")],
    check=True,
)
(out / "build.json").write_text(
    json.dumps(
        {
            "source_commit": subprocess.check_output(
                ["git", "-C", str(src.parent), "rev-parse", "HEAD"], text=True
            ).strip(),
            "compiler": subprocess.check_output([str(cpp), "--version"], text=True),
            "binary_sha256": hashlib.sha256(
                (out / "Relate.exe").read_bytes()
            ).hexdigest(),
            "shim": (out / "sys/resource.h").read_text(),
            "filesystem_cleanup_patch": "RmDir warns and returns on failure; leaves open scratch files intact on Windows",
        },
        indent=2,
    )
)

FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/gamma_smc_ts
COPY . .

# The same rootless installer used on HPC supplies pinned Python/SLiM plus the
# native compiler and libraries, then builds Gamma-SMC. Tests run in CI.
RUN chmod +x scripts/*.sh && bash scripts/bootstrap_uv.sh --skip-tests

ENV PATH="/opt/gamma_smc_ts/.venv/bin:/opt/gamma_smc_ts/.native/bin:${PATH}"
ENV SLIM_BIN="/opt/gamma_smc_ts/.native/bin/slim"
ENV GAMMA_SMC_BIN="/opt/gamma_smc_ts/bin/gamma_smc"
ENV LD_LIBRARY_PATH="/opt/gamma_smc_ts/.native/lib"

ENTRYPOINT ["/opt/gamma_smc_ts/scripts/aou.sh"]
CMD ["--help"]

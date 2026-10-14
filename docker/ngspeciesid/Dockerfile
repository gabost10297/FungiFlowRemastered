FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Medaka model matching the basecaller (dna_r10.4.1_e8.2_400bps_hac@v5.2.0)

RUN medaka tools download_models --models r1041_e82_400bps_hac_v5.2.0

WORKDIR /data
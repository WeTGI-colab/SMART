# ==============================================================================
# SMART — Somatic Mutation Annotation and Reporting Tool
# ==============================================================================
# Multi-stage build:
#   Stage 1: Pull VEP image
#   Stage 2: Pull GATK image
#   Stage 3: Final image with all tools combined
# ==============================================================================

ARG BUILD_PLATFORM=linux/amd64
ARG VEP_VERSION=114.0
ARG GATK_VERSION=4.6.0.0

FROM --platform=${BUILD_PLATFORM} ensemblorg/ensembl-vep:release_${VEP_VERSION} AS vep-base

FROM --platform=${BUILD_PLATFORM} broadinstitute/gatk:${GATK_VERSION} AS gatk-base

# ==============================================================================
# Final stage: combine everything into one image
# ==============================================================================
FROM --platform=${BUILD_PLATFORM} ubuntu:22.04

ARG VEP_VERSION=114.0
ARG GATK_VERSION=4.6.0.0
LABEL maintainer="SMART Pipeline Team"
LABEL description="SMART — Somatic Mutation Annotation and Reporting Tool"
LABEL smart.vep.version="${VEP_VERSION}"
LABEL smart.gatk.version="${GATK_VERSION}"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ---------- System dependencies ----------
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    wget \
    curl \
    gzip \
    bzip2 \
    tabix \
    bcftools \
    samtools \
    openjdk-17-jre-headless \
    python3 \
    python3-pip \
    python-is-python3 \
    perl \
    cpanminus \
    libdbi-perl \
    libdbd-mysql-perl \
    libwww-perl \
    libjson-perl \
    libmodule-build-perl \
    libarchive-zip-perl \
    libbio-perl-perl \
    libtry-tiny-perl \
    libset-intervaltree-perl \
    liblist-moreutils-perl \
    libbio-db-hts-perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---------- GATK (copy from gatk-base stage) ----------
COPY --from=gatk-base /gatk /opt/gatk
ENV PATH="/opt/gatk:${PATH}"

# ---------- VEP (copy from vep-base stage) ----------
COPY --from=vep-base /opt/vep /opt/vep
ENV PATH="/opt/vep/src/ensembl-vep:${PATH}"
ENV PERL5LIB=/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules

# ---------- Python dependencies ----------
# SpliceAI and a TensorFlow it works with. The released package uses tf.keras, which Keras 3
# (TF >= 2.16) removes, so the version is pinned rather than left to resolve. This is what makes
# '--spliceai compute' possible: the model scores against the transcripts we report, instead of
# looking up tables that carry no transcript and stop at +/-50 bp. See docs/spliceai_local.md.
RUN pip install --no-cache-dir "tensorflow-cpu==2.15.*" "spliceai==1.3.1"

COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# ---------- OncoKB Annotator ----------
RUN git clone https://github.com/oncokb/oncokb-annotator.git /opt/oncokb-annotator && \
    pip3 install --no-cache-dir -r /opt/oncokb-annotator/requirements/common.txt || true

# ---------- Pipeline scripts ----------
COPY scripts/ /opt/smart/scripts/
COPY reference_lists/svig_uk_canonical_variants.tsv reference_lists/svig_uk_vcep_population_thresholds.tsv /opt/smart/reference_lists/
RUN chmod +x /opt/smart/scripts/*.sh /opt/smart/scripts/*.py 2>/dev/null || true

# ---------- Working directory ----------
WORKDIR /output

# ---------- Entrypoint ----------
COPY VERSION /opt/smart/VERSION
COPY entrypoint.sh /opt/smart/entrypoint.sh
RUN chmod +x /opt/smart/entrypoint.sh

ENTRYPOINT ["/opt/smart/entrypoint.sh"]
CMD ["--help"]
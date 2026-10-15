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
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# SpliceAI and a TensorFlow it works with, installed last so that its constraints are the ones
# pip resolves against. The released package uses tf.keras, which Keras 3 (TF >= 2.16) removes,
# so the version is pinned. TF 2.15 is built against the NumPy 1.x ABI, which is why
# requirements.txt caps numpy below 2 -- installing this first and letting a later step upgrade
# numpy produces an image where 'import tensorflow' fails with 'numpy.core.umath failed to
# import', and the pipeline only discovers it at the point of scoring. This is what makes
# '--spliceai compute' possible. See docs/spliceai_local.md.
RUN pip3 install --no-cache-dir "tensorflow-cpu==2.15.*" "spliceai==1.3.1" && \
    python3 -c "import tensorflow, spliceai; print('tensorflow', tensorflow.__version__, 'ok')"

# ---------- OncoKB Annotator ----------
# Pinned. An unpinned clone made the annotator's version a function of the build date
# rather than of SMART's version, so two images carrying the same SMART version could
# annotate differently and neither was reproducible. Measured across the published
# images: smart:1.0.7 shipped annotator 4bfe720 (2026-06-25) while 1.0.8, 2.0.0 and
# 2.1.0 all shipped 77521cd (2026-08-25) -- a difference nothing recorded. It also makes
# a collaborator's bug report unanswerable, because the code they ran cannot be
# identified from the version they name (issue 32).
#
# Moving the pin is a deliberate act: the annotator decides how every variant is queried
# (see resolve_query_type_for_row, whose fallback to Genomic_Change issue 32 turns on),
# so a change there changes clinical output. Bump it with a measured run, not silently.
ARG ONCOKB_ANNOTATOR_COMMIT=77521cd2fc568422d87d84b73c24ca65b670a8e4
RUN git clone https://github.com/oncokb/oncokb-annotator.git /opt/oncokb-annotator && \
    git -C /opt/oncokb-annotator checkout --quiet "$ONCOKB_ANNOTATOR_COMMIT" && \
    [ "$(git -C /opt/oncokb-annotator rev-parse HEAD)" = "$ONCOKB_ANNOTATOR_COMMIT" ] && \
    echo "oncokb-annotator pinned at $ONCOKB_ANNOTATOR_COMMIT"
# Separate layer, because the `|| true` here is for the requirements file (optional
# upstream) and must not extend to the clone, the checkout or the pin assertion above:
# on one chain it would leave an image built from an unknown commit reporting success.
RUN pip3 install --no-cache-dir -r /opt/oncokb-annotator/requirements/common.txt || true

# Remove the annotator's ASCII fold. `append_annotation_to_file` ends with
#
#     rowstr = '\t'.join(row)
#     rowstr = rowstr.encode('ascii', 'ignore').decode('ascii')
#     outf.write(rowstr + "\n")
#
# so every row it writes loses any character outside ASCII, with errors='ignore' -- deleted,
# not replaced, and not reported. It rewrites every row of the MAF, so the loss is not limited
# to the columns OncoKB adds: ClinVar disease names, OncoKB summaries and drug names are all
# affected. `Diffuse_low-grade_glioma,_MAPK_pathway-altered` reaches the MAF as
# `MAPK_pathwayaltered`, which is not a term that exists and no longer matches its source.
# See issue 27.
#
# Safe to remove: the fold is the last operation before the write and feeds nothing. The API
# calls happen earlier, on the parsed row values, so nothing downstream in the annotator sees
# `rowstr`. The container's Python reports UTF-8 as its preferred encoding, so the write
# succeeds.
#
# Asserted in three steps, because a silent no-op would put us back to deleting characters
# without knowing: the target must be present before (so an upstream removal is noticed rather
# than assumed), absent after, and the file must still parse.
RUN f=/opt/oncokb-annotator/AnnotatorCore.py && \
    grep -q "rowstr.encode('ascii', 'ignore')" "$f" && \
    sed -i "s|^\( *\)rowstr = rowstr.encode('ascii', 'ignore').decode('ascii')|\1# SMART: ASCII fold removed, see issue 27 -- it deleted every non-ASCII character|" "$f" && \
    ! grep -q "rowstr.encode('ascii', 'ignore')" "$f" && \
    python3 -c "import ast; ast.parse(open('$f').read())" && \
    echo "AnnotatorCore.py: ASCII fold removed"

# ---------- Pipeline scripts ----------
COPY scripts/ /opt/smart/scripts/
# Every reference list the Tier 4 classifier resolves a default path for. Two were
# missing, and the consequence was silent rather than fatal: without
# svig_uk_functional_assays.tsv the O10/B6 assay route cannot load and the rule falls
# back to the OncoKB mutation-effect proxy capped at supporting (1 point instead of the
# assay table's 2), and without svig_uk_gene_specific_splice.tsv the TP53 exon 4/6/9
# synonymous provision the guideline requires never applies. Nothing in the output said
# the route had been unavailable, so a container run and a local run gave different
# classifications for the same variant with no visible reason.
COPY reference_lists/svig_uk_canonical_variants.tsv \
     reference_lists/svig_uk_vcep_population_thresholds.tsv \
     reference_lists/svig_uk_functional_assays.tsv \
     reference_lists/svig_uk_functional_assay_scores.tsv \
     reference_lists/svig_uk_gene_specific_splice.tsv \
     /opt/smart/reference_lists/
RUN chmod +x /opt/smart/scripts/*.sh /opt/smart/scripts/*.py 2>/dev/null || true

# ---------- Working directory ----------
WORKDIR /output

# ---------- Entrypoint ----------
COPY VERSION /opt/smart/VERSION
COPY entrypoint.sh /opt/smart/entrypoint.sh
RUN chmod +x /opt/smart/entrypoint.sh

ENTRYPOINT ["/opt/smart/entrypoint.sh"]
CMD ["--help"]
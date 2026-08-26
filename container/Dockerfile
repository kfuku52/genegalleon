# syntax=docker/dockerfile:1.7@sha256:a57df69d0ea827fb7266491f2813635de6f17269be881f696fbfdf2d83dda33e

FROM mambaorg/micromamba:noble@sha256:a0a4da2d5315f0d5b9d98196f2e9ecf50152c750bed4023bb127baade7e8faca AS system

USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/bin:${PATH} \
    R_LIBS_SITE=/opt/conda/lib/R/library

RUN rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /var/cache/apt/archives/partial/* \
 && apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y --no-install-recommends \
    augustus \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    default-jre-headless \
    fastp \
    gawk \
    gfortran \
    git \
    less \
    libarchive-tools \
    libboost-dev \
    libbpp-phyl-dev \
    libcurl4-openssl-dev \
    libeigen3-dev \
    libssl-dev \
    libxml2-dev \
    locales \
    parallel \
    pigz \
    perl \
    pkg-config \
    procps \
    rsync \
    tini \
    util-linux \
    wget \
    zip \
 && sed -i 's/^# *en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
 && locale-gen en_US.UTF-8 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Isolate the Notung download so Git source revisions in the runtime stage do
# not invalidate it.
FROM system AS notung

ARG NOTUNG_DOWNLOAD_PAGE="https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip"
ARG NOTUNG_DOWNLOAD_HOST_IP="128.2.205.60"
ARG NOTUNG_ZIP_SHA256="81cbff670ab4d2416c01eba503f81c454aa5a724b0982373dd17510113882ae6"

COPY container/scripts/resolve_notung_latest.sh /opt/pg/scripts/resolve_notung_latest.sh
COPY container/scripts/download_url.sh /opt/pg/scripts/download_url.sh
COPY container/scripts/download_url_with_fallback.sh /opt/pg/scripts/download_url_with_fallback.sh
COPY container/scripts/verify_sha256.sh /opt/pg/scripts/verify_sha256.sh
COPY container/scripts/extract_notung_jar.py /opt/pg/scripts/extract_notung_jar.py

RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 \
 && chmod +x /opt/pg/scripts/*.sh \
 && notung_url="$(/opt/pg/scripts/resolve_notung_latest.sh "${NOTUNG_DOWNLOAD_PAGE}")" \
 && tmp_zip="$(mktemp /tmp/notung.XXXXXX.zip)" \
 && NOTUNG_DOWNLOAD_HOST_IP="${NOTUNG_DOWNLOAD_HOST_IP}" /opt/pg/scripts/download_url_with_fallback.sh "${notung_url}" "${tmp_zip}" \
 && /opt/pg/scripts/verify_sha256.sh "${tmp_zip}" "${NOTUNG_ZIP_SHA256}" \
 && python3 /opt/pg/scripts/extract_notung_jar.py --archive "${tmp_zip}" --destination /Notung.jar \
 && rm -f "${tmp_zip}" \
 && apt-get purge -y python3 \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

FROM system AS runtime

ARG TARGETARCH
ARG TARGETPLATFORM

# Resolve the complete conda environment in one transaction. Keeping source
# revisions and validation files below this layer preserves the expensive
# solver/install cache when those fast-changing inputs are updated.
COPY container/env /opt/pg/env
COPY container/scripts/install_env.sh /opt/pg/scripts/install_env.sh

RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.mamba/pkgs,sharing=locked \
 chmod +x /opt/pg/scripts/install_env.sh \
 && export USE_DEFAULTS_CHANNEL=0 STRICT_CHANNEL_PRIORITY=1 MAMBA_NO_LOW_SPEED_LIMIT=1 \
 && base_required="/opt/pg/env/base.required.txt" \
 && optional_file="/opt/pg/env/base.${TARGETARCH}.optional.txt" \
 && reconciled_required="/tmp/pg.base.required.txt" \
 && if [[ -f "/opt/pg/env/base.${TARGETARCH}.required.txt" ]]; then base_required="/opt/pg/env/base.${TARGETARCH}.required.txt"; fi \
 && if [[ ! -f "${optional_file}" ]]; then optional_file="/dev/null"; fi \
 && sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "${base_required}" /opt/pg/env/base.r.required.txt | awk '!seen[$0]++' > "${reconciled_required}" \
 && /opt/pg/scripts/install_env.sh base "${reconciled_required}" "${optional_file}" python=3.12 \
 && rm -f "${reconciled_required}"

COPY container/testdata /opt/pg/testdata
COPY container/scripts/ensure_generax_stable.sh /opt/pg/scripts/ensure_generax_stable.sh

RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.mamba/pkgs,sharing=locked \
 chmod +x /opt/pg/scripts/ensure_generax_stable.sh \
 && /opt/pg/scripts/ensure_generax_stable.sh base

ARG PAML_REPO_URL="https://github.com/iqtree/paml.git"
ARG PAML_REPO_REF=""
ARG PAML_REPO_SHA=""

COPY container/source_branches.env /opt/pg/source_branches.env
COPY container/scripts/ensure_iq2mc_toolchain.sh /opt/pg/scripts/ensure_iq2mc_toolchain.sh
COPY container/scripts/download_url.sh /opt/pg/scripts/download_url.sh
COPY container/scripts/verify_sha256.sh /opt/pg/scripts/verify_sha256.sh
COPY container/scripts/fetch_git_repo.sh /opt/pg/scripts/fetch_git_repo.sh

RUN chmod +x /opt/pg/scripts/*.sh \
 && source /opt/pg/source_branches.env \
 && paml_repo_ref="${PAML_REPO_REF:-${GG_SOURCE_PAML_REPO_REF}}" \
 && PAML_REPO_URL="${PAML_REPO_URL}" \
    PAML_REPO_REF="${paml_repo_ref}" \
    PAML_REPO_SHA="${PAML_REPO_SHA}" \
    /opt/pg/scripts/ensure_iq2mc_toolchain.sh base

# Install Python-only compatibility/security updates once, after the conda
# solver has finished, so wheels cannot pull the scientific stack away from
# the NumPy range required by both PyMOL and ClipKIT.
COPY container/pip-compatibility.requirements.txt /opt/pg/pip-compatibility.requirements.txt

RUN --mount=type=cache,target=/root/.cache/pip \
 micromamba run -n base pip install \
      -r /opt/pg/pip-compatibility.requirements.txt \
 && micromamba run -n base python -c "import ete4, genomepy, lxml, ortools, PIL" \
 && micromamba run -n base python -m pip check

# Native tools and CRAN fallbacks change much less frequently than the Git
# sources below, so give them an independent reusable layer.
ARG TESTNH_TARBALL_SHA256="598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2"
ARG CAFE5_TARBALL_SHA256="71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8"

COPY container/scripts/install_nonconda_fallbacks.sh /opt/pg/scripts/install_nonconda_fallbacks.sh

RUN --mount=type=cache,target=/root/.cache/pip \
 chmod +x /opt/pg/scripts/install_nonconda_fallbacks.sh \
 && TESTNH_TARBALL_SHA256="${TESTNH_TARBALL_SHA256}" \
 CAFE5_TARBALL_SHA256="${CAFE5_TARBALL_SHA256}" \
 /opt/pg/scripts/install_nonconda_fallbacks.sh cafe5 mapnh astral r-cran verify-r

COPY --from=notung /Notung.jar /usr/local/bin/Notung.jar

RUN mkdir -p /home /usr/local/db/Pfam_LE /usr/local/db/jaspar /opt/pg/logs /opt/conda/envs \
 && ln -sfn /opt/conda /opt/conda/envs/base \
 && printf '%s\n' \
   'export MAMBA_ROOT_PREFIX=/opt/conda' \
   'eval "$(micromamba shell hook --shell bash)"' \
   'conda() { micromamba "$@"; }' \
   'export PATH=/usr/local/bin:${PATH}' \
   > /home/.bashrc \
 && cp /home/.bashrc /root/.bashrc

# Declare mutable source inputs only immediately before their install layer.
# BuildKit includes in-scope ARGs in RUN cache keys even when a command does
# not reference them, so early declarations would invalidate the conda layer.
ARG KFU52_REPO_REF=""
ARG KFU52_AMALGKIT_AUTO_SELECT_REF="0"
ARG KFU52_AMALGKIT_BRANCH_CANDIDATES="master,kfdevel,devel"
ARG KFU52_AMALGKIT_REPO_REF=""
ARG KFU52_AMALGKIT_REPO_SHA=""
ARG KFU52_CDSKIT_REPO_REF=""
ARG KFU52_CDSKIT_REPO_SHA=""
ARG KFU52_CSUBST_REPO_REF=""
ARG KFU52_CSUBST_REPO_SHA=""
ARG KFU52_NWKIT_REPO_REF=""
ARG KFU52_NWKIT_REPO_SHA=""
ARG BUSCO_REPO_URL="https://gitlab.com/ezlab/busco.git"
ARG BUSCO_MIRROR_REPO_URL=""
ARG BUSCO_REPO_REF=""
ARG BUSCO_REPO_SHA=""
ARG KFL1OU_REPO_URL="https://github.com/kfuku52/kfl1ou.git"
ARG KFL1OU_REPO_REF=""
ARG KFL1OU_REPO_SHA=""
ARG KFFRACTBIAS_REPO_URL="https://github.com/kfuku52/kfFractBias.git"
ARG KFFRACTBIAS_REPO_REF=""
ARG KFFRACTBIAS_REPO_SHA=""
ARG KFTOOLS_REPO_URL="https://github.com/kfuku52/kftools.git"
ARG RKFTOOLS_REPO_URL="https://github.com/kfuku52/rkftools.git"
ARG RADTE_REPO_URL="https://github.com/kfuku52/RADTE.git"
ARG KFTOOLS_REPO_REF=""
ARG RKFTOOLS_REPO_REF=""
ARG RADTE_REPO_REF=""
ARG KFTOOLS_REPO_SHA=""
ARG RKFTOOLS_REPO_SHA=""
ARG RADTE_REPO_SHA=""

COPY container/scripts/ensure_kfuku52_tools_latest.sh /opt/pg/scripts/ensure_kfuku52_tools_latest.sh
COPY container/scripts/install_helvetica_compat_font.sh /opt/pg/scripts/install_helvetica_compat_font.sh

RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.mamba/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.cache/pip \
 chmod +x /opt/pg/scripts/*.sh \
 && export USE_DEFAULTS_CHANNEL=0 STRICT_CHANNEL_PRIORITY=1 MAMBA_NO_LOW_SPEED_LIMIT=1 \
 && source /opt/pg/source_branches.env \
 && mkdir -p /opt/pg/logs \
 && kfu52_repo_ref="${KFU52_REPO_REF:-${GG_SOURCE_NWKIT_REPO_REF}}" \
 && amalgkit_repo_ref="${KFU52_AMALGKIT_REPO_REF:-${GG_SOURCE_AMALGKIT_REPO_REF}}" \
 && cdskit_repo_ref="${KFU52_CDSKIT_REPO_REF:-${GG_SOURCE_CDSKIT_REPO_REF}}" \
 && csubst_repo_ref="${KFU52_CSUBST_REPO_REF:-${GG_SOURCE_CSUBST_REPO_REF}}" \
 && nwkit_repo_ref="${KFU52_NWKIT_REPO_REF:-${GG_SOURCE_NWKIT_REPO_REF}}" \
 && busco_repo_ref="${BUSCO_REPO_REF:-${GG_SOURCE_BUSCO_REPO_REF}}" \
 && paml_repo_ref="${PAML_REPO_REF:-${GG_SOURCE_PAML_REPO_REF}}" \
 && kfl1ou_repo_ref="${KFL1OU_REPO_REF:-${GG_SOURCE_KFL1OU_REPO_REF}}" \
 && kffractbias_repo_ref="${KFFRACTBIAS_REPO_REF:-${GG_SOURCE_KFFRACTBIAS_REPO_REF}}" \
 && kftools_repo_ref="${KFTOOLS_REPO_REF:-${GG_SOURCE_KFTOOLS_REPO_REF}}" \
 && rkftools_repo_ref="${RKFTOOLS_REPO_REF:-${GG_SOURCE_RKFTOOLS_REPO_REF}}" \
 && radte_repo_ref="${RADTE_REPO_REF:-${GG_SOURCE_RADTE_REPO_REF}}" \
 && amalgkit_repo_sha="${KFU52_AMALGKIT_REPO_SHA:-}" \
 && cdskit_repo_sha="${KFU52_CDSKIT_REPO_SHA:-}" \
 && csubst_repo_sha="${KFU52_CSUBST_REPO_SHA:-}" \
 && nwkit_repo_sha="${KFU52_NWKIT_REPO_SHA:-}" \
 && kfl1ou_repo_sha="${KFL1OU_REPO_SHA:-}" \
 && kffractbias_repo_sha="${KFFRACTBIAS_REPO_SHA:-}" \
 && kftools_repo_sha="${KFTOOLS_REPO_SHA:-}" \
 && rkftools_repo_sha="${RKFTOOLS_REPO_SHA:-}" \
 && radte_repo_sha="${RADTE_REPO_SHA:-}" \
 && amalgkit_source_ref="${amalgkit_repo_sha:-${amalgkit_repo_ref}}" \
 && cdskit_source_ref="${cdskit_repo_sha:-${cdskit_repo_ref}}" \
 && csubst_source_ref="${csubst_repo_sha:-${csubst_repo_ref}}" \
 && nwkit_source_ref="${nwkit_repo_sha:-${nwkit_repo_ref}}" \
 && kfl1ou_source_ref="${kfl1ou_repo_sha:-${kfl1ou_repo_ref}}" \
 && kffractbias_source_ref="${kffractbias_repo_sha:-${kffractbias_repo_ref}}" \
 && busco_source_ref="${BUSCO_REPO_SHA:-${busco_repo_ref}}" \
 && paml_source_ref="${PAML_REPO_SHA:-${paml_repo_ref}}" \
 && kftools_source_ref="${kftools_repo_sha:-${kftools_repo_ref}}" \
 && rkftools_source_ref="${rkftools_repo_sha:-${rkftools_repo_ref}}" \
 && radte_source_ref="${radte_repo_sha:-${radte_repo_ref}}" \
 && KFU52_REPO_REF="${kfu52_repo_ref}" \
    KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    KFU52_AMALGKIT_REPO_REF="${amalgkit_repo_ref}" \
    KFU52_AMALGKIT_REPO_SHA="${amalgkit_repo_sha}" \
    KFU52_CDSKIT_REPO_REF="${cdskit_repo_ref}" \
    KFU52_CDSKIT_REPO_SHA="${cdskit_repo_sha}" \
    KFU52_CSUBST_REPO_REF="${csubst_repo_ref}" \
    KFU52_CSUBST_REPO_SHA="${csubst_repo_sha}" \
    KFU52_NWKIT_REPO_REF="${nwkit_repo_ref}" \
    KFU52_NWKIT_REPO_SHA="${nwkit_repo_sha}" \
    /opt/pg/scripts/ensure_kfuku52_tools_latest.sh base \
 && KFU52_REPO_REF="${KFU52_REPO_REF}" \
    KFL1OU_REPO_URL="${KFL1OU_REPO_URL}" \
    KFL1OU_REPO_REF="${kfl1ou_repo_ref}" \
    KFL1OU_REPO_SHA="${kfl1ou_repo_sha}" \
    /opt/pg/scripts/install_nonconda_fallbacks.sh kfl1ou \
 && rm -rf /tmp/busco \
 && /opt/pg/scripts/fetch_git_repo.sh "${BUSCO_REPO_URL}" "${busco_source_ref}" /tmp/busco "${BUSCO_MIRROR_REPO_URL}" \
 && micromamba run -n base pip install --upgrade --force-reinstall --no-deps --no-build-isolation /tmp/busco \
 && rm -rf /tmp/busco \
 && rm -rf /tmp/kftools \
 && /opt/pg/scripts/fetch_git_repo.sh "${KFTOOLS_REPO_URL}" "${kftools_source_ref}" /tmp/kftools \
 && micromamba run -n base pip install --upgrade --force-reinstall --no-deps --no-build-isolation /tmp/kftools \
 && rm -rf /tmp/kftools \
 && rm -rf /tmp/kffractbias \
 && /opt/pg/scripts/fetch_git_repo.sh "${KFFRACTBIAS_REPO_URL}" "${kffractbias_source_ref}" /tmp/kffractbias \
 && micromamba run -n base pip install --upgrade --force-reinstall --no-deps --no-build-isolation /tmp/kffractbias \
 && rm -rf /tmp/kffractbias \
 && /opt/pg/scripts/install_helvetica_compat_font.sh base \
 && rm -rf /tmp/rkftools \
 && /opt/pg/scripts/fetch_git_repo.sh "${RKFTOOLS_REPO_URL}" "${rkftools_source_ref}" /tmp/rkftools \
 && micromamba run -n base R CMD INSTALL /tmp/rkftools \
 && rm -rf /tmp/rkftools \
 && rm -rf /tmp/radte \
 && /opt/pg/scripts/fetch_git_repo.sh "${RADTE_REPO_URL}" "${radte_source_ref}" /tmp/radte \
 && if [[ ! -s /tmp/radte/radte.r ]]; then \
      echo "RADTE script was not found in ${RADTE_REPO_URL}"; \
      exit 1; \
    fi \
 && install -m 0755 /tmp/radte/radte.r /usr/local/bin/radte.r \
 && rm -rf /tmp/radte \
 && printf '%s\t%s\n' \
      source revision \
      amalgkit "${amalgkit_source_ref}" \
      cdskit "${cdskit_source_ref}" \
      csubst "${csubst_source_ref}" \
      nwkit "${nwkit_source_ref}" \
      BUSCO "${busco_source_ref}" \
      paml "${paml_source_ref}" \
      kfl1ou "${kfl1ou_source_ref}" \
      kfFractBias "${kffractbias_source_ref}" \
      kftools "${kftools_source_ref}" \
      rkftools "${rkftools_source_ref}" \
      RADTE "${radte_source_ref}" \
      > /opt/pg/logs/source_revisions.tsv

# Install GeneGalleon's repository-owned R package after its external R
# dependencies. Keeping this in its own layer makes treevis-only edits cheap to
# validate while ensuring every workflow consumer sees the packaged API.
COPY workflow/support/treevis /opt/pg/src/genegalleon.treevis
RUN micromamba run -n base R CMD INSTALL /opt/pg/src/genegalleon.treevis \
 && micromamba run -n base Rscript -e \
      "stopifnot(requireNamespace('genegalleon.treevis', quietly = TRUE))"

# Validation is intentionally last and isolated: editing checks must not force
# package solving or source builds to run again.
ARG SECURITY_REFRESH_EPOCH=""

# The system-package layer near the start of this image is intentionally highly
# cacheable. Refresh the Ubuntu security index again in a small late layer whose
# daily epoch participates in the shared build-input fingerprint.
RUN [[ "${SECURITY_REFRESH_EPOCH}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] \
 || { echo "SECURITY_REFRESH_EPOCH must use YYYY-MM-DD format" >&2; exit 1; } \
 && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /var/cache/apt/archives/partial/* \
 && apt-get update \
 && apt-get upgrade --with-new-pkgs -y \
 && apt-get --simulate dist-upgrade > /tmp/gg-apt-upgrade-simulation.txt \
 && if awk '/^Inst / { found = 1 } END { exit found ? 0 : 1 }' /tmp/gg-apt-upgrade-simulation.txt; then \
      echo "Upgradeable system packages remain after the security refresh." >&2; \
      exit 1; \
    fi \
 && rm -f /tmp/gg-apt-upgrade-simulation.txt \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /var/cache/apt/archives/partial/* \
 && printf '%s\n' "${SECURITY_REFRESH_EPOCH}" > /opt/pg/logs/security_refresh_epoch.txt

COPY container/spec /opt/pg/spec
COPY container/scripts/validate_runtime.sh /opt/pg/scripts/validate_runtime.sh
COPY container/scripts/post_install_notes.sh /opt/pg/scripts/post_install_notes.sh

RUN chmod +x /opt/pg/scripts/validate_runtime.sh \
 && runtime_required="/opt/pg/spec/required_commands.tsv" \
 && runtime_optional="/opt/pg/spec/optional_commands.${TARGETARCH}.tsv" \
 && if [[ -f "/opt/pg/spec/required_commands.${TARGETARCH}.tsv" ]]; then runtime_required="/opt/pg/spec/required_commands.${TARGETARCH}.tsv"; fi \
 && /opt/pg/scripts/validate_runtime.sh \
      --required "${runtime_required}" \
      --optional "${runtime_optional}" \
      --report /opt/pg/logs/runtime_validation_${TARGETARCH}.tsv \
 && micromamba run -n base python -m pip check \
 && micromamba clean -i -l -y \
 && bash /opt/pg/scripts/post_install_notes.sh "${TARGETARCH}"

ARG BUILD_DATE=""
ARG VCS_REF="unknown"
ARG GG_VERSION="unknown"
ARG BUILD_INPUT_HASH=""

LABEL org.opencontainers.image.title="GeneGalleon workflow container" \
      org.opencontainers.image.description="Reproducible multi-architecture runtime for GeneGalleon (amd64 + arm64)" \
      org.opencontainers.image.source="https://github.com/kfuku52/genegalleon" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.version="${GG_VERSION}" \
      io.genegalleon.build-input="${BUILD_INPUT_HASH}" \
      io.genegalleon.security-refresh-epoch="${SECURITY_REFRESH_EPOCH}"

WORKDIR /workspace
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]

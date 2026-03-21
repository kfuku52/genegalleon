# syntax=docker/dockerfile:1.7

FROM mambaorg/micromamba:noble

ARG TARGETARCH
ARG TARGETPLATFORM
ARG NOTUNG_DOWNLOAD_PAGE="https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip"
ARG NOTUNG_DOWNLOAD_HOST_IP="128.2.205.60"
ARG KFU52_REPO_REF="master"
ARG KFU52_AMALGKIT_AUTO_SELECT_REF="1"
ARG KFU52_AMALGKIT_BRANCH_CANDIDATES="master,kfdevel,devel"
ARG KFU52_AMALGKIT_REPO_REF=""
ARG KFU52_AMALGKIT_REPO_SHA=""
ARG KFU52_CDSKIT_REPO_SHA=""
ARG KFU52_CSUBST_REPO_SHA=""
ARG KFU52_NWKIT_REPO_SHA=""
ARG BUSCO_REPO_URL="https://gitlab.com/ezlab/busco.git"
ARG BUSCO_MIRROR_REPO_URL=""
ARG BUSCO_REPO_REF="6.0.0"
ARG BUSCO_REPO_SHA="6278721a1916f6da310e03ec9674099028c927a4"
ARG PAML_REPO_URL="https://github.com/iqtree/paml.git"
ARG PAML_REPO_REF="master"
ARG PAML_REPO_SHA="8daeead6b55523f375d9ac56dcfac38373ef8a2e"
ARG KFL1OU_REPO_URL="https://github.com/kfuku52/kfl1ou.git"
ARG KFL1OU_REPO_REF=""
ARG KFL1OU_REPO_SHA=""
ARG KFTOOLS_REPO_URL="https://github.com/kfuku52/kftools.git"
ARG RKFTOOLS_REPO_URL="https://github.com/kfuku52/rkftools.git"
ARG RADTE_REPO_URL="https://github.com/kfuku52/RADTE.git"
ARG KFTOOLS_REPO_REF=""
ARG RKFTOOLS_REPO_REF=""
ARG RADTE_REPO_REF=""
ARG KFTOOLS_REPO_SHA=""
ARG RKFTOOLS_REPO_SHA=""
ARG RADTE_REPO_SHA=""
ARG TESTNH_TARBALL_SHA256="598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2"
ARG CAFE5_TARBALL_SHA256="71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8"

USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/bin:${PATH} \
    R_LIBS_SITE=/opt/conda/lib/R/library

RUN rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /var/cache/apt/archives/partial/* \
 && apt-get update \
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
    libboost-dev \
    libbpp-phyl-dev \
    libcurl4-openssl-dev \
    libeigen3-dev \
    libssl-dev \
    libxml2-dev \
    parallel \
    pigz \
    perl \
    pkg-config \
    procps \
    rsync \
    tini \
    unzip \
    wget \
    zip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY container/env /opt/pg/env
COPY container/spec /opt/pg/spec
COPY container/testdata /opt/pg/testdata
COPY container/scripts/install_env.sh /opt/pg/scripts/install_env.sh
COPY container/scripts/ensure_generax_stable.sh /opt/pg/scripts/ensure_generax_stable.sh
COPY container/scripts/ensure_kfuku52_tools_latest.sh /opt/pg/scripts/ensure_kfuku52_tools_latest.sh
COPY container/scripts/ensure_iq2mc_toolchain.sh /opt/pg/scripts/ensure_iq2mc_toolchain.sh
COPY container/scripts/resolve_notung_latest.sh /opt/pg/scripts/resolve_notung_latest.sh
COPY container/scripts/download_url.sh /opt/pg/scripts/download_url.sh
COPY container/scripts/download_url_with_fallback.sh /opt/pg/scripts/download_url_with_fallback.sh
COPY container/scripts/verify_sha256.sh /opt/pg/scripts/verify_sha256.sh
COPY container/scripts/fetch_git_repo.sh /opt/pg/scripts/fetch_git_repo.sh
COPY container/scripts/install_helvetica_compat_font.sh /opt/pg/scripts/install_helvetica_compat_font.sh
COPY container/scripts/install_nonconda_fallbacks.sh /opt/pg/scripts/install_nonconda_fallbacks.sh
COPY container/scripts/validate_runtime.sh /opt/pg/scripts/validate_runtime.sh

# nwkit currently imports ete4 modules directly, so install ete4 explicitly
# via pip and assert importability during build.
# Keep this layer focused on architecture/base toolchain setup so it can be
# reused even when Git refs change.
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.mamba/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.cache/pip \
 chmod +x /opt/pg/scripts/*.sh \
 && export USE_DEFAULTS_CHANNEL=0 STRICT_CHANNEL_PRIORITY=1 MAMBA_NO_LOW_SPEED_LIMIT=1 \
 && base_required="/opt/pg/env/base.required.txt" \
 && reconciled_required="/tmp/pg.base.required.txt" \
 && if [[ -f "/opt/pg/env/base.${TARGETARCH}.required.txt" ]]; then base_required="/opt/pg/env/base.${TARGETARCH}.required.txt"; fi \
 && sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "${base_required}" /opt/pg/env/base.r.required.txt | awk '!seen[$0]++' > "${reconciled_required}" \
 && /opt/pg/scripts/install_env.sh base "${base_required}" /opt/pg/env/base.${TARGETARCH}.optional.txt python=3.12 \
 && /opt/pg/scripts/ensure_generax_stable.sh base \
 && micromamba install -y --retry-clean-cache --repodata-ttl 3600 -n base -c conda-forge -c bioconda pymol-open-source \
 && micromamba run -n base pip install --no-cache-dir ete4 \
 && micromamba run -n base python -c "import ete4" \
 && PAML_REPO_URL="${PAML_REPO_URL}" \
    PAML_REPO_REF="${PAML_REPO_REF}" \
    PAML_REPO_SHA="${PAML_REPO_SHA}" \
    /opt/pg/scripts/ensure_iq2mc_toolchain.sh base \
 && /opt/pg/scripts/install_env.sh base "${reconciled_required}" /opt/pg/env/base.r.optional.txt r-base \
 && kfl1ou_repo_ref="${KFL1OU_REPO_REF:-${KFU52_REPO_REF}}" \
 && TESTNH_TARBALL_SHA256="${TESTNH_TARBALL_SHA256}" \
    CAFE5_TARBALL_SHA256="${CAFE5_TARBALL_SHA256}" \
    KFU52_REPO_REF="${KFU52_REPO_REF}" \
    KFL1OU_REPO_URL="${KFL1OU_REPO_URL}" \
    KFL1OU_REPO_REF="${kfl1ou_repo_ref}" \
    KFL1OU_REPO_SHA="${KFL1OU_REPO_SHA}" \
    /opt/pg/scripts/install_nonconda_fallbacks.sh \
 && rm -f "${reconciled_required}"

# Keep fast-changing GitHub-derived installs separate to reduce rebuild scope.
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.mamba/pkgs,sharing=locked \
 --mount=type=cache,target=/root/.cache/pip \
 export USE_DEFAULTS_CHANNEL=0 STRICT_CHANNEL_PRIORITY=1 MAMBA_NO_LOW_SPEED_LIMIT=1 \
 && runtime_required="/opt/pg/spec/required_commands.tsv" \
 && runtime_optional="/opt/pg/spec/optional_commands.${TARGETARCH}.tsv" \
 && if [[ -f "/opt/pg/spec/required_commands.${TARGETARCH}.tsv" ]]; then runtime_required="/opt/pg/spec/required_commands.${TARGETARCH}.tsv"; fi \
 && kftools_repo_ref="${KFTOOLS_REPO_REF:-${KFU52_REPO_REF}}" \
 && rkftools_repo_ref="${RKFTOOLS_REPO_REF:-${KFU52_REPO_REF}}" \
 && radte_repo_ref="${RADTE_REPO_REF:-${KFU52_REPO_REF}}" \
 && busco_source_ref="${BUSCO_REPO_SHA:-${BUSCO_REPO_REF}}" \
 && kftools_source_ref="${KFTOOLS_REPO_SHA:-${kftools_repo_ref}}" \
 && rkftools_source_ref="${RKFTOOLS_REPO_SHA:-${rkftools_repo_ref}}" \
 && radte_source_ref="${RADTE_REPO_SHA:-${radte_repo_ref}}" \
 && KFU52_REPO_REF="${KFU52_REPO_REF}" \
    KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    KFU52_AMALGKIT_REPO_SHA="${KFU52_AMALGKIT_REPO_SHA}" \
    KFU52_CDSKIT_REPO_SHA="${KFU52_CDSKIT_REPO_SHA}" \
    KFU52_CSUBST_REPO_SHA="${KFU52_CSUBST_REPO_SHA}" \
    KFU52_NWKIT_REPO_SHA="${KFU52_NWKIT_REPO_SHA}" \
    /opt/pg/scripts/ensure_kfuku52_tools_latest.sh base \
 && rm -rf /tmp/busco \
 && /opt/pg/scripts/fetch_git_repo.sh "${BUSCO_REPO_URL}" "${busco_source_ref}" /tmp/busco "${BUSCO_MIRROR_REPO_URL}" \
 && micromamba run -n base pip install --no-cache-dir --upgrade --force-reinstall --no-deps --no-build-isolation /tmp/busco \
 && rm -rf /tmp/busco \
 && rm -rf /tmp/kftools \
 && /opt/pg/scripts/fetch_git_repo.sh "${KFTOOLS_REPO_URL}" "${kftools_source_ref}" /tmp/kftools \
 && micromamba run -n base pip install --no-cache-dir --upgrade --force-reinstall --no-deps --no-build-isolation /tmp/kftools \
 && rm -rf /tmp/kftools \
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
 && /opt/pg/scripts/validate_runtime.sh \
      --required "${runtime_required}" \
      --optional "${runtime_optional}" \
      --report /opt/pg/logs/runtime_validation_${TARGETARCH}.tsv \
 && micromamba clean -i -l -y

COPY container/scripts/post_install_notes.sh /opt/pg/scripts/post_install_notes.sh

RUN mkdir -p /home /usr/local/db/Pfam_LE /usr/local/db/jaspar /opt/pg/logs /opt/conda/envs \
 && ln -sfn /opt/conda /opt/conda/envs/base \
 && printf '%s\n' \
   'export MAMBA_ROOT_PREFIX=/opt/conda' \
   'eval "$(micromamba shell hook --shell bash)"' \
   'conda() { micromamba "$@"; }' \
   'export PATH=/usr/local/bin:${PATH}' \
   > /home/.bashrc \
 && cp /home/.bashrc /root/.bashrc

# Notung is Java-based, so this works on both amd64/arm64.
RUN notung_url="$(/opt/pg/scripts/resolve_notung_latest.sh "${NOTUNG_DOWNLOAD_PAGE}")" \
 && tmp_zip="$(mktemp /tmp/notung.XXXXXX.zip)" \
 && tmp_dir="$(mktemp -d /tmp/notung.XXXXXX)" \
 && NOTUNG_DOWNLOAD_HOST_IP="${NOTUNG_DOWNLOAD_HOST_IP}" /opt/pg/scripts/download_url_with_fallback.sh "${notung_url}" "${tmp_zip}" \
 && unzip -q "${tmp_zip}" -d "${tmp_dir}" \
 && notung_jar="$(find "${tmp_dir}" -type f -name 'Notung-2.9*.jar' | sort -V | tail -n 1)" \
 && if [[ -z "${notung_jar}" ]]; then echo "Notung jar was not found in ${notung_url}"; exit 1; fi \
 && cp "${notung_jar}" /usr/local/bin/Notung.jar \
 && chmod 0755 /usr/local/bin/Notung.jar \
 && rm -f "${tmp_zip}" \
 && rm -rf "${tmp_dir}" \
 && bash /opt/pg/scripts/post_install_notes.sh "${TARGETARCH}"

ARG BUILD_DATE=""
ARG VCS_REF="unknown"
ARG GG_VERSION="unknown"

LABEL org.opencontainers.image.title="GeneGalleon multi-arch candidate container" \
      org.opencontainers.image.description="Reproducible multi-arch build scaffold for GeneGalleon (amd64 + arm64)" \
      org.opencontainers.image.source="https://github.com/kfuku52/genegalleon" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.version="${GG_VERSION}"

WORKDIR /workspace
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]

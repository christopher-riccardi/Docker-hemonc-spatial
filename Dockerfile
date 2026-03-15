# Dockerfile for hemonc-spatial
FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL maintainer="christopher-riccardi"
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# --- system packages required by the requested R packages ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    gcc \
    gfortran \
    cmake \
    pkg-config \
    wget \
    curl \
    git \
    ca-certificates \
    locales \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libsqlite3-dev libhdf5-serial-dev \
    libgdal-dev libgeos-dev \
    libopenblas-dev liblapack-dev libatlas-base-dev libgsl-dev \
    default-jdk \
    python3-dev python3-pip \
    && rm -rf /var/lib/apt/lists/*

# locale
RUN locale-gen en_US.UTF-8 || true
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# copy and run package installer
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm -f /tmp/install_packages.R

EXPOSE 8787
CMD ["/init"]

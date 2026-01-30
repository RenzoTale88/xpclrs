# Use the official NVIDIA CUDA image as the base image
FROM ubuntu:18.04

# Copy source code
COPY . /xpclrs

# Set up the environment
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libopenblas-dev curl build-essential clang libclang-dev ca-certificates && \
    rm -rf /var/lib/apt/lists/* && \
    apt autoclean -y && apt autoremove -y && \
    update-ca-certificates

# Install Rust and Cargo
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Add Cargo to PATH
ENV PATH="/root/.cargo/bin:${PATH}"

# Compile the tool
WORKDIR /xpclrs
RUN /root/.cargo/bin/cargo build --release &&\
    cp target/release/xpclrs /usr/local/bin/xpclrs &&\
    rm -rf /xpclrs
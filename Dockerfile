FROM fedora:34
LABEL maintainer="Rodrigo Neumann Barros Ferreira <rneumann@br.ibm.com>"
EXPOSE 9080

# Create environment variables and folders
ARG USER_ID
ENV ROOT_DIR /root/simulator
ENV PYTHONPATH ${ROOT_DIR}/api
RUN mkdir -p ${ROOT_DIR}/results
WORKDIR ${ROOT_DIR}
VOLUME ${ROOT_DIR}/results

# Installing and upgrading OS dependencies
RUN dnf --best --allowerasing -y install \
  armadillo-devel \
  cmake \
  findutils \
  gcc-c++ \
  glog-devel \
  gtest-devel \
  lcov \
  libffi-devel \
  ninja-build \
  pipenv \
  pkgconfig \
  rapidjson \
  redhat-rpm-config \
  sundials-devel \
  tclap \
  gperftools \
  && dnf clean all

# Installing and upgrading API dependencies
COPY Pipfile ${ROOT_DIR}/
RUN pipenv install --system --skip-lock --clear

# Copying the source code (in increasing order of update frequency)
COPY build          ${ROOT_DIR}/build
COPY bin            ${ROOT_DIR}/bin
COPY docs           ${ROOT_DIR}/docs
COPY cmake          ${ROOT_DIR}/cmake
COPY util           ${ROOT_DIR}/util
COPY CMakeLists.txt ${ROOT_DIR}/
COPY test           ${ROOT_DIR}/test
COPY src            ${ROOT_DIR}/src
COPY gmm            ${ROOT_DIR}/gmm

# Build flow-simulator.x in RELEASE mode
RUN cmake -GNinja -DCMAKE_BUILD_TYPE=Release -B${ROOT_DIR}/build -H${ROOT_DIR} && \
    cmake --build ${ROOT_DIR}/build --target flow-simulator.x

# Build flow-simulator-test.x in DEBUG mode
RUN cmake -GNinja -DCMAKE_BUILD_TYPE=Debug -B${ROOT_DIR}/build -H${ROOT_DIR} && \
    cmake --build ${ROOT_DIR}/build --target flow-simulator-test.x

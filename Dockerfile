FROM chemiskyy/smartplus-stack

# install smartplus
RUN export SMART_URL=https://github.com/smartplus-team/smartplus.git && \
    export SMART_NAME=smartplus && \
    export SMART_SOURCE_DIR=${PREFIX}/source/smartplus && \
    export SMART_BUILD_DIR=${PREFIX}/build/smartplus && \
    export SMART_INSTALL_DIR=/usr && \
    cd ${PREFIX}/source && \
    git clone ${SMART_URL} ${SMART_NAME} && \
    mkdir -p ${SMART_BUILD_DIR} && \
    cd ${SMART_BUILD_DIR} && \
    cmake \
        -D CMAKE_INSTALL_PREFIX=${SMART_INSTALL_DIR} \
        ${SMART_SOURCE_DIR} && \
    make && \
    make test && \
    make install && \
    cd ${PREFIX} && \
    rm -rf ${SMART_BUILD_DIR} && \
    rm -rf ${SMART_SOURCE_DIR}
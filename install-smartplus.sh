apt-get update && apt-get install -y build-essential gfortran wget git cmake libssl-dev libpng-dev libfreetype6-dev npm nodejs nodejs-legacy libblas-dev liblapack-dev libatlas-dev python2.7-dev python-pip libboost-all-dev && apt-get clean && rm -rf /var/lib/apt/lists/*

pip install --no-cache-dir --upgrade pip numpy matplotlib scipy pandas jupyter && npm install -g configurable-http-proxy && rm -rf ~/.npm

PREFIX=/scratch
mkdir -p ${PREFIX} && cd ${PREFIX} && mkdir archive && mkdir source && mkdir build && mkdir install

BOOST_DIR=/usr

# install Boost.Numpy
export BNUMPY_URL=https://github.com/ndarray/Boost.NumPy.git
export BNUMPY_NAME=Boost.Numpy
export BNUMPY_SOURCE_DIR=${PREFIX}/source/Boost.Numpy
export BNUMPY_BUILD_DIR=${PREFIX}/build/Boost.Numpy 
cd ${PREFIX}/source && git clone ${BNUMPY_URL} ${BNUMPY_NAME} && mkdir -p ${BNUMPY_BUILD_DIR} && cd ${BNUMPY_BUILD_DIR} && cmake -D Boost_NO_BOOST_CMAKE=ON -D BOOST_ROOT=${BOOST_DIR} -D CMAKE_INSTALL_PREFIX=${BOOST_DIR} -D LIB_SUFFIX=""  ${BNUMPY_SOURCE_DIR} 
make -j4 install && cd ${PREFIX} && rm -rf ${BNUMPY_SOURCE_DIR} && rm -rf ${BNUMPY_BUILD_DIR}

ARMA_DIR=/usr
#ENV LD_LIBRARY_PATH=${ARMA_DIR}/lib:${LD_LIBRARY_PATH}

# Install armadillo
export ARMA_URL=http://sourceforge.net/projects/arma/files/armadillo-7.400.1.tar.xz 
export ARMA_SHA1=c69bdd341dce731e8b7d35b9605e0e4b27e475e6 
export ARMA_ARCHIVE=${PREFIX}/archive/armadillo-7.400.1.tar.xz 
export ARMA_SOURCE_DIR=${PREFIX}/source/armadillo/7.400.1 
export ARMA_BUILD_DIR=${PREFIX}/build/armadillo/7.400.1 
export ARMA_INSTALL_DIR=/usr
wget --quiet ${ARMA_URL} --output-document=${ARMA_ARCHIVE} && echo "${ARMA_SHA1} ${ARMA_ARCHIVE}" | sha1sum -c && mkdir -p ${ARMA_SOURCE_DIR} && tar -xf ${ARMA_ARCHIVE} -C ${ARMA_SOURCE_DIR} --strip-components=1 && mkdir -p ${ARMA_BUILD_DIR} && cd ${ARMA_BUILD_DIR} 
cmake -D CMAKE_INSTALL_PREFIX=${ARMA_INSTALL_DIR} ${ARMA_SOURCE_DIR} && make -j4 install 
cd ${PREFIX} && rm -rf ${ARMA_ARCHIVE} && rm -rf ${PREFIX}/source/armadillo && rm -rf ${PREFIX}/build/armadillo

# install tini
export TINI_VERSION=0.10.0 && wget --quiet https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini && chmod +x tini && mv tini /usr/local/bin/

cd
mkdir smartplus-stack
cd smartplus-stack
# install smartplus
export SMART_URL=https://github.com/smartplus-team/smartplus.git 
export SMART_NAME=smartplus 
export SMART_SOURCE_DIR=${PREFIX}/source/smartplus 
export SMART_BUILD_DIR=${PREFIX}/build/smartplus 
export SMART_INSTALL_DIR=/usr 
cd ${PREFIX}/source && git clone ${SMART_URL} ${SMART_NAME} && mkdir -p ${SMART_BUILD_DIR} && cd ${SMART_BUILD_DIR} && cmake -D CMAKE_INSTALL_PREFIX=${SMART_INSTALL_DIR} ${SMART_SOURCE_DIR} && make && make test && make install && cd ${PREFIX} && rm -rf ${SMART_BUILD_DIR} && rm -rf ${SMART_SOURCE_DIR}

#!/bin/ksh
set -au

version='Release'
#version='Debug'
#version='Production'
#version=Bit
#version='RelWithDebInfo'

build_dir=${AQ_JEDI_BLD}
install_dir=${AQ_JEDI_DIR}
mkdir -p ${build_dir}
mkdir -p ${install_dir}
cd ${build_dir}
ecbuild --prefix=${build_dir} \
        --build=${version} \
        -DCMAKE_INSTALL_PREFIX=${install_dir} \
	-DENABLE_LORENZ95_MODEL=OFF \
	-DENABLE_OOPS_DOC="OFF" \
        ${AQ_JEDI_SRC}/bundle \
|| return $?
make update || return $?
make -j 24 || return $?
make install || return $?

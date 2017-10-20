#!/bin/sh
# Copyright 2016 Station X, Inc. or its affiliates. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

usage() { echo "Usage: $0 [-m] [-n <package_name>] [-d <destination_s3_bucket_folder>]" 1>&2; exit 1; }
MOSAIC=0
while getopts ":n:d:s" opt; do
  case $opt in
    m)
      echo "will do a mosaic build" >&2
      MOSAIC=1
      ;;
    n)
      echo "package name: $OPTARG" >&2
      PACKAGE_NAME=$OPTARG
      ;;
    d)
      echo "destination s3 bucket/folder: $OPTARG" >&2
      S3_PATH=$OPTARG
      # remove trailing slash
      S3_PATH=${S3_PATH%/}
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${PACKAGE_NAME}" ]; then
    usage
fi

if [ -z "${S3_PATH}" ]; then
    usage
fi

# keep track of current directory since it has handler and we will need to package it
INITIAL_DIR=`pwd`

# to install extra packages, we need to enable epel repository
sudo yum-config-manager --enable epel

# First, make sure everything is up-to-date:
sudo yum -y update
sudo yum -y upgrade

# install everything
# readline is needed for rpy2, and fortran is needed for R
sudo yum install -y python27-devel python27-pip gcc gcc-c++ readline-devel libgfortran.x86_64 R.x86_64
# epel packages:
sudo yum install -y gdal-devel.x86_64 proj-devel.x86_64

# build survival R function if requested
if [ $MOSAIC == 1 ]; then
    # make sure R can find local directory
    echo "export R_LIBS=~/HOME/lambda" > .Renviron
    cd /tmp
    sudo R CMD BATCH ~/install_R_packages.R
fi

# setup virtualenv and install rpy2. Need version <2.9 for python2.7 support
virtualenv ~/env && source ~/env/bin/activate
pip install "rpy2<2.9"

# create a directory called lambda for our package
mkdir $HOME/lambda && cd $HOME/lambda
# copy R 
cp -rL /usr/lib64/R/* $HOME/lambda/
# Use ldd on R executable to find shared libraries, and copy all of the ones that were not already on the box
cp /usr/lib64/R/lib/libR.so $HOME/lambda/lib/
cp /usr/lib64/libgomp.so.1 $HOME/lambda/lib/  
#missing in this location  cp /usr/lib64/libblas.so.3 $HOME/lambda/lib/
cp /usr/lib64/libgfortran.so.3 $HOME/lambda/lib/
cp /usr/lib64/libquadmath.so.0 $HOME/lambda/lib/

# we also need to grab this one (as we learned from trial and error)
#cp /usr/lib64/liblapack.so.3 $HOME/lambda/lib/
cp /usr/lib64/R/modules/lapack.so $HOME/lambda/lib/

# newly added
cp /usr/lib64/R/lib/libRblas.so $HOME/lambda/lib
cp /usr/lib64/R/lib/libRlapack.so $HOME/lambda/lib
cp /usr/lib64/R/lib/libRrefblas.so $HOME/lambda/lib
cp /usr/lib64/libtre.so.5 $HOME/lambda/lib
cp /usr/lib64/liblzma.so.5 $HOME/lambda/lib
cp /usr/lib64/libicuuc.so.50 $HOME/lambda/lib
cp /usr/lib64/libicui18n.so.50 $HOME/lambda/lib
cp /usr/lib64/libicudata.so.50 $HOME/lambda/lib
cp /usr/lib64/libstdc++.so.6 $HOME/lambda/lib
#cp /lib64/libpthread.so.0 $HOME/lambda/lib
#cp /lib64/libc.so.6 $HOME/lambda/lib
#cp /lib64/libm.so.6 $HOME/lambda/lib
#cp /lib64/libreadline.so.6 $HOME/lambda/lib
#cp /lib64/libpcre.so.0 $HOME/lambda/lib
#cp /lib64/libbz2.so.1 $HOME/lambda/lib
#cp /lib64/libz.so.1 $HOME/lambda/lib
#cp /lib64/librt.so.1 $HOME/lambda/lib
#cp /lib64/libdl.so.2 $HOME/lambda/lib
#cp /lib64/libgcc_s.so.1 $HOME/lambda/lib
#cp /lib64/libtinfo.so.5 $HOME/lambda/lib

# gdal,proj.4 etc
cp /usr/lib64/libgdal.so.1 $HOME/lambda/lib/
cp /usr/lib64/libproj.so.0.6.6 $HOME/lambda/lib/
#./usr/lib64/libproj.so.0.6.6
# /usr/lib64/libgdal.so.1

# copy R executable to root of package
cp $HOME/lambda/bin/exec/R $HOME/lambda/

#Add the libraries from the activated Python virtual environment
cp -r $VIRTUAL_ENV/lib64/python2.7/site-packages/* $HOME/lambda
# we could copy all of $VIRTUAL_ENV/lib/python2.7/site-packages/, but let's grab the esseentials only
cp -r $VIRTUAL_ENV/lib/python2.7/site-packages/singledispatch* $HOME/lambda
# added by adriaan because above missing:
cp -r $VIRTUAL_ENV/lib/python2.7/dist-packages/singledispatch* $HOME/lambda

cd $HOME/lambda
echo "zipping $HOME/${PACKAGE_NAME}"
zip -r9 $HOME/${PACKAGE_NAME} *
echo "changing to directory $INITIAL_DIR"
cd $INITIAL_DIR
echo "adding handler.py to zip"
zip -r9 $HOME/${PACKAGE_NAME} handler.py

# copy to S3
aws s3 cp $HOME/${PACKAGE_NAME} ${S3_PATH}/${PACKAGE_NAME}

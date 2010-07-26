#!/bin/bash

mkdir -p $ARCH
mkdir -p include


mkdir tmp
cd tmp
tar -xzf ../archives/libsvm-2.85.tar.gz 
cd libsvm-2.85/
make clean
make
cp svm.o ../../${ARCH}/
cp svm.h ../../include/
cd ../..
rm -rf tmp

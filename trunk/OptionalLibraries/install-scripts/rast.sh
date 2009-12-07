#!/bin/bash

mkdir -p $ARCH
mkdir -p include

mkdir tmp
cd tmp
tar -xzf ../archives/rast-feature.tgz 
cd rast-feature/
make
cp rast.o ../../$ARCH
cp h/*.h ../../include/
cd ../..
rm -rf tmp

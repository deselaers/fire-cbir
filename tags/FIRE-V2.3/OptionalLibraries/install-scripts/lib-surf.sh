#!/bin/bash

mkdir -p $ARCH
mkdir -p include

mkdir tmp
cd tmp
tar -xzf ../archives/SURF-V1.0.9.tar.gz 
cd SURF-V1.0.9/
#cp libSurf.a libSurf.so ../../$ARCH
cp *.h ../../include/
cd ../..
rm -rf tmp

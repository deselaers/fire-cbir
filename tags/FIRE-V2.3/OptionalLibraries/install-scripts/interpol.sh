#!/bin/bash

cat <<EOF
This script tries to put all files required for the interpolation function in place

Download interpol.tar from    http://bigwww.epfl.ch/thevenaz/interpolation/

Please put the archive to $PWD/archives/interpol.tar and press enter to continue
EOF
read

mkdir -p $ARCH
mkdir -p include
mkdir tmp
cd tmp
tar -xf ../archives/interpol.tar 
cd interpol/
g++ -c coeff.c 
g++ -c interpol.c 
cp interpol.o coeff.o ../../$ARCH/
cp coeff.h  interpol.h  ../../include/
cd ../../
rm -rf tmp/

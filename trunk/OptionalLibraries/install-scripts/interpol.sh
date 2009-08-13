#!/bin/bash
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

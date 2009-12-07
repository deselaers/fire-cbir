#!/bin/bash
mkdir -p $ARCH
mkdir -p include


cat <<EOF
This script tries to put all files required for KD trees to the right location.

Please put the archive to $PWD/archives/libknn.tgz and press enter to continue
EOF

read

OLDP=$PWD
mkdir -p tmp include $ARCH
cd tmp
tar -xzf ../archives/libknn.tgz 
cd src/
cd libknn/
make clean
make linux
cp include/knn_api.h ../../../include/
cp lib/libknn.a ../../../$ARCH/
cd $OLDP
rm -rf tmp

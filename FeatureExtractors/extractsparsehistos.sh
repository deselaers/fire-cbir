#!/bin/bash

export QSUBOPT="-m a -m e -m b -cwd -j y -w e -S /bin/bash -v LD_LIBRARY_PATH=/home/thomasd/fire-lib -l h_vmem=1G -l h_rt=24:00:00 "

FIREPATH=$HOME/work/src/fire/bin/$ARCH
exsh=${FIREPATH}/extractsparsepatchhistogram

DB=$1

TRAIN=train
TEST=test

OPTS="-sc 0.25 --gray -st 4 -pos -r 6"
SUFF="s0.25-gray-r6"

#OPTS="-sc 0.25 --gray -st 4 -pos -r 4"
#SUFF="s0.25-gray-r4"


trainjid=1
trainjid=`qsubmit -n exsh-${DB}-${TRAIN}-${SUFF} $HOME/work/src/fire/bin/$ARCH/extractsparsepatchhistogram $OPTS -mo 3 -spca pcas/sparse-medaat-${SUFF}-mo3-mp3:5:10:15-ws7.pca.gz -mp 3 5 10 15 -ws 7 -su gray-${SUFF}-mo3-mp3:5:10:15-ws7.medaat.sparsehisto.gz --filelist lists/process/${DB}-${TRAIN}  | grep Job: | cut -f 2 -d ' '`
testjid=`qsubmit -w $trainjid -n exsh-${DB}-${TEST}-${SUFF} $HOME/work/src/fire/bin/$ARCH/extractsparsepatchhistogram $OPTS -lpca pcas/sparse-medaat-${SUFF}-mo3-mp3:5:10:15-ws7.pca.gz -mp 3 5 10 15 -ws 7 -su gray-${SUFF}-mo3-mp3:5:10:15-ws7.medaat.sparsehisto.gz --filelist lists/process/${DB}-${TEST} | grep Job: | cut -f 2 -d ' '`


echo "$trainjid $testjid"
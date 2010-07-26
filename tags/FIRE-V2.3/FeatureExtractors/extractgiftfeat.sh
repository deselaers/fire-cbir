#!/bin/bash

GIFTEX=~gass/gift-0.1.15b/bin/gift-extract-features
UNPACKFTS=$HOME/work/src/fire/tools/unpackfts.pl
FTS2SH=$HOME/work/src/fire/tools/ftsdescr2sh.py

SUFFIX=gift.sparsehistogram


function giftex() {
 fn=`mktemp`.ppm
 convert -geometry 256x256! $1 ${fn}
 ${GIFTEX} $fn
 rm ${fn} ${fn//.ppm}
 perl ${UNPACKFTS} ${fn//.ppm}.fts > $1.giftfeat
 ${FTS2SH} $1.giftfeat 
}


if [[ $# -eq 0 || "$1" = "-h" ]] ; then
  cat <<EOF
USAGE: `basename $0` [options]  (--filelist <filelist>|| --images image1 [image2 [image3 [image4 .. ] ] ] )
  
EOF
exit
fi

i=1


cat <<EOF
 fts0: local color feature
 fts1: global color feature
 fts2: local texture feature
 fts3: global texture feature
EOF

while [ $i -le $# ] ; do
  case $1 in
    --filelist)
        FILELIST=$2
	echo filelist= $FILELIST
        shift 2
        ;;
    --images)
        FILELIST=/tmp/tmpfilelist
        shift 1
        rm -f  $FILELIST
        while [ $i -le $# ] ; do
          eval echo \$$i >> $FILELIST
          i=$((i+1))
        done
        ;;
    *)
        echo "unknown option: $1"
        exit
  esac
done

echo opt= $OPT suffix= $SUFFIX filelist=$FILELIST
cat $FILELIST | (while read filename ; do echo $filename; giftex $filename ; done)

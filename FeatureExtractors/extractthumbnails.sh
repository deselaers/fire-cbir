#!/bin/bash


if [[ $# -eq 0 || "$1" = "-h" ]] ; then
  cat <<EOF
USAGE: `basename $0` [options]  (--filelist <filelist>|| --images image1 [image2 [image3 [image4 .. ] ] ] )

OPTIONS: 
  --32x32
  --Xx32
  --suffix <suffix>
  
EOF
exit
fi

i=1

while [ $i -le $# ] ; do
  case $1 in
    --suffix)
       SUFFIX=$2
       shift 2
       ;;
    --32x32)
        OPT="${OPT} -resize 32x32!"
        SUFFIX=${SUFFIX:-tn.32x32.png}
        shift
        ;;
    --Xx32)
        OPT="${OPT} -resize x32"
        SUFFIX=${SUFFIX:-tn.Xx32.png}
        shift
        ;;
    --colors)
        OPT="${OPT} -colors $2"
        shift 2
        ;;
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
cat $FILELIST | (while read filename ; do echo $filename; convert $OPT $filename $filename.$SUFFIX ; done)

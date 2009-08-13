#!/bin/bash

#set -v

### ----------------------------------------------------------------------
# SETTINGS
#
export LD_LIBRARY_PATH=$HOME/work/mpeg7/LIB/lib/
XMMAIN=$HOME/work/mpeg7/newsrc/XMMain_linux6.1.exe


if [[ $# -lt 2 || "$1" = "-h" ]] ; then
  cat <<EOF
USAGE: `basename $0` [options] <featuretype> <filelist>
  <featuretype> is one of 
      scd: scalable color descriptor
      csd: color structure descriptor
      dcd: dominant color descriptor
      cld: color layout descriptor
      tbd: texture browsing descriptor
      htd: homogeneous texture descriptor
      leh: local edge histogram
      3ss: 3d shape spectrum
      con: contrast based shape descriptor
      rsd: region shape descriptor
  <filelist> is a filelist containing one image basename per line
      the basename of an image is the string describing it in the "file"  line
      of the filelist used in FIRE later
 
  Available options:
      -spath <path>  to specify the path to be prepended to the basename
                     to find the image from the current working directory.
                     this parameter is not necessary when you run this script from
                     the place where the basenames in the filelist are correct to find
                     the images
      -dpath <path>  to specify where the data is stored. Usually the current working 
                     directory is used for this.
      -fname <filename> to specify the prefix of the filenames which are generated
EOF
exit 0
fi

spath=`pwd`
dpath=`pwd`
fname=`pwd`

while [ "${1:0:1}" = "-" ] ; do
    case $1 in
        "-spath")
                spath=$2
                shift 2
                ;;
        "-dpath")
                dpath=$2
                shift 2
                ;;
        "-fname")
                fname=$2
                shift 2
                ;;
        *)
            echo "Option $i unknown"
            shift 1
            ;;
    esac
done

featuretype=$1
filelist=$2

# no filename was given: take default
if [ -n $fname ] ; then
  fname=$featuretype
fi

# make filenames for the files to be generated
exlist=$filelist.ex
lstfile=$fname.lst
parfile=$fname.par
bitstream=$fname.mp7
# check whether the files are located in a specified directory
if [ -n $dpath ] ; then
    mkdir -p $dpath
    lstfile=$dpath/$fname.lst
    exfile=$dpath/`basename $filelist.ex` 
    parfile=$dpath/$fname.par
    bitstream=$dpath/$fname.mp7
fi

# create the appropriate list for feature extraction
if [ -n $spath ] ; then
    cat $filelist | sed 's:^:'$spath'/:' > $exfile
else
    cp $filelist $exfile
fi


#cp the listfile to the right location
cp $filelist $lstfile

type=""
case $featuretype in 
    scd) # scalable color descriptor
    type=ScalableColor
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
NumberOfBitplanesDiscarded 3
NumberOfCoefficients 64
CodingMode 0
EOF
    ;;
    csd) #  color structure
    type=ColorStructure
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
ColorQuantSize 256
CodingMode 0
EOF
    ;;
    dcd) # dominant color
    type=DominantColor
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
ColorSpacePresent			0
ColorQuantizationPresent	0
VariancePresent				0
SpatialCoherency			0
EOF
        ;;
    cld) # color layout
    type=ColorLayout
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
EOF
        ;;
    tbd) # Texture Browsing Descriptor ------------------------------
    type=TextureBrowsing
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
/* layer 0: basic layer (3 components): 1: full layer (5 components) */
layer 0
EOF
        ;;
    htd) # Homogeneous Texture Descriptor ------------------------------
    type=HomogeneousTexture
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 
/* layer: (0:base-layer 32-components  1:full-layer 62components) */
layer 1
/* default matching option for Client  n: normal matching* r:rotation invariant matching s:scale invariant matching rs:rotation and scale invariant matchin */
option  n
EOF
        ;;
    leh) # Local Edge Histogram ------------------------------
    type=EdgeHistogram
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
EOF
        ;;
    3ss) # 3d shape spectrum ------------------------------
    type=3DShapeSpectrum
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
NoOfBins 100
NoOfBits 12
Metric 0
EOF
        ;;
    con) # contour shape descriptor ------------------------------
    type=ContourShape
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0

EOF
        ;;
    rsd) # region shape descriptor ------------------------------
    type=RegionShape
    cat > $parfile <<EOF
Application	${type}Client
ListFile $lstfile
Bitstream $bitstream
NoOfMatches `wc -l $lstfile`
CodingMode 0
EOF
        ;;
    *) # unknown feature type
        echo "Unknown feature type: $featuretype"
        ;;
esac

# parameter files generated

# starting extraction
echo $XMMAIN -p$parfile -a${type}Server -b$bitstream -l$exfile

$XMMAIN -p$parfile -a${type}Server -b$bitstream -l$exfile

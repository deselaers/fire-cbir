#!/bin/bash

if [ $# -ne 2 ] ; then
  echo "USAGE: install-fire.sh <destinationpath> <imagepath>"
  echo "  where destinationpath is the path where FIRE will be installed" 
  echo "  and in imagepath some images (.jpg) can be found."
  echo "IMPORTANT: both paths should be absolute paths, i.e. starting with /"
  exit
fi 

DSTPATH=$1
IMGPATH=$2

# get architecture
if ARCH=`uname -p ` ; then echo "ARCH set" ; else ARCH=`uname -m` ; fi
export ARCH
# -------------------- settings --------------------
FIREFILE=firev2.1.0.tgz
APACHEFILE=httpd-2.0.55.tar.gz
APACHEURL=http://apache.mirrorplus.org/httpd/$APACHEFILE
INTERPOLURL=http://bigwww.epfl.ch/thevenaz/interpolation/interpol.tar

# -------------------- DESTINATION --------------------
mkdir $DSTPATH

cd $DSTPATH
wget 'http://www-i6.informatik.rwth-aachen.de/~deselaers/files/firev2.1.0.tgz'
tar -xvzf $FIREFILE

wget $APACHEURL
tar -xvzf $APACHEFILE

# -------------------- BUILD APACHE --------------------

cd httpd-2.0.54
./configure --prefix=$DSTPATH/apache
make 
make install

cd $DSTPATH/apache/conf

# configure apache to listen on port 8123
cp httpd.conf httpd.conf.bak
cat httpd.conf.bak | sed 's/Listen 80/Listen 8123/' > httpd.conf

cd $DSTPATH


# -------------------- BUILD FIRE --------------------
cd fire
make directories
cd lib/$ARCH

# ---- get interpolation library and compile ----
wget $INTERPOLURL
tar -xf interpol.tar
cd interpol
g++ -c coeff.c
g++ -c interpol.c


cd $DSTPATH/fire
make all

# ---- copy cgi into apache installation ----
cp $DSTPATH/fire/cgi/fire.py $DSTPATH/apache/cgi-bin/fire.cgi
cp $DSTPATH/fire/cgi/img.py $DSTPATH/apache/cgi-bin/img.py
cp $DSTPATH/fire/python/firesocket.py $DSTPATH/apache/cgi-bin/firesocket.py
cp $DSTPATH/fire/cgi/fire-td.ico $DSTPATH/apache/htdocs
cat $DSTPATH/fire/cgi/config.py | sed "s:DESTPATH:$DSTPATH:g" > $DSTPATH/apache/cgi-bin/config.py
mkdir $DSTPATH/apache/htdocs/images
cp $DSTPATH/fire/cgi/positive.png  $DSTPATH/fire/cgi/negativ.png $DSTPATH/fire/cgi/neutral.png $DSTPATH/apache/htdocs/images


# -------------------- EXTRACT FEATURES FROM IMAGES --------------------
cd $IMGPATH
find . | grep -i jpg > purefiles
$DSTPATH/fire/bin/$ARCH/extractcolorhistogram --color --suffix color.histo.gz --filelist purefiles
cat > $DSTPATH/filelist <<EOF
FIRE_filelist
suffix color.histo.gz
classes no
path $IMGPATH
EOF
cat $IMGPATH/purefiles | sed 's/^/file /' >> $DSTPATH/filelist

# -------------------- START APACHE --------------------

$DSTPATH/apache/bin/httpd

# -------------------- START FIRE --------------------
echo "Now you can access your images via http://localhost:8123/cgi-bin/fire.cgi"

$DSTPATH/fire/bin/$ARCH/fire -f $DSTPATH/filelist -D -r 10



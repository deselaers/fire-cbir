#!/usr/bin/python

import sys
import cgi
import re
import Image

sys.stderr=sys.stdout
try:
    form=cgi.FieldStorage()
    imagefilename=form["image"].value
    outWidth=-1
    outHeight=-1
    maxDim=-1

    if(form.has_key("max")):
        maxDim=int(form["max"].value)
    if(form.has_key("width")):
        outWidth=int(form["width"].value)
    if(form.has_key("height")):
        outHeight=int(form["height"].value)
    inImg=Image.open(imagefilename)
    inWidth=inImg.size[0]
    inHeight=inImg.size[1]

    if maxDim!=-1:
        if inWidth>inHeight:
            outWidth=maxDim
        else:
            outHeight=maxDim
    
    if (outWidth==-1 and outHeight==-1):
        outImg=inImg
    elif (outWidth==-1):
        scaleFactor=float(outHeight)/float(inHeight)
        outWidth=int(scaleFactor*inWidth)
        outImg=inImg.resize((outWidth,outHeight),Image.BILINEAR)
    elif (outHeight==-1):
        scaleFactor=float(outWidth)/float(inWidth)
        outHeight=int(scaleFactor*inHeight)
        outImg=inImg.resize((outWidth,outHeight),Image.BILINEAR)
    else:
        outImg=inImg.resize((outWidth,outHeight),Image.BILINEAR)
    contenttype="image/jpg"
    print "Content-Type: "+contenttype+"\n"
    outImg=outImg.convert("RGB")
    outImg.save(sys.stdout,"jpeg")
except:
    contenttype="text/html"
    print "Content-Type: "+contenttype+"\n"
    print "Access not permitted"


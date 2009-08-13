#!/usr/bin/env python

import sys,gzip

class Feature:
    def __init__(self):
        self.posx=0
        self.posy=0
        self.scale=0
        self.vec=[]

class LocalFeatures:
    def __init__(self):
        self.winsize=0
        self.dim=0
        self.subsampling=0
        self.padding=0
        self.numberOfFeatures=0
        self.varthreshold=0
        self.zsize=0
        self.filename=""
        self.imagesize=(0,0)
        self.features=[]
        self.selffilename=""

    def load(self,filename):
        f=gzip.GzipFile(filename)
        self.selffilename=filename

        for l in f:
            toks=l.split()
            if toks[0]=="FIRE_localfeatures":
                fine=True
            elif toks[0]=="winsize":
                self.winsize=int(toks[1])
            elif toks[0]=="dim":
                self.dim=int(toks[1])
            elif toks[0]=="subsampling":
                self.subsampling=int(toks[1])
            elif toks[0]=="padding":
                self.padding=int(toks[1])
            elif toks[0]=="numberOfFeatures":
                self.numberOfFeatures=int(toks[1])
            elif toks[0]=="varthreshold":
                self.varthreshold=float(toks[1])
            elif toks[0]=="zsize":
                self.zsize=int(toks[1])
            elif toks[0]=="filename":
                self.filename=toks[1]
            elif toks[0]=="imagesize":
                self.imagesize=(int(toks[1]),int(toks[2]))
            elif toks[0]=="features":
                if self.numberOfFeatures!=int(toks[1]):
                    print "Weird... inconsistent number of Features and features expected in the file"
            elif toks[0]=="feature":
                f=Feature()
                toks.pop(0)
                f.posx=float(toks.pop(0))
                f.posy=float(toks.pop(0))
                f.scale=float(toks.pop(0))
                f.vec=map(lambda x: float(x), toks)
                self.features+=[f]
                
            else:
                print "'%s' is an unknown keyword... ignoring and continuiing nonetheless."%(toks[0])

    




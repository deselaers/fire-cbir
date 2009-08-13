#!/usr/bin/env python
import gzip

# class for FIRE histogramfeatures
class SparseHistogramFeature:
    def __init__(self):
        self.bins={}
        self.data={}
        self.steps=[]
        self.stepSize=[]
        self.min=[]
        self.max=[]
        self.dim=0
        self.numberofBins=0
        self.counter=0
        self.filename=""

    def save(self,filename):
        f=gzip.GzipFile(filename,"w")
        print >>f,"FIRE_sparse_histogram"
        print >>f,"# sparse histogram file for FireV2"
        print >>f,"dim",self.dim
        print >>f,"counter",self.counter
        print >>f,"steps",reduce(lambda x,y: str(x)+" "+str(y),self.steps)
        print >>f,"min",reduce(lambda x,y: str(x)+" "+str(y),self.min)
        print >>f,"max",reduce(lambda x,y: str(x)+" "+str(y),self.max)
        print >>f,"bins",len(self.bins)
        for i in self.bins:
            print >>f,"data",i,self.bins[i]
        f.close()


    def load(self,filename):
        f=gzip.GzipFile(filename)
        self.filename=filename

        for line in f:
            toks=line.split()

            if toks[0]=="FIRE_sparse_histogram":
                fine=True
            elif line[0]=="#":
                fine=True                
            elif toks[0]=="dim":
                self.dim=int(toks[1])
            elif toks[0]=="counter":
                self.counter=int(toks[1])
            elif toks[0]=="steps":
                self.steps=map(lambda x: int(x), toks[1:])
            elif toks[0]=="min":
                self.min=map(lambda x: float(x),toks[1:])
            elif toks[0]=="max":
                self.max=map(lambda x: float(x),toks[1:])
            elif toks[0]=="bins":
                self.numberofBins=int(toks[1])
            elif toks[0]=="data":
                self.bins[toks[1]]=int(toks[2])
                self.data[toks[1]]=float(toks[2])/float(self.counter)
            else:
                print "Unparsed line: '"+line+"'."
                
            
                           

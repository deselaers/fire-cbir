#!/usr/bin/env python
import gzip

# class for FIRE histogramfeatures
class HistogramFeature:
    def __init__(self):
        self.steps=[]
        self.bins=[]
        self.min=[]
        self.max=[]
        self.stepsize=[]
        self.dim=0
        self.counter=0
        self.filename=""
        self.data=[]

    def save(self,filename):
        f=gzip.GzipFile(filename,"w")
        print >>f,"FIRE_histogram"
        print >>f,"# created with histogramfeature.py"
        print >>f,"dim",self.dim
        print >>f,"counter",self.counter
        print >>f,"steps"," ".join(map(lambda x: str(x),self.steps))
        print >>f,"min"," ".join(map(lambda x: str(x), self.min))
        print >>f,"max"," ".join(map(lambda x: str(x), self.max))
        print >>f,"data"," ".join(map(lambda x: str(x), self.data))
        f.close()

    def load(self,filename):
        f=gzip.GzipFile(filename)
        self.filename=filename

        for line in f:
            toks=line.split()

            if toks[0]=="FIRE_histogram":
                fine=True
            elif toks[0]=="#":
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
            elif toks[0]=="data":
                self.bins=map(lambda x: int(x),toks[1:])
                try:
                    self.data=map(lambda x: float(x)/float(self.counter), self.bins)
                except:
                    print "Error reading file:", filename,"probably division by zero, continuing anyway."

                   
            else:
                print "Unparsed line: '"+line+"'."
                
            
                           

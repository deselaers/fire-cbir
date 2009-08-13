#!/usr/bin/python

import sys

#this script merges fire distance files.
#
# as command line parameters, arbitrary many distance files are given and 
# the merged distance file is printed.

infiles=[]
line=[]
for i in sys.argv[1:]:
    infiles.append(open(i,"r"))
    line.append("")

line[0]=infiles[0].readline()
count=0
while(line[0]):
    for i in range(1,len(infiles)):
        line[i]=infiles[i].readline()

    if line[0].startswith("nofdistances"):

        nDist=0
        nImg=0
        
        for l in line:
           
            toks=l.split()
            nDist+=int(toks[2])
            nImg=int(toks[1])
        print "nofdistances",nImg,nDist

    elif not line[0].startswith("#"):
       print count,
       for l in line:
           toks=l.split()
           for i in range(1,len(toks)):
               print float(toks[i]),
       print
       count+=1
    line[0]=infiles[0].readline()

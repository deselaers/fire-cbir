#!/usr/bin/python

import sys,re

import filelist    

filelist=filelist.FileList()
filelist.load(sys.argv[1])

for i in filelist:
    cls=i[1]
    fname=i[0]
    print fname,
    for j in filelist:
        if j[1]==cls and j[0]!=fname:
            print j[0],
    print


#print filelist

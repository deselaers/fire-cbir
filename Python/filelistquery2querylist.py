#!/usr/bin/python

import sys,re,filelist

if len(sys.argv)<3:
    print """USAGE filelistquery2querylist.py <querylist> <dblist>"""
    sys.exit(5)

qlist=filelist.FileList()
qlist.load(sys.argv[1])

dblist=filelist.FileList()
dblist.load(sys.argv[2])

for i in qlist:
    cls=i[1]
    fname=i[0]
    print fname,
    for j in dblist:
        if j[1]==cls:
            print j[0],
    print

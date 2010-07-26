#!/usr/bin/python

import sys,re
import filelist,porterstemmer
querylist=filelist.FileList()
filelist=filelist.FileList()

filelist.load(sys.argv[1])
querylist.load(sys.argv[2])

stemmer=porterstemmer.PorterStemmer()

if not filelist.descriptions:
    print "no descriptions, this program is not appropriate"
    sys.exit(10)

for i in querylist:
    cls=i[1]
    desc=[]
    #print "Before stemming: ",i[2]
    for w in i[2]:
        w=stemmer.stem(w,0,len(w)-1)
        desc+=[w]
    #print "After stemming:",desc
    rels={}
    for j in filelist:
        desc2=[]
        for w in j[2]:
            w=stemmer.stem(w,0,len(w)-1)
            if w in desc:
                rels[j[0]]=1

    print i[0],
    for f in rels:
        print f,
    print
            
#print filelist

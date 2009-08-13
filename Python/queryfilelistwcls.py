#!/usr/bin/python

# easy script for automatic batches with fire
# takes database files for queries
# calculates error rate only

import sys, socket, re, os, string, time, traceback,numarray
import firesocket,filelist
# forget own name in commandline options
sys.argv.pop(0)

server="localhost"
port=12960
database=""
querybase=""
l1o=False
quit=False

while len(sys.argv)>0:
    argv=sys.argv.pop(0)
    if argv=="-s":
        server=sys.argv.pop(0)
    elif argv=="-p":
        port=int(sys.argv.pop(0))
    elif argv=="-f":
        database=sys.argv.pop(0)
    elif argv=="-q":
        querybase=sys.argv.pop(0)
    elif argv=="-l1o":
        l1o=True
    elif argv=="-x":
        quit=True
    else:
        print "Unknown option:",argv
        print """
USAGE: queryfilelistwcls.py <options>
      -h  show this help
      -s  server (default: localhost)
      -p  port (default: 12960)
      -f  database 
      -q  database for query
      -l1o make leaving one out with database
      -x  exit after having finished
"""
        sys.exit(10)

print "SETTINGS: server=",server,"port=",port,"database=",database,"querybase=",querybase,"l1o=",l1o

if database=="" or querybase=="" and not l1o:
    print """
USAGE: queryfilelistwcls.py <options>
      -h  show this help
      -s  server (default: localhost)
      -p  port (default: 12960)
      -f  database 
      -q  database for query
      -l1o make leaving one out with database
      -x  exit after having finished
"""
    sys.exit(10)
    

print "SETTINGS: server=",server,"port=",port,"database=",database,"querybase=",querybase,"l1o=",l1o
sys.stdout.flush()
f=filelist.FileList()
f.load(database)

q=filelist.FileList()
if not l1o:
    q.load(querybase)
else:
    q=f

if not f.classes:
    print "Need classes in database file"
    sys.exit(10)

if not l1o and not q.classes:
    print "Need classes in querybase file"
    sys.exit(10)
    

s=firesocket.FIRESocket()
s.connect(server,port)
try:
    s.sendcmd("info")
    res=s.getline()
    status=re.split(" ",res)
    keyword=status.pop(0)
    dbSize=int(status.pop(0))
    if dbSize!=len(f.files):
        print "database in retriever and database loaded are of different size:", dbSize,"!=",len(f.files)
        s.sendcmd("bye")
        time.sleep(1)
        sys.exit(10)

    print "database has",dbSize,"images."
    if not l1o:
        s.sendcmd("setresults 1"); res=s.getline()
    else:
        s.sendcmd("setresults 2"); res=s.getline()
        
    classified=0; errors=0; correct=0


    no_classes=0
    for qf in q.files:
        if qf[1]>no_classes:
            no_classes=qf[1]
    for df in f.files:
        if df[1]>no_classes:
            no_classes=df[1]
    no_classes+=1

    print no_classes,"classes."
    confusionmatrix=numarray.zeros((no_classes,no_classes),numarray.Int)
    counter=0
    for qf in q.files:
        counter+=1
        qname=qf[0]
        qcls=qf[1]
        s.sendcmd("retrieve "+qname)
        res=s.getline()
        res=re.sub('[ ]*$','',res) # and parse the results
        results=re.split(" ",res)
        if not l1o:
            if len(results)!=2:
                print "expected 2 tokens, got",len(results)
                sys.exit(10)
        else:
            if len(results)!=4:
                print "expected 4 tokens, got",len(results)
                sys.exit(10)
            results.pop(0)
            results.pop(0)
                
        returnedimage=results.pop(0)
        score=float(results.pop(0))

        for df in f.files:
            if not l1o:
                if df[0]==returnedimage:
                    dcls=df[1]
                    break
            else:
                if df[0]==returnedimage and df[0]!=qname:
                    dcls=df[1]
                    break
                    
        classified+=1
        if dcls==qcls:
            correct+=1
        else:
            errors+=1

        confusionmatrix[dcls][qcls]+=1
        print "Query "+str(counter)+"/"+str(len(q.files))+":",qname,"("+str(qcls)+") NN:",returnedimage,"("+str(dcls)+")","ER:",float(errors)/float(classified)*100
        sys.stdout.flush()



    print "RESULT: ER:",float(errors)/float(classified)*100
    print confusionmatrix
    time.sleep(1)
    if quit:
        s.sendcmd("quit")
    else:
        s.sendcmd("bye")
    time.sleep(1) #don't kill the server by exiting tooo fast 

except KeyboardInterrupt, e:
    s.sendcmd("bye")
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast 

except Exception, e:
    if quit:
        s.sendcmd("quit")
    else:
        s.sendcmd("bye")
    print e
    time.sleep(1)

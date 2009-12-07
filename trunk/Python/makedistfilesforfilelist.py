#!/usr/bin/python

# easy script for automatic batches with fire
# takes database files for queries
# calculates error rate only

import sys, socket, re, os, string, time, traceback
#sys.path.append("/u/deselaers/work/src/fire/python")
import firesocket,filelist
# forget own name in commandline options
sys.argv.pop(0)

server="localhost"
port=12960
database=""
querybase=""
l1o=False
quit=False

start=-1
stop=-1

while len(sys.argv)>0:
    argv=sys.argv.pop(0)
    if argv=="-s":
        server=sys.argv.pop(0)
    elif argv=="-p":
        port=int(sys.argv.pop(0))
    elif argv=="-q":
        querybase=sys.argv.pop(0)
    elif argv=="-x":
        quit=True
    elif argv=="-S":
        suffix=sys.argv.pop(0)
    elif argv=="-start":
        start=int(sys.argv.pop(0))
    elif argv=="-stop":
        stop=int(sys.argv.pop(0))
    else:
        print "Unknown option:",argv
        print """
USAGE: makedistfilesforfilelist.py <options>
      -h  show this help
      -s  server (default: localhost)
      -p  port (default: 12960)
      -q  database for query
      -S  suffix
      -start <n>, -stop <n> to process a subset of the querylist, -1 for full database
      -x  exit after having finished
"""
        sys.exit(10)

print "SETTINGS: server=",server,"port=",port,"querybase=",querybase

if querybase=="":
    print """
USAGE: makedistfilesforfilelist.py <options>
      -h  show this help
      -s  server (default: localhost)
      -p  port (default: 12960)
      -q  database for query
      -S  suffix
      -start <n>, -stop <n> to process a subset of the querylist
      -x  exit after having finished
"""
    sys.exit(10)
    

sys.stdout.flush()

q=filelist.FileList()
q.load(querybase)

s=firesocket.FIRESocket()
s.connect(server,port)

if start==-1:
    start=0
if  stop==-1:
    stop=len(q.files)
    
try:
    for qf in q.files[start:stop]:
        cmd="savedistances "+qf[0]+" "+qf[0]+"."+suffix
        s.sendcmd(cmd)
        res=s.getline()
        print res

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

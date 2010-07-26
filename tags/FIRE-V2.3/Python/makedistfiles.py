#!/usr/bin/python
import sys, socket, re, os, string, time
sys.path.append("/u/deselaers/work/src/fire/python")
import firesocket,filelist

#-----------------------------------------------------------------
s=firesocket.FIRESocket()

try:
    if len(sys.argv) < 5:
        print """USAGE:
        makedistfiles.py <querylist> <server> <port> <suffix> [-q]
        """
    else:
        ql=open(sys.argv[1],"r")
        host=sys.argv[2]
        port=int(sys.argv[3])
        suffix=sys.argv[4]
        try:
            s.connect(host,port)
        except:
            print "Connecting to server "+host+" on port "+str(port)+" failed."

        lines=ql.readlines()
        lines=map(lambda line:
                  re.sub("\n","",line), lines)
        ql.close()

        i=0
        for line in lines:
            print i,"/",len(lines)
            i+=1
            tokens=re.split(" ",line)
            file=tokens[0]
            cmd="savedistances "+file+" "+file+"."+suffix
            print "SENDING: ",cmd
            s.sendcmd(cmd)
            res=s.getline()
            print res
            sys.stdout.flush()
            

    sys.stdout.flush()
    if "-q" in sys.argv:
        s.sendcmd("quit")
    else:
        s.sendcmd("bye")
        
except Exception, e:
    s.sendcmd("bye")
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast

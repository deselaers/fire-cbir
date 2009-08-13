#!/usr/bin/python
import sys, socket, re, os, string, time
class mysocket:
    '''Socket stolen somewhere from the net'''
    def __init__(self, sock=None):
        if sock is None:
            self.sock = socket.socket(
                socket.AF_INET, socket.SOCK_STREAM)
        else:
            self.sock = sock
    def connect(self,host, port):
        self.sock.connect((host, port))
    def send(self,msg):
        totalsent = 0
        MSGLEN=len(msg)
        while totalsent < MSGLEN:
            sent = self.sock.send(msg[totalsent:])
            if sent == 0:
                raise RuntimeError, "socket connection broken"
            totalsent = totalsent + sent
    def sendcmd(self,cmd):
        cmd = cmd + "\r\n"
        self.send(cmd)
    def getResult(self):
        result=""
        line=""
        while line[0:3] != "end":
            line=self.getline()
            if line[0:3] != "end":
                result=result+line
        return result
    def receive(self):
        msg = ''
        while len(msg) < MSGLEN:
            chunk = self.sock.recv(MSGLEN-len(msg))
            if chunk == '':
                raise RuntimeError,  "socket connection broken"
            msg = msg + chunk
        return msg
    def getline(self):
        msg=''
        chunk='a'
        while chunk != '\n':
            chunk=self.sock.recv(1)
            if chunk == '':
                raise RuntimeError, "socket connection broken"
            if chunk != '\n':
                msg=msg+chunk
        return msg

#-----------------------------------------------------------------
s=mysocket()
try:
    host=sys.argv[1]
    if len(sys.argv) >=3:
        port=int(sys.argv[2])
        try:
            s.connect(host,port)
        except:
            print "Connecting to server "+host+" on port "+str(port)+" failed."
            
        s.sendcmd("listfiles")
        fileline=s.getline()
        files=string.split(fileline," ")
        s.sendcmd("setresults 2")
        s.getline()
        
        classified=0
        correct=0
        error=0
        
        for file in files:
            if len(file)>0:
                cmd="retrieve +"+file
                s.sendcmd(cmd)
                res=s.getline()
                parts=string.split(res," ")
                result=parts[2]
                
                print res
                
                #get classes
                s.sendcmd("class "+file)
                qclass=s.getline()
                s.sendcmd("class "+result)
                rclass=s.getline()
                
                qc=int(qclass)
                rc=int(rclass)
                
                if(rc==qc):
                    correct=correct+1
                else:
                    error=error+1
                
                classified=classified+1
                
                er=float(error)/float(classified)*100.0
                print file+" wrong="+str(error)+" correct="+str(correct)+" classified="+str(classified)+" ER="+str(er)
            
        s.sendcmd("bye")
        
except Exception, e:
    s.sendcmd("bye")
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast
    



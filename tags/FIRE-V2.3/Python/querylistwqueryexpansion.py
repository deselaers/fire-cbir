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

    if len(sys.argv) != 4:
        print """USAGE:
        querylistwqueryexpansion.py <querylist> <server> <port>
        """
    else:
        ql=open(sys.argv[1],"r")
        host=sys.argv[2]
        port=int(sys.argv[3])
        try:
            s.connect(host,port)
        except:
            print "Connecting to server "+host+" on port "+str(port)+" failed."
            
            
        files=ql.readlines()
        ql.close()
        s.sendcmd("setresults 10")
        s.sendcmd("setextensions 2")
    
        for file in files:
            file=re.sub('\n','',file)
            cmd="retrieveandsaveranks 1030 "+file+".ranks "+file
            s.sendcmd(cmd)
            res=s.getline()
            print res
            
        s.sendcmd("bye")
            
except Exception, e:
    s.sendcmd("bye")
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast

import socket, os, sys,re 
# ----------------------------------------------------------------------
# socket implementation for comunication with the server system
# running on any computer which is reachable via network sockets
# here we need mainly the functions
# * connect
# * sendcmd
# * getline
# * close
# for the communication
# ----------------------------------------------------------------------
class FIRESocket:
    '''Socket stolen somewhere from the net'''
    def __init__(self, sock=None):
        if sock is None:
            self.sock = socket.socket(
                socket.AF_INET, socket.SOCK_STREAM)
        else:
            self.sock = sock
            #self.sock.settimeout(10.0)
    def connect(self,host, port):
        self.sock.connect((host, port))
    def send(self,msg):
        #print "PUT:",msg; sys.stdout.flush()
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
        #print "GOTR:",result; sys.stdout.flush()
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
            if chunk!= '\r' and chunk != '\n':
                msg=msg+chunk
        msg=re.sub("[ ]$","",msg)
        #print "GOTL:",msg; sys.stdout.flush()
        return msg
    def flush(self):
        chunk=self.sock.recv(5000)
        return

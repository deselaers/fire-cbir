#! /usr/bin/python
#-----------------------------------------------------------------------------
# This file is part of the FIRE -- Flexible Image Retrieval System
#
# FIRE is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# FIRE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FIRE; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#----------------------------------------------------------------------------

__author__ = "Jens Forster <jens.forster@rwth-aachen.de>"
__version__= "0.1" 

import sys
import re
import string
from fireproxyserver import *
import socket
from firesocket import *


class FIREProxy:
    
    def __init__(self):
        self.commandLineServers = False
        self.helpBool = False
        self.port = 12963
        self.retrievers = []
        self.logging = "default"
        self.servers = []
        self.commands = ["-c","-l","+l","-p","-s","-n"]
        error = self.parseCommandline()
        if not self.helpBool and not error:
            self.start()
        
        
    def help(self):
        print " -h                 yields this help"
        print " -n <number>        specifies how many retrieval servers will interact with"
        print "                    the proxy. this option is mandatory if -s is not used"
        print "                    this means that -s or -n must be present"
        print "                    if -s and -n are present the servers specified after -s"
        print "                    will be used"
        print " -p <port>          port on which proxyserver will listen"
        print "                    default is 12963"
#        print " -c <config location> 
        print " -l                                 deactivates proxyserver logging to a file"
        print " +l <file location>                 proxyserver will log to spezified file"
        print " -s < retrieverip,retrieverport > a list of retrieval ip,port tupels"
        print "                                  format: ip,port ip,port ip,port etc." 
        self.helpBool = True
        
    def setport(self,port):
        self.port = port
    
    def setlog(self,location=None):
        if location == None:
            self.logging = None
        elif location == "default":
            self.logging = "default"
        else:
            self.logging = location
    
    def setservers(self,servers):
        self.servers += servers
            
    def parseCommandline(self):
        error = False
        arg = sys.argv
        arg.pop(0)
        if arg != []:
            if "-h" in arg:
                arg = []
                self.help()
            while arg != []:
                if "-p" in arg:
                    idx = arg.index("-p")
                    port = int(arg[idx+1])
                    self.setport(port)
                    arg.pop(idx+1)
                    arg.pop(idx)
                    del idx
                elif "-l" in arg:
                    arg.remove("-l")
                    self.setlog()
                elif "+l" in arg:
                    idx = arg.index("+l")
                    location = arg[idx+1]
                    self.setlog(location)
                    arg.pop(idx+1) 
                    arg.pop(idx)
                    del idx
                    del location
                elif "-s" in arg:
                    idx = arg.index('-s')
                    jdx = idx +1
                    leng = len(arg)
                    servers = []
                    while jdx < leng and arg[jdx] not in self.commands :
                        valList = re.split(",",arg[jdx])
                        ip = valList[0]
                        port = int(valList[1])
                        servers.append((ip,port))
                        jdx += 1
                    jdx = jdx - 1
                    while jdx >= idx:
                        arg.pop(jdx)
                        jdx = jdx -1
                    self.commandLineServers = True
                    self.setservers(servers)
                elif "-n" in arg:
                    idx = arg.index('-n')
                    if not self.commandLineServers:
                        error = self.listenForServers(int(arg[idx+1]))
                    arg.pop(idx+1)
                    arg.pop(idx)    
            return error 
        # ende der groﬂen if bedingung
        # no command line commands were given, so the proxyport is default
        # and the retrieval servers have to anounce themselfes
        else:
            # no commandline arguments given but nethertheless we need to know where the
            # retrieval servers are located
            # note that all other arguments are set to default
            self.help()
            return error

    def listenForServers(self,count):
        # first create a temporary socket for listenig
        address_family = socket.AF_INET
        socket_type = socket.SOCK_STREAM
        request_queue_size = 50
        ans = ""
        localserver = []
        localIncorrectAddr = []
        rightservers = []
        error = False
        self.socket = socket.socket(address_family,socket_type)
        self.socket.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
        self.socket.bind(('0.0.0.0',self.port))
        self.socket.listen(request_queue_size)
        rbuffsize = -1
        for idx in range(count):
            print("WAITING FOR RETRIEVAL SERVER "+str(idx)+" TO CONNECT")
            server, server_address = self.socket.accept()
            rfile = server.makefile('rb',rbuffsize)
            ans = rfile.readline()
            localserver.append(ans)
            localIncorrectAddr.append(server_address)
            rfile.close()
            server.close()
        self.socket.close()  
        # we are only intrested in the ports
        if localIncorrectAddr == [] or localserver == []:
            print "NO RETRIEVALSERVER DECLARED THEMSELVES"
            print "ABORTING STARTING SEQUENCE"
            error = True
            return(error)
        localserver = map(lambda s: string.split(s),localserver)
        # the number 3 is due to the answer "Here I am <port>" from the retrieval servers
        localserver = map(lambda s: int(s[3]),localserver)
        #  now get the ips and combine them witch the right port
        for idx in range(len(localserver)):
            ip , falseport = localIncorrectAddr[idx]
            rightservers.append((ip,localserver[idx]))
        self.setservers(rightservers)
        return(error)    
        
    def start(self):
        proxyserver = FIREProxyServer(('0.0.0.0',self.port),self.servers,FIREProxyHandler,self.logging)        
        proxyserver.serve()
        
def main():
    proxy = FIREProxy()

if __name__ == "__main__":
    main()
    

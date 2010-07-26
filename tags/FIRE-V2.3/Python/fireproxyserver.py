#------------------------------------------------------------------------------
# proxyserver implementation for communication with the web frontend and fire 
# server system
# server must be reachable via network sockets in order to be of any use
# 
#------------------------------------------------------------------------------

__author__ = "Jens Forster <jens.forster@rwth-aachen.de>"
__version__= "0.1" 

from SocketServer import *
from firesocket import *

# Needed for stringhandling
# re provides regular expressions and methods
import re
import string

# Needed for suitable choosing of random images
# hence each retrieval server will only have access to unique
# part of the database
import random

# Needed for logging
import logging
#from logging.handlers import TimedRotatingFileHandler
import sys
import os
import os.path


__all__ =["FIREProxyServer", "FIREProxyHandler"]

class FIREProxyServer(TCPServer):
    
    request_queque_size = 1
    allow_reuse_address = True
    
    def __init__(self,server_address, retr_addrs, RequestHandlerClass,logDestination=None):
        TCPServer.__init__(self,server_address,RequestHandlerClass)
        if len(retr_addrs) == 0:
            raise RuntimeError, "no retrieval servers spezified"
        
        else:
            self.retr_addrs = retr_addrs
            self.retr_sock = []
        # initializing loggin system
        self.log = self._server_log(logDestination)
        self.log.log(logging.INFO,'Initializing proxyserver')
        self.log.log(logging.DEBUG,'Logging system and retrieval servers list ready')
        # initializing boolean values for request and server termination
        # and user authorization
        self.notBye = True # true -> web interface didn't quit yet
        self.notQuit = True # true -> proxy and retrieval servers still running
        self.notProxyquit = True # false -> proxy going down, retrieval servers will remain online
        self.auth = False
        self.sendErrorsFrom = [] # a list which will contain retrieval servers which crashed or aren't reachable anymore
        self.results = 0 # given results per retrieval server
        self.log.log(logging.INFO,"Proxyserver running on "+str(server_address))
    
    def sendToRetrievers(self,cmd,upperBound=None):
        if upperBound != None:
            if upperBound < 0:
                raise RuntimeError, "no retrieval servers specified for sending, Proxyserver: method sendToRetrievers: upperBound "+str(upperBound)
            relevantRetr = zip(self.retr_sock,range(upperBound))
            for sock,idx in relevantRetr:
                try:
                    sock.sendcmd(cmd)
                    self.log.log(logging.getLevelName("SEND "+str(self.retr_addrs[idx])),cmd)
                except:
                # iff no socket error occured until now
                    name,port = self.retr_addrs[idx]
                    if self.sendErrorsFrom == []:
                        self.sendErrorsFrom.append(str(name))
                        self.log.log(logging.CRITICAL,"socket connection to retrieval server "+str((name,port))+" did break" )
                        raise RuntimeError, "socket connection to retrieval server "+str((name,port))+" did break"
                    elif str(name) in self.sendErrorsFrom:
                        pass
                    else:
                        self.sendErrorsFrom.append(str(name))
                        self.log.log(logging.CRITICAL,"socket connection to retrieval server "+str((name,port))+" did break")
        else:
            for sock,addr in zip(self.retr_sock,self.retr_addrs):
                try:
                    sock.sendcmd(cmd)
                    self.log.log(logging.getLevelName("SEND "+str(addr)),cmd)
                except:
                    name,port = addr
                    self.log.log(logging.INFO,"sendErrorsFrom: "+str(self.sendErrorsFrom))
                    if self.sendErrorsFrom == []:
                        self.sendErrorsFrom.append(str(name))
                        self.log.log(logging.CRITICAL,"socket connection to retrieval server "+str(addr)+" did break")
                        raise RuntimeError, "socket connection to retrieval server "+str(addr)+" did break"
                    elif str(name) in self.sendErrorsFrom:
                        self.log.log(logging.DEBUG,"error in socket communication for retrieval server "+str(name)+" already known")
                    else:
                        self.sendErrorsFrom.append(str(name))
                        self.log.log(logging.CRITICAL,"socket connection to retrieval server "+str(addr)+" did break")
    
    def getFromRetrievers(self):
        answer = []
        for sock,addr in zip(self.retr_sock,self.retr_addrs):
            try:
                ans = sock.getline()
            except:
                self.log.log(logging.CRITICAL,"retrieval server "+str(addr)+" doesn't respond")
                name,port = addr
                self.sendErrorsFrom.append(str(name))
                raise RuntimeError, "retrieval server "+str(addr)+" didn't respond in method getFromRetrievers"
            answer.append(ans)
            self.log.log(logging.getLevelName("RECV "+str(addr)),ans)
        return answer
          
          
    def handle_request(self):
        """Handle one request in terms of one connection possibly containing
        several actions, possibly blocking."""
        errorCode = 0
        try:
            request, client_address = self.get_request()
            self.log.log(logging.INFO,"client "+str(client_address)+" connected")
        except:
            self.notQuit = False
            self.log.log(logging.ERROR,"an error occured during client's connection attempt")
            self.log.log(logging.ERROR,"maybe the socket/pipe did break or there are network problems")
            self.log.log(logging.ERROR,"proxyserver shutting down (retrieval server system may still be online)")
            try:
                request.close()
                self.socket.close()
            except:
                pass
            errorCode = 1
            return errorCode
        else:
            try:
                for serverNum in range(len(self.retr_addrs)):
                    serverNum_sock = FIRESocket()
                    self.retr_sock.append(serverNum_sock)
            except:
                self.notQuit = False
                self.log.log(logging.ERROR,"couldn't create necessary network sockets for communication with retrieval servers")
                self.log.log(logging.ERROR,"proxyserver shutting down (retrieval server system may still be online)")
                request.close()
                self.socket.close()
                while self.retr_sock != []:
                    firesock = self.retr_sock.pop(0)
                    firesock.sock.close()
                    del firesock
                errorCode = 2
                return errorCode
            else:
                self.log.log(logging.DEBUG,"network sockets for communication with retrieval servers created")
        if self.verify_request(request, client_address):
            try:
                self.process_request(request, client_address)
                return errorCode
            except RuntimeError, msg:
                self.handle_error(request,msg)
                errorCode = 3
                return errorCode
            except:
                self.handle_error(request,"an error involving networksockets or an unknown error occured")
                errorCode = 3
                return errorCode
                
    def process_request(self , request, client_addr):
        self.finish_request(request, client_addr)
        try:
            self.close_request(request)
        except:
            pass
       
    def close_request(self, client):
        client.close()
        self.log.log(logging.INFO,"client disconnected")
        self.auth = False
        # the network sockets designated for communication with retrieval servers
        # will only be closed if bye, quit or an error are encountered
        if self.notQuit:
            self.sendToRetrievers("bye")
            while self.retr_sock != []:
                firesock = self.retr_sock.pop(0)
                firesock.sock.close()
                del firesock
    
    def serve(self):
        self.log.log(logging.INFO,"proxy handling requests")
        # iff the error code is True no error occured while handling
        errorCode = 0
        while self.notQuit and self.notProxyquit :
            errorCode = self.handle_request()
            self.log.log(logging.DEBUG,"errorCode "+str(errorCode))
        if errorCode == 0:
            self.server_close()
        else:
            logging.shutdown()
    
    def server_close(self):
        # no auth necessary due to prior authentification
        if self.notProxyquit:
            self.log.log(logging.INFO,"sending shutdown signal to retrieval servers")
            self.sendToRetrievers("quit")
        self.socket.close()
        # closing and deleting network sockets
        while self.retr_sock != []:
            sock = self.retr_sock.pop(0)
            sock.sock.close()
            del sock
        if self.notProxyquit:
            self.log.log(logging.INFO,"proxyserver and retrieval servers stopped")
        else:
            self.log.log(logging.INFO,"proxyserver stopped, retrieval servers running")
        self.log.log(logging.DEBUG,"shutting logging system down")
        logging.shutdown()
   
    # printing the error message in to the logging file, closing connection to all retrieval servers
    # and shutting down
    def handle_error(self,request,msg):
       self.notQuit = False
       self.log.log(logging.ERROR,"an error occured while handling requests")
       self.log.log(logging.ERROR,"for further information view the log file")
       self.log.log(logging.DEBUG,msg)
       self.log.log(logging.ERROR,"shutting proxyserver down")
       request.close()
       self.sendToRetrievers("bye")
       self.socket.close()
       while self.retr_sock != []:
           sock = self.retr_sock.pop(0)
           sock.sock.close()
           del sock
       del self.retr_sock
       # logging.shutdown()
           
    # Note that _server_log will always return a logger which will always log global statements to stdout
    # and if intended will log everything to a file depending on the message's importance level 
    def _server_log(self,logDestination):
        if logDestination == None:
            location = None
        else:
            # exists is only needed if in future an other file logger should be used
            # that doesn't check fileexistance
            location, exists = self._logFilehandling(logDestination)
        # Setting up the logging System
        # The root logger shall log to stdout, all debug msg will be found in the files
        # Adding new message levels to the logging system
        # DEBUG (10) and NOTSET (0) needn't to be altered
        logging.addLevelName(1020,"CRITICAL")
        logging.addLevelName(1010,"ERROR")
        logging.addLevelName(1005,"WARNING")
        logging.addLevelName(1000,"INFO")
        logging.addLevelName(11,"RECV (Web/Bash)")
        logging.addLevelName(12,"SEND (Web/Bash)")
        logging.INFO = 1000
        logging.ERROR = 1010
        logging.CRITICAL = 1020
        logging.WARNING = 1005
        # so the levels 13 up to 999 are free for retrieval servers
        # because every server needs two numbers this logging system is capable of 493 retrieval servers
        # now levels for the retrieval servers are added
        # they are formated as SEND/RECV (IP-Address/Name) 
        level = 13
        for addrs in self.retr_addrs:
            logging.addLevelName(level,"SEND "+str(addrs))
            logging.addLevelName(level+1,"RECV "+str(addrs))
            level += 2
        # now configurating the root logger which will log only to stdout
        sys.stdout = sys.stderr
        root = logging.getLogger('')
        root.setLevel(logging.DEBUG)
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(levelname)-10s %(message)s','%a, %d-%m-%Y %H:%M:%S')
        console.setFormatter(formatter)
        root.addHandler(console)
        # setting up the main file logger
        if location != None:
#            when = 'midnight'
#            filesPerDay = 1
#            backUpCount = 7
#            file = TimedRotatingFileHandler(location,when,filesPerDay,backUpCount)
            if exists:
                file = logging.FileHandler(location,'a')
            else:
                file = logging.FileHandler(location,'w')
            file.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s %(levelname)-6s %(message)s','%a, %d-%m-%Y %H:%M:%S')
            file.setFormatter(formatter)
            root.addHandler(file)
        return root
        
    def _logFilehandling(self,logDestination):
        if logDestination == "default":
            cwd = os.getcwd()
            # determinating whether mac, unix/linux or windows
            if os.name in ["mac","nt","posix"]:                 
                if os.name == "nt":
                    # windows os
                    osPathTag = "\\"
                else:
                    # unix/linux os or macintosch
                    osPathTag = "/"
                try:
                    # checking whether dir exists
                    location = cwd+osPathTag+"fireproxylogs"
                    exist = os.path.isdir(location)
                    if not exist:
                        os.mkdir(location)
                    location += osPathTag+"proxylog"
                    del exist
                    del osPathTag
                except:
                    raise RuntimeError, "error initalizing logging system (default logging location)" 
            else:
                raise RuntimeError, "OS not supported for logging"        
            del cwd
        # iff the user specified a storing location
        # outer if-clause
        else:
            # if only a directory is given (logDestination must be terminated by the os's
            # path separator (e.g / on linux \ on windows os)) default logfile name will be used
            head,tail = os.path.split(logDestination)
            osType = os.name
            if osType in ["mac","posix"]:
                osPathTag = "/"
            elif osType == "nt":
                osPathTag = "\\"
            else:
                raise RuntimeError, "OS not supported for logging"
            # checking whether a path name is given
            if head != "":
                # checking whether the user set all path separators' right
                head = os.path.normcase(head)
                head = os.path.normpath(head)
                # figuring out whether the supplied directories exist and
                # if not which have to be created
                if not os.path.isdir(head):
                    dirToCreate = [head]
                    newHead, newTail = os.path.split(head)
                    while newTail != "" and (not os.path.isdir(newHead)):
                        dirToCreate.append(newHead)
                        newHead, newTail = os.path.split(newHead)
                    while dirToCreate != []:
                        os.mkdir(dirToCreate.pop())
                    del newHead
                    del newTail
                # path is valid now checking if a filename was supplied
                if tail != "":
                    location = logDestination
                else:
                    location = logDestination+osPathTag+"proxylog"
                del dirToCreate
            # no path given by the user
            # tail must be != "" in this case otherwise it would be the default case
            else:
                # assuming current workdirectory was intended
                cwd = os.getcwd()
                location = cwd+osPathTag+tail
            del osPathTag
            del osType
            del head
            del tail
        # Outer If Clause ends here
        # testing if the file exist or not
        exists = True
        try:
            file = open(location,'r+')
        except IOError:
            exists = False
        else:
            file.close()
        return (location,exists)
                
    # Overidden to be useless, because this proxy shall quit when all retrievers quit
    def serve_forever(self):
        pass
    
class FIREProxyHandler(StreamRequestHandler):
    
      
    def __init__(self, request, client_addr, server):
        StreamRequestHandler.__init__(self,request, client_addr, server)
        
    def sendToWeb(self,msg):
        self.wfile.write(msg+"\r\n")
        self.server.log.log(logging.getLevelName("SEND (Web/Bash)"),msg)
        
    def _validateSettings(self,settings):
        self.server.log.log(logging.DEBUG,"valSet Argu:"+repr(settings))
        leng = len(settings)
        self.server.log.log(logging.DEBUG,"len: "+repr(leng))
        if leng:
            if settings[0][0] != "filelist":
                retur = (False,"Settingsstream starts with wrong keyword")
            else:
                retur = (True, "")
            if leng > 1:
                lenSingle = map(lambda single: len(single),settings)
                boolList = map(lambda single: single == lenSingle[0],lenSingle)
                if not (reduce(lambda fst,snd: fst == snd,boolList)):
                    retur = (False,"At least one retrieval server didn't send complete settings")
                # Fine all retrieval servers did send there complete settings
                # or all made the same mistake :D
                lenSetting = lenSingle[0]     
                # deleting helper variables
                del boolList
                # checking rest
                for idx in range(lenSetting):
                    if (idx != 1):
                        for jdx in range(leng):
                            if (settings[0][idx] != settings[jdx][idx]):
                                # settings differ and singaling it
	         		return  (False,"settings of servers "+str(0)+" and "+str(jdx)+" differ in "+ settings[0][idx]+" "+settings[jdx][idx])
	         	    else:
	                        retur = (True,"")	
            return retur
        else:
             return (False,"no info information at all received")
    
    def setup(self):
        StreamRequestHandler.setup(self)
        error = False
        try:
            for retr, addr in zip(self.server.retr_sock,self.server.retr_addrs):
                   name, port = addr
                   try:
                       retr.connect(name,port)
                       self.server.log.log(logging.DEBUG,"connection to retrival server "+str(addr)+" established")
                   except:
                       # sending error message to web frontend
                       self.sendToWeb("failure "+str(name)+" down")
                       self.server.sendErrorsFrom.append(str(name))                     
                       # closing
                       error = True
                       self.server.log.log(logging.CRITICAL,"retrieval server "+str(addr)+" not reachable")
        except: 
            pass
        else:
            if error:
                self.finish()
                raise RuntimeError, "some of the retrieval servers arent't reachable via network"
            else:
                self.server.notBye = True
                # signal: connection to retrieval servers established
                # currently unix specific linefeed
                self.server.log.log(logging.INFO,"connection to retrieval server system established")
               
     
    
    def handle(self): 
        retrieveAndMetaCommands = ["retrieve","expand","retrieveandsaveranks","metaretrieve","metaexpand","textretrieve","textexpand"]
        # handling actions until "bye", "quit" or doomsday            
        while self.server.notBye and self.server.notQuit:
            
            webMsg = self.rfile.readline()
            self.server.log.log(logging.DEBUG,webMsg)
            webMsg = re.sub("[\r\n\t\v\f]","",webMsg)
            self.server.log.log(logging.getLevelName("RECV (Web/Bash)"),webMsg)
            webKeyS = string.split(webMsg)
            keyWord = webKeyS.pop(0) # grabbing the very first element
            if keyWord == "info":
                self.server.sendToRetrievers("info")
                answer = self.server.getFromRetrievers()
                # convert list of single strings in usable format
                # e.g ["hello world","how do you do?"] becomes
                # [['hello','world'],['how','do','you','do?']]
                settings = map(lambda strin: string.split(strin),answer)
                bool, comment = self._validateSettings(settings)
                if bool:
                    sum = 0
                    # summing up total databank size
                    for idx in range(len(settings)):
                        sum += int(settings[idx][1])
                    settings[0][1] = sum
                    self.server.results = int(settings[0][3])
                    self.sendToWeb(re.sub("[[,'\]]","",str(settings[0])))
                    # deleting helper variables
                    del bool
                    del answer
                    del comment
                    del settings
                else:
                    self.finish()
                    raise RuntimeError, "handle():"+comment
            elif keyWord == "password":
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                tokens = map(lambda strin: string.split(strin),answer)
                bool = True
                idx = 0
                while bool and idx < len(answer):
                    if tokens[idx][0] != "ok":
                        bool = False
                    idx += 1
                if not bool:
                    self.sendToWeb(answer[idx-1])
                    self.server.auth = False
                else:
                    self.sendToWeb("ok")
                    self.server.auth = True
                # deleting helper variables
                del answer
                del tokens
                del idx
                del bool
            elif keyWord == "bye":
                self.server.notBye = False
            elif keyWord == "quit":
                self.server.notBye = False
                if self.server.auth:
                    self.server.notQuit = False
                else:
                    self.server.sendToRetrievers("bye")
            elif keyWord in retrieveAndMetaCommands:
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                self.server.log.log(logging.DEBUG,"answer(handle()): "+repr(answer))
                # converting retrieval servers' answers from strings to list of
                # strings like above and afterwards merging all such lists into 
                # one
                if keyWord in ["metaretrieve","metaexpand","textretrieve","textexpand"]:
                    listed = map(lambda string: re.split("[ :]",string),answer)
                else:
                    listed = map(lambda strin: string.split(strin),answer)
                listed = reduce(lambda fst,snd: fst+snd,listed)
                if "nometainformation" in listed:
                    self.sendToWeb("nometainformation")
                elif "notextinformation" in listed:
                    self.sendToWeb("notextinformation")
                elif "" in listed:
                    self.server.log.log(logging.INFO,"at least one retrival server answered with an empty string")
                    self.server.log.log(logging.ERROR,"retrival aborted")
                elif not self.server.results == 0:
                    # for sorting purposes changing data structure to a list of tupels
                    # consistent of distance as a float and filename as a string
                    toSort = [(float(listed[i+1]),listed[i]) for i in range(0,len(listed),2)]
                    toSort.sort()
                    toSort.reverse()
                    retrievalString = ""
                    for tup in toSort[:self.server.results]:
                        dist, name = tup
                        if keyWord in ["metaretrieve","metaexpand","textretrieve","textexpand"]:
                            retrievalString += name+":"+str(dist)
                        else:
                            retrievalString += name+" "+str(dist)
                        retrievalString += " "
                    # removing trailing blank
                    retrievalString = re.sub(" $","",retrievalString)
                    self.sendToWeb(retrievalString)
                    #deleting helper variables
                    del toSort
                    del retrievalString
                else:
                    self.server.log.log(logging.INFO,"Proxyserver assumes that settings for results equals 0")
                    self.server.log.log(logging.INFO,"Please use setresults to set results to an usefull value")
                    self.server.log.log(logging.DEBUG,"possible mismath between retrieval servers settings and proxyserver")
                    self.server.log.log(logging.ERROR,"retrieval aborted")
                # deleting helper variables
                del answer
                del listed
            elif keyWord == "help":
                self.server.sendToRetrievers(webMsg,1)
                answer = self.server.retr_sock[0].getline()
                ans = string.split(answer)
                try:
                    ans.remove("filelist")
                except:
                    pass
                answer = "proxy only: quitproxy ||normal fire commands: "
                for cmd in ans:
                    answer += cmd
                    answer += " "
                answer += "||not supported in proxy mode: filelist"
                self.sendToWeb(answer)
                del answer
            elif keyWord == "savedistances":
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                fst = answer[0]
                bool = True
                for name in answer:
                    if name != fst:
                        bool = False
                if not bool:
                    self.finish()
                    raise RuntimeError,"distance filenames differ"
                self.sendToWeb(fst)
                del answer
                del bool
                del fst
            elif keyWord == "random":
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                tokens = map(lambda strin: string.split(strin),answer)
                tokenList = reduce(lambda fst,snd:fst+snd,tokens)
                # choosing randomly images from all random images
                # because every retrieval server only accesses an unique database
                # subset
                randomImg = random.sample(tokenList,self.server.results)
                self.sendToWeb(re.sub("[[',\]]","",str(randomImg)))
                # deleting helper variables
                del answer
                del tokens
                del tokenList
                del randomImg
            elif keyWord == "saverelevances":
                self.server.sendToRetrievers(webMsg)
            elif keyWord in ["setscoring","setresults","setextensions","setdist","setweight"]:
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                tokens = map(lambda strin: string.split(strin),answer)
                keyToken = re.sub("set","",keyWord)
                # checking if all retrieval servers accepted new settings
                # first checking keyword
                boolList = map(lambda list: list[0] == keyToken,tokens)
                boolAll = reduce(lambda fst,snd:fst == snd,boolList)
                # checking iff snd argument must be checked
                # in case it will be "=" testing is omitted
                if keyToken in ["weight","dist"]:
                    sndArg = tokens[0][1]
                    boolList = map(lambda list: list[1] == sndArg,tokens)
                    boolSnd = reduce(lambda fst,snd: fst == snd,boolList)
                    boolAll = boolAll and boolSnd
                    del boolSnd
                # checking third argument
                trdArg = tokens[0][2]
                boolList = map(lambda list: list[2] == trdArg,tokens)
                boolTrd = reduce(lambda fst,snd:fst == snd,boolList)
                boolAll = boolAll and boolTrd
                if boolAll:
                    self.sendToWeb(answer[0])
                    if keyToken == "results" and tokens[0][0] != "Invalid":
                        self.server.results = int(tokens[0][2])
                else:
                    self.server.sendToRetrievers(keyWord,1)
                    alterMsg = self.server.retr_sock[0].getline()
                    self.sendToWeb(alterMsg)
                    del alterMsg
                del boolAll
                del boolTrd
                del boolList
                del answer
                del tokens
            elif keyWord == "listfiles":
                self.server.sendToRetrievers(webMsg)
                answer = self.server.getFromRetrievers()
                result = ""
                for idx in range(len(answer)):
                    result += answer[idx]
                    result += " "
                # deleting trailing blank
                result = re.sub(" $","",result)
                self.sendToWeb(result)
                del answer
                del result
            elif keyWord == "class":
                self.server.sendToRetrievers(webMsg)
                self.server.log.log(logging.INFO,webMsg+" send")
                answer = self.server.getFromRetrievers()
                tokens = map(lambda strin: string.split(strin),answer)
                tokens = reduce(lambda fst,snd: fst+snd,tokens)
                # since retrieval servers will only return a single value 
                # it is only neccessary to destinguish the three possible modes of the class command
                # for error handling
                # do all retrieval servers have classes or don't they
                if "yes" in tokens:
                    for item in tokens:
                        if "yes" != item:
                            self.finish()
                            raise RuntimeError, "not all retrieval servers do have classes"
                    self.sendToWeb("yes")
                elif "no" in tokens:
                    for item in tokens:
                        if "no" != item:
                            self.finish()
                            raise RuntimeError, "not all retrieval servers don't have classes"
                    self.sendToWeb("no")
                else:
                    try:
                        intList = []
                        for item in tokens:
                            intItem = int(item)
                            if intItem >= 0:
                                intList.append(intItem)
                    except:
                        self.finish()
                        raise RuntimeError, "neiter yes/no nor numbers were returned for classes"
                    else:
                        boolList = map(lambda it: intList[0]==it,intList)
                        booL = reduce(lambda fst, snd: fst and snd,boolList)
                        if booL:
                            self.sendToWeb(str(intList[0]))
                        else:
                            self.finish()
                            raise RuntimeError, "retrieval servers classes differ"
                        del boolList
                        del booL
                        del tokens
                        del answer
            elif keyWord == "metafeatureinfo":
                self.server.sendToRetrievers(webMsg)
                self.server.log.log(logging.INFO,webMsg+" send")
                answer = self.server.getFromRetrievers()
                answer = map(lambda string: re.sub(" $","",string),answer)
                returnstring = ""
                for strings in answer:
                    returnstring += strings
                    returnstring += " "
                self.sendToWeb(re.sub(" $","",returnstring))
                del answer
                del returnstring
            # Note that the commands filelist is dumped on purpose
            # because every retrieval server accesses it's own databade subset
            # The interactor command is dumped too.
            elif keyWord in ["filelist","interactor"]:
                self.sendToWeb(keyWord+" is not avaible since fireproxy is in use")
            elif keyWord == "quitproxy":
                self.server.notBye = False
                self.server.notProxyquit = False
            else:
                self.sendToWeb("Unknown command: "+webMsg)
        del retrieveAndMetaCommands
        del webKeyS
        del webMsg

#!/usr/bin/python

import sys, socket, re, os, string, time, traceback
#sys.path.append("../python")
import firesocket

def retrieveData(param):
  data = {}
  s=firesocket.FIRESocket()
  try:
        ql=open(sys.argv[1],"r")
        host=sys.argv[2]
        port=int(sys.argv[3])
        useX="-x" in sys.argv
        verbose="-v" in sys.argv
        try:
            time.sleep(1) #don't kill the server by exiting tooo fast 
            s.connect(host,port)
        except:
            print "Connecting to server "+host+" on port "+str(port)+" failed."

        fileprefix=""
        if len(sys.argv)>4:
            if sys.argv[4]!="-q" and sys.argv[4]!="-x" and sys.argv[4]!="-v":
                fileprefix=sys.argv[4]

        lines=ql.readlines()

        # get number of files in retriever
        s.sendcmd("info")
        res=s.getline()
        status=re.split(" ",res)
        status.reverse()
        keyword=status.pop()
        dbSize=int(status.pop())
        print "database has",dbSize,"images."
        s.sendcmd("setresults "+str(dbSize)); res=s.getline(); 
        ql.close()
        
        classified=0

        ###query images and evaluate results
        for line in lines:
            line=re.sub('\n','',line)
            relevant=re.split(" ",line)
            query=fileprefix+relevant[0] # this is the name of the query image
            relevant=relevant[1:] # this are the names of the images relevant to the query image
            cmd="retrieve +"+query
            s.sendcmd(cmd) 
            res=s.getline() # get the results from the retriever
            res=re.sub('[ ]*$','',res) # and parse the results
            results=re.split(" ",res)
            results.reverse()
            returnedimages=[] 
            while len(results)>0:
                returnedimage=results.pop()
                score=results.pop()
                returnedimage=re.sub(fileprefix,"",returnedimage)
                returnedimages=returnedimages+[returnedimage]

            # --------------------------------------------------------------
            # performance evaluation
            # --------------------------------------------------------------
            
            classified=classified+1

            if returnedimages[0] != query:
                # if the query image is not the first returned, we assume that it is
                # not part of the database, and thus we have to start with the first
                # returned image for query evaluation. This is a workaround, usually we start
                # at position 1, to keep this valid we insert a new one on position 0
                # this makes the old position 0 to position 1 and everything hereafter
                # stays valid
                returnedimages.insert(0,"DUMMYIMAGE")

            
        
            key = classified
            data[key] = [query, relevant, returnedimages]
            #print "key", key, "data", data[key]
            

            
            print "query:", query, "image no.:", classified, "of", dbSize, "relevant:", relevant[:5], "returnedimages:", returnedimages[:5]
            if verbose:
                print "results:",returnedimages[0],"->0",
            sys.stdout.flush()
            
        s.sendcmd("info")
        res=s.getline()
        print "SETTINGS:",res
        
    
        sys.stdout.flush()
        time.sleep(1)
        if "-q" in sys.argv:
            s.sendcmd("quit")
        else:
            s.sendcmd("bye")
        time.sleep(1) #don't kill the server by exiting tooo fast 
  except KeyboardInterrupt, e:
    s.sendcmd("bye")
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast 
  except Exception, e:
    if "-q" in sys.argv:
        s.sendcmd("quit")
    else:
        s.sendcmd("bye")
    traceback.print_exc(file=sys.stdout)
#    t=sys.last_traceback
#    traceback.print_tb(t)
    print e

    time.sleep(1) #don't kill the server by exiting tooo fast 
    
  param["data"] = data
  param["dbSize"] = dbSize
  param["useX"] = useX
  param["verbose"] = verbose
 
      
def calculate(param):
  
  data = param["data"]
  dbSize = param["dbSize"]
  useX = param["useX"]
  verbose = param["verbose"]

  # init variables for quantitative comparison
  classified=0; correct=0; error=0; AMP=0.0; ARE=0.0; ARP05=0.0;APRP=0.0; MAP=0.0
  AP1=0; AP20=0; AP50=0; APNrel=0; AR100=0; ARank1=0; ARanktilde=0;
  APRarea=0.0; APRgraph=[(0.0,0),(0.1,0),(0.2,0),(0.3,0),(0.4,0),(0.5,0),(0.6,0),(0.7,0),(0.8,0),(0.9,0),(1.0,0)]
  if useX:
    g=Gnuplot.Gnuplot()
    g.title('PR graph')
    g('set data style linespoints')
    g('set xrange [0.0:1.0]')
    g('set yrange [0.0:1.0]')  
    
  data.keys().sort()
  for j in data.keys():
    #print "j", j
    relevant = data[j][1]
    returnedimages = data[j][2]
      
    # calculate performance measures: Error rate
    classified=classified+1
    if returnedimages[1] in relevant:
      correct=correct+1
    else:
      error=error+1

    # calculate performance measures: ARE, AMP, Rank1, Ranktilde
    MPQ=0.0; REQ=0.0; MPQ_div=0.0; REQ_div=0.0; i=1; RP05=-1.0; PRP=-1.0; AvPr=0.0;

    for Ri in range(1,len(returnedimages)):
      if returnedimages[Ri] in relevant:
        AvPr+=(float(i)/float(Ri))
        MPQ+=float(dbSize-Ri)/float(dbSize-(i))
        REQ+=float(Ri); REQ_div+=float(i)
        #print "i:", i
        if i==1:
          Rank1=Ri
        i+=1
        if verbose:
          print returnedimages[Ri],"->",Ri,
          
    print
    Ranktilde=float(REQ-((len(relevant)*(len(relevant)-1))/2))/float(dbSize*len(relevant))

    MPQ*=float(100)/float(len(relevant))
    if REQ_div!=0:
      REQ/=REQ_div
    AMP+=MPQ
    ARE+=REQ
    AvPr*=float(100)/float(len(relevant))      

    # calculate performance measures: P,R, PR-graph
    Nrelret=0;
    PR=[]
    for Ri in range(1,len(returnedimages)):
      if returnedimages[Ri] in relevant:
        Nrelret+=1
      P=float(Nrelret)/float(Ri)
      R=float(Nrelret)/float(len(relevant))
      PR+=[(P,R)]

    # smoothing
    Ris=range(1,len(returnedimages)-1)
    Ris.reverse()
    maxP=0
    for Ri in Ris:
      if PR[Ri][0]>maxP:
        maxP=PR[Ri][0]
      else:
        PR[Ri]=(maxP,PR[Ri][1])
      if PR[Ri][0]>0.5 and RP05<0:
        RP05=PR[Ri][1]
                    
      if PR[Ri][0]>PR[Ri][1] and PRP<0:
        PRP=PR[Ri][0]

      if RP05<0:
        RP05=0
      if PRP<0:
        PRP=0
 
    # pick the eleven values R=0.0,0.1,...,1.0 from the PR values for the PRgraph
    i=0
    PRgraph=[]
    for r in range(11):
        R=float(r)/float(10)
        while i<len(PR) and PR[i][1] < R:
            i+=1
        if i>len(PR)-1:
            i=len(PR)-1
        P=PR[i][0]
        PRgraph+=[(R,P)]
    if PRgraph[0][1] < PRgraph[1][1]:
        PRgraph[0]=(PRgraph[0][0],PRgraph[1][1])

    PRarea=0.0
    for i in PRgraph[1:-1]:
        PRarea+=i[1]
    PRarea+=0.5*(PRgraph[0][1]+PRgraph[-1][1])
    PRarea*=0.1        

    P1=PR[0][0]
    P20=0
    if len(PR)>20:
        P20=PR[20][0]
    P50=0
    if len(PR)>50:
        P50=PR[50][0]
    PNrel=PR[len(relevant)][0]
    if len(PR)>100:
        R100=PR[100][1]
    else:
        R100=1

    AP1   +=P1   
    AP20  +=P20  
    AP50  +=P50  
    APNrel+=PNrel
    AR100 +=R100 
    ARank1+=Rank1
    ARanktilde+=Ranktilde
    APRarea+=PRarea
    ARP05+=RP05
    APRP+=PRP
    MAP+=AvPr
    for i in range(len(APRgraph)):
        APRgraph[i]=(APRgraph[i][0],APRgraph[i][1]+PRgraph[i][1])

    tmpAPRgraph=[]
    for i in range(len(APRgraph)):
        tmpAPRgraph+=[(APRgraph[i][0],APRgraph[i][1]/float(classified))]

    if useX:
        g.plot(Gnuplot.Data(PRgraph,title="PRgraph current query"),Gnuplot.Data(tmpAPRgraph,title="PRgraph average"))
    print "   AMP: %1.2f" % (float(AMP)/float(classified)),
    print "ARE: %1.2f" % (float(ARE)/float(classified)),
    print "MAP: %1.2f" % (float(MAP)/float(classified)),
    print "ER: %1.2f:" % (float(error)*100/float(classified)),
    print "P1: %1.2f" % (float(AP1)/float(classified)),
    print "P20: %1.2f" % ( float(AP20)/float(classified)),
    print "P50: %1.2f" % ( float(AP50)/float(classified)),
    print "PNrel: %1.2f" % ( float(APNrel)/float(classified)),
    print "R100: %1.2f" % (float(AR100)/float(classified)),
    print "Rank1: %1.2f" % (float(ARank1)/float(classified)),
    print "Ranktilde: %1.2f" % (float(ARanktilde)/float(classified)),
    print "PRarea: %1.2f" % (float(APRarea)/float(classified)),
    print "R(P=0.5): %1.2f" % (float(ARP05)/float(classified)),
    print "P(R=P): %1.2f" % (float(APRP)/float(classified)),
    #print "PRgraph: ",PRgraph
    print "classfied: %i" % (classified)    
    
  for i in range(len(APRgraph)):
    APRgraph[i]=(APRgraph[i][0],APRgraph[i][1]/float(classified))

  print "RESULT:",
  print "   AMP: %1.2f" % (float(AMP)/float(classified)),
  print "ARE: %1.2f" % (float(ARE)/float(classified)),
  print "MAP: %1.2f" % (float(MAP)/float(classified)),
  print "ER: %1.2f:" % (float(error)*100/float(classified)),
  print "P1: %1.2f" % (float(AP1)/float(classified)),
  print "P20: %1.2f" % ( float(AP20)/float(classified)),
  print "P50: %1.2f" % ( float(AP50)/float(classified)),
  print "PNrel: %1.2f" % ( float(APNrel)/float(classified)),
  print "R100: %1.2f" % (float(AR100)/float(classified)),
  print "Rank1: %1.2f" % (float(ARank1)/float(classified)),
  print "Ranktilde: %1.2f" % (float(ARanktilde)/float(classified)),
  print "PRarea: %1.2f" % (float(APRarea)/float(classified)),
  print "R(P=0.5): %1.2f" % (float(ARP05)/float(classified)),
  print "P(R=P): %1.2f" % (float(APRP)/float(classified)),
  print "classfied: %i" % (classified)
  print "PRgraph: ",APRgraph
  
    
               
if __name__ == "__main__":

  if len(sys.argv) < 4:
    print """USAGE:
      querylist.py <querylist> <server> <port> [<fileprefix>] [-q] [-x]
         -q: quit the server after finished eval
         -x: show PR graph after each query on X display
         -v: verbose
         <fileprefix> should not be necessary
      """
  else:
    param = {}
    retrieveData(param)  
    calculate(param)
      
      
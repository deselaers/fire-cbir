#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import sys,cgi,re,Image,Gnuplot,string,os
sys.path+=[os.getcwd()] 
from types import *
from firesocket import *

try:
    form=cgi.FieldStorage()

    if form.has_key("server"):
        server=form["server"].value
    else:
        server="localhost"
    if form.has_key("port"):
        port=int(form["port"].value)
    else:
        port=12960

    if form.has_key("image"):
        image=form["image"].value
    else:
        image="100.jpg"

    if form.has_key("feature"):
        feature=form["feature"].value
    else:
        feature="1"

    s=FIRESocket()

    s.connect(server,port)
    s.sendcmd("feature "+image+" "+feature)
    line=s.getline()
    suffix=string.split(line)[1]
    lines=[]
    line=s.getline()
    while line!="end":
        lines+=[line]
        line=s.getline()

    s.sendcmd("bye")

    if suffix.endswith(".histo") or suffix.endswith(".histo.gz") or suffix.endswith("oldhisto") or suffix.endswith("oldhisto.gz"):
        for l in lines:
            if l.startswith("data"):
                tokens=string.split(l)
                tokens.pop(0)
                tokens=map(lambda t:
                           float(t),tokens)
                bins=len(tokens)
                maximum=max(tokens)
                g=Gnuplot.Gnuplot()
                g('set terminal png size 300,200')
                g('unset xtics')
                g('unset ytics')
                g('set data style boxes')
                print "Content-Type: image/png\n"
                sys.stdout.flush()
                g.plot(Gnuplot.Data(tokens))
    elif suffix.endswith(".vec.gz") or suffix.endswith(".oldvec.gz") or suffix.endswith(".vec"):
        print "Content-Type: text/html\n"
        sys.stdout.flush()
        print "<html><body>\n"
        for l in lines:
            if l.startswith("data"):
                print l[5:]
        print "</body></html>\n"
    elif suffix.endswith(".textID"):
        print "Content-Type: text/html\n"
        sys.stdout.flush()
        print "<html><body>\n"
        for l in lines:
            print l
        print "</body>\n</html>"
    elif  suffix.endswith("txt") or suffix.endswith("txt.gz"):
        print "Content-Type: text/html\n"
        sys.stdout.flush()
        print "<html><body>\n"
        for l in lines:
            print l
        print "</body></html>\n"
    elif suffix.endswith(".png"):
        for l in lines:
            if l.startswith("sizes"):
                tok=string.split(l)
                xsize=int(tok[1])
                ysize=int(tok[2])
                zsize=int(tok[3])
            if l.startswith("data"):
                tok=string.split(l)
                tok.pop(0)
                if zsize==1:
                    i=Image.new("RGB",(xsize,ysize))
                    for x in range(xsize):
                        for y in range(ysize):
                            g=int(float(tok.pop(0))*255)
                            i.putpixel((x,y),(g,g,g))
                elif zsize==3:
                    i=Image.new("RGB",(xsize,ysize))
                    for x in range(xsize):
                        for y in range(ysize):
                            r=int(float(tok.pop(0))*255)
                            g=int(float(tok.pop(0))*255)
                            b=int(float(tok.pop(0))*255)
                            i.putpixel((x,y),(r,g,b))
        print "Content-Type: image/png\n"
        sys.stdout.flush()
        i.save(sys.stdout,"png")
                
    else:
        print "Content-Type: text/html\n"
        sys.stdout.flush()

        print suffix
        print "feature not supported: ",suffix
        for l in lines:
            print l+"<br>"
        print "</body></html>\n"
        
    
except IOError:
    contenttype="text/html"
    print "Content-Type: "+contenttype+"\n"
    print "Something failed"


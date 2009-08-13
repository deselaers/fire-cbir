#!/usr/bin/python

from string import split, lower, strip
from random import shuffle

import sys
import os
from os.path import join, getsize

ARCH=os.environ["ARCH"]


# returns the string after a certain option string 
def getStringAfter(name, array, default):
  for i in range(len(array)):
    if array[i] == name:
      if (len(array) == i+1):
        return default
      else:
        return array[i+1]
  return default
  

  

def usage(commandline):
  print "note: before you use this program make \"make all\" in the firedirectory. "
  print commandline[0], "(-np|-npfromdb|-ndb|-ndbfromp|-showf|-e|-h) [options]"
  print
  print "mandatory options (choose exactly one of them):"
  print "  -h        shows this help"
  print "  -showf    show all possible features that can be extracted"
  print "  -e        extracts features. options: -limages -ldb -ld -selectf "
  print "            -p -db -d -c -f -q -n"
  print "  -np       creates a new purefiles file using the filenames of the "
  print "            files in the directory. options: -d -p"
  print "  -npfromdb creates a new purefiles file from the data of the dbfile."
  print "            options: -d -db -p"
  print "  -ndb      creates a new database file using the filenames of the "
  print "            files in the directory. options: -d -db"
  print "  -ndbfromp creates a new database file from the data of the purefiles "
  print "            file. options: -d -p -db"
  print
  print "options:"
  print "  -d dir    directory of where the files are. default: ."
  print "  -p path   path of an [existing] purefiles file. default: purefiles"
  print "  -db path  path of an [existing] database file. default: list"
  print "  -c path   configfile that shall be used. "
  print "            default:", os.path.dirname(sys.argv[0]) + "/features.conf"
  print "  -f dir    directory where fire resides. "
  print "            default:", os.path.dirname(sys.argv[0]) + "../../bin/"+ARCH+"/"
  print "  -q (day|week)"
  print "            put all jobs into the queuingsystem. default: no queuing"
  print "  -n        name of the jobs in the queue. default: extractor"
  print "  -selectf feature1 [feature2 [... [featureN]]"
  print "            extracts only feature1 ... featureN of the config file"
  print "  -v        verbose"
  print "  -limages image1 [image2 [... [imageN]]"
  print "            extracts features of image1 ... imageN only"
  print "  -ldb      extracts the images from the dbfile"
  print "            (hint: it is smarter to create a db file first, using" 
  print "            -ndb or -ndbfromp)"
  print "  -ld       extracts the images of the directory"
  print "            (hint: it is smarter to create a purefiles file first,"
  print "            using -np or -npfromdb)"

  
  
  
class Extractor:

  def __init__(self):
    self.datadir = "."
    self.purefiles = "purefiles"
    self.dbfile = "list"
    self.progdir = os.path.dirname(sys.argv[0])  # path of the tool
    self.firedir = self.progdir + "/../../bin/"+ARCH+"/"
    self.config = "features.conf"
    self.extensionlist = ["jpg", "png", "gif", "jpeg", "tif", "tiff"]  
    
    self.images = []
    self.features = {}
    self.selectedFeatures = []
    
    self.queue = None
    self.name = "extractor"

    
  # gets for the attributes
  
  def setPurefiles(self, new):
    self.purefiles = new    
      
  def setDBFile(self, new):
    self.dbfile = new
    
  def setDatadir(self, new):
    self.datadir = new
   
  def setFiredir(self, new):
    self.firedir = new
    
  def setConfigfile(self, new):
    self.config = new
    
  def setExtensionlist(self, newl):
    self.extensionlist = newl
    
  def setQueue(self, new):
    self.queue = new
    
  def setName(self, new):
    self.name = new
    
    
  # printing      
    
  def printInfo(self):
    print "config info of extractor"
    print "datadir", self.datadir
    print "purefiles", self.purefiles
    print "dbfile", self.dbfile
    print "progdir", self.progdir
    print "firedir", self.firedir
    print "configfile", self.config
    print "extensionlist", self.extensionlist
  
  # prints available features
  def printFeatures(self):
    for feature in self.features.keys():
      print feature, self.features[feature]
      
  def printImages(self):
    if self.images == []:
      print "error: no images loaded."
    else:
      print "images:"
      for image in self.images:
        print image    
    
    
  # loading 
      
  # load images from the directory
  def loadImagesFromDir(self):
    print "loading images from directory."
    self.images = []
    for root, dirs, files in os.walk(self.datadir):
      for file in files:
        if self.checkExtension(file):
          if (root == self.datadir):
            self.images.append(file)
          else:
            self.images.append(root[len(self.datadir)+1:] + "/" + file)  
    if self.images == []:
      print "error: no images found." 

  # loads images from database file 
  # (does not consider path information of the db file yet)         
  def loadImagesFromDB(self):   
    print "loading images from dbfile." 
    f = open(self.dbfile, "r")  
    lines = f.readlines()
    f.close()
    for line in lines:
      if strip(split(line)[0]) == "file":
        self.images.append(strip(split(line)[1]))
        #print "file:", strip(split(line)[1])      
    if self.images == []:
      print "error: no images found."
      
  # loads features from the config file
  def loadFeatures(self):
    print "loading features file."
    try:
      f = open(self.config, "r")  
    except IOError: 
      print "warning: no specific config file found. I try the one in the tooldir."
      f = open(self.progdir + "/features.conf")
    lines = f.readlines()
    f.close()
    for line in lines:
      if line[0] == "#":
        #print "comment: ", line
        pass
      elif len(strip(line)) < 1:
        pass
      else:
        key = strip(split(line)[0])
        #print "key", key
        feature = strip(line)[len(key)+1:]
        #print "feature", feature
        self.features[key] = feature 
 
  # loads an existing purefiles list
  def loadImagesFromPurefiles(self):
    print "loading images from purefiles."
    self.images = []
    file = open(self.purefiles, 'r')
    lines = file.readlines()
    file.close()
    num = len(lines) 
    print "info: number of images =", num
    for line in lines:
      self.images.append(strip(line))
     
  # loads images from commandline
  def loadImagesFromCommandline(self, argv):
    self.images = []
    print "loading images from commandline."
    catch = False
    for command in argv:
      #print "command:", command
      if command == "-limages":
        #print "image tag found."
        catch = True
      elif command[0] == "-":
        #print "end tag found."
        catch = False
      else:
        if catch:
          self.images.append(command)
    if self.images == []:
      print "error: no images found"
          
  def loadSelectedFeaturesFromCommandline(self, argv):
    self.selectedFeatures = []
    print "loading selected features from commandline."
    catch = False
    for command in argv:
      #print "command:", command
      if command == "-selectf":
        #print "image tag found."
        catch = True
      elif command[0] == "-":
        #print "end tag found."
        catch = False
      else:
        if catch:
          #print "added", command
          self.selectedFeatures.append(command)
    if self.selectedFeatures == []:
      print "error: no features found" 
             
        
  # checking      
        
  # check if image file has a correct extension
  def checkExtension(self, filename):
    extension = split(filename, ".")[-1]  
    if lower(extension) in self.extensionlist:
      return True
    return False
    
  def pureFileExists(self):
    try:
      f = open(self.purefiles, "r")  
    except IOError: 
      return False
    else:
      return True
  
  # writing    
    
  # writes the purefilesfile
  def writePurefiles(self):
    print "writing purefiles file."
    if self.images == []:
      print "error: cannot write a purefiles file, because image list is empty."
    else:
      f = open(self.purefiles, "w")
      for image in self.images:
        f.write(image)
        f.write("\n")
      f.close()
    
  # creates the database file  
  def writeDBFile(self):
    # possible feature: dbfiles with feature directories?
    # possible feature: non recursive traversion
    print "writing dbfile."
    f = open(self.dbfile, "w")
    f.write("FIRE_filelist\n")
    f.write("classes no\n")
    f.write("descriptions no\n")
    f.write("featuredirectories no\n")
    f.write("path this\n")
    for file in self.images:
      f.write("file ")
      f.write(file)
      f.write("\n")
    f.close()
    
    
  # selection
    
  # selects features
  def selectFeatures(self):
    if self.selectedFeatures == []:
      print "error: there are no features to select"
      print "hint: maybe wrong spelling?"
    temp = {}
    for key in self.features.keys():
      #print "key", key
      if key in self.selectedFeatures:
        temp[key] = self.features[key]
    self.features = temp
    
  
  # extraction    
        
  # extracts all features for all images
  def extractFeatures(self):
    for key in self.features.keys():
      self.extractFeature(key)
    
  # extracts one feature for all images
  def extractFeature(self, feature):
    if self.images != []:
      commandline = self.firedir + self.features[feature] + " --images"
      for image in self.images:
        commandline = commandline + " " + image
    else:  
      if not self.pureFileExists():
        print "error: purefiles file not existing"
        print "hint: create with -np a new file or use -limage, -ldb or -ld to load images"
        return False
      else:
        #self.printInfo()
        #print "self.features[feature]", self.features[feature]
        #self.printFeatures()
        commandline = self.firedir + self.features[feature] + " --filelist " + self.purefiles
    print "commandline:", commandline
    if self.queue == "day":
      queueline = "qsubmit -n " + self.name + "." + feature + " -m 0.4 -t 8:00:00 -bash  "
      print "job", queueline + commandline
      os.system(queueline + commandline)
    elif self.queue == "week": 
      queueline = "qsubmit -n " + self.name + "." + feature + " -m 0.4 -t 24:00:00 -bash  "
      print "job", queueline + commandline
      os.system(queueline + commandline)
    else:
      os.system(commandline)
      
      
      
if __name__ == "__main__":

  if "-h" in sys.argv:
    
    usage(sys.argv)
  
  else:
  
    e = Extractor()

    # adjusting file locations
    if "-p" in sys.argv:
      e.setPurefiles(getStringAfter("-p", sys.argv, None))
    if "-db" in sys.argv:
      e.setDBFile(getStringAfter("-db", sys.argv, None))
    if "-d" in sys.argv:
      e.setDatadir(getStringAfter("-d", sys.argv, None))
    if "-c" in sys.argv:
      e.setConfigfile(getStringAfter("-c", sys.argv, None))
    if "-f" in sys.argv:
      e.setFiredir(getStringAfter("-f", sys.argv, None))
    if "-q" in sys.argv:
      e.setQueue(getStringAfter("-q", sys.argv, None))
    if "-n" in sys.argv:
      e.setName(getStringAfter("-n", sys.argv, None))
    ## maybe later
    # if "-el" in sys.argv
    # e.setExtensionlist(...)  
  
    if "-v" in sys.argv:  
      e.printInfo()
    if "-np" in sys.argv: 
      e.loadImagesFromDir()
      e.writePurefiles()
    elif "-npfromdb" in sys.argv:
      e.loadImagesFromDB()
      e.writePurefiles()
    elif "-ndb" in sys.argv:
      e.loadImagesFromDir()
      e.writeDBFile()
    elif "-ndbfromp" in sys.argv:
      e.loadImagesFromPurefiles()
      e.writeDBFile()
    elif "-showf" in sys.argv:
      e.loadFeatures()
      e.printFeatures()
    elif "-e" in sys.argv:
      # loading the images (if no purefiles file is used)
      if "-limages" in sys.argv:
        e.loadImagesFromCommandline(sys.argv)
        if "-v" in sys.argv:
          e.printImages()
      elif "-ldb" in sys.argv:
        e.loadImagesFromDBFile()
      elif "-ld" in sys.argv:
        e.loadImagesFromDatadir()
      else:
        pass
      e.loadFeatures()
      if "-selectf" in sys.argv:
        e.loadSelectedFeaturesFromCommandline(sys.argv)
        e.selectFeatures()
      if "-v" in sys.argv:
        print "selected features:"
        e.printFeatures()
      e.extractFeatures()
    else:
      print "wrong usage. use -h for help."
  

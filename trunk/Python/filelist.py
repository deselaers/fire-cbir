import re,sys,string,os
class FileList:
    def __init__(self):
        self.files=[]
        self.classes=False
        self.featuredirectories=False
        self.descriptions=False
        self.path=""
        self.rawpath=""
        self.suffices=[]
        self.filename2idx={}

    def __str__(self):
        result=""
        result+="Suffices: "+str(self.suffices)
        result+="Files: "+str(self.files)
        return result
        
    def __getitem__(self,i):
        return self.files[i]

    def __getslice__(self,i):
        return self.files[i]

    def index(self,filename):
        if self.filename2idx.has_key(filename):
            return self.filename2idx[filename]
        else:
            print >>sys.stderr,filename,"not in database, returning -1"
            #print >>sys.stderr,self.filename2idx
            #print >>sys.stderr,self.files
            return -1
        #result=-1
        #for i in range(len(self.files)):
        #    if self.files[i][0]==filename:
        #        result=i
        #if result==-1:
        #    print "filelist.py: file",filename,"not in database"
        #return result
    def load(self, filename):
        fl=open(filename,"r")
        lines=fl.readlines()

        for l in lines:
            l=re.sub("\n","",l)
            tokens=string.split(l)             
            if tokens[0] == "file":
                f=tokens[1]
                if self.classes==True:
                    cls=int(tokens[2])
                    if self.descriptions:
                        desc=tokens[3:]
                        self.files+=[[f,cls,desc]]
                    else:
                        self.files+=[[f,cls]]
                elif self.descriptions:
                    desc=tokens[2:]
                    self.files+=[[f,0,desc]]
                else:
                    self.files+=[[f,]]
            elif tokens[0] == "classes":
                if tokens[1]=="yes":
                    self.classes=True
            elif tokens[0] == "descriptions" or tokens[0] == "description":
                if tokens[1]=="yes":
                    self.descriptions=True
            elif tokens[0] == "featuredirectories":
                if tokens[1]=="yes":
                    self.featuredirectories=True
            elif tokens[0] == "path":
                self.rawpath=tokens[1]
                self.path=self.rawpath
                
                if self.path.startswith("this"):
                    self.path=self.path.replace("this+","")
                    offs=os.path.dirname(filename)
                    cwd=os.path.realpath(os.path.curdir)
                    self.path=offs+"/"+self.path+"/"

                
            elif tokens[0] == "suffix":
                self.suffices+=[tokens[1]]
        self.initMap()

    def add(self, filename, class_=0, description=""):
        if description == "":
            assert(self.descriptions == False)
            self.files.append([filename, class_, description])
        else:
            assert(self.descriptions == True)
            self.files.append([filename, class_])

    def save(self, filename):
        if filename=="-":
            f=sys.stdout
        else:
            f=open(filename,"w")

        print >>f,"FIRE_filelist"
        if self.classes:
            print >>f,"classes yes"
        if self.descriptions:
            print >>f,"descriptions yes"
        if self.featuredirectories:
            print >>f,"featuredirectories yes"
        print >>f,"path",self.rawpath
        for s in self.suffices:
            print >>f,"suffix",s

        for fi in self.files:
            print >>f, "file %s" % " ".join(map(lambda x: str(x), fi))

    def initMap(self):
        self.filename2idx={}
        for i in range(len(self.files)):
            self.filename2idx[self.files[i][0]]=i

    def basenameify(self):
        for i in range(len(self.files)):
            if len(self.files[i])==3:
                self.files[i]=(os.path.basename(self.files[i][0]),self.files[i][1], self.files[i][2])
            if len(self.files[i])==2:
                self.files[i]=(os.path.basename(self.files[i][0]),self.files[i][1])
            if len(self.files[i])==1:
                self.files[i]=(os.path.basename(self.files[i][0]),)

    def nameWithSuffix(self, item, suffix):
        assert item < len(self.files)
        assert suffix < len(self.suffices)
        return self.files[item][0] + "." + self.suffices[suffix]

    # generator for images with full path
    def fullpaths(self):
        for image in self.files:
            yield os.path.join(self.rawpath, image[0])
    
    def get_full_path_for_item(self, item, suffixidx=-1):
        name = os.path.join(self.rawpath, self.files[item][0])
        if suffixidx > -1 and len(self.suffices) > 0:
            name += "." + self.suffices[suffixidx]
        return name
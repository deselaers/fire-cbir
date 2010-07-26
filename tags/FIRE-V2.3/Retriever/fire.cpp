/*
  This file is part of the FIRE -- Flexible Image Retrieval System
  
  FIRE is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your
  option) any later version. 
  
  FIRE is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.
  
  You should have received a copy of the GNU General Public License
  along with FIRE; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

/** This file contains the main (or init) class for the FIRE image
    retrieval engine. This file does not contain any image retrieval
    functionality but mere initialization of the necessary code.
    
    The image retrieval engine can be found mainly in the files
    
    - server.hpp, server.cpp -- processing of options, network
    communication initialization of retriever and database.
  
    - retriever.hpp, retriever.cpp -- logic for image retrieval. Here
    the requests that are given to the server are processed, image
    comparision operations are processed and here the settings for how
    images are processed are stored.
 
    - database.hpp, database.cpp -- here the features representing the
    images are managed. All the features are contained in classes and
    here the management is done.
*/

#include <iostream>
#include <csignal>
#include "diag.hpp"
#include "getpot.hpp"
#include "server.hpp"

using namespace std;

/// short method to print help on command line options of the
/// retrieval framework
void USAGE() {
  cout << "USAGE: " << endl
       << "  -h,--help                   give help" << endl
       << "  -c,--config <configfile>    read the specified config file" << endl
       << "  -s,--server <port>          start server on given port" << endl
       << "  -f,--filelist <filelist>    load filelist" << endl
       << "  -d,--dist <nr> <distname>   set dist nr to distname" << endl
       << "  -D,--defaultdists           get default distance for each feature" << endl
       << "  -w,--weight <nr> <weight>   set weight nr to weight" << endl
       << "  -r,--results <nr>           set number of results" << endl
       << "  -e,--expansion <nr>         set the number of images used for query expansion" << endl
       << "                              note: 0(default)=disabled, 1=use the first result" << endl
       << "                              (often the query itself), 2=use the first two" << endl
       << "                              results, 3=the first three results, etc." << endl
       << "  -R,--relevanceFile <file>   set the file, into which relevances are saved" << endl
       << "  -l,--logfile <file>         set the logfile" << endl
       << "  -S,--startAfterInit <cmd>   run the specified command after the database " << endl
       << "                              is completely loaded (in the background) " << endl
       << "                              in cmd, %P is substituted with the port number" << endl
       << "  -P,--proxy <proxy> <port>   announce your host and port to proxy after init" << endl
       << "  -C,--scoring <scorer>       to set the algorithm that is used to combine distances" << endl
       << "                              vectors into scores. default: linear, maxent, maxent2nd, maxent1st2nd" << endl
       << "  -p,--password <pwd>         set the password which is necessary for changing settings."<< endl
       << "                              if not set: everybody is allowed to change everything."<< endl
       << "  -I,--interactor <interactorspecification> to specify an interactor." << endl
       << "                              example: -I 'ACTION=min(0,1):ACTION=max(2,3,4):ACTION=mul(3,5)'" << endl
       << "                              actions specified are min, max, and mul, each of them takes" << endl
       << "                              arbitrary many position parameters" << endl
       << "  -q,--queryCombiner <combiner> to specify the query combination technique." << endl
       << "                              relevance feedback techniques: adding, nn, to be continued..." <<endl
       << "  --reRanker <reranker>   to specify the reranking technique."<<endl
       << "                              reranking: none, cluster,greedy,visavis,...."<<endl
       << "  -B,--batch <file>           take the commands from the file instead of listening to a socket" << endl
       << "                              by default the complete ranking is returned for each query" << endl
       << "  -F,--filter <filtersequence> set filter for filtered retrieval; the first distance will be" << endl
       << "                              be computed on the whole database and the amount of" << endl
       << "                              amountToBeUsed best resulting images is used for computing" << endl
       << "                              the second distance and so on. these 'best' result from a combined" << endl
       << "                              scoring of all previous computations; syntax:" << endl
       << "                              distanceID:amountToBeUsed-distanceID:amountToBeUsed" << endl
       << "                              example: 0:1000-3:500-2:20" << endl
       << "  -u,--dontload <featureindexsequence> only useable if also the -F/--filter option is used." << endl
       << "                              the first feature specified in the filtersequence following" << endl
       << "                              -F/--filter will always be loaded at startup. all features specified" << endl
       << "                              in the featureindexsequence won't be loaded into memory except they" << endl
       << "                              are present in the filtersequence. if so these features are only loaded" << endl
       << "                              when needed for the computation." << endl
       << "                              syntax: idx:idx:idx...  the idx is zero based" << endl
       << "                              example: 0:3:2 will force FIRE not to load features 0 2 3 at startup" << endl
       << "  -U,--defdontload            forces FIRE only to load the features specified in the filter sequence" << endl
       << "                              at startup. during runtime there will be no more feature loading." << endl
       << "                              only usable when also -F/--filter is used otherwise ignored." << endl
       << " --cache <filename>           use sqlite cache from that file" << endl
       << "  -t,--type2bin <file>        override the type2bin-path set in the filelist" << endl
       << endl;
  exit(20);
}

void  INThandler(int)
{
  ERR << "SIG_SEGV. Stacktrace" << endl;
  stackTrace();
  exit(20);
}

int main(int argc, char **argv)
{
  GetPot cl(argc,argv);

  signal(SIGSEGV,INThandler);

  Server server;

  vector<string> ufos=cl.unidentified_options(46,
                      "-h", "--help", "-c", "--config", "-s",//5
                      "--server", "-f", "--filelist", "-d", "--dist", //10
                      "-D", "--defaultdists", "-w", "--weight", "-r",//15
                      "--results", "-e", "--expansion", "-R", "--relevanceFile",//20
                      "-l", "--logfile", "-S", "--startAfterInit", "-C",//25
                      "--scoring", "-p", "--password","-I","--interactor",//30
                      "-P","--proxy","-B","--batch","-F", //35
                      "--filter","-u","--dontload","-U","--defdontload",//40
                                              "-t", "--type2bin","--cache","-q","--queryCombiner", //45
                                              "--reRanker"); //46

  if(ufos.size()!=0)
  {
    for(vector<string>::const_iterator i=ufos.begin();i!=ufos.end();++i)
    {
      cout << "Unknown option detected: " << *i << endl;
    }
    USAGE();
  }

  if (cl.search(2,"-h", "--help"))
  {
    USAGE();
  }
  else if(cl.search(2,"-c","--config"))
  {
    DBG(10) << "Reading configfile" << endl;
    GetPot cf(cl.follow("config.pot",2,"-c","--config"));
    cl.absorb(cf);
  }
  server.parseConfig(cl);
  server.initialize();

  if(cl.search(2,"-B","--batch"))
  {
    server.batch();
  }
  else
  {
    server.start();
  }

  DBG(1) << "Commandline was:" ; for(int i=0;i<argc;++i) { BLINK(1)<< " " << argv[i]; }
  DBG(1) << "Resulting in this config:" << endl;
  cl.print();
}

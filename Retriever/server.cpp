/*
  This file is part of the FIRE -- Flexible Image Retrieval System
 
  FIRE is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  FIRE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with FIRE; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <pthread.h>
#include <cstring>
#include <errno.h>
#include "server.hpp"

#include "net.hpp"
#include "diag.hpp"
#include "logfile.hpp"
#include "runprogram.hpp"
#include "ScopeTimer.h"

using namespace std;

typedef struct DataForThread_ {
  Server* server;
  Socket* client;
} DataForThread;


void Server::initialize() {
  retriever_.initialize();
}

void* threadProcess(void *data)
{
#ifdef __USE_PTHREADS_FOR_SERVER__
  pthread_detach(pthread_self());
#endif

  
  DataForThread& dfs=*((DataForThread*)data);
  
  Socket &client=*(dfs.client);
  Server &serv=(*dfs.server);
  bool authorized=(serv.password()=="");
  ServerStatus processCommandReturned=ServerStatusGOOD;

  DBG(2) << "Client connected" << endl;
  serv.log().log("client connected");
  bool notBye=true;
  while(notBye)
  {
    string toClient;
    // get one line
    string commandLine=client.getline();
    serv.log().log("RECV: "+commandLine);

    // process the line
    processCommandReturned=serv.processCommand(commandLine,toClient,authorized);

    // everything fine?
    switch(processCommandReturned)
    {
    case ServerStatusGOOD:
      break;
    case ServerStatusBYE:
      notBye=false;
      authorized=false;
      break;
    case ServerStatusQUIT:
      serv.notQuit()=false;
      notBye=false;
      authorized=false;
      break;
    case ServerStatusBAD:
      break;
    }
    serv.log().log("SEND: "+toClient);
    toClient+="\r\n";
    client << toClient;
  } // while (notBye)
  // the client has said goodbye or quit
  client.close();
  serv.log().log("client connection closed");
  
  //delete dfs.client;
  
  delete data;
  
#ifdef __USE_PTHREADS_FOR_SERVER__
  pthread_exit(0);
#endif
  return NULL;
}

static const CommandType CMD_INFO=10000;
static const CommandType CMD_RETRIEVE=10001;
static const CommandType CMD_SETDIST=10002;
static const CommandType CMD_SETWEIGHT=10003;
static const CommandType CMD_FILELIST=10004;
static const CommandType CMD_SETRESULTS=10005;
static const CommandType CMD_RANDOM=10006;
static const CommandType CMD_LISTFILES=10007;
static const CommandType CMD_SAVEDISTANCES=10009;
static const CommandType CMD_HELP=10010;
static const CommandType CMD_QUIT=10011;
static const CommandType CMD_BYE=10012;
static const CommandType CMD_CLASS=10013;
static const CommandType CMD_EXPAND=10014;
static const CommandType CMD_SAVERELEVANCES=10015;
static const CommandType CMD_RETRIEVEANDSAVERANKS=10016;
static const CommandType CMD_SETEXTENSION=10017;
static const CommandType CMD_METARETRIEVE=10018;
static const CommandType CMD_SETSCORING=10019;
static const CommandType CMD_PASSWORD=10020;
static const CommandType CMD_TEXTRETRIEVE=10021;
static const CommandType CMD_METAEXPAND=10022;
static const CommandType CMD_TEXTEXPAND=10023;
static const CommandType CMD_METAFEATUREINFO=10024;
static const CommandType CMD_INTERACTOR=10025;
static const CommandType CMD_IMAGE=10026;
static const CommandType CMD_FEATURE=10027;
static const CommandType CMD_SETFILTER=10028;
static const CommandType CMD_NEWFILE=10029;

Server::Server()  :  port_(12960), retriever_(),batchfile_(""), notQuit_(true)
{
  map_["info"]=CMD_INFO;
  map_["retrieve"]=CMD_RETRIEVE;
  map_["setdist"]=CMD_SETDIST;
  map_["filelist"]=CMD_FILELIST;
  map_["setweight"]=CMD_SETWEIGHT;
  map_["setresults"]=CMD_SETRESULTS;
  map_["random"]=CMD_RANDOM;
  map_["listfiles"]=CMD_LISTFILES;
  map_["savedistances"]=CMD_SAVEDISTANCES;
  map_["help"]=CMD_HELP;
  map_["quit"]=CMD_QUIT;
  map_["bye"]=CMD_BYE;
  map_["class"]=CMD_CLASS;
  map_["expand"]=CMD_EXPAND;
  map_["metaexpand"]=CMD_METAEXPAND;
  map_["metafeatureinfo"]=CMD_METAFEATUREINFO;
  map_["textexpand"]=CMD_TEXTEXPAND;
  map_["saverelevances"]=CMD_SAVERELEVANCES;
  map_["retrieveandsaveranks"]=CMD_RETRIEVEANDSAVERANKS;
  map_["setextensions"]=CMD_SETEXTENSION;
  map_["metaretrieve"]=CMD_METARETRIEVE;
  map_["textretrieve"]=CMD_TEXTRETRIEVE;
  map_["setscoring"]=CMD_SETSCORING;
  map_["password"]=CMD_PASSWORD;
  map_["interactor"]=CMD_INTERACTOR;
  map_["image"]=CMD_IMAGE;
  map_["feature"]=CMD_FEATURE;
  map_["setfilter"]=CMD_SETFILTER;
  map_["newfile"]=CMD_NEWFILE;

}

bool Server::parseFilter(const char* str)
{
  bool filterOK = false;
  string fst;
  string snd;
  vector< pair<uint,uint> > filter;
  while (*str != '\0')
  {
    fst.clear();
    snd.clear();
    while(*str != ':' && *str != '\0')
    {
      if (isalpha(*str))
      {
        filterOK = false;
        return(filterOK);
      }
      fst.append(str,1);
      str++;
    }
    if (*str == '\0')
    {
      filterOK = false;
      return(filterOK);
    }
    str++;
    while(*str != '-' && *str != '\0')
    {
      if (isalpha(*str) || *str==':')
      {
        filterOK = false;
        return(filterOK);
      }
      snd.append(str,1);
      str++;
    }
    if (*str == '-')
    {
      str++;
    }
    /// converting numbers in string format to uint, since atoi is only capable
    /// of c-style strings fst and snd must be converted.
    filter.push_back(make_pair(atoi(fst.c_str()),atoi(snd.c_str())));
  }
  retriever_.setFilter(filter);
  filterOK = true;

  DBG(105) << "Filter parsed: " ; for(uint i=0;i<filter.size();++i) { BLINK(105) << "   " << filter[i].first << ":" << filter[i].second; }
  return(filterOK);
}

void Server::parseConfig(GetPot &config)
{
  // this variable is needed for setting up the partial loading part when the features
  // are loaded for the very first time at fire start up
  string partialLoadingString="empty";

  if(config.search(2,"-F","--filter"))
  {
    string unparsed = config.follow("empty",2,"-F","--filter");
    if (parseFilter(unparsed.c_str()))
    {
      retriever_.setFilterApply(true);
      DBG(10) << "Using filter " << unparsed << " for retrieval" << endl;
    }
    else
    {
      retriever_.setFilterApply(false);
      retriever_.clearFilter();
      DBG(10) << "Error in filter String " << unparsed << endl;
      DBG(10) << "Aborting filtering and using normal retrieval settings" << endl;
    }
  }

  if(config.search(2,"-q","--queryCombiner")) {
    retriever_.setQueryCombiner(config.follow("adding:POS=1.0:NEG=0.833",2,"-q","--queryCombiner"));
  }

  if(config.search("--reRanker")) {
    retriever_.setReranking(config.follow("cluster:CONS=100:RR=20:CLUSTERS=5","--reRanker"));
  } 
  if(config.search(2,"-U","--defdontload"))
  {
    partialLoadingString="default";
    retriever_.setPartialLoadingApply(true);
  }

  if(config.search(2,"-u","--dontload"))
  {
    string unparsed=config.follow("empty",2,"-u","--dontload");
    if(unparsed!="empty")
    {
      // the parsing is not possible here, because the amount of different features
      // is not known yet
      // but a sanity check is possible
      const char* str = unparsed.c_str();
      bool errformat = false;
      // note it is not checked if there are two or more : present in a row
      while(*str != '\0' && !errformat)
      {
        if(!isdigit(*str) && *str!=':')
        {
          errformat=true;
          DBG(10) << "Error in format string " << unparsed << " following -u/--dontload option" << endl;
          DBG(10) << "format string must start with a digit; aborting partial feature loading" << endl;
        }
        str++;
      }
      if(!errformat && retriever_.setPartialLoadingApply(true))
      {
        partialLoadingString=unparsed;
      }
      else
      {
        ERR << "-u/--dontload was used without -F/--filter option or there are illegal characters in the format string" << endl;
        ERR << "ignoring -u/--dontload option, loading all features" << endl;
      }
    }
    else
    {
      retriever_.setPartialLoadingApply(false);
      ERR << "Error after -u/--dontload option: No featureindexsequence given" << endl;
      //DBG(10) << "Aborting partial loading for filtered retrieval and loading all features" << endl;
    }
  }

  if(config.search(2,"-f","--filelist"))
  {
    string filelistname=config.follow("list.txt",2,"-f","--filelist");
    string result=retriever_.filelist(filelistname,partialLoadingString);
    if(result=="filelist FAILURE")
    {
      ERR << "Error in filelist. Please correct it and start again." << endl;
      exit(1);
    }
    DBG(10) << result << endl;
  }
  if(config.search(2,"-R","--relevancefile"))
  {
    relevanceFile_=config.follow("",2,"-R","--relevancefile");
  }

  retriever_.setScoring("linear");

  
  log_=LogFile(config.follow("",2,"-l","--logfile"));

  if(config.search(2,"-D","--defaultdists"))
  {
    for(uint i=0;i<retriever_.numberOfSuffices();++i)
    {
      retriever_.dist(i,distanceMaker_.getDefaultDistance(retriever_.featureType(i)));
      retriever_.weight(i,1.0);
      DBG(10) << "default dist[" << i << "]=" << retriever_.dist(i)->name() << endl;
    }
  }

  config.init_multiple_occurrence(); // can have multiple distances
  while(config.search(2,"-d","--dist"))
  {
    uint idx=config.follow(0,2,"-d","--dist");
    string distname=config.next("basedist");  // default distance measure (dummy!)
    retriever_.weight(idx,1.0); // set default weight for configured distances
    retriever_.dist(idx,distanceMaker_.makeDistance(distname));
    DBG(10) << "dist[" << idx << "]=" << retriever_.dist(idx)->name() << endl;
  }

  config.init_multiple_occurrence(); // can have multiple weights
  while(config.search(2,"-w","--weight"))
  {
    uint idx=config.follow(0,2,"-w","--weight");
    double w=config.next(1.0);
    retriever_.weight(idx,w);
    DBG(10) << "weight[" << idx << "]=" << w << endl;
  }

  retriever_.setScoring(config.follow("linear",2,"-C","--scoring"));
  
  if(config.search(2,"-I","--interactor"))
  {
    retriever_.setInteractor(config.follow("",2,"-I","--interactor"));
  }

  config.enable_loop(); // enable returning to the beginning of command line (was disabled for multiplicity before)
  if(config.search(2,"-S","--startAfterInit"))
  {
    startAfterInit_=config.follow("echo \"No program defined\"",2,"-S","--startAfterInit");
  }

  if(config.search(2,"-P","--proxy"))
  {
    proxyhost_=config.follow("localhost",2,"-P","--proxy");
    proxyport_=config.next(12963);
    DBG(10) << "Will notify proxy at " << proxyhost_ << ":"<<proxyport_ << " about my existence." << endl;
  }

  if(config.search(2,"-r","--results"))
  {
    retriever_.results(config.follow(9,2,"-r","--results"));
    DBG(10) << "results=" << retriever_.results() << endl;
  }

  if(config.search(2,"-l","--log"))
  {
    logfile_=config.follow("log.txt",2,"-l","--log");
    DBG(10) << "logfile="<< logfile_ << endl;
  }

  if(config.search(2,"-s","--server"))
  {
    port_=config.follow(12961,2,"-s","--server");
    DBG(10) << "port="<< port_ << endl;
  }

  if(config.search(2,"-B","--batch"))
  {
    batchfile_=config.follow("batch",2,"-B","--batch");
    DBG(10) << "batchfile="<< batchfile_ << endl;
  }
  
  if(config.search("--cache")) {
    retriever_.setCache(config.follow("cache.sqlite3.db","--cache"));
  }


  if(config.search(2,"-e","--expansion"))
  {
    retriever_.extensions(config.follow(2,2,"-e","--expansion"));
    DBG(10) << "expansions=" << config.follow(2,2,"-e","--expansion") << endl;
  }

  if(config.search(2,"-p","--password"))
  {
    password_=config.follow("",2,"-p","--password");
    DBG(10) << "password necessary for changing settings: " << password_ << endl;
  }

  if(config.search(2,"-t","--type2bin"))
  {
    t2bpath=config.follow("",2,"-t","--type2bin");
    DBG(10) << "type2bin-file has been reset to: " << t2bpath << endl;
  }
}


void Server::tokenize(const string &cmdline, vector<string> &tokens)
{
  istringstream iss(cmdline);
  while(!iss.eof())
  {
    string t;
    iss >> t;
    tokens.push_back(t);
  }
}

CommandType Server::mapCommand(const string &cmd)
{
  return map_[cmd];
}



ServerStatus Server::processCommand(const ::std::string& commandline, ::std::string& toclient, bool & authorized)
{

  ScopeTimer st1("Server::processCommand");
  toclient=string("");
  ServerStatus result=ServerStatusGOOD;

  ostringstream os;

  vector<string> tokens(0);
  tokenize(commandline,tokens);

  CommandType command=mapCommand(tokens[0]);

  switch(command) {
  case CMD_PASSWORD: {
    if(tokens.size()>1) {
      if(tokens[1]==password_) {
        authorized=true;
      }
    }
    if(authorized) {
      os << "ok";
    } else {
      os << "permission denied";
    }
    break;
  }
  case CMD_QUIT: {
    if(authorized) { // only authorized users may quit the server
      result=ServerStatusQUIT;
    } else {
      result=ServerStatusBYE;
    }
    break;
  }
  case CMD_BYE: { 
    result=ServerStatusBYE;
    break;
  }
  case CMD_INFO: {
    os << retriever_.info();
    break;
  }
  case CMD_HELP: {
    for(map<const string, CommandType>::const_iterator i=map_.begin();i!=map_.end();++i) {
      os << i->first << " ";
    }
    break; 
  }
  case CMD_METARETRIEVE: {
    string querystring(commandline,commandline.find(" ")+2,commandline.size());
    vector<ResultPair> results=retriever_.metaretrieve(querystring);
    if(results.size()==0) {
      os <<"nometainformation";
    } else {
      sort(results.rbegin(), results.rend());
      for(uint i=0;i<results.size() && i<retriever_.results() && i<retriever_.numberOfFilelistEntries();++i) {
        os << retriever_.filelist(results[i].second) << " " << results[i].first << " ";
      }
    }
    break;
  }
  case CMD_TEXTRETRIEVE: {
    string querystring(commandline,commandline.find(" ")+2,commandline.size());
    DBG(10) << "calling textretrieve(" << querystring << ")" << endl;
    vector<ResultPair> results=retriever_.textretrieve(querystring);
    if(results.size()==0) {
      os <<"notextinformation";
    } else {
      sort(results.rbegin(), results.rend());
      for(uint i=0;i<results.size() && i<retriever_.results() && i<retriever_.numberOfFilelistEntries();++i) {
        os << retriever_.filelist(results[i].second) << " " << results[i].first << " ";
      }
    }
    break;
  }
  case CMD_RETRIEVE:     // normal retrieval
  case CMD_RETRIEVEANDSAVERANKS:   //normal retrieval + save ranks file
  case CMD_EXPAND: {    // the n-th step of normal retrieval (that is show more results)
    ScopeTimer st1("Server::processCommand -> CMD_RETRIEVE");
    uint queriesStartFrom, // where do the names of the queries
      // start, there are some more
      // parameters at the start of the
      // commandline for CMD_EXPAND and
      // CMD_RETRIEVEANDSAVERANKS
      resultsStep=0;       // which step of showing results are
    // we. if command was CMD_RETRIEVE, we
    // are in step 0, in expand this is
    // given as first parameter.
    if(command==CMD_RETRIEVE) {
      queriesStartFrom=1;
    } else if(command==CMD_RETRIEVEANDSAVERANKS) {
      queriesStartFrom=3;
    } else if(command==CMD_EXPAND) {
      queriesStartFrom=2;
      istringstream iss(tokens[1]);
      iss >> resultsStep;
    } else {
      queriesStartFrom=1;ERR << "unrecognized command, setting queriesStartFrom=1" << endl;
    }
    
    vector<string> posQueriesNames;
    vector<string> negQueriesNames;
    vector<ResultPair> results;
    string tmp;
    for(uint i=queriesStartFrom;i<tokens.size();++i) {
      if(tokens[i][0]=='+') {
        tokens[i].erase(0,1);
        posQueriesNames.push_back(tokens[i]);
      }
      else if(tokens[i][0]=='-') {
        tokens[i].erase(0,1);
        negQueriesNames.push_back(tokens[i]);
      } else {
        //no sign given: assuming +
        posQueriesNames.push_back(tokens[i]);
      }
    }
    retriever_.retrieve(posQueriesNames, negQueriesNames,results);
        
    // output!
    sort(results.rbegin(), results.rend());
    for(uint i=resultsStep*retriever_.results();i<(resultsStep+1)*retriever_.results() and i<retriever_.numberOfFilelistEntries();++i) {
      os << retriever_.filelist(results[i].second) << " " << results[i].first << " ";
    }
    
    if(command==CMD_RETRIEVEANDSAVERANKS) { //now save ranks:
      istringstream iss(tokens[1]);
      uint nOfRanks;
      iss >> nOfRanks;
      string filename=tokens[2];
      
      DBG(10) << "saving " << nOfRanks << " ranks to " << filename << endl;
      
      ofstream os(filename.c_str());
      if(!os) { ERR << "Error opening logfile:" << filename << endl;}
      os << "# " << commandline << endl;
      for(uint i=0;i<nOfRanks;++i) {
        os << i << " " << retriever_.filelist(results[i].second) << " " << results[i].first << endl;
      }
      os.close();
    }
    break;
  }
  case CMD_METAEXPAND: {
    DBG(10) << "metaexpand" << endl;
    uint resultsStep;
    istringstream iss(tokens[1]);
    iss >> resultsStep;
    
    string querystring(commandline,commandline.find(" ")+2,commandline.size());
    DBG(10) << "calling textretrieve(" << querystring << ")" <<endl;
    vector<ResultPair> results=retriever_.textretrieve(querystring);
    if(results.size()==0) {
      os <<"nometainformation";
    } else {
      sort(results.rbegin(), results.rend());
      for(uint i=resultsStep*retriever_.results(); i<results.size() && i<(resultsStep+1)*retriever_.results() && i<retriever_.numberOfFilelistEntries();++i) {
        os << retriever_.filelist(results[i].second) << " " << results[i].first << " ";
      }
    }
    break;
  }
  case CMD_TEXTEXPAND: {
    DBG(10) << "textexpand" << endl;
    uint resultsStep;
    istringstream iss(tokens[1]);
    iss >> resultsStep;
    
    string qstring(commandline,commandline.find(" ")+2,commandline.size());
    DBG(10) << "calling textretrieve("<<qstring<<")"<<endl;
    vector<ResultPair> results=retriever_.textretrieve(qstring);
    if(results.size()==0) {
      os <<"notextinformation";
    } else {
      sort(results.rbegin(), results.rend());
      for(uint i=resultsStep*retriever_.results(); i<results.size() && i<(resultsStep+1)*retriever_.results() && i<retriever_.numberOfFilelistEntries();++i) {
        os << retriever_.filelist(results[i].second) << " " << results[i].first << " ";
      }
    }
    break;
  }

  case CMD_SAVEDISTANCES: { // save distance vectors for ONE query iamge to a file
    string imgname=tokens[1];  // query image
    string filename=tokens[2]; // distance file
    retriever_.saveDistances(imgname, filename);
    os << filename;
    break;
  }
  case CMD_SAVERELEVANCES: {
    if(relevanceFile_ != "") {
      ofstream os(relevanceFile_.c_str(),::std::ios::out|::std::ios::app);
      os << commandline << endl;
      os.close();
    } else {
      ERR << "Requested to save relevances, but no relevance file given." << endl;
    }
  }
  case CMD_SETSCORING: {
    if(tokens.size()==2 && authorized) { // only authorized users may change this
      retriever_.setScoring(tokens[1]);
      os << "scoring = " << retriever_.scorer()->type();
    } else {
      os << retriever_.availableScorings();
    }
    break;
  }
  case CMD_SETEXTENSION: { // set value for query EXPANSION (called EXTENSION due to legacy code)
    if(tokens.size()==2 && authorized) {
      istringstream iss(tokens[1]);
      uint i;
      iss >> i;
      retriever_.extensions(i);
      os << "extensions = " << i ;
    } else {
      os << "wrong syntax or not authorized";
    }
    break;
  }
  case CMD_SETDIST: {
    if(tokens.size()==3 && authorized) { // only authorized users may change this
      istringstream iss(tokens[1]);
      uint i;
      iss >> i;
      os << retriever_.dist(i,distanceMaker_.makeDistance(tokens[2]));
    } else {
      ostringstream oss("");
      os << distanceMaker_.availableDistances();
    }
    break;
  }
  case CMD_SETWEIGHT: {
    if(tokens.size()==3 && authorized) { // only authorized users may change this
      istringstream iss(tokens[1]);
      uint i;
      iss >> i;
      
      istringstream iss2(tokens[2]);
      double w=0.0;
      iss2 >> w;
      os << retriever_.weight(i,w);
    } else {
      os << "Invalid syntax: setweight idx weight";
    }
    break;
  }
  case CMD_FILELIST: { 
    if(tokens.size()==2 && authorized) {
        os << retriever_.filelist(tokens[1]);
    } else {
      os << "Invalid syntax: filelist filename";
    }
    break;
  }
  case CMD_SETRESULTS: {
    if(tokens.size()==2 && authorized) {
      istringstream iss(tokens[1]);
      uint i;
      iss >> i;
      os << retriever_.results(i);
    } else {
      os << "Invalid syntax: setresults nr";
    }
    break;
  }
  case CMD_RANDOM: {
    if(tokens.size() == 2) {
      uint i;
      istringstream iss(tokens[1]);
      iss >> i;
      os << retriever_.random(i);
    } else if(tokens.size()==1) {
      os << retriever_.random(retriever_.results());
    } else {
      os << "Invalid syntax: random [nr]";
    }
    break;
  }
  case CMD_LISTFILES: {
    os << retriever_.filelistEntries();
    break;
  }
  case CMD_CLASS: {
    if(tokens.size()==1) {
      if(retriever_.haveClasses()) {
        os << "yes";
      } else {
        os << "no";
      } 
    } else {
      if(retriever_.haveClasses()) {
        os << retriever_.clas(tokens[1]);
      } else {
        os << "-1";
      }
    }
    break;
  }
  case CMD_METAFEATUREINFO: {
    // Get and output metafeature info
    // This looks scary, but the mfi structure is just a pair of
    // vectors of strings. The first vector contains the keys,
    // the second vector contains example values for the keys.
    
    ::std::string res;
    ::std::pair< ::std::vector< ::std::string >,::std::vector< ::std::string > > mfi;
    ::std::vector< ::std::string > keys, vals;
    ::std::vector< ::std::string >::iterator key_it, val_it;
    mfi = retriever_.getMetaFeatureInfo();
    keys = mfi.first;
    vals = mfi.second;
    
    for(key_it=keys.begin(),val_it=vals.begin(); key_it!=keys.end(), val_it!=vals.end(); ++key_it,++val_it) {
      os << *key_it << ":" << *val_it << " ";
    }
    break;
  }
  case CMD_INTERACTOR: {
    if(tokens.size()!=2) {
      os <<"Unknown Syntax! USAGE: interactor <interactorspecification>" ;
    } else {
      retriever_.setInteractor(tokens[1]);
      os << "interactor set";
    }
    break;
  }
  case CMD_IMAGE: {
    if(tokens.size()!=2) {
      os << "Unknown syntax! USAGE: image <imagename>";
    } else {
      retriever_.printImageInformation(tokens[1],os);
    }
    break;
  }
  case CMD_FEATURE: {
    if(tokens.size()!=3) {
      os << "Unknown syntax! USAGE: feature <imagename> <featureno>";
    } else {
      retriever_.printFeatureInformation(tokens[1],tokens[2],os);
    }
    break;
  }
  case CMD_SETFILTER: {
    if(tokens.size()!= 2) {
      os << "Invalid syntax: setfilter <filtersequence>";
    } else if (!tokens[1].empty() && parseFilter(tokens[1].c_str())) {
      retriever_.setFilterApply(true);
      os << "filter = " << tokens[1];
    } else {
      retriever_.setFilterApply(false);
      retriever_.clearFilter();
      os << "filter = empty";
    }
    break;
  }
  case CMD_NEWFILE: {
    
    DBG(50) << "newfile 1" << endl;
    if(tokens.size()!=3) {
      os << "Unknown Syntax! USAGE: newfile <mode> <absolute_path_to_image>" << endl;
      os << "Options for <mode>:" << endl << "  0 : delete file after retrieval" << endl << "  1 : keep file, but don't add to database" << endl;
      os << "  2 : keep file and add to database" << endl << "  3 : keep file, add to database and add to filelist & purefiles";
      break;
    }

    DBG(50) << "newfile 2" << endl;
     
    if(tokens[1] != "0" && tokens[1] != "1" && tokens[1] != "2" && tokens[1] != "3") {
      os << "Unknown Syntax! USAGE: newfile <mode> <absolute_path_to_image>" << endl;
      os << "Options for <mode>:" << endl << "  0 : delete file after retrieval" << endl << "  1 : keep file, but don't add to database" << endl;
      os << "  2 : keep file and add to database" << endl << "  3 : keep file, add to database and add to filelist & purefiles";
      break;
    }
    
    DBG(50) << "newfile 3" << endl;

    char *end;      // unnecessary but needed since this is a stupid function ;-)
    uint mode=strtoul(tokens[1].c_str(), &end, 10);
    if(t2bpath.empty()) {t2bpath=retriever_.getT2bPath(); } 

    ::std::string newfile=tokens[2];    // absolute path to the image
    ::std::string imagename=newfile.substr(newfile.rfind('/')+1, newfile.length());    // filename of the image
    ::std::string imagepath=retriever_.getPath()+"/";

    DBG(50) << "newfile 4" << endl;

    // checks if an image with the same name already exists and renames the new image if neccesary
    uint count=1;
    ::std::string imagename_backup=imagename;
    while(fileExists(imagepath+imagename)) {
      imagename = imagename_backup;
      ostringstream number;
      number << "_" << count;
      imagename.insert(imagename.rfind('.'), number.str());
      count++;
    }

    DBG(50) << "newfile 5" << endl;

    // copies the image into the imagepath
    if(FileCopy(newfile.c_str(), (imagepath+imagename).c_str()) != 0) {
      os << "Imagefile doesn't exist!";
      break;
    }

    DBG(50) << "newfile 6" << endl;

    // gets all used features and puts them into a vector
    vector<string> features(0);
    tokenize(retriever_.getFeatures(),features);
    vector<string> used_features(0);    // stores all actually used features for deletion

    DBG(50) << "newfile 7" << endl;

    // opens file type2bin, gets the path to the feature-extractors and executes them on the image
    bool process_successful=true;
    ifstream is; is.open(t2bpath.c_str());
    if(!is.good()) {
      ERR << "Unable to open file type2bin. Corrupted t2bpath!" << endl;
      os << "Server configuration error.";
      break;
    } else {
      ::std::string line;
      ::std::string bin_command;
      
      while(!is.eof()) {
        getline(is,line);
        if(line=="") break;
        bin_command=line.substr(line.find(' ')+1, line.length());   // extractor-executable
        line.erase(line.find(' '), line.length()+1);    // feature-suffix
        
        for(uint i=0;i<features.size();++i) {
          if(features[i]==line) {
            if(system((bin_command+" --suffix "+line+" --images \""+imagepath+imagename+"\"").c_str())!=0) {
              os << "Server configuration error." << endl;
              ERR << "Unable to process image file with" << bin_command << endl;
              process_successful=false;
              break;
            }
            used_features.push_back(line);
          }
        }
      }
    }
    is.close();
    
    DBG(50) << "newfile 8" << endl;
    if(process_successful) {
      // save image in database (mode 2,3)
      if(mode > 1) { retriever_.loadQuery(imagename); }
      
      // uses the existing command "retrieve" and shows the output
      ::std::string retrieve_output;

      processCommand("retrieve "+imagename, retrieve_output, authorized);
    DBG(50) << "newfile 9" << endl;

      //      retriever_.retrieve(posQueriesNames, negQueriesNames,results);      
      cout << retrieve_output << endl;
      os << retrieve_output;
    DBG(50) << "newfile 10" << endl;

    }
    
    // delete generated files (mode 0)
    if(mode == 0 || process_successful==false) {
      if(remove((imagepath+imagename).c_str())!=0) {
        ERR << "Couldn't delete file " << imagepath+imagename << ". Please delete manually." << endl;
      }
      for(uint i=0;i<used_features.size();++i) {
        if(remove((imagepath+imagename+"."+used_features[i]).c_str())!=0) {
          ERR << "Couldn't delete file " << imagepath+imagename+"."+used_features[i] << ". Please delete manually." << endl;
        }
      }
    }
    DBG(50) << "newfile 11" << endl;
    
    // add file to filelist & purefiles (mode 3)
    if(mode == 3 && process_successful==true) {
      ofstream app;
      if(fileExists(filelistname.c_str())) {
        app.open(filelistname.c_str(), ofstream::out | ofstream::app);
        app << "file ./" << imagename << endl;
        app.close();
      } else {
        ERR << "Couldn't open filelist. Please add the image \"" << imagename << "\" manually." << endl;
      }
      if(fileExists((imagepath+"purefiles").c_str())) {
        app.open((imagepath+"purefiles").c_str(), ofstream::out | ofstream::app);
        app << "./" << imagename << endl;
        app.close();
      } else {
        ERR << "Couldn't open purefiles. Please add the image \"" << imagename << "\" manually." << endl;
      }
    }
    DBG(50) << "newfile 12" << endl;
    break;
  }
  default:
    {
      DBG(2) << "Received unknown: " << commandline << endl;
      os << "Unknown command: " << commandline;
      break;
    }
  }
  toclient=os.str();
  return result;
}


void Server::start() {
  
  /// socket stuff
  ServerSocket server(port_);
  uint p=port_;
  while(!server.listening() && p<port_+10) {
    ++p;
    DBG(20) << "Trying next port: " << p << endl;
    server=ServerSocket(p);
  }
  port_=p;
  
  if(!server.listening()) {
    uint p;
    for(p=port_+1;p<port_+10&&!server.listening();++p) {
      DBG(20) << "Trying next port: " << p << endl;
      server=ServerSocket(p);
    }
    port_=p;
    if(!server.listening()) {
      ERR << "Unable to open network socket" << endl;
      exit(20);
    }
  }
  
  /// if proxy was configured tell him where I amn
  if(proxyhost_.size()!=0) {
    DBG(10) << "Proxy: here I am" << endl;
    announceYourselfToProxy();
  } else {
    DBG(10) << "No proxy" << endl;
  }
  
  /// if i have to start another program right after init, do it here ("-S")
  if(startAfterInit_.size()!=0) {
    uint pos=startAfterInit_.find("%P");
    if(pos < startAfterInit_.size()) {
      ostringstream oss;
      string start, end;
      start.assign(startAfterInit_,0,pos);
      end.assign(startAfterInit_,pos+2,startAfterInit_.size());
      oss << start << port_ << end;
      
      startAfterInit_=oss.str();
    }
    DBG(10) << "Starting " << startAfterInit_ << endl;
    runInBackground(startAfterInit_);
  }

  /// and now we are in process commands from network mode
  string commandLine;

  string toClient;
  DataForThread *dfs;
  while (notQuit_) {
    DBG(2) << "Waiting for connections on port " << port_ << endl;
    Socket * client=server.acceptPointer();

    dfs=new DataForThread(); dfs->server=this; dfs->client=client;

#ifdef __USE_PTHREADS_FOR_SERVER__
    pthread_t thr;
    int errcode=pthread_create(&thr, NULL , threadProcess, (void*)dfs);
    if(errcode!=0)
    {
      ERR << "Problems with starting thread!: "<< strerror(errno) << endl;
      exit(20);
    }
#else
    threadProcess( (void*)dfs);
#endif


  } // while (notQuit_)

  // that's it, were are going to terminate
  server.close();

  if(startAfterInit_.size()!=0)
  {
    waitForAllChildren();
  }
}


// tell a proxy where I am
void Server::announceYourselfToProxy() const
{
  Socket s(proxyhost_,proxyport_);
  if(s.connected())
  {
    ostringstream os;
    os << "listening " << port_ ;
    s.send(os.str());
    s.close();
  }
  else
  {
    ERR << "Cannot reach proxy. Exiting" << endl;
  }
}

// batch mode processing: read a file containing commandos and process these
void Server::batch()
{
  ifstream commandfile(batchfile_.c_str());
  ::std::string cmdline;
  ::std::string output;
  bool auth=true;

  //in batch mode it is probably wanted to get the complete ranking but this can be reduced
  retriever_.results(retriever_.numberOfFilelistEntries());

  //have status before really starting
  cout << "STATUS: " << retriever_.info() << endl;

  if( commandfile.good())
  {
    getline(commandfile,cmdline);

    while( not commandfile.eof())
    {
      cout << "RECV: " << cmdline << endl;
      processCommand(cmdline,output,auth);
      cout << "SEND: " << output << endl;
      getline(commandfile,cmdline);
    }
  }
}

int Server::FileCopy ( const char *src, const char *dst ) {
#define BUFSZ 16000
#define COPY_ERROR      -1
#define COPY_OK          0

  char            *buf;
  FILE            *fi;
  FILE            *fo;
  unsigned        amount;
  unsigned        written;
  int             result;

  buf = new char[BUFSZ];

  fi = fopen( src, "rb" );

  result = COPY_OK;
  if(fi == NULL) {
    result = COPY_ERROR;
    cerr << "Problem opening input file" << std::endl;
  } else {
    fo = fopen( dst, "wb" );
    if(fo == NULL) {
      result = COPY_ERROR;
      cerr << "Problem opening output file" << std::endl;
    }
  }
  
  if (result == COPY_OK) {
    do {
      amount = fread( buf, sizeof(char), BUFSZ, fi );
      if (amount) {
        written = fwrite( buf, sizeof(char), amount, fo );
        if (written != amount) {
          result = COPY_ERROR; // out of disk space or some other disk err?
          cerr << "Problem while writing" << std::endl;
        }
      }
    } // when amount read is < BUFSZ, copy is done
    while ((result == COPY_OK) && (amount == BUFSZ));
    fclose(fi);
    fclose(fo);
  }
  delete [] buf;
  return(result);
}

bool Server::fileExists(::std::string file) {
  ifstream is;
  is.open(file.c_str());
  if(is.good()) {
    is.close();
    return true;
  } else {
    is.close();
    return false;
  }
}

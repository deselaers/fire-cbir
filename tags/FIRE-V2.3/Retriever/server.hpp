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
#ifndef __server_hpp__
#define __server_hpp__

#include <map>
#include <string>
#include "distancemaker.hpp"
#include "getpot.hpp"
#include "retriever.hpp"
#include "logfile.hpp"
#include "net.hpp"

typedef uint CommandType;

typedef uint ServerStatus;
static const ServerStatus ServerStatusGOOD=0;
static const ServerStatus ServerStatusBYE=1;
static const ServerStatus ServerStatusQUIT=2;
static const ServerStatus ServerStatusBAD=3;
static const ServerStatus ServerStatusAUTH=4;

/**
 * Server class: Handling of network connectivity, parsing of
 * commands, initialization of retrieval engine
 */
class Server
{

private:
  /// the port on which this server listens
  uint port_;

  /// the retrieval engine
  Retriever retriever_;

  /// to set up the distances
  DistanceMaker distanceMaker_;

  /// the name of the log file if one is written
  ::std::string logfile_;

  /// the name of the relevance file if one is written, this is a file
  /// where marked relevances from the web-interface are saved, when
  /// the "save relevances" button was pressed
  ::std::string relevanceFile_;

  /// the command line of a program that is started after the system is up and running
  ::std::string startAfterInit_;

  /// the file that is read instead of opening a server socket and
  /// waiting for connections if batch mode is specified
  ::std::string batchfile_;

  /// the log file object
  LogFile log_;

  /// the map where the commands are matched with the unique identifiers
  /// this map is initialized in the constructor
  ::std::map<const ::std::string,CommandType> map_;


  /// here the password is stored which is necessary for changing
  /// settings. When this string is empty, everybody can change all
  /// settings.
  ::std::string password_;

  /// if there is a proxy running, it should be on the host specified here
  ::std::string proxyhost_;

  /// port where the proxy is running
  uint proxyport_;

  /// function for parsing the filter commandline and pushing value pairs to
  /// retriever member variable filter
  bool parseFilter(const char* str);

  /// path to the filelist
  ::std::string filelistname;

  /// path to type2bin-file
  ::std::string t2bpath;

  bool notQuit_;

public:

  /// constructor, only initialization of variables
  Server() ;

  /// read the configuration from a GetPot object which can be read
  /// from command line or from a file
  /// this may include reading filelists and files (e.g. images, features)
  void parseConfig(GetPot& config);

  /// start the server, wait for incoming connections, parse the
  /// commands and pass them to the retriever
  void start();

  /// process a specified batch file
  void batch();

  /// map a command to a command identifier
  CommandType mapCommand(const ::std::string &cmd);

  /// tokenize a line by white spaces
  void tokenize(const ::std::string &cmdline, ::std::vector< ::std::string > &tokens);

  /// process a command.
  ServerStatus processCommand(const ::std::string& commandline, ::std::string& toclient, bool & authorized);

  /// tell a proxy where the server is running
  void announceYourselfToProxy() const;

  /// File copy function prototype by M. Ratcliff with little changes by Fabian Schwahn
  int FileCopy(const char *src, const char *dst);
  /// checks whether a file exists or not
  bool fileExists(::std::string file);

  bool& notQuit() {return notQuit_;}

  const string& password() {return password_;}

  LogFile& log() {return log_;}

  /// initialize all components
  void initialize();


};

#endif


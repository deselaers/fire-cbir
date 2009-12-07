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
#ifndef __net_hpp__
#define __net_hpp__

#include <string>
#include "diag.hpp"
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream> 
#include "diag.hpp"


/** 
 * Basic network socket class.
 * @author Thomas Deselaers <thomas@deselaers.de>
 * $Date: 2008-08-15 16:35:49 +0200 (Fri, 15 Aug 2008) $
 * $Id: net.hpp 744 2008-08-15 14:35:49Z deselaers $
 * 
 * This class implements a network socket. There are several operators
 * for reading and writing form and to this socket. It should be quite
 * easy to handle and is basicly inspired by the Socket class from JAVA.
 * 
 */
class Socket {
public:

  ///default constructor.
  Socket();


  ///default destructor
  ~Socket();
    
  /**standard constructor.
   *
   * @param host connect to this host
   * @param port connect to this port
   */
  Socket(::std::string host, int port);
    
  /** copy constructor */
  Socket(const Socket& s);

  /** make a connection.
   *
   * @param host connect to this host
   * @param port connect to this port
   */
  void connect(::std::string host, int port);
    
  /// close a connection.
  void close();

  ///send the ::std::string to the other side of the socket.
  void send(::std::string text);

  int getline(char *line,int len=1024, char delim='\r');
    
  /** get a whole line of text from the socket.
   * @param delim the character which marks the linebreak
   */
  ::std::string getline(char delim='\r');
    
  ///receive a ::std::string from the other side of the socket.
  ::std::string receive();

  ///return true, if the socket is connected.
  const bool connected() const ;
    
  ///return a reference to the file descriptor of the socket.
  int & socketPointer();

  ///return a reference to the pointer to the Address element of the socket.
  struct sockaddr_in* & socketAddress();
    
  ///return the size of the address element of the socket.
  int& len();
    
  ///return a reference to the file descriptor of the socket, const version,
  const int & socketPointer() const;
    
  ///return a reference to the pointer to the Address element of the socket, const version.
  const struct sockaddr_in*  socketAddress() const;
    
  ///return the size of the address element of the socket, const version.
  const int& len() const ;

  ///send a ::std::string to the socket.
  Socket & operator<<(const ::std::string& src);

  ///receive a ::std::string from the socket.
  Socket & operator>>(::std::string &dest);
    
  ///send a double to the socket (will be passed as ::std::string).
  Socket & operator<<(const double &src);
    
  ///receive a double from the socket (a ::std::string is received and interpreted as double).
  Socket & operator>>(double &dest);
    
    
  ///send a int to the socket (will be passed as ::std::string).
  Socket & operator<<(const int &src);
    
  ///receive a int from the socket (a ::std::string is received and interpreted as int).
  Socket & operator>>(int &dest);
    
private:
  struct sockaddr_in *socketAddress_;
  int len_;
  int socketPointer_;
  bool connected_;
};

/** ServerSocket class. 
 * @author Thomas Deselaers <thomas@deselaers.de>
 * $Date: 2008-08-15 16:35:49 +0200 (Fri, 15 Aug 2008) $
 * $Id: net.hpp 744 2008-08-15 14:35:49Z deselaers $
 * 
 * This class implements a ServerSocket, which is able to listen for
 * clients, accepts connections and return Sockets for incoming clients.
 */
class ServerSocket {
public:
  ///create an unbound ServerSocket.
  ServerSocket();
    
  ///create a ServerSocket listening on port.
  ServerSocket(int port);

  ///wait for an incoming connection and return it.
  Socket accept() const;
  Socket* acceptPointer() const;

  ///make the ServerSocket listening on port. 
  void open(int port);
    
  ///close the listener.
  void close();
    
  ///return the port where the server is listening, const version.
  const int &port() const;
    
  ///return the port where the server is listening, const version.
  int &port();
    
  /** Return true, if the server socket is listening for connections. */
  const bool listening() const;
    
private:
  int listenerFilePointer_;
  struct sockaddr_in serverAddress_;
  int port_;
  bool listening_;
};


#endif

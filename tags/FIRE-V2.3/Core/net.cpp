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

#include <cerrno>
#include <errno.h>
#include "net.hpp"
using namespace std;

//////////////////////////////////////////////////////////////////////
// ServerSocket

ServerSocket::ServerSocket() {
  listening_=false;
}

ServerSocket::ServerSocket(int port) {
  this->open(port);
}
  

Socket * ServerSocket::acceptPointer() const {
  int clientPointer;
  int clientLen=0;
  
  struct sockaddr* clientAddr=NULL;
  
  Socket* result=new Socket();
  DBG(15) << "Waiting for connection" << ::std::endl;
  clientPointer=::accept(listenerFilePointer_, (struct sockaddr *)0, (socklen_t *)clientLen);
  if(clientPointer <0) {
    ERR << "Error accepting connection: " << strerror(errno) << ::std::endl;
  } else {
    DBG(15) << "Got connection:" <<clientPointer << ::std::endl;
  }
  
  result->socketAddress()=(struct sockaddr_in*)clientAddr;
  result->socketPointer()=clientPointer;
  result->len()=clientLen;
  return result;
}

Socket ServerSocket::accept() const {
  int clientPointer;
  int clientLen=0;

  struct sockaddr* clientAddr=NULL;
    
  Socket result;
  DBG(15) << "Waiting for connection" << ::std::endl;
  clientPointer=::accept(listenerFilePointer_, (struct sockaddr *)0, (socklen_t *)clientLen);
  if(clientPointer <0) {
    ERR << "Error accepting connection: " << strerror(errno) << ::std::endl;
  } else {
    DBG(15) << "Got connection:" <<clientPointer << ::std::endl;
  }

  result.socketAddress()=(struct sockaddr_in*)clientAddr;
  result.socketPointer()=clientPointer;
  result.len()=clientLen;
  return result;
}
  
const bool ServerSocket::listening() const {
  return listening_;
}

void ServerSocket::open(int port) {
  DBG(15) << "Opening port " <<port<< "." << ::std::endl;
  int serverLen;
  int res;
  port_=port;
  listenerFilePointer_=socket(AF_INET, SOCK_STREAM, 0);
    
  serverAddress_.sin_family=AF_INET;
  serverAddress_.sin_addr.s_addr=htonl(INADDR_ANY);
  serverAddress_.sin_port=htons(port_);
    
  serverLen=sizeof(serverAddress_);
  res=bind(listenerFilePointer_, (struct sockaddr *) &serverAddress_, serverLen);
  listening_=true;    
  if(res<0) { 
    ERR << "Problem binding the server socket on port "<< port << ": "<< strerror(errno) << ::std::endl;
    listening_=false;
  }
  DBG(15) << "Listening on port " << port <<" now." << ::std::endl;
  res=listen(listenerFilePointer_,5);
    
  if(res<0) { ERR << "Problem with activating the listener" << ::std::endl;}

}
    
void ServerSocket::close() {
  ::close(listenerFilePointer_);
  listening_=false;
}
    
const int &ServerSocket::port() const {
  return port_;
}

int &ServerSocket::port() {
  return port_;
}
  

//////////////////////////////////////////////////////////////////////
// Socket

Socket::Socket() {
  connected_=false;
  //socketAddress_=(struct sockaddr_in *) malloc( sizeof(struct sockaddr_in) );
}

Socket::~Socket() {
  free(socketAddress_);
}


Socket::Socket(const Socket& s) {
  socketAddress_=(struct sockaddr_in *)malloc(sizeof(struct sockaddr_in));
  memcpy(socketAddress_,s.socketAddress_,sizeof(struct sockaddr_in));
  connected_=s.connected_;
}


Socket::Socket(string host, int port) {
  connected_=false;
  socketAddress_=(struct sockaddr_in *)malloc(sizeof(struct sockaddr_in));
  connect(host.c_str(), port);
}
  
void Socket::connect(string host, int port) {
  struct hostent *hostlookup;
  connected_=true;

  hostlookup=gethostbyname(host.c_str());
  if(hostlookup==NULL) {
    ERR << "Cannot lookup host '"<< host<<"'." << ::std::endl;
    connected_=false;
  }
    
  socketPointer_=socket(AF_INET,SOCK_STREAM,0);
  if(socketPointer_<0) {
    ERR << "Problem allocation socket" << ::std::endl;
    connected_=false;
  }
  socketAddress_->sin_family=AF_INET;
  memcpy(&(socketAddress_->sin_addr), hostlookup->h_addr, hostlookup->h_length);
  socketAddress_->sin_port=htons(port);
    
  len_=sizeof(struct sockaddr_in);
  int result=::connect(socketPointer_, (struct sockaddr *) socketAddress_, len_);
  if(result == -1) {
    ERR << "Cannot connect" << ::std::endl;
    connected_=false;
  }
}
  
const bool Socket::connected() const {
  return connected_;
}
  
void Socket::close() {
  ::close(socketPointer_);
  connected_=false;
}
  
int & Socket::socketPointer() {
  return socketPointer_;
}
const int & Socket::socketPointer() const {
  return socketPointer_;
}
  
struct sockaddr_in*& Socket::socketAddress() {
  return this->socketAddress_;
}
const struct sockaddr_in* Socket::socketAddress() const {
  return this->socketAddress_;
}
  
int & Socket::len() {
  return len_;
}
const int & Socket::len() const {
  return len_;
}

Socket & Socket::operator<<(const string& src) {
  send(src);
  return *this;
}
  
Socket & Socket::operator>>(string &dest) {
  dest=receive();
  return *this;
}

Socket & Socket::operator<<(const double &src) {
  ostringstream out;
  out << src;
  this->send(out.str());
  return *this;
}
  
Socket & Socket::operator>>(double &dest) {
  istringstream in(this->receive());
  in >> dest;
  return *this;
}

Socket & Socket::operator<<(const int &src) {
  ostringstream out;
  out << src;
  this->send(out.str());
  return *this;
}
  
Socket & Socket::operator>>(int &dest) {
  istringstream in(this->receive());
  in >> dest;
  return *this;
}

int Socket::getline(char *line, int len, char delim) {
  char a[1];
  int count=0;
  int res;
  len=0;
  do {
    res=read(socketPointer_, &a,1);
    if(res<0) ERR << "error in getline" << ::std::endl;
    if(*a != delim) {
      line[count]=*a;
      ++count;
    }
  } while(*a != delim || count < len-1);
  res=read(socketPointer_, &a,1); //eat second part of the \r\n    
  line[count]=(char)0;
  return count;
}
  
string Socket::getline(char delim) {
  char a[1];
  string result;
  int res;
  do {
    res=read(socketPointer_, &a,1);
    if(res<0) { ERR << "error in getline: problem reading character" << ::std::endl; break;}
    if(res==0) { ERR << "EOF in getline" << endl; break;}
    if(*a != delim) {result += *a;}
  } while (*a != delim);
  if(res<=0) { // had error or eof 
    // force disconnect of client, hoping that fire survives this.
    result="bye";
  } else {
    res=read(socketPointer_, &a,1); //eat second part of the \r\n
    if(res<=0) ERR << "error in getline: problem discarding CR" << ::std::endl;
  }
  return result;
}
  
void Socket::send(string text) {
  int res;

  char * tosend = new char[text.length()+1];
  strcpy(tosend, text.c_str());
  res=write(socketPointer_,tosend,text.length());
  delete[] tosend;
  if(res<0) {
    ERR << "Error in send: " << strerror(errno) << ::std::endl;
  }
}
  
string Socket::receive() {
  char received[1024];
  string result;
  int res=read(socketPointer_,&received,1024);
  if(res<0) {
    ERR << "Error in receive: " << strerror(errno) <<::std::endl;
    result=string("");
  } else {
    received[res]=(char)0;
    result=string(received);
  }
  return result;
}





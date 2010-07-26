#include "dist_textfeature.hpp"
#include "textfeature.hpp"
#include "net.hpp"

using namespace std;


void TextFeatureDistance::start(const BaseFeature* queryFeature) {
  const TextFeature* query=dynamic_cast<const TextFeature*>(queryFeature);

  DBG(15) << "Doing text-based retrieval from file: '" << query->value() <<"'."<< endl;
  query_wmir(":qfile "+query->value());
  queryfile_ = query->value();
 
}


double TextFeatureDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  
  const TextFeature* db=dynamic_cast<const TextFeature*>(databaseFeature);
  const TextFeature* query=dynamic_cast<const TextFeature*>(queryFeature);

  // There are two situations in which this can be called:
  // A normal image-based retrieval and a text-based retrieval
  // In the latter case, query_wmir is called before this 
  // function from textretrieve() to fill the rsv_table. Then,
  // this function is called with an empty textfeature as query.
  
  // In the former case, however, the TextFeatureDistance::start
  // function has filled the rsv_table
  
  if(rsv_table_.size() > 0) {
    if(db && query) {
      if(rsv_table_.find(db->value()) != rsv_table_.end()) {
        double dist = max_rsv_ - rsv_table_[db->value()];
        DBG(20) << VAR(db->value()) << " " << VAR(dist)  << endl;

        /*
          DBG(10) << "found one with filename "<<db->value()
          << " and RSV " << rsv_table_[db->value()]
          << " => dist: " << dist
          << endl;
        */
        return dist;
      } else {
        DBG(15) << "no matching  entry in RSV table found for '" <<db->value()<< "'." << endl;
        return max_rsv_; // Just to have a value that is much higher than the value for the found documents
      }
    } else {
      ERR << "Features not comparable" << ::std::endl;
      if(!db) {
        ERR << "no db" << ::std::endl;
      } else {
        ERR << "no query" << ::std::endl;
      }
      return -1.0;
    }

  } else {
    return 10000.0;
  }
}

// Fill the rsv table with results from wmir with the query-string
void TextFeatureDistance::query_wmir(const string& _query) {
  rsv_table_.clear();

  // Get RSVs from WMIR
  DBG(15) << "trying to contact WMIR ...";
  Socket sock(server_, port_);
  if(!sock.connected()) {
    DBG(10) << "not OK" << endl;
    ERR << "Could not connect! Make sure WMIR is started." << endl;
  } else {
    DBG(15) << "OK" << endl;
  }
  
  // Query WMIR and store the result in the rsv_table
  DBG(15) << "Sending query" << endl;
  sock << _query + "\r\n";
  string line;
  DBG(15) << "getting count" << endl;
  line = sock.getline();
  unsigned int count = atoi(line.c_str());
  DBG(15) << "It's " << count << endl;
    
  max_rsv_ = 0.0;
  if(count == 0) {
    DBG(10) << "No results!" << endl;
    max_rsv_ = 10000.0;
  } else {
    for(unsigned int i = 0; i < count; ++i) {
      line = sock.getline();
      istringstream is(line);
      string filename, rsv_s;
      is >> filename;
      is >> rsv_s;
      double rsv = atof(rsv_s.c_str());
      if(rsv > max_rsv_) {
        max_rsv_ = rsv;
      }
      rsv_table_[filename] = rsv;
      DBG(15) << filename << " " << rsv << endl;
    }
  }
  
  DBG(15) << "max_rsv="<<max_rsv_ << endl;
  
  DBG(15) << "Disconnecting" << endl;
  sock << ":bye\r\n";
  sock.close();
}

void TextFeatureDistance::getServerSettings(::std::string &server, unsigned &port, ::std::string &language)
{
	server = server_;
	port = port_;
	language = language_;
	return;
}


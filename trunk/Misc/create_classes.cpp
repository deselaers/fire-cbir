#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include "getpot.hpp"

using namespace std;

int main(int argc, char **argv) {
    ::std::string filelist, dir;
    GetPot cl(argc,argv);
    vector<string> ufos=cl.unidentified_options(4, "-f", "--filelist", "-d", "--directory");
    if(ufos.size()!=0) {
        cout << "USAGE:" << endl
             << "create_classes (-f,--filelist <filelist> | -d,--directory <dir>)" << endl
             << endl;
        exit(20);
    }

    if(cl.search(2,"-f","--filelist")) {
        filelist=cl.follow("",2,"-f","--filelist");
    } else if(cl.search(2,"-d","--directory")) {
        dir=cl.follow("",2,"-d","--directory");
    } else {
        cout << "USAGE:" << endl
             << "create_classes (-f,--filelist <filelist> | -d,--directory <dir>)" << endl
             << endl;
        exit(20);
    }
    
    if(!dir.empty()) {  // directory method
        chdir(dir.c_str());
        system("find . > allfiles");
        system("cat allfiles | sed 's/^/file /' > tmpdata");
        system("bash tmpdata | grep image > purefiles");
		ifstream is; is.open("purefiles");
        ofstream os;
        if(!is.good()) {
            cout << "ERROR!" << endl
                 << "  Internal error. Could not open tmpdata." << endl
                 << endl;
            exit(20);
        } else {
            os.open("filelist", ofstream::out | ofstream::trunc);
            os << "FIRE_filelist" << endl
               << "# add some suffixes" << endl
               << "classes yes" << endl
               << "path " << dir << endl
               << "# add path to type2bin file" << endl;
            
            ::std::string line;
            while(!is.eof()) {
                getline(is,line);
                if(line=="") break;
                istringstream iss(line);
                ::std::string file;
                iss >> file;
                file = file.substr(0, file.length()-1);
                os << "file " << file << endl;
            }
        }
        is.close();
        os.close();
        remove("allfiles");
        remove("tmpdata");
        remove("purefiles");
        filelist = dir+"filelist";
    }
    
    if(!filelist.empty()) {     // filelist method
        ifstream is; is.open(filelist.c_str());
		ofstream os;
        if(!is.good()) {
            cout << "ERROR!" << endl
                 << "  Could not open filelist." << endl
                 << endl;
            exit(20);
        } else {
            os.open((filelist+"_with_classes").c_str(), ofstream::out | ofstream::trunc);
            ::std::string line;
            vector<string> classes(0);
            while(!is.eof()) {
                    getline(is,line);
                    if(!is.eof()) {
                        istringstream iss(line);
                        string keyword;
                        iss >> keyword;
                        if("file"==keyword) {
                            uint file_class=0;
                            string filename;
                            iss >> filename;
                            ::std::string filepath=filename.substr(0, filename.rfind('/'));
                            bool class_exists=false;
                            for(uint i=0; i<classes.size(); i++) {
                                if(classes[i] == filepath) {
                                    class_exists=true;
                                    file_class=i;
                                    break;
                                }
                            }
                            if(class_exists==false) {
                                classes.push_back(filepath);
                                file_class=classes.size()-1;
                            }
                            os << "file " << filename << "     " << file_class << endl;
                        } else if("classes"==keyword) {
                            os << "classes yes" << endl;
                        } else {
                            os << line << endl;
                        }
                }
            }
        }
        is.close();
        os.close();
    }
    
    return 0;
}

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "getpot.hpp"
#include "diag.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

void USAGE() {
  cout << "USAGE: randomfilelistcreator [options] --ifilelist input filelist --ofilelist output filelists" << endl
       << " Options: " << endl
       << "  -h, --help:          shows this help" << endl
       << "  -tr, --trainsamples: number of images for training file" << endl
       << "  -te, --testsamples:  number of images for test file" << endl
       << "  -ctr, --checktrain:  check if there are enough train files" << endl
       << "  -cte, --checktest:   check if there are enough test files" << endl;
  exit(1);
}


int main(int argc, char** argv) {

  GetPot cl(argc, argv);

  if (cl.search(2, "-h", "--help")) {
    USAGE();
  }

  int trainSamples = 15;
  if (cl.search(2, "-tr", "--trainsamples")) {
    trainSamples = cl.next(15);
  }

  int testSamples = 50;
  if (cl.search(2, "-te", "--testsamples")) {
    testSamples = cl.next(50);
  }

  bool checkTest = false, checkTrain = false;
  if (cl.search(2, "-cte", "--checktest")) {
    checkTest = true;
  }
  if (cl.search(2, "-ctr", "--checktrain")) {
    checkTrain = true;
  }

  vector<string> imageFiles;
  if (!cl.search(1, "--ifilelist")) {
    USAGE();
  } else {
    string filename = cl.next(" ");
    ifstream ifs; ifs.open(filename.c_str());
    if (!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " << filename << ". Aborting." << endl;
      exit(20);
    }
    string fname;
    while (!ifs.eof()) {
      getline(ifs, fname);      
      if (fname != "") {
	imageFiles.push_back(fname);
      }
    }
    ifs.close();
  }

  string ofilelist = "";
  if (!cl.search(1, "--ofilelist")) {
    USAGE();
  } else {
    ofilelist = cl.next(" ");
  }


  ofstream ofs; 
  srand(time(NULL));

  if (trainSamples > 0) {
    string trainFile = ofilelist;
    trainFile.append("-train");
    ofs.open(trainFile.c_str());
    for (int i = 0; i < trainSamples; i++) {
      if (imageFiles.empty()) {
	if (checkTrain) {
	  ERR << "Base filelist contains too less entries !" << endl;
	  exit(2);
	}
	break;
      }
      int nextImage = rand()%imageFiles.size();
      ofs << imageFiles[nextImage].c_str() << endl;
      vector<string>::iterator it = imageFiles.begin();
      for (int j = 0; j < nextImage; j++) {
	it++;
      }
      imageFiles.erase(it);
    }
    ofs.close();
  }

  if (testSamples > 0) {
    string testFile = ofilelist;
    testFile.append("-test");
    
    srand(time(NULL));
    ofs.open(testFile.c_str());
    for (int i = 0; i < testSamples; i++) {
      if (imageFiles.empty()) {
	if (checkTest) {
	  ERR << "Base filelist contains too less entries !" << endl;
	  exit(2);
	}
	break;
      }
      int nextImage = rand()%imageFiles.size();
      ofs << imageFiles[nextImage].c_str() << endl;
      vector<string>::iterator it = imageFiles.begin();
      for (int j = 0; j < nextImage; j++) {
	it++;
      }
      imageFiles.erase(it);
    }
    ofs.close();
  }

}

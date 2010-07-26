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
#include <string>
#include <vector>
#include "sift.hpp"

#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"
#include "imagefeature.hpp"
#include "localfeatures.hpp"
#include "imagelib.hpp"
#include "pca.hpp"
#include "database.hpp"

using namespace std;


class SIFTExtractor {
private:

  string suffix;
  string pcasuffix;
  uint   levels;
  int    octaves;
  int    first;
  float  thresh;
  float  edgeThreshold;
  int n_strongest;
  bool use_color;
  bool do_pca;
  bool separate_rgb_pca;
  int pca_dim_in;
  int pca_dim_out;
  bool savePCA;
  bool loadPCA;
  string savepca_name;
  string loadpca_name;

public:

  SIFTExtractor():
    suffix("sift"),
    pcasuffix("pca.sift"),
    levels(3),
    octaves(-1),
    first(-1),
    thresh(0.04f / levels / 2.0f),
    edgeThreshold(10.0f),
    n_strongest(-1),
    use_color(false),
    do_pca(false),
    separate_rgb_pca(false),
    pca_dim_in(0),
    pca_dim_out(0),
    savePCA(false),
    loadPCA(false),
    savepca_name(""),
    loadpca_name("")
  {}

  void USAGE()
  {

    cout << endl << "USAGE:" << endl << endl
        << "extractsift [options] (--images filename1 [filename2 ... ]|--(fire)filelist <filelist>)" << endl << endl
        << "   Options:" << endl
        << "    -h, --help                   show this help" << endl
        << "    -s, --suffix <str>           suffix of output files (default: " << suffix << ")" << endl
        << "    -O, --octaves <int>          set the number of octaves (default: all)" << endl
        << "    -S, --level <int>            set the number of levels per octave (default: " << levels << ")" << endl
        << "    -F, --first-octave <int>     set the starting octave (default: " << first << ")" << endl
        << "    -t, --threshold <float>      set the SIFT threshold (default: " << thresh << ")" << endl
        << "    -e, --edge-threshold <float> set the edge rejection threshold (default: " << edgeThreshold << ")" << endl
        << "    -n, --n-strongest <int>      output only the n strongest keypoints (default: disabled)" << endl
        << "    -c, --color                  extract SIFT info for the three RGB channels" << endl
        << "    -g, --gray                   extract SIFT info on luminance of the image (default)" << endl
        << "    -p, --pca <int>              reduce dimensionality to the given value with PCA" << endl
				<< "    -P, --rgbpca                 treat the RGB channels separately in the PCA (requires -p and -c)" << endl
        << "    --savePCA <filename>         save the calculated PCA to the given file" << endl
        << "    --loadPCA <filename>         don't calculate pca on the local features created, but load the given one" << endl
        << endl
        << "Note on PCA: If PCA is done, both the untransformed and the transformed feature files will be created." << endl
        << "The transformed file will have 'pcaN.' prepended to the suffix, where N is the number of dimensions." << endl << endl
				<< "If -P is used, the dimensionality given with -p is the desired total dimensionality of the three color channels."
        << endl;
    exit(20);
  }

  void run(int argc, char** argv) {
    GetPot cl(argc,argv);

    //
    // command line parsing
    //
    if(cl.search(2,"--help","-h")) USAGE();

    if(cl.size()<2) USAGE();
    suffix        = cl.follow(suffix.c_str(), 2, "-s", "--suffix");
    octaves       = cl.follow(octaves, 2,        "-O", "--octaves");
    levels        = cl.follow((int)levels, 2,    "-S", "--level");
    first         = cl.follow(first, 2,          "-F", "--first-octave");
    thresh        = cl.follow(thresh, 2,         "-t", "--threshold");
    edgeThreshold = cl.follow(edgeThreshold, 2,  "-e", "--edge-threshold");
    n_strongest   = cl.follow(n_strongest, 2,    "-n", "--n-strongest");

    if(cl.search(2,"-c","--color")) {
      use_color = true;
    }

    if(cl.search(2,"-g","--gray")){
      use_color = false;
    }

    PCA *pca = NULL, *pca_r = NULL, *pca_g = NULL, *pca_b = NULL;

    if(cl.search(2,"-P","--rgbpca")) {
      separate_rgb_pca = true;
    } else {
      separate_rgb_pca = false;
    }

    if(cl.search(2,"-p","--pca")) {
      do_pca = true;

      if( (use_color) && (!separate_rgb_pca) ) {
        pca_dim_in = 384;
      } else {
        pca_dim_in = 128;
      }
    }

    if(cl.search("--savePCA")) {
      savePCA = true;
      savepca_name = cl.follow("sift.pca.gz","--savePCA");
    }

    if(cl.search("--loadPCA")) {
      loadPCA = true;
      loadpca_name = cl.follow("sift.pca.gz","--loadPCA");
    }

    // Make PCA objects
    if( do_pca || loadPCA ) {
      if(!separate_rgb_pca) {
        pca = new PCA(pca_dim_in);

        pca_dim_out = cl.follow(pca_dim_out, 2, "-p", "--pca");
      } else {
        pca_r = new PCA(pca_dim_in);
        pca_g = new PCA(pca_dim_in);
        pca_b = new PCA(pca_dim_in);

        pca_dim_out = cl.follow(pca_dim_out, 2, "-p", "--pca") / 3;
      }

      if(pca_dim_out == 0) {
        cout << "Please specify the PCA output dimensionality!" << endl;
        exit(1);
      }

      char pcasuffix_c[20] = {'\0'};
      sprintf(pcasuffix_c, "pca%u.%s", pca_dim_out, suffix.c_str());

      pcasuffix = pcasuffix_c;
    }

    //
    // get list of files to be processed
    //
    vector<string> infiles;

    if(cl.search("--images")) {
      string filename=cl.next(" ");;
      while(filename!=" ") {
        infiles.push_back(filename);
        filename=cl.next(" ");
      }
    } else if (cl.search("--filelist")) {
      string filename="test";
      igzstream ifs; ifs.open(cl.follow("list","--filelist"));
      if(!ifs.good() || !ifs) {
        ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
        exit(20);
      }
      while(!ifs.eof() && filename!="") {
        getline(ifs,filename);
        if(filename!="") {
          infiles.push_back(filename);
        }
      }
      ifs.close();
    } else if (cl.search("--firefilelist")) {
    	Database db;
    	db.loadFileList(cl.follow("list","--firefilelist"));
    	for(uint i=0; i<db.size(); ++i) {
    		infiles.push_back(db.path()+"/"+db.filename(i));
    	}
    } else {
      USAGE();
      exit(20);
    }

    //
    // processing the files
    //
    ImageFeature im;
    for(uint i=0;i<infiles.size();++i) {
        string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;

      if( !im.load(filename) ) {
        ERR << "Image load error! Exiting." << endl;
        exit(1);
      }


      float *raw_image_r = NULL, *raw_image_g = NULL, *raw_image_b = NULL;

      if( use_color ) {
        // Get color channels
        raw_image_r = (float*)malloc(sizeof(float) * im.xsize() * im.ysize());
        raw_image_g = (float*)malloc(sizeof(float) * im.xsize() * im.ysize());
        raw_image_b = (float*)malloc(sizeof(float) * im.xsize() * im.ysize());

        for(unsigned long y = 0; y < im.ysize(); ++y) {
          for(unsigned long x = 0; x < im.xsize(); ++x) {
            raw_image_r[im.xsize() * y + x] = im(x,y,0);
            raw_image_g[im.xsize() * y + x] = im(x,y,1);
            raw_image_b[im.xsize() * y + x] = im(x,y,2);
          }
        }
      }

      // Convert the image to gray (we need this for the keypoint extraction,
      // so it it done regardless if we use color or not
      im = makeGray(im, 2);

      // Get raw image
      float *raw_image = (float*)malloc(sizeof(float) * im.xsize() * im.ysize());
      for(unsigned long x = 0; x < im.xsize(); ++x) {
        for(unsigned long y = 0; y < im.ysize(); ++y) {
          raw_image[im.xsize() * y + x] = im(x,y,0);
        }
      }

      //
      // Gaussian scale space
      //
      DBG(10) << "Computing Gaussian scale space ... " << endl;

      int         O      = octaves;
      int const   S      = levels;
      int const   omin   = first;
      float const sigman = .5;
      float const sigma0 = 1.6 * powf(2.0f, 1.0f / S);

      // optionally autoselect the number number of octaves
      // we downsample up to 8x8 patches
      if(O < 1)
        O = std::max(int (std::floor(log2(std::min(im.xsize(),im.ysize()))) - omin -3), 1);

      // initialize scalespace
      VL::Sift *sift = new VL::Sift(raw_image, im.xsize(), im.ysize(), sigman, sigma0, O, S, omin, -1, S+1);

      VL::Sift *sift_r = NULL, *sift_g = NULL, *sift_b = NULL;

      if( use_color ) {
        sift_r = new VL::Sift(raw_image_r, im.xsize(), im.ysize(), sigman, sigma0, O, S, omin, -1, S+1);
        sift_g = new VL::Sift(raw_image_g, im.xsize(), im.ysize(), sigman, sigma0, O, S, omin, -1, S+1);
        sift_b = new VL::Sift(raw_image_b, im.xsize(), im.ysize(), sigman, sigma0, O, S, omin, -1, S+1);
      }

      DBG(10) << "(" << O << " octaves from " << omin << ", " << S << " levels per octave) done." << endl;


      //
      // SIFT detector
      //
      DBG(10) << "Running SIFT detector (threshold:" << thresh << ", edge-threshold: " << edgeThreshold << ") ..." << flush;

      sift->detectKeypoints(thresh, edgeThreshold);

      DBG(10) << " done (" << sift->keypointsEnd() - sift->keypointsBegin() << " keypoints)." << endl;

      if(n_strongest > 0) {
        DBG(10) << "Selecting " << n_strongest << " strongest ... " << flush;
        sift->selectNStrongest(n_strongest);
        DBG(10) << "OK" << endl;
      }

      //
      // SIFT orientations and descriptors
      //
      DBG(10) << "Extracting orientations and descriptors." << endl;

      vector<pair<double,double> > positions;
      vector<double> sigmas;
      vector<double> orientations;
      vector<vector<double> > descriptors;

      for(VL::Sift::KeypointsConstIter iter = sift->keypointsBegin(); iter != sift->keypointsEnd(); ++iter )
      {
        // detect orientations
        VL::float_t angles[8];
        int nangles = sift->computeKeypointOrientations(angles, *iter);

        // compute descriptors
        for(int a = 0 ; a < nangles ; ++a) {
          positions.push_back(make_pair(iter->x, iter->y));

          sigmas.push_back(iter->sigma);
          orientations.push_back(angles[a]);

          VL::float_t descr[128], descr_r[128], descr_g[128], descr_b[128];

          if( !use_color ) {
            sift->computeKeypointDescriptor(descr, *iter, angles[a]);
          } else {
            sift_r->computeKeypointDescriptor(descr_r, *iter, angles[a]);
            sift_g->computeKeypointDescriptor(descr_g, *iter, angles[a]);
            sift_b->computeKeypointDescriptor(descr_b, *iter, angles[a]);
          }

          vector<double> descrvec;

          if( !use_color ) {
            for(int i = 0; i < 128; ++i) {
              descrvec.push_back(descr[i]);
            }

            descriptors.push_back(descrvec);
          } else {
            // 128 * red + 128 * green + 128 * blue = 384

            for(int i = 0; i < 128; ++i) {
              descrvec.push_back(descr_r[i]);
            }

            for(int i = 0; i < 128; ++i) {
              descrvec.push_back(descr_g[i]);
            }

            for(int i = 0; i < 128; ++i) {
              descrvec.push_back(descr_b[i]);
            }

            descriptors.push_back(descrvec);
          }

          if(do_pca && !loadPCA) {
            if(use_color && separate_rgb_pca) {
              vector<double> descrvec_r, descrvec_g, descrvec_b;
              int index = 0;

              for(int i = 0; i < 128; ++i) {
                descrvec_r.push_back(descrvec.at(index++));
              }

              for(int i = 0; i < 128; ++i) {
                descrvec_g.push_back(descrvec.at(index++));
              }

              for(int i = 0; i < 128; ++i) {
                descrvec_b.push_back(descrvec.at(index++));
              }

              pca_r->putData(descrvec_r);
              pca_g->putData(descrvec_g);
              pca_b->putData(descrvec_b);
            } else {
              pca->putData(descrvec);
            }
          }
        } // next angle
      } // next keypoint


      //
      // Save as local feature
      //
      LocalFeatures lf;
      lf.filename_ = filename;
      if(use_color) { lf.dim_ = 384; }
      else { lf.dim_ = 128; }
      lf.imageSizeX_ = im.xsize();
      lf.imageSizeY_ = im.ysize();

      for(unsigned int i=0; i<positions.size(); ++i) {
        lf.addLocalFeature(descriptors[i], make_pair((uint)positions[i].first,(uint)positions[i].second), (uint)sigmas[i]);
      }

      lf.save(filename+"."+suffix);

      free(raw_image);

      delete sift;

      if(use_color) {
        delete sift_r;
        delete sift_g;
        delete sift_b;
        free(raw_image_r);
        free(raw_image_g);
        free(raw_image_b);
      }

      DBG(20) << "Finished with '" << filename << "'." << endl;
    }

    if(do_pca && !loadPCA) {
      if(use_color && separate_rgb_pca) {
        pca_r->dataEnd();
        pca_g->dataEnd();
        pca_b->dataEnd();
        DBG(10) << "Calculating PCA" << endl;
        pca_r->calcPCA();
        pca_g->calcPCA();
        pca_b->calcPCA();
        DBG(10) << "OK" << endl;
      } else {
        pca->dataEnd();
        DBG(10) << "Calculating PCA" << endl;
        pca->calcPCA();
        DBG(10) << "OK" << endl;
      }

      if(savePCA) {
        pca->save(savepca_name);

        DBG(10) << "PCA saved as " << savepca_name << endl;
      }
    }

    if( loadPCA ) {
      DBG(10) << "Loading PCA from " << loadpca_name << endl;
      pca->load(loadpca_name);
    }

    if( do_pca || loadPCA ) {
      DBG(10) << "PCA transforming features" << endl;

      // Do a second iteration over the data for PCA transformation

      for(uint i=0;i<infiles.size();++i) {
        LocalFeatures lf;

        string lffilename=infiles[i]+"."+suffix;

        lf.load(lffilename);

        if(use_color && separate_rgb_pca) {
          for(uint j = 0; j < lf.numberOfVectors(); ++j) {
            // Split the vector into RGB first
            vector<double> descrvec_r, descrvec_g, descrvec_b;
            int index = 0;

            for(int i = 0; i < 128; ++i) {
              descrvec_r.push_back(lf[j].at(index++));
            }

            for(int i = 0; i < 128; ++i) {
              descrvec_g.push_back(lf[j].at(index++));
            }

            for(int i = 0; i < 128; ++i) {
              descrvec_b.push_back(lf[j].at(index++));
            }

            // Peform the PCA
            vector<double> transformed_r = vector<double>(pca_r->transform(descrvec_r, pca_dim_out));
            vector<double> transformed_g = vector<double>(pca_g->transform(descrvec_g, pca_dim_out));
            vector<double> transformed_b = vector<double>(pca_b->transform(descrvec_b, pca_dim_out));

            // Put the data back into the local feature
            lf[j].resize(pca_dim_out * 3);

            index = 0;

            for(int k = 0; k < pca_dim_out; ++k) {
              lf[j][index++]=transformed_r[k];
            }

            for(int k = 0; k < pca_dim_out; ++k) {
              lf[j][index++]=transformed_g[k];
            }

            for(int k = 0; k < pca_dim_out; ++k) {
              lf[j][index++]=transformed_b[k];
            }
          }

          lf.dim() = pca_dim_out * 3;
        } else {
          for(uint j=0; j<lf.numberOfVectors(); ++j) {
            vector<double> transformed = vector<double>(pca->transform(lf[j], pca_dim_out));
            lf[j].resize(pca_dim_out);
            for(int k = 0; k < pca_dim_out; ++k) {
              lf[j][k]=transformed[k];
            }
          }

          lf.dim() = pca_dim_out;
        }
        lf.save(infiles[i]+"."+pcasuffix);
      }
      DBG(10) << "OK" << endl;
    }

  }
};

int main(int argc, char** argv)
{
  SIFTExtractor se;
  se.run(argc,argv);
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}

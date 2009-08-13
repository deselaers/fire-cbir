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

#include <iostream>
#include <sstream>
#include "float.h"
#include "gzstream.hpp"
#include "gaussiandensity.hpp"
#include "localfeatures.hpp"
#include "diag.hpp"
#include "getpot.hpp"
#include "database.hpp"
#include "positionclusterer.hpp"

using namespace std;

#define CL_APPEARANCE  1
#define CL_POSITION 2
#define CL_RELATIVE_POSITION 4

// represents a classified image, with a number and a probability
// used to sort images according to the classification accuracy
struct ClassifiedImage {
  int nr;
  double p;
  ClassifiedImage() {
    nr = 0;
    p = 0.0;
  }
};

// sort function for 'ClassifiedImage' struct (sort according to probability)
struct sortClassifications : public binary_function<ClassifiedImage, ClassifiedImage, bool> {
  bool operator()(ClassifiedImage x, ClassifiedImage y) { return x.p < y.p; }
};

// data structure to represent (untied) clusters
struct ClusterData {
  GaussianDensity density;
  int clazz;
  vector<ClusterPosition> positions;
  vector< vector<ClusterPosition > > positionDifferences;
  vector<bool> hasPositionDifference;
  double weight;
  // initialize with default values
  ClusterData() {
    weight = 0.0;
  }
};

// options affecting the emission probability calculation
struct EPOptions {
  vector<bool> useDim;
  bool square, dimensionPooling;
};

// prints the usage
void USAGE() {
  cout << "USAGE: trainclusterweightsdiscriminative [options] --models <filelist> --train <filelists> [classification]" << endl
       << " Classification: one of 'appearance', 'position', or both" << endl
       << " Options: " << endl
       << "  -e, --epsilon <epsilon>            : update parameter for iteration" << endl
       << "  -eq, --eqapriori                   : equal a priori probabilities" << endl
       << "  -h, --help                         : shows this help" << endl
       << "  -i, --iterations <numiter>         : iterate the given number of times" << endl
       << "  -l, --loadfile <filename>          : load clusterweights from given file" << endl
       << "  -p, --productrule                  : use product rule instead of sum rule in training" << endl
       << "  -s, --savefile <filename>          : save clusterweights in given file" << endl
       << "  -w, --worstimages <n>              : use only the worst n% (0.0 < n <= 100.0) images for updating" << endl
       << "                                       clusterweights in each iteration" << endl
       << "  -ei, --eachIter <n>                : save clusterweights each n-th iteration (default:1 = every)" << endl
       << "  --test <filelists>                 : for classification test data" << endl
       << "  -var, --variance <value>           : set variance of all clusters to 'value/PI'" << endl
       << "  -mvar, --multiplyvariance <factor> : multiply variances of all clusters by the given factor" << endl;
  exit(1);
}

// loads a cluster model stored in the file with the given filename
void loadModel(const string filename, int classNr, vector<ClusterData>& clusters) {
  igzstream is; is.open(filename.c_str());
  vector<GaussianDensity> densities;
  DBG(15) << "loading cluster model for class " << classNr << " from file " << filename << endl;
  if (!is || !is.good()) {
    ERR << "Unable to read file '" << filename << "'." << endl;
  } else {
    string line;
    getline(is, line);
    if (!is.good()) {
      ERR << "Error reading file '" << filename << "'." << endl;
    } else {
      if (line != "# EM/LBG clustering model file") {
        ERR << "This is probably not an EM-Model file" << endl
            << "Continuing anyway" << endl;
      }
      while (!is.eof()) {
        getline(is,line);
        if (!is.eof()) {
          istringstream iss(line);
          string keyword;
          iss >> keyword;
          if (keyword == "gaussians" ) {
            uint size;
            iss >> size;
            densities.resize(size);
          } else if (keyword == "density") {
            uint no;
            iss >> no;
            iss >> keyword;
            if (keyword == "elements") {
              iss >> densities[no].elements;
            } else if (keyword == "dim") {
              iss >> densities[no].dim;
            } else if (keyword == "mean") {
              GaussianDensity& c = densities[no];
              c.mean = vector<double>(c.dim);
              for (uint i = 0; i < c.dim; ++i) {
                iss >> c.mean[i];
              }
            } else if (keyword == "variance") {
              GaussianDensity& c = densities[no];
              c.sigma = vector<double>(c.dim);
              for (uint i = 0; i < c.dim; ++i) {
                iss >> c.sigma[i];
              }
            } else {
              ERR << "Reading density received unknown keyword in position 3: '" << keyword << "'." << endl;
            }
          } else if ((keyword == "splitmode") || (keyword == "maxsplits") || (keyword == "stopWithNClusters") || (keyword == "disturbMode") ||
                     (keyword == "poolMode") || (keyword == "dontSplitBelow") || (keyword == "iterationsBetweenSplits") || (keyword == "epsilon") ||
                     (keyword == "minObservationsPerCluster") || (keyword == "distance")) {
            // ignore these keywords
          } else {
            ERR << "Unknown keyword '" << keyword << "'." << endl;
          }
        }
      }
    }
    is.close();
  }
  for (int c = 0; c < (int) densities.size(); c++) {
    ClusterData newCluster;
    newCluster.clazz = classNr;
    newCluster.density = densities[c];
    clusters.push_back(newCluster);
  }
  DBG(20) << "read " << densities.size() << " clusters for class " << classNr << endl;
}

// returns the number of classes found in the given database
// the method does not look for distinct class labels, but
// returns the highest class label + 1
// therefore, classes are assumed to be labeled starting from
// 0 to n-1 
int getNumClasses(const Database& db) {
  int maxClassNr = -1;
  for (int n = 0; n < (int) db.size(); n++) {
    int clazz = db[n]->clas();
    if (clazz < 0) {
      ERR << "found image with class " << clazz << " in database!" << endl;
      exit(1);
    } else {
      if (clazz > maxClassNr) {
      	maxClassNr = clazz;
      }
    }
  }
  return maxClassNr + 1;
}


void calcInitialClusterweights(vector<ClusterData>& clusters, int numClasses) {
  DBG(15) << "calculating initial weights" << endl;
  vector<int> localFeaturesPerClass(numClasses, 0);
  for (int c = 0; c < (int) clusters.size(); c++) {
    localFeaturesPerClass[clusters[c].clazz] += clusters[c].density.elements;
  }
  vector<double> weightSum(numClasses, 0.0);
  for (int c = 0; c < (int) clusters.size(); c++) {
    clusters[c].weight = double(clusters[c].density.elements) / double(localFeaturesPerClass[clusters[c].clazz]);
    weightSum[clusters[c].clazz] += clusters[c].weight;
    DBG(20) << "for cluster " << c << " (class " << clusters[c].clazz << "): " << clusters[c].weight << endl;
  }
  for (int k = 0; k < numClasses; k++) {
    if (fabs(weightSum[k] - 1.0) > 0.00001) {
      ERR << "sum of all clusterweights for class " << k << " should be 1.0, but is " << weightSum[k] << endl;
      exit(1);
    }
  }
}

void loadInitialClusterweights(vector<ClusterData>& clusters, const string& filename) {
  ifstream ifs; ifs.open(filename.c_str());
  if (!ifs.good()) {
    ERR << "clusterweights file '" << filename << "' could not be opened!" << endl;
    exit(1);
  }
  int numClusters;
  ifs >> numClusters;
  for (int c = 0; c < numClusters; c++) {
    int clazz;
    ifs >> clazz;
    if (clusters[c].clazz != clazz) {
      ERR << "clusterweights file and clusterfilelist do not match !" << endl;
      exit(1);
    }
    ifs >> clusters[c].weight;
  }
  int end;
  ifs >> end;
  if (end != -1) {
    ERR << "clusterweights file did not end with '-1' !" << endl;
    exit(1);
  }
  ifs.close();
}

void saveClusterWeights(const vector<ClusterData>& clusters, const string& filename) {
  ofstream ofs; ofs.open(filename.c_str());
  if (!ofs || !ofs.good()) {
    ERR << "cannot save clusterweights to file '" << filename << "'" << endl;
    exit(1);
  }
  ofs << clusters.size() << endl;
  for (int c = 0; c < (int) clusters.size(); c++) {
    ofs << clusters[c].clazz << endl;
    ofs << clusters[c].weight << endl;
  }
  ofs << "-1" << endl;
  ofs.close();
}

// modifies the clusters' variances, if desired, and sets flags for classification accordingly
void adaptClusterVariances(vector<ClusterData>& clusters, double varianceValue, double varianceFactor, bool& square, bool& dimensionPooling) {
  if (varianceValue != 0.0) {
    // a fixed pooled variance is desired
    if (varianceFactor != 0.0) {
      ERR << "Options -var and -mvar are mutually exclusive !" << endl;
      exit(1);
    } else {
      square = false;
      dimensionPooling = true;
      DBG(10) << "setting all clusters' variances to " << (varianceValue / M_PI) << "/M_PI" << endl;
      for (int c = 0; c < (int) clusters.size(); c++) {
        for (int d = 0; d < (int) clusters[c].density.dim; d++) {
          clusters[c].density.sigma[d] = varianceValue / M_PI;
        }
      }
    }
  } else {
    square = true;
    if (varianceFactor != 0.0) {
      DBG(10) << "multiplying all clusters' variances by " << varianceFactor << endl;
      for (int c = 0; c < (int) clusters.size(); c++) {
        for (int d = 0; d < (int) clusters[c].density.dim; d++) {
          clusters[c].density.sigma[d] *= varianceFactor;
        }
      }
    }
    // test if variances are pooled by the clustering process
    dimensionPooling = true;
    if (dimensionPooling) {
      for (int c = 0; c < (int) clusters.size(); c++) {
        if (dimensionPooling) {
          double curVariance = clusters[c].density.sigma[0];
          for (int d = 1; d < (int) clusters[c].density.dim; d++) {
            if (clusters[c].density.sigma[d] != curVariance) {
              dimensionPooling = false;
              break;
            }
          }
        }
      }
    }
    if (dimensionPooling) {
      DBG(10) << "cluster variances are dimension-pooled" << endl;
    } else {
      DBG(10) << "cluster variances are not dimension-pooled" << endl;
    }
  }
}

// precalculates the normalization values for the cluster emission probabilities
// this is possible since the normalization values are independent from the local features
// in case of no dimension pooling, this helps speedup the emission probability calculation
void calcClusterNormalizationValues(const vector<ClusterData>& clusters, const EPOptions& options, vector<double>& clusterNormalizationValues) {
  clusterNormalizationValues.resize(clusters.size());
  for (int c = 0; c < (int) clusters.size(); c++) {
    clusterNormalizationValues[c] = 0.0;
    for (int d = 0; d < (int) clusters[c].density.dim; d++) {
      if (options.useDim[d]) {
        clusterNormalizationValues[c] += log(2.0 * M_PI * clusters[c].density.sigma[d]);
      }
    }
  }
}

double calc_position_probability(const pair<double, double> mean, const PositionCovarianceMatrix& cov, const pair<double, double>& actualPosition) {
  double xDiff = actualPosition.first - mean.first;
  double yDiff = actualPosition.second - mean.second;
  PositionCovarianceMatrix inverted = cov.invert();
  return (1.0 / (2.0 * M_PI * sqrt(cov.getDeterminant()))) *
    exp(-0.5 * (xDiff * xDiff * inverted.xx + xDiff * yDiff * (inverted.xy + inverted.yx) + yDiff * yDiff * inverted.yy));
}


// calculates the (logarithmic) cluster emission probability p(x|c,k)
// 
double calc_log_x_given_c_k(const GaussianDensity& clustermodel, double clusterNormalizationValue, 
                            const vector<double>& x_l, const EPOptions& options) {
  double pixelNorm = 0.0, pixelDist = 0.0;

  // TODO: remove this consistency check for the sake of efficieny !
  if (x_l.size() != clustermodel.dim) {
    ERR << "size of local feature (" << x_l.size() << ") and of clustermodel (" << clustermodel.mean.size() << ") do not match !" << endl;
    exit(3);
  }

  int dims = (int) x_l.size();

  for (int d = 0; d < dims; d++) {
    if (options.useDim[d]) {
      double tmp = 0.0;
      if (options.square) {
        tmp = x_l[d] - clustermodel.mean[d];
        tmp *= tmp;
      } else {
        tmp = fabs(x_l[d] - clustermodel.mean[d]);
      }
      if (!options.dimensionPooling) {
        pixelDist += (tmp / clustermodel.sigma[d]);
      } else {
        pixelDist += tmp;
      }
    }
  }

  if (options.dimensionPooling) {
    pixelDist /= clustermodel.sigma[0];
  }
  pixelNorm = clusterNormalizationValue;

  //  return -0.5 * ((1.0 - alpha) * pixelDist + alpha * posDist + (1.0 - alpha) * pixelNorm + alpha * posNorm);
  return -0.5 * (pixelDist + pixelNorm);
}


double calc_relative_position_probability(const LocalFeatures& localFeatures, int l, int k, const vector< vector<int> >& nearestClusterSaved,
                                          const pair<double, double>& position, const vector<ClusterData>& clusters) {

  double relpos_given_c_k = 0.0;
  int c = nearestClusterSaved[l][k];

  int count = 0;
  for (int l_ = 0; l_ < (int) localFeatures.numberOfFeatures(); l_++) {
    if (l_ != l) {
      int c_ = nearestClusterSaved[l_][k];
      if (clusters[c].hasPositionDifference[c_]) {
        count++;
        pair<double, double> relativePosition = pair<double, double>(double(localFeatures.position(l_).x) / double(localFeatures.imageSizeX()),
                                                                     double(localFeatures.position(l_).y) / double(localFeatures.imageSizeY()));
        relativePosition.first -= position.first;
        relativePosition.second -= position.second;

        int numElements = 0;
        double relpos = 0.0;
        for (int i = 0; i < (int) clusters[c].positionDifferences[c_].size(); i++) {
          relpos += double(clusters[c].positionDifferences[c_][i].elements) * 
            calc_position_probability(clusters[c].positionDifferences[c_][i].mean,
                                      clusters[c].positionDifferences[c_][i].covariance, relativePosition);
          if (isnan(relpos)) {
            cout << calc_position_probability(clusters[c].positionDifferences[c_][i].mean,
                                              clusters[c].positionDifferences[c_][i].covariance, relativePosition) << endl;
            cout << c << "," << c_ << "," << i << endl;
            cout << clusters[c].positionDifferences[c_][i].mean.first << "," << clusters[c].positionDifferences[c_][i].mean.second << endl;
            cout << clusters[c].positionDifferences[c_][i].covariance.xx << "," << clusters[c].positionDifferences[c_][i].covariance.xy << "," << clusters[c].positionDifferences[c_][i].covariance.yx << "," << clusters[c].positionDifferences[c_][i].covariance.yy << endl;
            cout << relativePosition.first << "," << relativePosition.second << endl;
            cout << "# elements:" << endl;
            cout << clusters[c].positionDifferences[c_][i].elements << endl;
            exit(0);
          }
          numElements += clusters[c].positionDifferences[c_][i].elements;
        }
        relpos_given_c_k += relpos / double(numElements);
      }
    }
  }

  relpos_given_c_k /= (double) count;
  return relpos_given_c_k;
}

// calculates for a single local feature x_l the probability p(x_l|k) for each k
// and stores them in the 'p_x_given_k' vector
// if desired, the individual (weighted) emission probabilities that
// sum up to p(x|k) are stored in the given vector 'saved_weighted_emission_probability'
void calc_x_given_k(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, const LocalFeatures& localFeatures,
                    int l, const vector< vector<int> >& nearestClusterSaved, const pair<double, double>& position, const EPOptions& options, 
                    int classificationMode, vector<double>& saved_weighted_emission_probability, vector<double>& p_x_given_k) {

  int numClusters = (int) clusters.size();
  int numClasses = (int) p_x_given_k.size();
  bool save = (saved_weighted_emission_probability.size() > 0);
  if (save) {
    if (saved_weighted_emission_probability.size() != clusters.size()) {
      ERR << "size of 'save_weighted_emission_probability' and clustermodel does not match !" << endl;
      ERR << saved_weighted_emission_probability.size() << endl;
      ERR << clusters.size() << endl;
      exit(1);
    }
  }

  vector<double> relpos_given_c_k(numClasses, 0.0);
  if ((classificationMode & CL_RELATIVE_POSITION) > 0) {
    for (int k = 0; k < numClasses; k++) {
      relpos_given_c_k[k] = calc_relative_position_probability(localFeatures, l, k, nearestClusterSaved, position, clusters);
    }
  }

  for (int c = 0; c < numClusters; c++) {
    double weighted_emission_probability;
    if ((classificationMode & CL_POSITION) > 0) {
      double pos_given_c_k = 0.0;
      int numElements = 0;
      for (int c_ = 0; c_ < (int) clusters[c].positions.size(); c_++) {
        double cur_pos_given_c_k = calc_position_probability(clusters[c].positions[c_].mean, clusters[c].positions[c_].covariance, position);
        pos_given_c_k += cur_pos_given_c_k * (double) clusters[c].positions[c_].elements;
        numElements += clusters[c].positions[c_].elements;
      }
      pos_given_c_k /= (double) numElements;
      if ((classificationMode & CL_RELATIVE_POSITION) > 0) {
        weighted_emission_probability = clusters[c].weight * exp(calc_log_x_given_c_k(clusters[c].density, clusterNormalizationValues[c], localFeatures[l], options)) * pos_given_c_k * relpos_given_c_k[clusters[c].clazz];
      } else {
        weighted_emission_probability = clusters[c].weight * exp(calc_log_x_given_c_k(clusters[c].density, clusterNormalizationValues[c], localFeatures[l], options)) * pos_given_c_k;
      }
    } else {
      weighted_emission_probability = clusters[c].weight * exp(calc_log_x_given_c_k(clusters[c].density, clusterNormalizationValues[c], localFeatures[l], options));
    }
    if (save) {
      saved_weighted_emission_probability[c] = weighted_emission_probability;
    }
    if (clusters[c].clazz >= numClasses) {
      ERR << "cluster " << c << " belongs to class " << clusters[c].clazz << ", but p(x|k) is dimensioned only for " << numClasses << " classes!" << endl;
      exit(1);
    }
    p_x_given_k[clusters[c].clazz] += weighted_emission_probability;
  }
}

// calculates for a single local feature x_l the probability p(k|x_l) for each k
// and stores them in the p_k_given_x vector
// if desired, the probabilities p(x_l|k) that are used to calculate p(k|x_l) using Bayes
// are stored in the given vector 'saved_x_given_k'
void calc_k_given_x(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, const vector<double>& aPriori, 
                    const LocalFeatures& localFeatures, int l, const vector< vector<int> >& nearestClusterSaved, 
                    const pair<double, double>& position, const EPOptions& options, int classificationMode, 
                    vector<double>& saved_x_given_k, vector<double>& saved_weighted_emission_probability, vector<double>& p_k_given_x) {

  int numClasses = (int) aPriori.size();
  bool save = (saved_x_given_k.size() > 0);
  if (save) {
    if ((int) saved_x_given_k.size() != numClasses) {
      ERR << "saved_x_given_k incorrectly dimensionized !" << endl;
      exit(1);
    }
  }
  if ((int) p_k_given_x.size() != numClasses) {
    ERR << "saved_x_given_k incorrectly dimensionized !" << endl;
    exit(1);
  }

  vector<double> p_x_given_k(numClasses, 0.0);
  calc_x_given_k(clusters, clusterNormalizationValues, localFeatures, l, nearestClusterSaved,
                 position, options, classificationMode, saved_weighted_emission_probability, p_x_given_k);
  double denominator = 0.0;

  for (int k = 0; k < numClasses; k++) {
    if (save) {
      saved_x_given_k[k] = p_x_given_k[k];
    }
    double product = aPriori[k] * p_x_given_k[k];
    p_k_given_x[k] = product;
    denominator += product;
  }

  for (int k = 0; k < numClasses; k++) {
    // TODO: remove this check for the sake of efficiency !
    if ((p_k_given_x[k] / denominator < 0.0) || (p_k_given_x[k] / denominator > 1.0) || isnan(p_k_given_x[k] / denominator)) {
      ERR << "not a valid probability value in calc_k_given_x: " << (p_k_given_x[k] / denominator) << endl;
      ERR << p_k_given_x[k] << endl;
      ERR << denominator << endl;
      exit(1);
    }
    p_k_given_x[k] /= denominator;
  }
}

int getNearestCluster(const vector<double>& x_l, const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues,
                      int k, const EPOptions& options) {
  int bestCluster = -1;
  bool first = true;
  double maxProb = 0.0;
  for (int c = 0; c < (int) clusters.size(); c++) {
    if (clusters[c].clazz == k) {
      double prob = calc_log_x_given_c_k(clusters[c].density, clusterNormalizationValues[c], x_l, options);
      if (first || (prob >= maxProb)) {
      	maxProb = prob;
        bestCluster = c;
        first = false;
      }
    }
  }
  if ((bestCluster < 0) || (bestCluster >= (int) clusters.size())) {
    ERR << "wrong nearest cluster calculated: " << bestCluster << endl;
    exit(1);
  }
  return bestCluster;
}

double calcFderivation_sumrule(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, 
                               const vector< pair<LocalFeatures, uint> >& localfeatures, const vector< vector< vector<int> > >& nearestClusterSaved,
                               const vector<double>& aPriori, const EPOptions& options, 
                               int classificationMode, int k, int c,  vector< vector< vector<double> > >& saved_k_given_x, 
                               vector< vector< vector<double> > >& saved_x_given_k, vector< vector< vector<double> > >& saved_weighted_emission_probability,
                               vector<ClassifiedImage>& save_classification_order, double& probUpdatedImagesFormerIteration,
                               double& probUpdatedImages, int numImagesForUpdating, bool isFirstWeight) {

  int numTrainingSamples = (int) localfeatures.size();
  int numClasses = (int) aPriori.size();

  int ok = 0, error = 0;
  double sum = 0.0;
  vector<double> average_correct_class_probability(numClasses, 0.0);
  vector<int> class_items(numClasses, 0);
  vector<double> histogram(20, 0);
  vector<double> updates(numTrainingSamples);
  
  vector<int> fomerUpdatedImages(numImagesForUpdating);
  if (isFirstWeight) {
    for (int i = 0; i < numImagesForUpdating; i++) {
      fomerUpdatedImages[i] = save_classification_order[i].nr;
    }
  }
  
  for (int n = 0; n < numTrainingSamples; n++) {
    int k_n = (int) localfeatures[n].second;
    int numLocalFeatures = localfeatures[n].first.numberOfFeatures();
    double denominator = 0.0;
    double correct_class_probability = 0.0;
    for (int l = 0; l < numLocalFeatures; l++) {
      if (isFirstWeight) {
        pair<double, double> position = pair<double, double>(double(localfeatures[n].first.position(l).x) / double(localfeatures[n].first.imageSizeX()),
                                                             double(localfeatures[n].first.position(l).y) / double(localfeatures[n].first.imageSizeY()));
        calc_k_given_x(clusters, clusterNormalizationValues, aPriori, localfeatures[n].first, l, 
                       nearestClusterSaved[n], position, options, classificationMode,
                       saved_x_given_k[n][l], saved_weighted_emission_probability[n][l], saved_k_given_x[n][l]);
      }
      denominator += saved_k_given_x[n][l][k_n];
    }

    correct_class_probability = denominator / double(numLocalFeatures);
    if (isFirstWeight) {
      average_correct_class_probability[k_n] += correct_class_probability;
      class_items[k_n]++;
      if (correct_class_probability >= 0.5) {
        ok++;
      } else {
        error++;
        DBG(20) << "p(" << k << "|{x_1_l}_" << n << ")=" << correct_class_probability << endl;
      }
      int bin = (int) (correct_class_probability * 20);
      if (bin > 19) bin = 19;
      if (bin < 0) exit(1);
      histogram[bin]++;
    }

    if (isFirstWeight) {
      save_classification_order[n].nr = n;
      save_classification_order[n].p = correct_class_probability;
    }

    double product = 1.0 / denominator;
    double inner_sum = 0.0;
    for (int l = 0; l < numLocalFeatures; l++) {
      double u = aPriori[k_n] * saved_x_given_k[n][l][k_n];
      double u_;
      if (k == k_n) {
        u_ = aPriori[k_n] * saved_weighted_emission_probability[n][l][c];
      } else {
        u_ = 0.0;
      }
      double v = 0.0;
      for (int k_ = 0; k_ < numClasses; k_++) {
        v += aPriori[k_] * saved_x_given_k[n][l][k_];
      }
      double v_ = aPriori[k] * saved_weighted_emission_probability[n][l][c];
      inner_sum  += (u_ * v - u * v_) / (v * v);
    }
    product *= inner_sum;
    updates[n] = product;
  }

  if (isFirstWeight) {
    // calculate the new correct class probability for those images
    // which were formerly (i.e. in the last iteration) classified worst
    // this probability must have been increased, otherwise the step size (epsilon) is too high
    probUpdatedImagesFormerIteration = 0.0;
    for (int i = 0; i < numImagesForUpdating; i++) {
      int index = fomerUpdatedImages[i];
      if ((index < 0) || (index >= numTrainingSamples)) {
        ERR << "unexpected index for improvement calculation: " << index << endl;
        exit(1);
      }
      probUpdatedImagesFormerIteration += save_classification_order[index].p;
    }

    sort(save_classification_order.begin(), save_classification_order.end(), sortClassifications());

    probUpdatedImages = 0.0;
    for (int i = 0; i < numImagesForUpdating; i++) {
      probUpdatedImages += save_classification_order[i].p;
    }
    probUpdatedImagesFormerIteration /= double(numImagesForUpdating);
    probUpdatedImages /= double(numImagesForUpdating);
  }
  for (int n = 0; n < numImagesForUpdating; n++) {
    int index = save_classification_order[n].nr;
    if ((index < 0) || (index >= numTrainingSamples)) {
      ERR << "unexpected index for updating: " << index << endl;
      exit(1);
    }
    sum += updates[index];
  }

  if (isFirstWeight) {
    DBG(10) << "current error rate: " << (double(error) / double(ok + error)) << endl;
    double overall_average_correct_class_probability = 0.0;
    for (int k = 0; k < numClasses; k++) {
      DBG(10) << "avg. correct class probability for class " << k << ": " << (average_correct_class_probability[k] / double(class_items[k])) << endl;
      overall_average_correct_class_probability += (average_correct_class_probability[k] / double(class_items[k]));
    }
    overall_average_correct_class_probability /= double(numClasses);
    DBG(10) << "overall avg. correct class probability: " << overall_average_correct_class_probability << endl;
    DBG(15) << "using " << numImagesForUpdating << " images for update of clusterweights" << endl;
    DBG(15) << "histogram: " << endl;
    for (int b = 0; b < 20; b++) {
      DBG(15) << (double(b) / 20.0) << "-" << (double(b + 1) / 20.0) << "%: " << histogram[b] << endl;
    }
  }
  return sum;
}

double calcFderivation_productrule(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, 
                                   const vector< pair<LocalFeatures, uint> >& localfeatures, const vector< vector< vector<int> > >& nearestClusterSaved,
                                   const vector<double>& aPriori, const EPOptions& options, 
                                   int classificationMode, int k, int c, vector< vector< vector<double> > >& saved_k_given_x,
                                   vector< vector< vector<double> > >& saved_x_given_k, vector< vector< vector<double> > >& saved_weighted_emission_probability,
                                   vector<ClassifiedImage>& save_classification_order, double& probUpdatedImagesFormerIteration,
                                   double& probUpdatedImages, int numImagesForUpdating, bool isFirstWeight) {

  int numTrainingSamples = (int) localfeatures.size();
  int numClasses = (int) aPriori.size();

  int ok = 0, error = 0;
  double sum = 0.0;
  vector<double> average_correct_class_probability(numClasses, 0.0);
  vector<int> class_items(numClasses, 0);
  vector<double> histogram(20, 0);
  vector<double> updates(numTrainingSamples);
  
  vector<int> fomerUpdatedImages(numImagesForUpdating);
  if (isFirstWeight) {
    for (int i = 0; i < numImagesForUpdating; i++) {
      fomerUpdatedImages[i] = save_classification_order[i].nr;
    }
  }
  
  for (int n = 0; n < numTrainingSamples; n++) {
    int k_n = (int) localfeatures[n].second;
    int numLocalFeatures = localfeatures[n].first.numberOfFeatures();
    vector<double> correct_class_probability(numClasses, 0.0);
    for (int l = 0; l < numLocalFeatures; l++) {
      if (isFirstWeight) {
        pair<double, double> position = pair<double, double>(double(localfeatures[n].first.position(l).x) / double(localfeatures[n].first.imageSizeX()),
                                                             double(localfeatures[n].first.position(l).y) / double(localfeatures[n].first.imageSizeY()));
        calc_k_given_x(clusters, clusterNormalizationValues, aPriori, localfeatures[n].first, l, nearestClusterSaved[n],
                       position, options, classificationMode, saved_x_given_k[n][l], 
                       saved_weighted_emission_probability[n][l], saved_k_given_x[n][l]);
        for (int k_ = 0; k_ < numClasses; k_++) {
          correct_class_probability[k_] += log(saved_k_given_x[n][l][k_]);
        }
      }
    }

    if (isFirstWeight) {
      double ccp_sum = 0.0;
      for (int k_ = 0; k_ < numClasses; k_++) {
        correct_class_probability[k_] /= double(numLocalFeatures);
        correct_class_probability[k_] = exp(correct_class_probability[k_]);
        ccp_sum += correct_class_probability[k_];
      }
      correct_class_probability[k_n] /= ccp_sum;
      average_correct_class_probability[k_n] += correct_class_probability[k_n];
      class_items[k_n]++;
      if (correct_class_probability[k_n] >= 0.5) {
        ok++;
      } else {
        error++;
        DBG(20) << "p(" << k << "|{x_1_l}_" << n << ")=" << correct_class_probability[k_n] << endl;
      }
      int bin = (int) (correct_class_probability[k_n] * 20);
      if (bin > 19) bin = 19;
      if (bin < 0) exit(1);
      histogram[bin]++;
    }

    if (isFirstWeight) {
      save_classification_order[n].nr = n;
      save_classification_order[n].p = correct_class_probability[k_n];
    }

    double inner_sum = 0.0;
    for (int l = 0; l < numLocalFeatures; l++) {
      double u1;
      if (k == k_n) {
        u1 = aPriori[k_n] * saved_weighted_emission_probability[n][l][c];
      } else {
        u1 = 0.0;
      }
      double v1 = saved_x_given_k[n][l][k_n];
      double u2 = aPriori[k] * saved_weighted_emission_probability[n][l][c];
      double v2 = 0.0;
      for (int k_ = 0; k_ < numClasses; k_++) {
        v2 += aPriori[k_] * saved_x_given_k[n][l][k_];
      }
      inner_sum  += (u1 / v1) - (u2 / v2);
    }
    updates[n] = inner_sum / double(numLocalFeatures);
  }

  if (isFirstWeight) {
    // calculate the new correct class probability for those images
    // which were formerly (i.e. in the last iteration) classified worst
    // this probability must have been increased, otherwise the step size (epsilon) is too high
    probUpdatedImagesFormerIteration = 0.0;
    for (int i = 0; i < numImagesForUpdating; i++) {
      int index = fomerUpdatedImages[i];
      if ((index < 0) || (index >= numTrainingSamples)) {
        ERR << "unexpected index for improvement calculation: " << index << endl;
        exit(1);
      }
      probUpdatedImagesFormerIteration += save_classification_order[index].p;
    }
    
    sort(save_classification_order.begin(), save_classification_order.end(), sortClassifications());
    
    probUpdatedImages = 0.0;
    for (int i = 0; i < numImagesForUpdating; i++) {
      probUpdatedImages += save_classification_order[i].p;
    }
    probUpdatedImagesFormerIteration /= double(numImagesForUpdating);
    probUpdatedImages /= double(numImagesForUpdating);
  }
  for (int n = 0; n < numImagesForUpdating; n++) {
    int index = save_classification_order[n].nr;
    if ((index < 0) || (index >= numTrainingSamples)) {
      ERR << "unexpected index for updating: " << index << endl;
      exit(1);
    }
    sum += updates[index];
  }

  if (isFirstWeight) {
    DBG(10) << "current error rate: " << (double(error) / double(ok + error)) << endl;
    double overall_average_correct_class_probability = 0.0;
    for (int k = 0; k < numClasses; k++) {
      DBG(10) << "avg. correct class probability for class " << k << ": " << (average_correct_class_probability[k] / double(class_items[k])) << endl;
      overall_average_correct_class_probability += (average_correct_class_probability[k] / double(class_items[k]));
    }
    overall_average_correct_class_probability /= double(numClasses);
    DBG(10) << "overall avg. correct class probability: " << overall_average_correct_class_probability << endl;
    DBG(15) << "using " << numImagesForUpdating << " images for update of clusterweights" << endl;
    DBG(15) << "histogram: " << endl;
    for (int b = 0; b < 20; b++) {
      DBG(15) << (double(b) / 20.0) << "-" << (double(b + 1) / 20.0) << "%: " << histogram[b] << endl;
    }
  }
  return sum;
}

void calculateTestER_sumrule(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, 
                             const vector< pair<LocalFeatures, uint> >& localfeatures_test, const vector< vector< vector<int> > >& nearestClusterSaved,
                             const vector<double>& aPriori, const EPOptions& options, 
                             int classificationMode) {

  int numTestSamples = (int) localfeatures_test.size();
  int numClasses = (int) aPriori.size();

  DBG(10) << "calculating error rate for " << numTestSamples << " test samples..." << endl;
  double average_correct_class_probability = 0.0;
  vector<double> empty_x_given_k;
  vector<double> empty_weighted_emission_probability;
  int ok = 0, error = 0;

  for (int n = 0; n < numTestSamples; n++) {
    int k_n = (int) localfeatures_test[n].second;
    int numLocalFeatures = localfeatures_test[n].first.numberOfFeatures();
    double correct_class_probability = 0.0;

    for (int l = 0; l < numLocalFeatures; l++) {
      vector<double> p_k_given_x(numClasses, 0.0);
      pair<double, double> position = pair<double, double>(double(localfeatures_test[n].first.position(l).x) / double(localfeatures_test[n].first.imageSizeX()),
                                                           double(localfeatures_test[n].first.position(l).y) / double(localfeatures_test[n].first.imageSizeY()));
      calc_k_given_x(clusters, clusterNormalizationValues, aPriori, localfeatures_test[n].first, l, 
                     nearestClusterSaved[n], position, options, classificationMode,
                     empty_x_given_k, empty_weighted_emission_probability, p_k_given_x);
      correct_class_probability += p_k_given_x[k_n];
    }

    correct_class_probability /= double(numLocalFeatures);
    average_correct_class_probability += correct_class_probability;
    if (correct_class_probability >= 0.5) {
      ok++;
    } else {
      error++;
    }
  }
  average_correct_class_probability /= double(numTestSamples);
  DBG(15) << "ER is " << (double(error) / double(ok + error)) << ", avg. correct class probability: " << average_correct_class_probability << endl;
}

void calculateTestER_productrule(const vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, 
                                 const vector< pair<LocalFeatures, uint> >& localfeatures_test, const vector< vector< vector<int> > >& nearestClusterSaved,
                                 const vector<double>& aPriori, const EPOptions& options, int classificationMode) {

  int numTestSamples = (int) localfeatures_test.size();
  int numClasses = (int) aPriori.size();

  DBG(10) << "calculating error rate for " << numTestSamples << " test samples..." << endl;
  double average_correct_class_probability = 0.0;
  vector<double> empty_x_given_k;
  vector<double> empty_weighted_emission_probability;
  int ok = 0, error = 0;

  for (int n = 0; n < numTestSamples; n++) {
    int k_n = (int) localfeatures_test[n].second;
    int numLocalFeatures = localfeatures_test[n].first.numberOfFeatures();
    vector<double> correct_class_probability(numClasses, 0.0);

    for (int l = 0; l < numLocalFeatures; l++) {
      vector<double> p_k_given_x(numClasses, 0.0);
      pair<double, double> position = pair<double, double>(double(localfeatures_test[n].first.position(l).x) / double(localfeatures_test[n].first.imageSizeX()),
                                                           double(localfeatures_test[n].first.position(l).y) / double(localfeatures_test[n].first.imageSizeY()));
      calc_k_given_x(clusters, clusterNormalizationValues, aPriori, localfeatures_test[n].first, l, 
                     nearestClusterSaved[n], position, options, classificationMode,
                     empty_x_given_k, empty_weighted_emission_probability, p_k_given_x);
      for (int k = 0; k < numClasses; k++) {
        correct_class_probability[k] += log(p_k_given_x[k]);
      }
    }

    double norm = 0.0;
    for (int k = 0; k < numClasses; k++) {
      correct_class_probability[k] = exp(correct_class_probability[k] / double(numLocalFeatures));
      norm += correct_class_probability[k];
    }
    for (int k = 0; k < numClasses; k++) {
      correct_class_probability[k] /= norm;
    }

    average_correct_class_probability += correct_class_probability[k_n];
    if (correct_class_probability[k_n] >= 0.5) {
      ok++;
    } else {
      error++;
    }
  }
  average_correct_class_probability /= double(numTestSamples);
  DBG(15) << "ER is " << (double(error) / double(ok + error)) << ", avg. correct class probability: " << average_correct_class_probability << endl;
}

void trainWeights(vector<ClusterData>& clusters, const vector<double>& clusterNormalizationValues, 
                  const vector< pair<LocalFeatures, uint> >& localfeatures, const vector<double>& aPriori, const EPOptions& options,  
                  const vector< pair<LocalFeatures, uint> >& localfeatures_test, bool productRule, int classificationMode, int numIterations,
                  double epsilon, double percentWorstImages, const string outputfile, const int eachIter) {

  int numClusters = (int) clusters.size();
  int numClasses = (int) aPriori.size();
  int numImages = (int) localfeatures.size();
  vector< vector<double> > saved_k_given_X(localfeatures.size());
  vector< vector< vector<double> > > saved_k_given_x(localfeatures.size());
  vector< vector< vector<double> > > saved_x_given_k(localfeatures.size());
  vector< vector< vector<double> > > saved_weighted_emission_probability(localfeatures.size());
  for (int n = 0; n < (int) localfeatures.size(); n++) {
    saved_k_given_X[n].resize(numClasses);
    saved_k_given_x[n].resize(localfeatures[n].first.numberOfFeatures());
    saved_x_given_k[n].resize(localfeatures[n].first.numberOfFeatures());
    saved_weighted_emission_probability[n].resize(localfeatures[n].first.numberOfFeatures());
    for (int l = 0; l < (int) localfeatures[n].first.numberOfFeatures(); l++) {
      saved_k_given_x[n][l].resize(numClasses);
      saved_x_given_k[n][l].resize(numClasses);
      saved_weighted_emission_probability[n][l].resize(numClusters);
    }
  }

  vector< vector< vector<int> > > nearestClusterSaved(localfeatures.size());
  vector< vector< vector<int> > > nearestClusterSaved_test(localfeatures_test.size());
  if ((classificationMode & CL_RELATIVE_POSITION) > 0) {
    DBG(15) << "calculating nearest clusters for all training patches" << endl;
    for (uint n = 0; n < localfeatures.size(); n++) {
      nearestClusterSaved[n].resize(localfeatures[n].first.numberOfFeatures(), vector<int>(numClasses));
      for (uint l = 0; l < localfeatures[n].first.numberOfFeatures(); l++) {
        for (int k = 0; k < numClasses; k++) {
          nearestClusterSaved[n][l][k] = getNearestCluster(localfeatures[n].first[l], clusters, clusterNormalizationValues, k, options);
        }
      }
    }

    if (!localfeatures_test.empty()) {
      DBG(15) << "calculating nearest clusters for all test patches" << endl;
      for (uint n = 0; n < localfeatures_test.size(); n++) {
        nearestClusterSaved_test[n].resize(localfeatures_test[n].first.numberOfFeatures(), vector<int>(numClasses));
        for (uint l = 0; l < localfeatures_test[n].first.numberOfFeatures(); l++) {
          for (int k = 0; k < numClasses; k++) {
            nearestClusterSaved_test[n][l][k] = getNearestCluster(localfeatures_test[n].first[l], clusters, clusterNormalizationValues, k, options);
          }
        }
      }
    }
  }

  vector<double> updates(numClusters);
  vector<ClassifiedImage> save_classification_order(localfeatures.size());
  double probUpdatedImagesInThisIteration = 0.0;
  double probUpdatedImagesInPreviousIteration = 0.0;
  double saveProbUpdatedImagesInPreviousIteration = 0.0;
  int stopTraining = 0;

  if (!localfeatures_test.empty()) {
    DBG(15) << "calculating initial error rate on testing data" << endl;
    if (productRule) {
      calculateTestER_productrule(clusters, clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
    } else {
      calculateTestER_sumrule(clusters, clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
    }
  }
  
  for (int iteration = 0; iteration < numIterations; iteration++) {
    
    BLINK(15) << endl;
    DBG(15) <<  "starting iteration " << (iteration + 1) << endl;
    
    bool firstCalculation = true;
    
    for (int c = 0; c < numClusters; c++) {
      DBG(20) << "current weight for cluster " << c << " of class " << clusters[c].clazz << " is: " << clusters[c].weight << endl;
      double f_deriv;
      int numImagesForUpdate = (int) round((double) numImages * percentWorstImages / 100.0);
      if (productRule) {
        f_deriv = calcFderivation_productrule(clusters, clusterNormalizationValues, localfeatures, nearestClusterSaved, aPriori, options, classificationMode, 
                                              clusters[c].clazz, c, saved_k_given_x, saved_x_given_k, saved_weighted_emission_probability, 
                                              save_classification_order, probUpdatedImagesInPreviousIteration, probUpdatedImagesInThisIteration, 
                                              numImagesForUpdate, firstCalculation);
      } else {
        f_deriv = calcFderivation_sumrule(clusters, clusterNormalizationValues, localfeatures, nearestClusterSaved, aPriori, options, classificationMode, 
                                          clusters[c].clazz, c, saved_k_given_x, saved_x_given_k, saved_weighted_emission_probability, 
                                          save_classification_order,  probUpdatedImagesInPreviousIteration, probUpdatedImagesInThisIteration, 
                                          numImagesForUpdate, firstCalculation);
      }
      firstCalculation = false;
      updates[c] = f_deriv;
      DBG(20) << "update cluster " << c << " class " << clusters[c].clazz << " with: " << updates[c] << endl;
    }

    if (iteration > 0) {
      DBG(15) << "prob. of images updated in last iteration is now " << probUpdatedImagesInPreviousIteration << endl;
      DBG(15) << "was " << saveProbUpdatedImagesInPreviousIteration << " before" << endl;
      if (probUpdatedImagesInPreviousIteration < saveProbUpdatedImagesInPreviousIteration) {
        epsilon /= 2.0;
        DBG(10) << "result got worse, decreasing epsilon to " << epsilon << endl;
      }
    }
    DBG(15) << "prob. of images updated in this iteration " << probUpdatedImagesInThisIteration << endl;

    if (probUpdatedImagesInThisIteration < saveProbUpdatedImagesInPreviousIteration) {
      stopTraining++;
      if (stopTraining >= 3) {
        DBG(10) << "training is finished!" << endl;
        if (productRule) {
          calculateTestER_productrule(clusters, clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
        } else {
          calculateTestER_sumrule(clusters, clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
        }
        if (outputfile != "") {
          saveClusterWeights(clusters, outputfile);
        }
        break;
      }
    }

    saveProbUpdatedImagesInPreviousIteration = probUpdatedImagesInThisIteration;

    vector<double> clusterWeightSum(numClasses, 0.0);
    int negative = 0;
    double maxRelativeUpdate = 0.0, maxAbsoluteUpdate = 0.0;
    for (int c = 0; c < numClusters; c++) {
      updates[c] *= epsilon;
      if (fabs(updates[c]) / clusters[c].weight > maxRelativeUpdate) {
        maxRelativeUpdate = fabs(updates[c]) / clusters[c].weight;
      }
      if (fabs(updates[c]) > maxAbsoluteUpdate) {
        maxAbsoluteUpdate = fabs(updates[c]);
      }
      clusters[c].weight += updates[c];
      if (clusters[c].weight < 0.0) {
        clusters[c].weight = 0.0;
        negative++;
      }
      clusterWeightSum[clusters[c].clazz] += clusters[c].weight;
    }
    DBG(15) << negative << " negative clusterweights were set to 0.0" << endl;
    DBG(15) << "max. rel. update: " << maxRelativeUpdate << endl;
    DBG(15) << "max. abs. update: " << maxAbsoluteUpdate << endl;
    for (int k = 0; k < numClasses; k++) {
      DBG(15) << "sum of all clusterweights of class " << k << " is: " << clusterWeightSum[k] << endl;
    }

    if (((iteration % 10) == 9) && (!localfeatures_test.empty())) {
      if (productRule) {
        calculateTestER_productrule(clusters,  clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
      } else {
        calculateTestER_sumrule(clusters, clusterNormalizationValues, localfeatures_test, nearestClusterSaved_test, aPriori, options, classificationMode);
      }
      if (outputfile != "") {
        saveClusterWeights(clusters, outputfile);
      }
    }
    
    if( (iteration+1) % eachIter == 0 ) {
      ostringstream sname;
      sname << outputfile << "-" << iteration << ".clusterweights";
      saveClusterWeights(clusters,sname.str());
    }


  }
}

void calculateClusterPositions(const Database& trainingData, vector<ClusterData>& clusters, int maxSplitsPerCluster, 
                               bool positionDifferences, int maxSplitsPerClusterDifferences, const vector<double>& clusterNormalizationValues, EPOptions& options) {
  PositionClusterer pc;
  pc.maxSplit() = maxSplitsPerCluster;
  pc.iterationsBetweenSplits() = 10;
  pc.stopWithNClusters() = 8;
  pc.minObservationsPerCluster() = 5;
  pc.dontSplitBelow() = 8;
  pc.epsilon() = 0.01;
  pc.diagonalize() = true;
  pc.minVariance()=0.01;
  PositionClusterer pcDifferences;
  pcDifferences.maxSplit() = maxSplitsPerClusterDifferences;
  pcDifferences.iterationsBetweenSplits() = 10;
  pcDifferences.stopWithNClusters() = 8;
  pcDifferences.minObservationsPerCluster() = 5;
  pcDifferences.dontSplitBelow() = 8;
  pcDifferences.epsilon() = 0.01;
  pcDifferences.diagonalize() = false;

  DBG(15) << "estimating cluster positions (max " << (int) pow(2.0, (double) maxSplitsPerCluster) << ")" << endl;
  if (positionDifferences) {
    DBG(15) << "estimating cluster position differences (max " << (int) pow(2.0, (double) maxSplitsPerClusterDifferences) << ")" << endl;
  }
  int numFiles = (int) trainingData.size();

  vector< vector< vector<double> > > inputdataDifferences(clusters.size(), vector< vector<double> >(clusters.size()));

  vector< vector<DoubleVector> > inputdata(clusters.size());
  vector< vector< uint > > numPositionDifferences(clusters.size(), vector<uint>(clusters.size(), 0));

  for (int n = 0; n < numFiles; n++) {
    int k = trainingData[n]->clas();
    string fname = trainingData.path()+"/"+trainingData.filename(n) + "." + trainingData.suffix(0);
    LocalFeatures localFeatures;
    localFeatures.load(fname);
    double imgSizeX = double(localFeatures.imageSizeX());
    double imgSizeY = double(localFeatures.imageSizeY());
    DBG(20) << "processing " << localFeatures.numberOfFeatures() << " local features from file " << fname << endl;
    vector<int> nearestClusterNum(localFeatures.numberOfFeatures());
    for (int l = 0; l < (int) localFeatures.numberOfFeatures(); l++) {
      nearestClusterNum[l] = getNearestCluster(localFeatures[l], clusters, clusterNormalizationValues, k, options);
      if ((nearestClusterNum[l] < 0) || (nearestClusterNum[l] > (int) clusters.size())) {
        ERR << "nearest cluster should not be " << nearestClusterNum[l] << endl;
        exit(1);
      }
      DoubleVector nextPositionItem(2);
      nextPositionItem[0] = double(localFeatures.position(l).x) / imgSizeX;
      nextPositionItem[1] = double(localFeatures.position(l).y) / imgSizeY;
      inputdata[nearestClusterNum[l]].push_back(nextPositionItem);
    }
    if (positionDifferences) {
      for (int l = 0; l < (int) localFeatures.numberOfFeatures(); l++) {
        for (int l2 = 0; l2 < (int) localFeatures.numberOfFeatures(); l2++) {
          if (l != l2) {
            if (inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].size() <= numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]]) {
              inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].push_back(((double(localFeatures.position(l2).x) / imgSizeX) - (double(localFeatures.position(l).x) / imgSizeX)));
              inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].push_back(((double(localFeatures.position(l2).y) / imgSizeY) - (double(localFeatures.position(l).y) / imgSizeY)));
            } else {
              uint position = numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]];
              inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]][position] = ((double(localFeatures.position(l2).x) / imgSizeX) - (double(localFeatures.position(l).x) / imgSizeX));
              inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]][position + 1] = ((double(localFeatures.position(l2).y) / imgSizeY) - (double(localFeatures.position(l).y) / imgSizeY));
            }
            numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]] += 2;
          }
        }
      }
    }
  }

  int totalrelativepositions = 0;
  for (uint c = 0; c < clusters.size(); c++) {
    vector<ClusterPosition> cl;
    DBG(15) << "  cluster " << c << " (" << inputdata[c].size() << " positions)" << endl;
    if(inputdata[c].size()==0) {
      inputdata[c].push_back(DoubleVector(2,1.0));
    }
    
    pc.run(inputdata[c], cl);
    clusters[c].positions.resize(cl.size());
    DBG(20) << "got " << cl.size() << " position clusters" << endl;
    for (int p = 0; p < (int) cl.size(); p++) {
      clusters[c].positions[p] = cl[p];
      DBG(25) << " " << p << "(" << cl[p].elements << "): mean=" << cl[p].mean.first << "," << cl[p].mean.second << endl;
      DBG(25) << "    covariance=" << cl[p].covariance.xx << "," << cl[p].covariance.xy << "," << cl[p].covariance.yx << "," << cl[p].covariance.yy << endl;
    }
    if (positionDifferences) {
      clusters[c].positionDifferences.resize(clusters.size());
      clusters[c].hasPositionDifference.resize(clusters.size());
      int relativeclusters = 0;
      for (uint c2 = 0; c2 < clusters.size(); c2++) {
        totalrelativepositions += (numPositionDifferences[c][c2] / 2);
        if (numPositionDifferences[c][c2] / 2 >= 10) {
          vector<ClusterPosition> result;
          DBG(20) << "clustering position differences for clusters " << c << " and " << c2 << " (" << (inputdataDifferences[c][c2].size() / 2) << " position differences)" << endl;
          vector< vector<double> > differences(numPositionDifferences[c][c2] / 2, vector<double>(2));
          int count = 0;
          for (uint i = 0; i < numPositionDifferences[c][c2]; i++) {
            differences[count][0] = inputdataDifferences[c][c2][i];
            i++;
            differences[count][1] = inputdataDifferences[c][c2][i];
            count++;
          }
          pcDifferences.run(differences, result);
          clusters[c].positionDifferences[c2].resize(result.size());
          DBG(20) << "got " << result.size() << " position difference clusters" << endl;
          for (int p = 0; p < (int) result.size(); p++) {
            DBG(25) << " " << p << "(" << result[p].elements << "): mean=" << result[p].mean.first << "," << result[p].mean.second << endl;
            DBG(25) << "    covariance=" << result[p].covariance.xx << "," << result[p].covariance.xy << "," << result[p].covariance.yx << "," << result[p].covariance.yy << endl;
            clusters[c].positionDifferences[c2][p] = result[p];
            if (result[p].covariance.getDeterminant() < 0.0) {
              ERR << "determinant is < 0 !" << endl;
              ERR << "(" << result[p].elements << "): mean=" << result[p].mean.first << "," << result[p].mean.second << endl;
              ERR << "      covariance=" << result[p].covariance.xx << "," << result[p].covariance.xy << "," << result[p].covariance.yx << "," << result[p].covariance.yy << endl;
              for (uint i = 0; i < differences.size(); i++) {
                ERR << differences[i][0] << "," << differences[i][1] << endl;
              }
              exit(1);
            }
          }
          clusters[c].hasPositionDifference[c2] = true;
          relativeclusters++;
        } else {
          DBG(20) << "not clustering position differences for clusters " << c << " and " << c2 <<
            ", only (" << inputdataDifferences[c][c2].size() << " position differences)" << endl;
          clusters[c].hasPositionDifference[c2] = false;
        }
      }
      DBG(15) << "  relative position clusters to " << relativeclusters << " clusters" << endl;
    }
  }
  if (positionDifferences) {
    DBG(15) << "  collected " << totalrelativepositions << " relative positions in total" << endl;
  }
}


int main(int argc, char** argv) {

  GetPot cl(argc, argv);

  if (cl.search(2, "-h", "--help")) {
    USAGE();
  }

  double epsilon = 0.001;
  if (cl.search(2, "-e", "--epsilon")) {
    epsilon = cl.next(0.001);
  }
  DBG(10) << "updating weights with epsilon=" << epsilon << endl;
  double percentWorstImages = 100.0;
  if (cl.search(2, "-w", "--worstimages")) {
    percentWorstImages = cl.next(100.0);
    if ((percentWorstImages <= 0.0) || (percentWorstImages > 100.0)) {
      USAGE();
    }
  }
  DBG(10) << "using " << percentWorstImages << "% of all images for update in each iteration";
  if (percentWorstImages < 100.0) {
    BLINK(10) << " (the worst ones)" << endl;
  } else {
    BLINK(10) << endl;
  }

  int numIterations = 100;
  if (cl.search(2, "-i", "--iterations")) {
    numIterations = cl.next(100);
  }
  DBG(10) << "iterate " << numIterations << " times" << endl;

  string outputfile = "";
  if (cl.search(2, "-s", "--savefile")) {
    outputfile = cl.next("");
    DBG(15) << "saving clusterweights to " << outputfile << endl;
  }
  string inputfile = "";
  if (cl.search(2, "-l", "--loadfile")) {
    inputfile = cl.next("");
    DBG(15) << "loading clusterweights from " << inputfile << endl;
  }

  Database modelDatabase;
  if (!cl.search(1, "--models")) {
    USAGE();
  } else {
    modelDatabase.loadFileList(cl.next(" "));
    if (modelDatabase.numberOfSuffices() > 1) {
      ERR << "found " << modelDatabase.numberOfSuffices() << " suffices in model database, only 1 is supported." << endl;
      exit(1);
    }
  }

  bool productRule = cl.search(2, "-p", "--productrule");
  if (productRule) {
    DBG(10) << "using product rule in classification and training" << endl;
  } else {
    DBG(10) << "using sum rule in classification and training" << endl;
  }

  int numClasses = getNumClasses(modelDatabase);
  vector<double> aPriori(numClasses);

  // load further options
  EPOptions options;
  int classificationMode = 0;
  if (cl.search(1, "appearance")) {
    classificationMode += CL_APPEARANCE;
  }
  int numPositionClustersPerCluster = 0, numPositionDifferencesClustersPerCluster = 0;
  if (cl.search(1, "position")) {
    classificationMode += CL_POSITION;
    numPositionClustersPerCluster = cl.next(0);
    if ((numPositionDifferencesClustersPerCluster < 0) || (numPositionDifferencesClustersPerCluster > 4)) {
      USAGE();
    }
    DBG(10) << "including positions in classification" << endl;
  }
  if (cl.search(1, "relativeposition")) {
    classificationMode += CL_RELATIVE_POSITION;
    numPositionDifferencesClustersPerCluster = cl.next(0);
    if ((numPositionDifferencesClustersPerCluster < 0) || (numPositionDifferencesClustersPerCluster > 4)) {
      USAGE();
    }
    DBG(10) << "including relative positions in classification" << endl;
  }
  if (classificationMode == 0) {
    USAGE();
  }

  // load models
  DBG(10) << "loading models..." << endl;
  vector<ClusterData> clusters;
  for (int m = 0; m < (int) modelDatabase.size(); m++) {
    string fname = modelDatabase.filename(m) + "." + modelDatabase.suffix(0);
    loadModel(modelDatabase.path()+"/"+fname, modelDatabase[m]->clas(), clusters);
  }
  DBG(15) << "loaded " << clusters.size() << " clustermodels for " << numClasses << " classes" << endl;
  
  // read local feature training files
  vector< pair<LocalFeatures, uint> > localfeatures;
  Database trainDb;
  if (!cl.search(1, "--train")) {
    USAGE();
  } else {
    trainDb.loadFileList(cl.next(" "));
    if (trainDb.numberOfSuffices() > 1) {
      ERR << "found " << trainDb.numberOfSuffices() << " suffices in training database, only 1 is supported." << endl;
      exit(1);
    }
    int totalFeatures = 0;
    for (int n = 0; n < (int) trainDb.size(); n++) {
      LocalFeatures lf;
      string fname = trainDb.path()+"/"+trainDb.filename(n) + "." + trainDb.suffix(0);
      lf.load(fname);
      totalFeatures += lf.numberOfFeatures();
      localfeatures.push_back(pair<LocalFeatures, uint>(lf, (uint) trainDb[n]->clas()));
    }
    DBG(15) << "read " << trainDb.size() << " training files with " << totalFeatures << " local features in total" << endl;
  }

  // read local feature training files
  vector< pair<LocalFeatures, uint> > localfeatures_test;
  Database testDb;
  if (cl.search(1, "--test")) {
    testDb.loadFileList(cl.next(" "));
    if (testDb.numberOfSuffices() > 1) {
      ERR << "found " << testDb.numberOfSuffices() << " suffices in testing database, only 1 is supported." << endl;
      exit(1);
    }
    int totalFeatures = 0;
    for (int n = 0; n < (int) testDb.size(); n++) {
      LocalFeatures lf;
      string fname = testDb.filename(n) + "." + testDb.suffix(0);
      lf.load(fname);
      totalFeatures += lf.numberOfFeatures();
      localfeatures_test.push_back(pair<LocalFeatures, uint>(lf, (uint) testDb[n]->clas()));
    }
    DBG(15) << "read " << testDb.size() << " test files with " << totalFeatures << " local features in total" << endl;
  }

  // calculate a priori probabilities for all classes
  if (cl.search(2, "-eq", "--eqapriori")) {
    for (int k = 0; k < numClasses; k++) {
      aPriori[k] = 1.0 / double(numClasses);
      DBG(15) << "a priori probability for class " << k << ": " << aPriori[k] << endl;
    }
  } else {
    for (int k = 0; k < numClasses; k++) {
      int filesK = 0;
      for (int n = 0; n < (int) trainDb.size(); n++) {
        if ((int) trainDb[n]->clas() == k) {
          filesK++;
        }
      }
      aPriori[k] = double(filesK) / double(trainDb.size());
      DBG(15) << "a priori probability for class " << k << ": " << aPriori[k] << endl;
    }
  }
  
  // determine if there are dimensions which are to be discarded for classification
  map<int, bool> dims2discard;
  int numDims2discard = 0;
  if (cl.search(2, "-dd", "--discardDimensions")) {
    string dim = cl.next(" ");;
    while ((dim != " ") && (dim.c_str()[0] != '-')) {
      int dim2discard = atoi(dim.c_str());
      DBG(10) << "discarding dimension " << dim2discard << " in classification" << endl;
      dims2discard[dim2discard] = true;
      dim = cl.next(" ");
      numDims2discard++;
    }
  }
  vector<bool> useDim(clusters[0].density.dim, true);
  DBG(10) << "local features have " << clusters[0].density.dim << " dimensions, ";
  if (numDims2discard == 0) {
    BLINK(10) << "using all of them" << endl;
  } else {
    BLINK(10) << "discarding the following:" << endl;
    for (int d = 0; d < (int) clusters[0].density.dim; d++) {
      if (dims2discard.find(d) != dims2discard.end()) {
      	DBG(10) << " dimension " << d << endl;
        useDim[d] = false;
      }
    }
  }
  options.useDim = useDim;

  // if the variances are pooled or the variance is in another way manipulated, this is done here:
  double varianceValue = 0.0, varianceFactor = 0.0;
  if (cl.search(2, "-var", "--variance")) {
    varianceValue = cl.next(0.0);
  }
  if (cl.search(2, "-mvar", "--multiplyvariance")) {
    varianceFactor = cl.next(0.0);
  }
  adaptClusterVariances(clusters, varianceValue, varianceFactor, options.square, options.dimensionPooling);
  // pre-calculate emission probability normalization values for the clusters based on their variances
  // this speeds calculation up, if clusters are not dimension-pooled
  vector<double> clusterNormalizationValues;
  calcClusterNormalizationValues(clusters, options, clusterNormalizationValues);

  if (inputfile == "") {
    calcInitialClusterweights(clusters, numClasses);
  } else {
    loadInitialClusterweights(clusters, inputfile);
  }

  if (((classificationMode & CL_POSITION) > 0) || ((classificationMode & CL_RELATIVE_POSITION) > 0)) {
    calculateClusterPositions(trainDb, clusters, numPositionClustersPerCluster, (classificationMode & CL_RELATIVE_POSITION) > 0,
                              numPositionDifferencesClustersPerCluster, clusterNormalizationValues, options);
  }

  int eachIter=cl.follow(1,2,"-eI","--eachIter");

  trainWeights(clusters, clusterNormalizationValues, localfeatures, aPriori, options, localfeatures_test, productRule, classificationMode, 
               numIterations, epsilon, percentWorstImages, outputfile,eachIter);
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);

}

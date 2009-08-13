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

using namespace std;

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
  double weight;
  double normalizationValue;
  // initialize with default values
  ClusterData() {
    weight = 0.0;
  }
};

// options affecting the emission probability calculation
struct EPOptions {
  vector<bool> useDim;
};

// prints the usage
void USAGE() {
  cout << "USAGE: trainclusterweightsdiscriminative [options] --models <filelist> --train <filelists>" << endl
       << " Options: " << endl
       << "  -e, --epsilon <epsilon>            : update parameter for iteration" << endl
       << "  -eq, --eqapriori                   : equal a priori probabilities" << endl
       << "  -h, --help                         : shows this help" << endl
       << "  -i, --iterations <numiter>         : iterate the given number of times" << endl
       << "  -l, --loadfile <filename>          : load clusterweights from given file" << endl
       << "  -s, --savemodels <filelist>        : save modified clustermodels to the files in the filelist" << endl
       << "  -w, --worstimages <n>              : use only the worst n% (0.0 < n <= 100.0) images for updating" << endl
       << "                                       clusterweights in each iteration" << endl
       << "  --test <filelists>                 : for classification test data" << endl;
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

// saves a cluster model in the file with the given filename
void saveModel(const string filename, int classNr, const vector<ClusterData>& clusters) {
  ERR << "saveModel(...) not yet implemented!" << endl;
  exit(1);
}

// returns the number of classes found in the given database
// the method does not look for distinct class labels, but
// returns the highest class label + 1
// therefore, classes are assumed to be labeled starting from
// 0 to n-1
uint getNumClasses(const Database& db) {
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
  if (maxClassNr < 0) {
    ERR << "unexpected: 'maxClassNr' in 'getNumClasses(...)' is " << maxClassNr << endl;
    exit(1);
  }
  return (uint) (maxClassNr + 1);
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

// precalculates the normalization values for the cluster emission probabilities
// this is possible since the normalization values are independent from the local features
// in case of no dimension pooling, this helps speedup the emission probability calculation
void calcClusterNormalizationValues(vector<ClusterData>& clusters, const EPOptions& options) {
  for (int c = 0; c < (int) clusters.size(); c++) {
    clusters[c].normalizationValue = 0.0;
    for (int d = 0; d < (int) clusters[c].density.dim; d++) {
      if (options.useDim[d]) {
        clusters[c].normalizationValue += log(2.0 * M_PI * clusters[c].density.sigma[d]);
      }
    }
  }
}


// calculates the (logarithmic) cluster emission probability p(x|c,k)
//
double calc_log_x_given_c_k(const GaussianDensity& clustermodel, double clusterNormalizationValue,
                            const vector<double>& x_l, const EPOptions& options) {
  double pixelNorm = 0.0, pixelDist = 0.0;
  // g
  if (x_l.size() != clustermodel.dim) {
    ERR << "size of local feature (" << x_l.size() << ") and of clustermodel (" << clustermodel.mean.size() << ") do not match !" << endl;
    exit(3);
  }
  int dims = (int) x_l.size();
  for (int d = 0; d < dims; d++) {
    if (options.useDim[d]) {
      double tmp = 0.0;
      tmp = x_l[d] - clustermodel.mean[d];
      tmp *= tmp;
      pixelDist += (tmp / clustermodel.sigma[d]);
    }
  }
  pixelNorm = clusterNormalizationValue;
  return -0.5 * (pixelDist + pixelNorm);
}

// calculates for a single local feature x_l the probability p(x_l|k) for each k
// and stores them in the 'p_x_given_k' vector
// if desired, the individual (weighted) emission probabilities that
// sum up to p(x|k) are stored in the given vector 'saved_weighted_emission_probability'
void calc_x_given_k(const vector<ClusterData>& clusters, const LocalFeatures& localFeatures, int l, const EPOptions& options,
                    vector<double>& p_x_given_k, vector<double>& saved_emission_probability) {

  int numClusters = (int) clusters.size();
  int numClasses = (int) p_x_given_k.size();

  for (int c = 0; c < numClusters; c++) {
    // TODO: remove this check for the sake of efficiency !
    if (clusters[c].clazz >= numClasses) {
      ERR << "cluster " << c << " belongs to class " << clusters[c].clazz << ", but p(x|k) is dimensioned only for " << numClasses << " classes!" << endl;
      exit(1);
    }
    double emission_probability = exp(calc_log_x_given_c_k(clusters[c].density, clusters[c].normalizationValue, localFeatures[l], options));
    if (!saved_emission_probability.empty()) {
      saved_emission_probability[c] = emission_probability;
    }
    p_x_given_k[clusters[c].clazz] += clusters[c].weight * emission_probability;
  }
}

// calculates for a single local feature x_l the probability p(k|x_l) for each k
// and stores them in the p_k_given_x vector
// if desired, the probabilities p(x_l|k) that are used to calculate p(k|x_l) using Bayes
// are stored in the given vector 'saved_x_given_k'
void calc_k_given_x(const vector<ClusterData>& clusters, const vector<double>& aPriori,
                    const LocalFeatures& localFeatures, int l, const EPOptions& options, vector<double>& p_k_given_x,
                    double& save_v, vector<double>& saved_emission_probability) {

  uint numClasses = aPriori.size();
  if (p_k_given_x.size() != numClasses) {
    ERR << "saved_x_given_k incorrectly dimensionized !" << endl;
    exit(1);
  }

  vector<double> p_x_given_k(numClasses, 0.0);
  calc_x_given_k(clusters, localFeatures, l, options, p_x_given_k, saved_emission_probability);
  double denominator = 0.0;

  for (uint k = 0; k < numClasses; k++) {
    double product = aPriori[k] * p_x_given_k[k];
    p_k_given_x[k] = product;
    denominator += product;
  }

  for (uint k = 0; k < numClasses; k++) {
    // TODO: remove this check for the sake of efficiency !
    if ((p_k_given_x[k] / denominator < 0.0) || (p_k_given_x[k] / denominator > 1.0) || isnan(p_k_given_x[k] / denominator)) {
      ERR << "not a valid probability value in calc_k_given_x: " << (p_k_given_x[k] / denominator) << endl;
      ERR << p_k_given_x[k] << endl;
      ERR << denominator << endl;
      exit(1);
    }
    p_k_given_x[k] /= denominator;
  }

  save_v = denominator;
}

double calcFderiv(const vector<ClusterData>& clusters, const vector< pair<LocalFeatures, uint> >& localfeatures,
                  const vector<double>& aPriori, const EPOptions& options, int k, int c, int d,
                  vector< vector<double> >& saved_u, vector< vector<double> >& saved_v,
                  vector< vector< vector<double> > >& saved_emission_probability, bool firstCalculation) {

  uint numTrainingSamples = localfeatures.size();
  uint numClasses = aPriori.size();

  int ok = 0, error = 0;
  double sum = 0.0;
  vector<double> average_correct_class_probability(numClasses, 0.0);
  vector<int> class_items(numClasses, 0);
  vector<double> histogram(20, 0);
  //vector<double> updates(numTrainingSamples);

  //vector<int> fomerUpdatedImages(numImagesForUpdating);
  //if (firstCalculation) {
  //  for (int i = 0; i < numImagesForUpdating; i++) {
  //    fomerUpdatedImages[i] = save_classification_order[i].nr;
  //  }
  //}

  for (uint n = 0; n < numTrainingSamples; n++) {
    uint k_n = localfeatures[n].second;
    uint numLocalFeatures = localfeatures[n].first.numberOfFeatures();
    double denominator = 0.0, numerator = 0.0;
    double correct_class_probability = 0.0;
    for (uint l = 0; l < numLocalFeatures; l++) {
      if (firstCalculation) {
        vector<double> p_k_given_x(numClasses, 0.0);
        calc_k_given_x(clusters, aPriori, localfeatures[n].first, l, options, p_k_given_x, saved_v[n][l], saved_emission_probability[n][l]);
        saved_u[n][l] = p_k_given_x[k_n] * saved_v[n][l];
        denominator += p_k_given_x[k_n];
      } else {
        denominator += (saved_u[n][l] / saved_v[n][l]);
      }
      double u_ = 0.0, v_ = 0.0;
      double tmp = localfeatures[n].first[l][d] - clusters[c].density.mean[d];
      tmp *= tmp;
      tmp = -tmp;
      tmp /= (clusters[c].density.sigma[d] * clusters[c].density.sigma[d]);
      tmp += (1.0 / clusters[c].density.sigma[d]);
      v_ = clusters[c].weight * saved_emission_probability[n][l][c] * (-0.5) * tmp;
      if (k == (int) k_n) {
        u_ = v_;
      }
      numerator += ((u_ * saved_v[n][l] - v_ * saved_u[n][l]) / (saved_v[n][l] * saved_v[n][l]));
    }

    correct_class_probability = denominator / double(numLocalFeatures);
    if (firstCalculation) {
      average_correct_class_probability[k_n] += correct_class_probability;
      class_items[k_n]++;
      if (correct_class_probability >= 0.5) {
        ok++;
      } else {
        error++;
      }
      int bin = (int) (correct_class_probability * 20);
      if (bin > 19) bin = 19;
      if (bin < 0) exit(1);
      histogram[bin]++;
    }

    //if (isFirstWeight) {
    //  save_classification_order[n].nr = n;
    //  save_classification_order[n].p = correct_class_probability;
    //}

    //updates[n] = numerator / denominator;
    sum += numerator / denominator;

  }

  /*
  if (firstClassification) {
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

  for (uint n = 0; n < numImagesForUpdating; n++) {
    int index = save_classification_order[n].nr;
    if ((index < 0) || (index >= numTrainingSamples)) {
      ERR << "unexpected index for updating: " << index << endl;
      exit(1);
    }
    sum += updates[index];
  }
  */

  if (firstCalculation) {
    DBG(10) << "current error rate: " << (double(error) / double(ok + error)) << endl;
    double overall_average_correct_class_probability = 0.0;
    for (uint k = 0; k < numClasses; k++) {
      DBG(10) << "avg. correct class probability for class " << k << ": " << (average_correct_class_probability[k] / double(class_items[k])) << endl;
      overall_average_correct_class_probability += (average_correct_class_probability[k] / double(class_items[k]));
    }
    overall_average_correct_class_probability /= double(numClasses);
    DBG(10) << "overall avg. correct class probability: " << overall_average_correct_class_probability << endl;
    //DBG(15) << "using " << numImagesForUpdating << " images for update of clusterweights" << endl;
    DBG(15) << "histogram: " << endl;
    for (int b = 0; b < 20; b++) {
      DBG(15) << (double(b) / 20.0) << "-" << (double(b + 1) / 20.0) << "%: " << histogram[b] << endl;
    }
  }
  return sum;
}


void calculateTestER(const vector<ClusterData>& clusters, const vector< pair<LocalFeatures, uint> >& localfeatures_test,
                     const vector<double>& aPriori, const EPOptions& options) {
  uint numTestSamples = localfeatures_test.size();
  uint numClasses = aPriori.size();

  DBG(10) << "calculating error rate for " << numTestSamples << " test samples..." << endl;
  double average_correct_class_probability = 0.0;
  int ok = 0, error = 0;

  for (uint n = 0; n < numTestSamples; n++) {
    uint k_n = localfeatures_test[n].second;
    uint numLocalFeatures = localfeatures_test[n].first.numberOfFeatures();
    double correct_class_probability = 0.0;
    double dummy_v;
    vector<double> dummy_emission_probabilities;

    for (uint l = 0; l < numLocalFeatures; l++) {
      vector<double> p_k_given_x(numClasses, 0.0);
      calc_k_given_x(clusters, aPriori, localfeatures_test[n].first, l, options, p_k_given_x, dummy_v, dummy_emission_probabilities);
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

void trainVariances(vector<ClusterData>& clusters, const vector< pair<LocalFeatures, uint> >& localfeatures,
                    const vector<double>& aPriori, const EPOptions& options, const vector< pair<LocalFeatures, uint> >& localfeatures_test,
                    int numIterations, double epsilon, double percentWorstImages, const vector<string>& outputfile) {

  uint numClusters = clusters.size();
  uint numClasses = aPriori.size();
  uint numDims = clusters[0].density.dim;

  vector< vector<double> > updates(numClusters, vector<double>(clusters[0].density.dim, 0.0));
  vector< vector<double> > saved_u(localfeatures.size());
  vector< vector<double> > saved_v(localfeatures.size());
  vector< vector< vector<double> > > saved_emission_probability(localfeatures.size());
  for (uint n = 0; n < localfeatures.size(); n++) {
    saved_u[n].resize(localfeatures[n].first.numberOfFeatures());
    saved_v[n].resize(localfeatures[n].first.numberOfFeatures());
    saved_emission_probability[n].resize(localfeatures[n].first.numberOfFeatures(), vector<double>(clusters.size()));
  }

  //vector<ClassifiedImage> save_classification_order(localfeatures.size());
  //double probUpdatedImagesInThisIteration = 0.0;
  //double probUpdatedImagesInPreviousIteration = 0.0;
  //double saveProbUpdatedImagesInPreviousIteration = 0.0;
  //int stopTraining = 0;
  //int numImagesForUpdate = (int) round((double) numImages * percentWorstImages / 100.0);

  if (!localfeatures_test.empty()) {
    DBG(15) << "calculating initial error rate on testing data" << endl;
    calculateTestER(clusters, localfeatures_test, aPriori, options);
  }


  for (int iteration = 0; iteration < numIterations; iteration++) {

    BLINK(15) << endl;
    DBG(15) <<  "starting iteration " << (iteration + 1) << endl;

    bool firstCalculation = true;
    for (uint c = 0; c < numClusters; c++) {
      for (uint d = 0; d < numDims; d++) {
        double f_deriv = calcFderiv(clusters, localfeatures, aPriori, options, clusters[c].clazz, c, d, saved_u, saved_v,
                                    saved_emission_probability, firstCalculation);
        firstCalculation = false;
        updates[c][d] = f_deriv;
      }
    }

    int numMinimal = 0;
    double minVarianceValue = 1E-08;
    double maxRelativeUpdate = 0.0, maxAbsoluteUpdate = 0.0;
    double maxRelativeUpdatePreviousValue = 0.0;
    double maxRelativeUpdateUpdateValue = 0.0;
    for (uint c = 0; c < numClusters; c++) {
      for (uint d = 0; d < numDims; d++) {
        updates[c][d] *= epsilon;
        if (fabs(updates[c][d]) / clusters[c].density.sigma[d] > maxRelativeUpdate) {
          maxRelativeUpdate = fabs(updates[c][d]) / clusters[c].density.sigma[d];
          maxRelativeUpdatePreviousValue = clusters[c].density.sigma[d];
          maxRelativeUpdateUpdateValue = updates[c][d];
        }
        if (fabs(updates[c][d]) > maxAbsoluteUpdate) {
          maxAbsoluteUpdate = fabs(updates[c][d]);
        }
        clusters[c].density.sigma[d] += updates[c][d];
        if (clusters[c].density.sigma[d] <= minVarianceValue) {
          numMinimal++;
          clusters[c].density.sigma[d] = minVarianceValue;
        }
      }
    }
    DBG(15) << numMinimal << " variance values have minimal value " << minVarianceValue << endl;
    DBG(15) << "max. rel. update: " << maxRelativeUpdate << endl;
    DBG(15) << "value " << maxRelativeUpdatePreviousValue << " was updated by " << maxRelativeUpdateUpdateValue << endl;
    DBG(15) << "max. abs. update: " << maxAbsoluteUpdate << endl;
    calcClusterNormalizationValues(clusters, options);

    if (((iteration % 10) == 9) && (!localfeatures_test.empty())) {
      calculateTestER(clusters, localfeatures_test, aPriori, options);
    }
    if (outputfile.size() > 0) {
      for (uint k = 0; k < numClasses; k++) {
        saveModel(outputfile[k], k, clusters);
      }
    }

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
  DBG(10) << "updating variances with epsilon=" << epsilon << endl;
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
    if (modelDatabase.numberOfSuffices() != 1) {
      ERR << "found " << modelDatabase.numberOfSuffices() << " suffices in model database, exactly 1 required!" << endl;
      exit(1);
    }
  }

  uint numClasses = getNumClasses(modelDatabase);
  vector<double> aPriori(numClasses);

  // load further options
  EPOptions options;

  // load models
  DBG(10) << "loading models..." << endl;
  vector<ClusterData> clusters;
  for (uint k = 0; k < modelDatabase.size(); k++) {
    string fname = modelDatabase.filename(k) + "." + modelDatabase.suffix(0);
    loadModel(fname, modelDatabase[k]->clas(), clusters);
  }
  DBG(15) << "loaded " << clusters.size() << " clustermodels for " << numClasses << " classes" << endl;

  string outputfilelist = "";
  vector<string> outputfiles;
  if (cl.search(2, "-s", "--savefile")) {
    outputfilelist = cl.next("");
    DBG(15) << "saving clustermodels to files in the filelist " << outputfilelist << endl;
    Database saveModelDatabase;
    saveModelDatabase.loadFileList(outputfilelist);
    if (saveModelDatabase.numberOfSuffices() != 1) {
      ERR << "filelist for saving clustermodels must have exactly 1 suffix !" << endl;
      exit(1);
    }
    if (saveModelDatabase.size() != modelDatabase.size()) {
      ERR << "model database contains " << modelDatabase.size() << " cluster files, but database for saving clustermodel contains " << saveModelDatabase.size() << endl;
      exit(1);
    }
    if (getNumClasses(saveModelDatabase) != numClasses) {
      ERR << "model database contains " << numClasses << " classes, but database for saving clustermodel contains " << getNumClasses(saveModelDatabase) << endl;
      exit(1);
    }
    outputfiles.resize(saveModelDatabase.size());
    for (uint k = 0; k < saveModelDatabase.size(); k++) {
      outputfiles[saveModelDatabase[k]->clas()] = saveModelDatabase.filename(k) + "." + saveModelDatabase.suffix(0);
    }
  }

  // read local feature training files
  vector< pair<LocalFeatures, uint> > localfeatures;
  Database trainDb;
  if (!cl.search(1, "--train")) {
    USAGE();
  } else {
    trainDb.loadFileList(cl.next(" "));
    if (trainDb.numberOfSuffices() > 1) {
      ERR << "found " << trainDb.numberOfSuffices() << " suffices in training database, exactly 1 is required!" << endl;
      exit(1);
    }
    int totalFeatures = 0;
    for (uint n = 0; n < trainDb.size(); n++) {
      LocalFeatures lf;
      string fname = trainDb.filename(n) + "." + trainDb.suffix(0);
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
      ERR << "found " << testDb.numberOfSuffices() << " suffices in testing database, exactly 1 is required!" << endl;
      exit(1);
    }
    int totalFeatures = 0;
    for (uint n = 0; n < testDb.size(); n++) {
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
    for (uint k = 0; k < numClasses; k++) {
      aPriori[k] = 1.0 / double(numClasses);
      DBG(15) << "a priori probability for class " << k << ": " << aPriori[k] << endl;
    }
  } else {
    for (uint k = 0; k < numClasses; k++) {
      int filesK = 0;
      for (uint n = 0; n < trainDb.size(); n++) {
        if (trainDb[n]->clas() == k) {
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

  // pre-calculate emission probability normalization values for the clusters based on their variances
  // this speeds calculation up, if clusters are not dimension-pooled
  vector<double> clusterNormalizationValues;
  calcClusterNormalizationValues(clusters, options);

  if (inputfile == "") {
    calcInitialClusterweights(clusters, numClasses);
  } else {
    loadInitialClusterweights(clusters, inputfile);
  }

  trainVariances(clusters, localfeatures, aPriori, options, localfeatures_test, numIterations, epsilon, percentWorstImages, outputfiles);

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);

}

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
#include "diag.hpp"
#include "getpot.hpp"
#include "gzstream.hpp"
#include "localfeatures.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "gaussiandensity.hpp"
#include "database.hpp"
#include "positionclusterer.hpp"

using namespace std;

#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480
#define DEFAULT_NUM_SPLITS 3
#define DEFAULT_NUM_REESTIMATIONS 10
#define DEFAULT_CLUSTER_NR 0
#define DEFAULT_RELATIVE_CLUSTER_NR 1

// data structure to represent (untied) clusters
struct ClusterData {
	GaussianDensity density;
	vector<ClusterPosition> positions;
	vector< vector<ClusterPosition> > positionDifferences;
	vector<bool> hasPositionDifference;
	uint clazz;
};

void USAGE() {
	cout
			<< "USAGE: visualizeclusterpositions <options> --model <cluster filelist> --train <local feature filelist> --image <image file>"
			<< endl << " Options:"
			<< "   -h,--help                show this help" << endl
			<< "   --width <n>              width of newly created images, default is "
			<< DEFAULT_WIDTH << endl
			<< "   --height <n>             height of newly created images, default is "
			<< DEFAULT_HEIGHT << endl
			<< "   --splits <n>             number of splits per position cluster, default is "
			<< DEFAULT_NUM_SPLITS << endl
			<< "   --reestimations <n>      number of reestimations per split, default is "
			<< DEFAULT_NUM_REESTIMATIONS << endl
			<< "   --clusterNr <n>          visualize position probabilities for cluster n, default is "
			<< DEFAULT_CLUSTER_NR << endl
			<< "   --relativeClusterNr <n>  relative positions to cluster n, default is "
			<< DEFAULT_RELATIVE_CLUSTER_NR << endl
			<< "   --onlyobs                visualize only the observations, no probabilities"
			<< endl
			<< "   --relativePositions      visualize relative instead of absolute positions, default is false"
			<< endl
			<< "   --nodisplay              do not display image, default is false"
			<< endl << endl;
}

void loadModel(const string filename, int classNr, vector<ClusterData>& clusters) {
	igzstream is;
	is.open(filename.c_str());
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
				getline(is, line);
				if (!is.eof()) {
					istringstream iss(line);
					string keyword;
					iss >> keyword;
					if (keyword == "gaussians") {
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
					} else if ((keyword == "splitmode") || (keyword
							== "maxsplits") || (keyword == "stopWithNClusters")
							|| (keyword == "disturbMode") || (keyword
							== "poolMode") || (keyword == "dontSplitBelow")
							|| (keyword == "iterationsBetweenSplits")
							|| (keyword == "epsilon") || (keyword
							== "minObservationsPerCluster") || (keyword
							== "distance")) {
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
	DBG(15) << "read " << densities.size() << " clusters for class " << classNr << endl;
}

// calculated the cluster emission probability p(x|c,k)
double calc_log_x_given_c_k(const GaussianDensity& clustermodel,
		const vector<double>& x_l, double clusterNormalizationValue) {
	double pixelNorm = 0.0, pixelDist = 0.0;
	uint dims = x_l.size();
	for (uint d = 0; d < dims; d++) {
		double tmp = 0.0;
		tmp = x_l[d] - clustermodel.mean[d];
		tmp *= tmp;
		pixelDist += (tmp / clustermodel.sigma[d]);
	}
	pixelNorm = clusterNormalizationValue;
	return -0.5 * (pixelDist + pixelNorm);
}

int getNearestCluster(const vector<double>& x_l,
		const vector<ClusterData>& clusters, uint k,
		const vector<double>& clusterNormalizationValues) {
	int bestCluster = -1;
	bool first = true;
	double maxProb = 0.0;
	for (uint c = 0; c < clusters.size(); c++) {
		if (clusters[c].clazz == k) {
			double prob = calc_log_x_given_c_k(clusters[c].density, x_l,
					clusterNormalizationValues[c]);
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

/*
 @param clusterNr: The cluster for which absolute and relative position clusters are to be estimated
 @param progessAbsolutePositions: contains the absolute position clusters for each split and reestimation step
 @param progressRelativePosition: contains the relative position clusters for each split and reestimation step
 @param relativeClusterNr: Relative position clusters are estimated for the cluster pair <clusterNr, relativeClusterNr>.
 Must not be equal to clusterNr!

 */
void calculateClusterPositions(const Database& trainingData,
		vector<ClusterData>& clusters, uint maxSplitsPerCluster,
		uint reestimations, const vector<double>& clusterNormalizationValues,
		uint clusterNr, vector<DoubleVector>& positionObservations,
		vector< vector<double> >& positionDifferenceObservations,
		vector< vector< vector<ClusterPosition> > >& progressAbsolutePositions,
		vector< vector< vector<ClusterPosition> > >& progressRelativePositions,
		bool positionDifferences, uint relativeClusterNr) {
	PositionClusterer pc;
	pc.maxSplit() = maxSplitsPerCluster;
	pc.iterationsBetweenSplits() = reestimations;
	pc.stopWithNClusters() = 16;
	pc.minObservationsPerCluster() = 5;
	pc.dontSplitBelow() = 8;
	pc.epsilon() = 0.01;
	pc.diagonalize() = false;
	PositionClusterer pcDifferences;
	pcDifferences.maxSplit() = maxSplitsPerCluster;
	pcDifferences.iterationsBetweenSplits() = reestimations;
	pcDifferences.stopWithNClusters() = 16;
	pcDifferences.minObservationsPerCluster() = 5;
	pcDifferences.dontSplitBelow() = 8;
	pcDifferences.epsilon() = 0.01;
	pcDifferences.diagonalize() = false;

	if (clusterNr == relativeClusterNr) {
		ERR << "'clusterNr' and 'relativeClusterNr' must not be equal!" << endl;
		exit(1);
	}

	DBG(15) << "estimating cluster positions (max " << (int) pow(2.0, (double) maxSplitsPerCluster) << ")" << endl;
	uint numFiles = trainingData.size();

	vector< vector<DoubleVector> > inputdata(clusters.size());
	// get relative positions: maxNumClusters \times maxNumClusters \times 2\cdot numberOfPositionDifferences (contains data in this form x_1,y_1,x_2,y_2,....x_N,y_N
	vector< vector< vector<double> > > inputdataDifferences(clusters.size(),
			vector< vector<double> >(clusters.size()));
	// wieviele relativen position f�r jedes clusterpaar
	vector< vector< uint> > numPositionDifferences(clusters.size(),
			vector<uint>(clusters.size(), 0));

	for (uint n = 0; n < numFiles; n++) {
		string fname = trainingData.filename(n) + "." + trainingData.suffix(0);
		LocalFeatures localFeatures;
		localFeatures.load(fname);
		double imgSizeX = double(localFeatures.imageSizeX());
		double imgSizeY = double(localFeatures.imageSizeY());
		DBG(20) << "processing " << localFeatures.numberOfFeatures() << " local features from file " << fname << endl;
		vector<uint> nearestClusterNum(localFeatures.numberOfFeatures());
		for (uint l = 0; l < localFeatures.numberOfFeatures(); l++) {
			nearestClusterNum[l] = getNearestCluster(localFeatures[l],
					clusters, trainingData[n]->clas(),
					clusterNormalizationValues);
			if (nearestClusterNum[l] >= clusters.size()) {
				ERR << "nearest cluster should not be " << nearestClusterNum[l] << endl;
				exit(1);
			}
			DoubleVector nextPositionItem(2);
			nextPositionItem[0] = double(localFeatures.position(l).x) / imgSizeX;
			nextPositionItem[1] = double(localFeatures.position(l).y) / imgSizeY;
			inputdata[nearestClusterNum[l]].push_back(nextPositionItem);
		}

		//collect relative position information:
		if (positionDifferences) {
			for (uint l = 0; l < localFeatures.numberOfFeatures(); l++) {
				for (uint l2 = 0; l2 < localFeatures.numberOfFeatures(); l2++) {
					if (l != l2) {

						double deltax=((double(localFeatures.position(l2).x) / imgSizeX) - (double(localFeatures.position(l).x) / imgSizeX));
						double deltay=((double(localFeatures.position(l2).y) / imgSizeY) - (double(localFeatures.position(l).y) / imgSizeY));

						if (inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].size()
								<= numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]]) {
							inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].push_back(deltax);
							inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]].push_back(deltay);
						} else {
							uint
									position =
											numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]];
							inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]][position]
									= (deltax);
							inputdataDifferences[nearestClusterNum[l]][nearestClusterNum[l2]][position + 1]
									= (deltay);
						}
						numPositionDifferences[nearestClusterNum[l]][nearestClusterNum[l2]]
								+= 2;
					}
				}
			}
		}
	}

	DBG(15) << "clustering positions for cluster " << clusterNr << " and position differences for cluster pair <" << clusterNr << ", " << relativeClusterNr << ">" << endl;
	vector<ClusterPosition> cl;
	DBG(15) << "clustering " << inputdata[clusterNr].size() << " positions" << endl;
	pc.run(inputdata[clusterNr], cl, progressAbsolutePositions, true);
	clusters[clusterNr].positions.resize(cl.size());
	DBG(20) << "got " << cl.size() << " position clusters" << endl;
	for (uint p = 0; p < cl.size(); p++) {
		clusters[clusterNr].positions[p] = cl[p];
		clusters[clusterNr].positions[p].covariance.calcCachedData();
		DBG(25) << " " << p << "(" << cl[p].elements << "): mean=" << cl[p].mean.first << "," << cl[p].mean.second << endl;
		DBG(25) << "    covariance=" << cl[p].covariance.xx << "," << cl[p].covariance.xy << "," << cl[p].covariance.yx << "," << cl[p].covariance.yy << endl;
	}

	// now cluster relative positions
	if (positionDifferences) {
		clusters[clusterNr].positionDifferences.resize(clusters.size());
		clusters[clusterNr].hasPositionDifference.resize(clusters.size());
		int relativeclusters = 0;
		// iterate over all appearance clusters (for a second time, we are in a loop over these already (c)
		if (numPositionDifferences[clusterNr][relativeClusterNr] / 2 >= 10) {
			// prepare data structure for result of clustering
			vector<ClusterPosition> result;

			// copy data into data structure that can be used for the position clusterer
			vector< vector<double> > differences(
					numPositionDifferences[clusterNr][relativeClusterNr] / 2,
					vector<double>(2));
			int count = 0;
			uint i=0;
			while (i<numPositionDifferences[clusterNr][relativeClusterNr]) {
				differences[count][0]
						= inputdataDifferences[clusterNr][relativeClusterNr][i];
				i++;
				differences[count][1]
						= inputdataDifferences[clusterNr][relativeClusterNr][i];
				count++;
				i++;
			}

			// now start clustering -> result is stored in "result"
			DBG(15) << "clustering " << differences.size() << " position differences" << endl;
			pcDifferences.run(differences, result, progressRelativePositions,
					true);
			positionDifferenceObservations = differences;

			clusters[clusterNr].positionDifferences[relativeClusterNr].resize(result.size());
			DBG(15) << "got " << result.size() << " position difference clusters" << endl;
			uint numElements = 0;
			for (uint p = 0; p < result.size(); p++) {
				DBG(25) << " " << p << "(" << result[p].elements << "): mean=" << result[p].mean.first << "," << result[p].mean.second << endl;
				DBG(25) << "    covariance=" << result[p].covariance.xx << "," << result[p].covariance.xy << "," << result[p].covariance.yx << "," << result[p].covariance.yy << endl;
				clusters[clusterNr].positionDifferences[relativeClusterNr][p]
						= result[p];
				clusters[clusterNr].positionDifferences[relativeClusterNr][p].covariance.calcCachedData();
				numElements
						+= clusters[clusterNr].positionDifferences[relativeClusterNr][p].elements;
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
			for (uint p = 0; p < result.size(); p++) {
				clusters[clusterNr].positionDifferences[relativeClusterNr][p].elements_log
						= log((double) clusters[clusterNr].positionDifferences[relativeClusterNr][p].elements
								/ (double) numElements);
			}
			clusters[clusterNr].hasPositionDifference[relativeClusterNr] = true;
			relativeclusters++;
		} else {
			DBG(15) << "not clustering position differences for clusters " << clusterNr << " and " << relativeClusterNr <<
			", only (" << inputdataDifferences[clusterNr][relativeClusterNr].size() << " position differences)" << endl;
			clusters[clusterNr].hasPositionDifference[relativeClusterNr]
					= false;
		}
	}

	positionObservations = inputdata[clusterNr];
}

// precalculates the normalization values for the cluster emission probabilities
// this is possible since the normalization values are independent from the local features
// in case of no dimension pooling, this helps speeding up the emission probability calculation
void calcClusterNormalizationValues(const vector<ClusterData>& clusters,
		vector<double>& clusterNormalizationValues) {
	clusterNormalizationValues.resize(clusters.size());
	for (uint c = 0; c < clusters.size(); c++) {
		clusterNormalizationValues[c] = 0.0;
		for (uint d = 0; d < clusters[c].density.dim; d++) {
			clusterNormalizationValues[c] += log(2.0 * M_PI
					* clusters[c].density.sigma[d]);
		}
	}
}

double calc_position_probability(const pair<double, double>& mean,
		const PositionCovarianceMatrix& cov,
		const pair<double, double>& actualPosition) {
	double xDiff = actualPosition.first - mean.first;
	double yDiff = actualPosition.second - mean.second;
	return cov.normValue * exp(-0.5 * (xDiff * xDiff * cov.inverted_xx + xDiff
			* yDiff * (cov.inverted_xy + cov.inverted_yx) + yDiff * yDiff
			* cov.inverted_yy));
}

void drawPositionProbabilities(ImageFeature& img, const ClusterData& cluster) {
	uint width = img.xsize();
	uint height = img.ysize();
	for (uint x = 0; x < width; x++) {
		for (uint y = 0; y < height; y++) {
			uint totalElements = 0;
			double prob = 0.0;
			for (uint i = 0; i < cluster.positions.size(); i++) {
				double p = calc_position_probability(cluster.positions[i].mean,
						cluster.positions[i].covariance, pair<double, double>(
								double(x) / double(width), double(y) / double(height)));
				totalElements += cluster.positions[i].elements;
				prob += p * (double) cluster.positions[i].elements;
			}
			prob /= (double) totalElements;
			img(x, y, 0) = 1.0 - prob;
			img(x, y, 1) = 1.0 - prob;
			img(x, y, 2) = 1.0 - prob;
		}
	}
	normalize(img);
}

void drawRelativePositionProbabilities(ImageFeature& img,
		const ClusterData& cluster, uint relativeClusterNr) {
	int width = img.xsize();
	int height = img.ysize();
	for (int x = -width; x < width; x++) {
		for (int y = -height; y < height; y++) {
			uint totalElements = 0;
			double prob = 0.0;
			for (uint i = 0; i
					< cluster.positionDifferences[relativeClusterNr].size(); i++) {
				double
						p =
								calc_position_probability(
										cluster.positionDifferences[relativeClusterNr][i].mean,
										cluster.positionDifferences[relativeClusterNr][i].covariance,
										pair<double, double>(double(x) / double(width), double(y) / double(height)));
				totalElements
						+= cluster.positionDifferences[relativeClusterNr][i].elements;
				prob
						+= p
								* (double) cluster.positionDifferences[relativeClusterNr][i].elements;
			}
			prob /= (double) totalElements;
			img((uint) (((double) (x + width)) / 2.0), (uint) (((double) (y
					+ height)) / 2.0), 0) = 1.0 - prob;
			img((uint) (((double) (x + width)) / 2.0), (uint) (((double) (y
					+ height)) / 2.0), 1) = 1.0 - prob;
			img((uint) (((double) (x + width)) / 2.0), (uint) (((double) (y
					+ height)) / 2.0), 2) = 1.0 - prob;
		}
	}
	normalize(img);
}

void drawObservations(ImageFeature& img, const vector<DoubleVector>& assigns) {
	for (uint i = 0; i < assigns.size(); i++) {
		int x = (int) (assigns[i][0] * (double) img.xsize());
		int y = (int) (assigns[i][1] * (double) img.ysize());
		img(x, y, 0) = 0;
		img(x, y, 1) = 1;
		img(x, y, 2) = 0;
		img(x + 1, y, 0) = 1;
		img(x + 1, y, 1) = 0;
		img(x + 1, y, 2) = 0;
		img(x - 1, y, 0) = 1;
		img(x - 1, y, 1) = 0;
		img(x - 1, y, 2) = 0;
		img(x, y + 1, 0) = 1;
		img(x, y + 1, 1) = 0;
		img(x, y + 1, 2) = 0;
		img(x, y - 1, 0) = 1;
		img(x, y - 1, 1) = 0;
		img(x, y - 1, 2) = 0;
	}
}

void drawRelativeObservations(ImageFeature& img,
		const vector< vector<double> >& observations) {
	for (uint i = 0; i < observations.size(); i++) {
		int x = (int) (((observations[i][0] + 1.0) / 2.0)
				* (double) img.xsize());
		int y = (int) (((observations[i][1] + 1.0) / 2.0)
				* (double) img.ysize());
		img(x, y, 0) = 0;
		img(x, y, 1) = 1;
		img(x, y, 2) = 0;
		img(x + 1, y, 0) = 1;
		img(x + 1, y, 1) = 0;
		img(x + 1, y, 2) = 0;
		img(x - 1, y, 0) = 1;
		img(x - 1, y, 1) = 0;
		img(x - 1, y, 2) = 0;
		img(x, y + 1, 0) = 1;
		img(x, y + 1, 1) = 0;
		img(x, y + 1, 2) = 0;
		img(x, y - 1, 0) = 1;
		img(x, y - 1, 1) = 0;
		img(x, y - 1, 2) = 0;
	}
}

int main(int argc, char**argv) {
	GetPot cl(argc, argv);

	if (cl.search(2, "-h", "--help")) {
		USAGE();
		exit(0);
	}

	int width= DEFAULT_WIDTH;
	int height= DEFAULT_HEIGHT;

	Database modelDatabase;
	if (!cl.search(1, "--model")) {
		USAGE();
		exit(1);
	} else {
		modelDatabase.loadFileList(cl.next(" "));
	}

	Database trainDb;
	if (!cl.search(1, "--train")) {
		USAGE();
		exit(1);
	} else {
		trainDb.loadFileList(cl.next(""));
	}

	if ((trainDb.numberOfSuffices() != 1) || (modelDatabase.numberOfSuffices()
			!= 1)) {
		ERR << "only 1 suffix allowed in clustermodel database and training image database !" << endl;
		exit(1);
	}

	// load models
	vector<ClusterData> clusters;
	DBG(10) << "loading models..." << endl;
	for (int m = 0; m < (int) modelDatabase.size(); m++) {
		string fname = modelDatabase.filename(m) + "."
				+ modelDatabase.suffix(0);
		loadModel(fname, modelDatabase[m]->clas(), clusters);
	}

	string imagefile = "";
	if (!cl.search(1, "--image")) {
		USAGE();
		exit(1);
	} else {
		imagefile = cl.next("");
	}

	bool onlyobs = cl.search(1, "--onlyobs");
	bool relativePositions = cl.search(1, "--relativePositions");

	ImageFeature img;
	// test if image exists
	ifstream is;
	is.open(imagefile.c_str());
	bool emptyImage;
	if (!is || !is.good()) {
		if (cl.search(1, "--width")) {
			width = cl.next(DEFAULT_WIDTH);
		}
		if (cl.search(1, "--height")) {
			height = cl.next(DEFAULT_HEIGHT);
		}
		// image does not exist, create empty one
		DBG(15) << "creating new image " << imagefile << " (" << width << "x" << height << ")" << endl;
		emptyImage = true;
		img.resize(width, height, 3);
	} else {
		// image exists, load it
		DBG(15) << "drawing into image " << imagefile << endl;
		img.load(imagefile);
		emptyImage = false;
		if (cl.search(1, "--width")) {
			DBG(15) << "--width parameter ignored" << endl;
		}
		if (cl.search(1, "--height")) {
			DBG(15) << "--height parameter ignored" << endl;
		}
	}

	uint maxSplitsPerCluster= DEFAULT_NUM_SPLITS;
	if (cl.search(1, "--splits")) {
		maxSplitsPerCluster = cl.next(DEFAULT_NUM_SPLITS);
	}
	uint numReestimations= DEFAULT_NUM_REESTIMATIONS;
	if (cl.search(1, "--reestimations")) {
		numReestimations = cl.next(DEFAULT_NUM_REESTIMATIONS);
	}

	vector<double> clusterNormalizationValues;
	calcClusterNormalizationValues(clusters, clusterNormalizationValues);

	vector<DoubleVector> positionObservations;
	vector< vector<double> > positionDifferenceObservations;

	vector< vector< vector<ClusterPosition> > > progressAbsolutePositions;
	vector< vector< vector<ClusterPosition> > > progressRelativePositions;

	uint clusterNr= DEFAULT_CLUSTER_NR;
	if (cl.search(1, "--clusterNr")) {
		clusterNr = cl.next(DEFAULT_CLUSTER_NR);
	}
	uint relativeClusterNr= DEFAULT_RELATIVE_CLUSTER_NR;
	if (cl.search(1, "--relativeClusterNr")) {
		relativeClusterNr = cl.next(DEFAULT_RELATIVE_CLUSTER_NR);
	}

	calculateClusterPositions(trainDb, clusters, maxSplitsPerCluster,
			numReestimations, clusterNormalizationValues, clusterNr,
			positionObservations, positionDifferenceObservations,
			progressAbsolutePositions, progressRelativePositions,
			relativePositions, relativeClusterNr);

	if (!onlyobs) {
		if (relativePositions) {
			drawRelativePositionProbabilities(img, clusters[clusterNr],
					relativeClusterNr);
		} else {
			drawPositionProbabilities(img, clusters[clusterNr]);
		}
	} else {
		if (emptyImage) {
			shift(img, 1.0);
			vector<double> black(3);
			black[0] = 0.0;
			black[1] = 0.0;
			black[2] = 0.0;
			rect(img, 0, 0, width, height, black, false);
			rect(img, 1, 1, width - 2, height - 2, black, false);
			normalize(img);
		}
	}

	if (relativePositions) {
		drawRelativeObservations(img, positionDifferenceObservations);
	} else {
		drawObservations(img, positionObservations);
	}
	if (!cl.search(1, "--nodisplay")) {
		img.display();
	}

	if (emptyImage) {
		img.save(imagefile);
	}

	if (!onlyobs) {
		uint numSplits;
		if (relativePositions) {
			numSplits = progressRelativePositions.size();
		} else {
			numSplits = progressAbsolutePositions.size();
		}
		for (uint splits = 0; splits < numSplits; splits++) {
			uint reestimsInSplit;
			if (relativePositions) {
				reestimsInSplit = progressRelativePositions[splits].size();
			} else {
				reestimsInSplit = progressAbsolutePositions[splits].size();
			}
			for (uint reestimations = 0; reestimations < reestimsInSplit; reestimations++) {
				cout << "split=" << splits << ", reestimation="
						<< reestimations << endl;
				ImageFeature tmpImg;
				tmpImg.resize(width, height, 3);
				ClusterData tmpCluster;
				tmpCluster.density = clusters[clusterNr].density;

				DBG(20) << "after splits=" << splits << ", after reestimations=" << (reestimations + 1) << ":" << endl;
				if (relativePositions) {
					tmpCluster.positionDifferences.resize(relativeClusterNr + 1);
					tmpCluster.positionDifferences[relativeClusterNr]
							= progressRelativePositions[splits][reestimations];
					for (uint i = 0; i
							< tmpCluster.positionDifferences[relativeClusterNr].size(); i++) {
						tmpCluster.positionDifferences[relativeClusterNr][i].covariance.calcCachedData();
						DBG(20) << " position nr. " << i << endl;
						DBG(20) << "  " << tmpCluster.positionDifferences[relativeClusterNr][i].mean.first << " " << tmpCluster.positionDifferences[relativeClusterNr][i].mean.second << endl;
						DBG(20) << "  " << tmpCluster.positionDifferences[relativeClusterNr][i].covariance.xx << " " << tmpCluster.positionDifferences[relativeClusterNr][i].covariance.xy
						<< " " << tmpCluster.positionDifferences[relativeClusterNr][i].covariance.yx << " " << tmpCluster.positionDifferences[relativeClusterNr][i].covariance.yy << endl;
						DBG(20) << "  " << tmpCluster.positionDifferences[relativeClusterNr][i].elements << endl;
					}
				} else {
					tmpCluster.positions
							= progressAbsolutePositions[splits][reestimations];
					for (uint i = 0; i < tmpCluster.positions.size(); i++) {
						tmpCluster.positions[i].covariance.calcCachedData();
						DBG(20) << " position nr. " << i << endl;
						DBG(20) << "  " << tmpCluster.positions[i].mean.first << " " << tmpCluster.positions[i].mean.second << endl;
						DBG(20) << "  " << tmpCluster.positions[i].covariance.xx << " " << tmpCluster.positions[i].covariance.xy
						<< " " << tmpCluster.positions[i].covariance.yx << " " << tmpCluster.positions[i].covariance.yy << endl;
						DBG(20) << "  " << tmpCluster.positions[i].elements << endl;
					}
				}

				tmpCluster.clazz = clusters[clusterNr].clazz;

				if (relativePositions) {
					drawRelativePositionProbabilities(tmpImg, tmpCluster,
							relativeClusterNr);
					drawRelativeObservations(tmpImg,
							positionDifferenceObservations);
				} else {
					drawPositionProbabilities(tmpImg, tmpCluster);
					drawObservations(tmpImg, positionObservations);
				}

				string tmpimgname;
				ostringstream oss(ostringstream::out);
				oss << "vis/" << splits << "-" << (reestimations) << ".png";
				DBG(15) << "saving image '" << oss.str() << "'" << endl;
				tmpImg.save(oss.str());
			}
		}
	}

}

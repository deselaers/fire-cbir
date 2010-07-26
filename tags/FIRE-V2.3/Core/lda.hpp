#ifndef __lda_hpp__
#define __lda_hpp__

#include <vector>
#include <string>

class LDA {
public:
  
  /** constructor which also makes the object know the dimensionality which is to be used. */
  LDA(int dim=0, int noCls=2);

  /** load the transformation from a file. 
   */
  void load(const ::std::string &filename);
  
  /** save the transformation to a file.
   */
  void save(const ::std::string &filename) const;

  void putDataFirst(const ::std::vector<double> &in, int clasz);
  void putDataSecond(const ::std::vector<double> &in, int clasz);
  
  
  /** make the object know, that there will be no further data. After
      this it is not a good idea to call putData again, and also it is
      probably not a good idea to call calcLDA, mean, covariance
      before this. */
  void dataEndFirst();
  void dataEndSecond();
  /** start the svd decomposition. 
   *
   * Attention: Does not take care whether the data is in a valid
   * state for this. Will apply the svd to anything contained in the
   * covariance. Before this is started the data has to be put into
   * with putData, and the process of doing this has to be finished
   * with dataEnd();
   */
  void calcLDA();
  
  /** return the dimensionality of the data to be put into */
  const int dim() const;
  
  /** return the mean
   * Attention: does not take care whether the process of estimation was
   * already ended. Will return whatever is contained in the mean matrix.
   */
  const ::std::vector<double> & mean() const;

  /** return the covariance matrix 
   * 
   * Attention: does not take care whether the process of estimation was
   * already ended. Will return whatever is contained in the covariance matrix.
   */
  const ::std::vector< ::std::vector <double> > & covariance() const ;
  
  /** return the eigenvalues */
  const ::std::vector<double>& eigenvalues() const {return eigenVals_;}

  /** return the eigenvector i*/
  const ::std::vector<double>& eigenvector(const int i) const {return eigenVecs_[i];}
  
  /** how many data elements have been put into this lda. */
  const int counter() const;
  
  /**
   * transform a vector due to the calculated lda. 
   * @param in vector to be transformed
   * @param dim the dimensionality of the vector to be returned. If 0 no dimensionality reduction will be applied.
   * 
   * Attention: This method does not take care if the lda is in a
   * valid state for transformation. If the svd has not yet been
   * applied the vector will be multiplied by whatever is contained in
   * the covariance matrix.
   */
  const ::std::vector<double> transform(const ::std::vector<double> &in, int dim=0) const;


  /**
   * transform a vector from lda-transformed space back into the
   * untransformed space. If the input vector is shorter than the
   * dimensionality of the untransformed space, zero padding is
   * automatically done.
   *
   * @param in vector to be transformed
   * 
   * Attention: This method does not take care if the lda is in a
   * valid state for transformation. If the svd has not yet been
   * applied the vector will be multiplied by whatever is contained in
   * the covariance matrix.
   */
  const ::std::vector<double> backTransform(::std::vector<double> in) const;

  
private:
  int dim_, noCls_, counter_;
  ::std::vector< ::std::vector<double> > withinScatter_, betweenScatter_, eigenVecs_, classMeans_;
  ::std::vector< double > eigenVals_, mean_;
  ::std::vector<int> observationsInClasses_;

};


#endif

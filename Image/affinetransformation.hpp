#ifndef __AFFINETRANSFORMATION_HPP__
#define __AFFINETRANSFORMATION_HPP__

#include <iostream>

/**
 *  rotate by alpha: [ cos(alpha)   sin(alpha)  ]
 *                   [ -sin(alpha)  cos(alpha)  ]
 *
 */

class AffineTransformation {
public:
  /// the matrix representing the Linear-Part of this affine transformation
  double A[2][2];
  /// the vector representing the translational part of this affine transformation
  double b[2];
  
  void Display(::std::ostream &os=::std::cout) const {
    os
      << "A=/ " << this->A[0][0] << " " << this->A[0][1] << " \\  b=/" << this->b[0] << " \\"<< ::std::endl
      << "  \\ "<< this->A[1][0] << " " << this->A[1][1] << " /     \\"<< this->b[1] << "/"<< ::std::endl;
  }
  
  AffineTransformation(double a1=0.0,double a2=0.0, double a3=0.0, double a4=0.0, double b0=0.0, double b1=0.0){
    A[0][0]=a1;
    A[0][1]=a2;
    A[1][0]=a3;
    A[1][1]=a4;
    b[0]=b0;
    b[1]=b1;
  }

  void transform(double &x, double &y) const {
    y=int((A[0][0]*y) + (A[0][1]*x)+ b[1]);
    x=int((A[1][0]*y) + (A[1][1]*x)+ b[0]);
  }

  void reverse(double &x, double &y) const {
    double inverted[2][2];
    inverted[0][0]=A[1][1]/(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    inverted[0][1]=A[0][1]/(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    inverted[1][0]=A[1][0]/(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    inverted[1][1]=A[0][0]/(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    
    y=inverted[0][0]*(y-b[1])+inverted[0][1]*(x-b[0]);
    x=inverted[1][0]*(y-b[1])+inverted[1][1]*(x-b[0]);
  }

  
};

inline ::std::ostream& operator<<(::std::ostream& os, const AffineTransformation & src) {
  src.Display(os);
  return os;
}


#endif

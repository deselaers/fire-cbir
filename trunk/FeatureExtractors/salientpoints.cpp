// THIS SOFTWARE IS BASED ON THE PAPER "Wavelet-based Salient Points
// for Image Retrieval" by E.Loupias, N.Sebe (Leiden Institute of
// Advanced Computer Science, Leiden University, The Netherlands,
// nicu@wi.leidenuniv.nl and S.Bres, J.-M.Jolion (Laboratoire
// Reconnaisannce de Formes et Vision, INSA Lyon,
// {loupias,jolion}@rfv.insa-lyon.fr AND THE ORIGINAL JAVA SOFTWARE
// WAS WRITTEN BY E.LOPIAS
//
// The transfer of this software into the FIRE System (and translating
// to c++) was done by Andre Hegerath and Thomas Deselaers
 
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


#include "salientpoints.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "wavelet.hpp"
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/timeb.h>

using namespace std;

struct PointSorter: public binary_function< pair<float, Point>, pair<float, Point>, bool > {
  bool operator()(pair<float, Point> x, pair<float, Point> y) {
    return x.first > y.first;
  }
};

/// convert a 1 dimensional array into a 2 dimensional array given height and width
vector< vector<float> > array1Dto2D( vector<float>& in, uint height, uint width ) {
  vector< vector<float> > out(height,vector<float>(width));
  for(uint l=0;l<height;++l) {
    for(uint w=0;w<width;++w) {
      out[l][w]=in[l*width+w];
    }
  }
  return out;
}


/// initialize salient point extraction
SalientPoints::SalientPoints(const ImageFeature& img) {
  // initialization
  BORD = 15; // BORD = border
  nb_niveaux = -1;
  
  ImageFeature grayImage = makeGray(img);
  //  floatPixels=vector<float>(grayImage.size());
  Pixels=vector<int>(grayImage.size());
  for (uint i = 0; i < grayImage.size(); i++){
    Pixels[i] = (int) floor(grayImage[i] * 255);
  }
  
  width = grayImage.xsize();
  height = grayImage.ysize();
  
  pixelsSet = true;
  imageSet = true;
}

// get gradient for the position (x,y) in pixels
int SalientPoints::grad(int x, int y) {
  if (x == 0 || x == height-1 || y == 0 || y == width-1) {
    return 0;
  } else {
    return (4 * Pixels[x * width + y] 
            - Pixels[x * width + y - 1] 
            - Pixels[(x - 1) * width + y] 
            - Pixels[x * width + y + 1] 
            - Pixels[(x + 1) * width + y]);          
  }
}

void SalientPoints::getMultires() {
  // x : width; y : height
  int i, j;
  nivX  = (int) floor(log((double) width)/log((double) 2));
  nivY  = (int) floor(log((double) height)/log((double) 2));
  nb_niveaux = min(nivX, nivY);
  
  for (j = width, i = nivX; i > 0; i--, j /=2) {
    limX[i] = j;
  }
  limX[0]=1;
  for (j = height, i = nivY; i > 0; i--, j/=2) {
    limY[i] = j;
  }
  limY[0] = 1;  	
}

Point SalientPoints::getNiveaux(int l, int c) {
  Point temp(-1,-1);
  // niveau du coefficient waveletNS (l,c), l:ligne, c:colonne
  // en utilisant limX et limY
  if (nb_niveaux == -1) { getMultires(); }

  int nl=0, nc=0;
  while (l>=limY[nl+1]) { nl++; }
  while( c>=limX[nc+1] ) { nc++; } 										
  
  temp.x = nl; temp.y = nc;
  return temp;
}


vector<Point> SalientPoints::dernierCoeff2pixels(const string& wavelet, int la, int ca, int niv_l, int niv_c) {
  
  int l_debut_zone=0, c_debut_zone=0;
  int l_fin_zone=limY[nivY-1], c_fin_zone=limX[nivX-1];
  
  if( niv_l == nivY ) {
    l_debut_zone = limY[nivY-1];
    l_fin_zone = limY[nivY];
  }
  if( niv_c == nivX ) {
    c_debut_zone = limX[nivX-1];
    c_fin_zone = limX[nivX];
  }
  if( niv_l != nivY && niv_c != nivX) {
    return vector<Point>(0);
  }
		
  int lr = la - l_debut_zone;
  int cr = ca - c_debut_zone;
		
  Point premierPixel;
  premierPixel.x = 2*lr;
  premierPixel.y = 2*cr;
		
  int l=premierPixel.x, c=premierPixel.y;

  if (wavelet=="haar") {
    vector<Point> coord(4);
    
    coord[0] = Point(l, c);
    coord[1] = Point(l + 1, c + 1);
    coord[2] = Point(l, c + 1);								
    coord[3] = Point(l + 1, c);
    return coord;
    
  } else if (wavelet=="daub4") {
    vector<Point> coord(16);
    
    int l2, c2;
    l2 = l+2;
    c2 = c+2;
	
    if( l+4>height || c+4>width) {
      return vector<Point>(0);
    }
    
    coord[0]=Point(l,c); coord[1]=Point(l+1,c+1); coord[2]=Point(l,c+1); coord[3]=Point(l+1,c); 
    coord[4]=Point(l2,c2); coord[5]=Point(l2+1,c2+1); coord[6]=Point(l2,c2+1); coord[7]=Point(l2+1,c2); 
    coord[8]=Point(l,c2); coord[9]=Point(l+1,c2+1); coord[10]=Point(l,c2+1); coord[11]=Point(l+1,c2); 
    coord[12]=Point(l2,c); coord[13]=Point(l2+1,c+1); coord[14]=Point(l2,c+1); coord[15]=Point(l2+1,c); 
    return coord;

  } else {
    DBG(10) << "Unknown wavelet selected: " << wavelet << endl;
    return vector<Point>(0);	
  }
}						
 

void SalientPoints::coordEnfant(int la, int ca, int niv_l, int niv_c, Point& coordEnfAbs) {
  if(ca == 0 && la == 0 )
    return;

  int dniv = nivX - nivY; 
  int diff = niv_l - niv_c + dniv;
  // pour gerer les images avec 2 nombre de niveaux differents
		
  if (diff==0) { //l_niveau == c_niveau ) {
    coordEnfAbs.x = limY[niv_l] + 2*(la-limY[niv_l-1]);
    coordEnfAbs.y = limX[niv_c] + 2*(ca-limX[niv_c-1]);
  } else {
    if( diff > 0 ) { //c_niveau < l_niveau ) {
      coordEnfAbs.x = limY[niv_l] + 2*(la-limY[niv_l-1]);
      coordEnfAbs.y = 2*ca;
    } else {
      coordEnfAbs.x = 2*la;
      coordEnfAbs.y = limX[niv_c] + 2*(ca-limX[niv_c-1]);
    }
  }
  return;
}


vector<Point> SalientPoints::getBlocEnfant(const string& wavelet, int la, int ca, int niv_l, int niv_c) {
  vector<Point> coord;
  
  if (wavelet== "haar") {
    coord=vector<Point>(4);
    coord[0]=Point( la, ca); 
    coord[1]=Point( la+1, ca+1); 
    coord[2]=Point( la, ca+1); 
    coord[3]=Point( la+1, ca);
    return coord;
  } else if (wavelet=="daub4") {
    int l_debut_zone = limY[niv_l], c_debut_zone=limX[niv_c];
    int l_fin_zone=limY[niv_l+1], c_fin_zone=limX[niv_c+1];
	
    int dniv = nivX-nivY;	//niv_col - niv_ligne
    int diff = niv_l - niv_c + dniv;
    // pour gerer les images avec 2 nombre de niveaux differents
		
    if (diff > 0) { //c_niveau < l_niveau ) {
      c_debut_zone = 0;
      c_fin_zone = limX[niv_l+dniv];
    }
    if (diff < 0) { //c_niveau > l_niveau ) {
      l_debut_zone = 0;
      l_fin_zone = limY[niv_c-dniv];
    }
		
    if (la+4>l_fin_zone || ca+4>c_fin_zone) {
      return vector<Point>(0);
    } 
							
    coord = vector<Point>(16);
    int la2,ca2;
    la2 = la+2;
    ca2 = ca+2;
				
    coord[0] =  Point( la, ca);
    coord[1] =  Point( la+1, ca+1);
    coord[2] =  Point( la, ca+1);								
    coord[3] =  Point( la+1, ca);

    coord[4] =  Point( la2, ca2);
    coord[5] =  Point( la2+1, ca2+1);
    coord[6] =  Point( la2, ca2+1);								
    coord[7] =  Point( la2+1, ca2);
    
    coord[8] =  Point( la, ca2);
    coord[9] =  Point( la+1, ca2+1);
    coord[10] =  Point( la, ca2+1);								
    coord[11] =  Point( la+1, ca2);
			
    coord[12] =  Point( la2, ca);
    coord[13] =  Point( la2+1, ca+1);
    coord[14] =  Point( la2, ca+1);								
    coord[15] =  Point( la2+1, ca);
    return coord;

  } else {
    return vector<Point>(0);
  }

}

Point SalientPoints::followPt(const string &wavelet, int la_parent, int ca_parent, float& som_multires,  vector<vector<int> >&  gradient) {
  Point Pniv = getNiveaux(la_parent,ca_parent);
  DBG(70)  << VAR(la_parent) << " " << VAR(ca_parent) << " " << VAR(Pniv.x) << " " << VAR(Pniv.y) << endl;
  Point lePt;
  int i;
  
  if( Pniv.x >= (nivY-1) || Pniv.y >= (nivX-1) ) { // nb_niveaux ) {
    // on est au dernier niveau, il faut donc renvoyer un point; attention, on a encore 
    // plusieurs pixels a ce niveau, et on va renvoyer celui qui a le plus fort gradient
    Point unPt(-1,-1);
    
    vector<Point> bloc = dernierCoeff2pixels(wavelet, la_parent, ca_parent, Pniv.x + 1, Pniv.y + 1);
    
    int blocSize = bloc.size();
    
    bool no_point = false;
    if (blocSize == 0) { no_point = true; }
    
    int grad;
    int max_grad = 0;
    int i_max = 0;                                  // i avec le gradient max
    if (no_point == false) {
      for (i = 0; i < blocSize; i++ ) {
        grad = gradient[bloc[i].x][bloc[i].y];                           
        if ((grad == max_grad) && (max_grad != 0) )
          no_point = true;
        
        if (grad > max_grad) {
          max_grad = grad;
          i_max = i;
        }
      }
    }
    
    if (no_point)
      lePt = Point(-1,-1);
    else {
      lePt =  Point(bloc[i_max].x, bloc[i_max].y);
    }
    som_multires = fabs(floatPixels[la_parent * width + ca_parent]);
    return lePt;
  }
                
  Point coordA_enfant = Point(0,0);   
  coordEnfant(la_parent, ca_parent, Pniv.x+1, Pniv.y+1, coordA_enfant);
  
  int la = coordA_enfant.x;
  int ca = coordA_enfant.y;
  
  vector<Point> bloc = getBlocEnfant(wavelet, la, ca, Pniv.x+1, Pniv.y+1 );
  int blocSize = bloc.size();
  if (blocSize == 0) {
    return Point(-1,-1);
  }
  
  float max = 0;
  int i_max = 0;
  for (i=0; i < blocSize; i++) {
    if (fabs(floatPixels[bloc[i].x * width + bloc[i].y]) > max) {
      max = fabs(floatPixels[bloc[i].x * width + bloc[i].y]);
      i_max = i;
    }
  }
  lePt = followPt(wavelet, bloc[i_max].x, bloc[i_max].y, som_multires, gradient);
  som_multires += fabs(floatPixels[la_parent * width + ca_parent]);
  
  return lePt;
}


vector<Point> SalientPoints::sPoints(const string& wavelet, vector<float>& floatPixels, int bord) {
  int l, c, i;
  int nb = 0;
  int nb_total = 0;
  
  getMultires();
  vector<float> map(height*width,0);
  vector< vector<int> > gradient(height,vector<int>(width));
  for(l = 0; l < height; l++) {
    for(c = 0;c < width; c++) {
      gradient[l][c] = grad(l, c);
    }
  }
  
  for (l = 0; l < height; l++) {
    for (c = 0; c < width; c++) {
      if ((fabs(floatPixels[l * width + c]) > 0.1) && (!(l == 0 && c == 0))) {
        float som_multires=0.0;
        Point unPt = followPt(wavelet, l, c, som_multires, gradient);
        DBG(70) << VAR(l) << " " << VAR(c) <<" " << VAR(som_multires) << endl;
        if ((unPt.x!=-1 || unPt.y!=-1) && 
            (unPt.x >= bord) && 
            (unPt.x < height - bord) && 
            (unPt.y >= bord) && 
            (unPt.y < width - bord )) {
          DBG(70) << " " <<VAR(unPt.x) <<" " << VAR(unPt.y)<< endl;
          if (map[unPt.x * width + unPt.y] == 0 ) {
            nb_total++;
          }
          map[unPt.x * width + unPt.y] = max(som_multires, map[unPt.x * width + unPt.y]);
          // si plusieurs coeff tombent sur le meme point, on garde le plus grand coeff
        }
      }
    }
  }
  if (nb_total == 0) {
    return vector<Point>(0);
  }
  
  vector< vector<float> > temp=array1Dto2D(map, height, width);					
  floatPixels=map;
  
  // ce qui permet d'avoir la valeur du coeff correspondant en cliquant sur le point
  vector< pair<float, Point > > liste(nb_total);
  i = 0;
  for (l = 0; l < height; l++) {
    for (c = 0; c < width; c++) {
      if (temp[l][c] > 0) {
        DBG(70) << VAR(i) << " " << VAR(temp[l][c]) << endl;
	Point p;
	p.x = c;
	p.y = l;
        liste[i++] = pair<float, Point>(temp[l][c], p);
      }
    }
  }

  // on trie la liste pour avoir tous les coeffs dans l'ordre croissant

  sort(liste.begin(), liste.end(), PointSorter());
  vector<Point> points(liste.size());
  for (int i = 0; i < (int) liste.size(); i++) {
    points[i] = liste[i].second;
  }
  return points;
  
  /*
  if (seuil <= 1) {					// on ne garde qu'un pourcentage
    double somme = 0;
    for (i = 0; i < nb_total; i++)
      somme += liste[i];
    double reSomme = 0;
    i = nb_total;														
    // on recupere tous les plus grands coeff tels que
    // leur somme represente un certain pourcentage de la somme totale
    do {
      i --;                        
      reSomme = reSomme + liste[i];                           
    } while (reSomme < somme * (100 - (float)(seuil * 100)) / 100); 
    seuil = liste[i];
  } else { // on garde un nombre fixe de points
    if (seuil > nb_total)
      seuil = nb_total;
    seuil = liste[nb_total - (int) seuil];
  }
												

  vector<vector<int> > pts(height,vector<int>(width));
  for (l = 0; l < height; l++) {
    for(c = 0; c < width; c++) {
      if (temp[l][c] >= seuil ) {
        if (display) {
          Pixels[l * width + c] = 255; 
          Pixels[(l - 1) * width + c] = 255; 
          Pixels[ (l+1) *width + c ] = 255; 
          Pixels[ l *width + c-1 ] = 255; 
          Pixels[ l *width + c+1 ] = 255; 
          Pixels[ (l-2) *width + c ] = 255; 
          Pixels[ (l-2) *width + c+1 ] = 255; 
          Pixels[ (l-1) *width + c+1 ] = 255; 
          Pixels[ l *width + c-2 ] = 255; 
          Pixels[ l *width + c+2 ] = 255; 
          Pixels[ l *width + c+3 ] = 255; 
          Pixels[ (l+1) *width + c-2 ] = 255; 
          Pixels[ (l+1) *width + c-1 ] = 255; 
          Pixels[ (l+1) *width + c+1 ] = 255; 
          Pixels[ (l+1) *width + c+2 ] = 255; 
          Pixels[ (l+1) *width + c+3 ] = 255; 
          Pixels[ (l+2) *width + c ] = 255; 
          Pixels[ (l+2) *width + c+1 ] = 255; 
          Pixels[ (l+3) *width + c ] = 255; 
          Pixels[ (l+3) *width + c+1 ] = 255; 
        }
        nb++;
        pts[l][c] = 255;
      } else {
        pts[l][c] = 0;
      }
    }
  }
  
  vector< Point > lesPoints(nb);
  i = 0;
  for (l = 0; l < height; l++) {
    for (c = 0; c < width; c++) {
      if (pts[l][c] == 255) {
        lesPoints[i++]=Point(c, l);
      }
    }
  }
  
  return lesPoints;
  */
}

vector<Point> SalientPoints::getSalientPoints(int nrPoints) {
  vector<int> bidon(width*height);
  Wavelet wav("daub4", false, 1, *this, bidon);
  //  for(uint i=0;i<floatPixels.size();++i) {
  //    cout << floatPixels[i] << " ";
  //  }cout << endl;
  struct timeb timeNow, timeLastChanged;
  int elapsedTime, second, milliSecond;
  ftime (&timeLastChanged);
  
  vector<Point> pointVector=sPoints("daub4", floatPixels, BORD);
  if(pointVector.size()>nrPoints) {
    pointVector.resize(nrPoints);
  } else {
    DBG(10) << "Could not find " << nrPoints << " salient points. Only produced " << pointVector.size() << "." << endl;
  }

  ftime (&timeNow);
  
  second = (int) difftime (timeNow.time, timeLastChanged.time);
  milliSecond = timeNow.millitm - timeLastChanged.millitm;
  elapsedTime = second * 1000 + milliSecond + 1;
  DBG(20) << "Calculation took " << elapsedTime << " msecs." << endl;
  
  return pointVector;
}

vector<Point> SalientPoints::getSalientPoints() {
  vector<int> bidon(width*height);
  Wavelet wav("daub4", false, 1, *this, bidon);
  vector<Point> pointVector=sPoints("daub4", floatPixels, BORD);
  return pointVector;
}

vector<int>& SalientPoints::getPixels() {
  return Pixels;
}


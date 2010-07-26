// THIS SOFTWARE IS BASED ON THE PAPER
// "Wavelet-based Salient Points for Image Retrieval"
// by E.Loupias, N.Sebe (Leiden Institute of Advanced Computer Science, Leiden University, The Netherlands, nicu@wi.leidenuniv.nl
// and S.Bres, J.-M.Jolion (Laboratoire Reconnaisannce de Formes et Vision, INSA Lyon, {loupias,jolion}@rfv.insa-lyon.fr
// AND SOFTWARE WRITTEN BY E.LOPIASusing namespace std;

#include <vector>
#include "wavelet.hpp"
#include <cmath>

using namespace std;

vector<int> linearisation(int dim,  vector<int>& in) {
  vector<int> out(dim);
  
  int histogram[256];
  double histSom[256];
  int l,i;
  
  for(i=0;i<256;i++) histogram[i]=0;
  for(l=0;l<dim;l++) {
    i=in[l];
    histogram[i]++;
  }
  
  histSom[0]= histogram[0];
  for(i=1;i<256;i++)
    histSom[i] = histSom[i-1] + histogram[i];
  for(i=0;i<256;i++)
    histSom[i] = 255 * histSom[i] / dim;
  
  for(l=0;l<dim;l++) {
    i=in[l];
    out[l] = (int) histSom[i];
  }
  
  return out;
}

vector<int> recadrage(int dim,  vector<float>& in, int min, int max) {
  // recadre 'in' entre 'min' et 'max'
  vector<int> out(dim);
  
  int l;
  float i;
  float imagemin= (float) (1 << 30);
  float imagemax= (float) (- (1 << 30));
  for(l = 0; l < dim; l++) {
    i = in[l];
    if (imagemin>i) imagemin=i;
    if (imagemax<i) imagemax=i;
  }
  
  if( imagemax == imagemin ) return vector<int>(0);
  
  for(l = 0;l < dim; l++) {
    i = in[l];
    out[l]= (int) (min + (i - imagemin) * (max - min) / (imagemax-imagemin));
  }
  
  return out;
}

vector<int> recadrage(int dim,  vector<int>& in, int min, int max) {
  vector<int> out(dim);
  
  int l,i;
  int imagemin = (1 << 30);
  int imagemax = -(1 << 30);
  for (l = 0;l < dim; l++) {
    i = in[l];
    if (imagemin > i) imagemin = i;
    if (imagemax < i) imagemax = i;
  }
  
  if( imagemax == imagemin ) return vector<int>(0);
  
  for(l=0;l<dim;l++) {
    i=in[l];
    out[l]=(int) (min + (i-imagemin)*(max-min)/(imagemax-imagemin));
  }
  
  return out;
}

Wavelet::Wavelet(const string& m_Wavelet, bool std_, int isign, SalientPoints& m_Exemple, vector<int>& transfoAffichable) {
  wavelet = m_Wavelet;
  std = std_;

  if (isign > 0) {				// transfo directe, on recadre dans tableau d'entiers
    m_Exemple.floatPixels = transformee( m_Wavelet, std, m_Exemple, isign );
    m_Exemple.floatSet = true;
	
    int l,c;
    vector<float> ftemp(m_Exemple.height * m_Exemple.width);
    for(l=0;l<m_Exemple.height;l++) {
      for(c = 0; c < m_Exemple.width; c++) {
        ftemp[l * m_Exemple.width + c] = fabs(m_Exemple.floatPixels[l * m_Exemple.width + c]);
      }
    }
    
    //float moy = ftemp[0];
    ftemp[0] = 0;
    
    vector<int> temp2 = recadrage(m_Exemple.height * m_Exemple.width, ftemp, 0, 255);
    vector<int> temp3 =  linearisation(m_Exemple.height * m_Exemple.width, temp2);
    vector<int> temp = recadrage(m_Exemple.height * m_Exemple.width, temp3, 0, 255);
    
    for(int i=0;i<m_Exemple.height*m_Exemple.width;++i) {
      transfoAffichable[i]=temp[i];
    }
    
  } else {								// transfo inverse
    m_Exemple.floatPixels = transformee(m_Wavelet, std, m_Exemple.floatPixels, m_Exemple.width, m_Exemple.height, isign);
    m_Exemple.floatSet = false;
	
    // conversion en tableau d'int pour affichage
    for(int i=0; i<m_Exemple.width * m_Exemple.height; i++)
      transfoAffichable[i]= (int) m_Exemple.floatPixels[i];		
  }
  
  return;
}

vector<float> Wavelet::transformee(const string& wave, bool std, SalientPoints& m_Exemple, int isign) {
  // il faut d'abord lire les pixels, pour lire l'image, et avoir la dimension...
  vector<int> pixels = m_Exemple.getPixels();
  int DimX = m_Exemple.getWidth();
  int DimY = m_Exemple.getHeight();
  int l,c;
  vector<float> FloatImage(DimX*DimY);
		
  for(l=0;l<DimY;l++)
    for(c=0;c<DimX;c++)
      FloatImage[l*DimX + c]=(float) pixels[l*DimX+c];

  return transformee(wave, std, FloatImage, DimX, DimY, isign);
}


vector<float> Wavelet::transformee(const string & wave, bool std, vector<float>& floatPixels, int DimX, int DimY, int isign) {
  if (wave=="haar"){
    ERR << "haar tranformation not yet implemented!" << endl;
    exit(1);  
  } else if (wave=="daub4") {
    if (std) {
      return daub4Std(floatPixels, DimX, DimY, isign);
    } else {
      return daub4NS(floatPixels, DimX, DimY, isign);
    }
  } else {
    ERR << "Unknown transformaion: " << wave << endl;
    exit(1);
  }
}

void Wavelet::daub4_step(vector<float>& data, int nn, int isign) {
  int ns2,i,j;
  vector<float> temp(nn);

  double C0 = 0.4829629131445341;
  double C1 = 0.8365163037378079;
  double C2 = 0.2241438680420134;
  double C3 = -0.1294095225512604;
  
  for(i=0; i < nn;i++) 
    temp[i] = data[i];
  ns2 = nn / 2;
  
  if (isign > 0) {
    
    // Calcul des valeurs
    for(i = 0, j = 0; j < nn-3; j+=2, i++) {
      data[i] = (float) (C0*temp[j] + C1*temp[j+1] + C2*temp[j+2] + C3*temp[j+3]);
      data[ns2+i] = (float) (C3*temp[j] - C2*temp[j+1] + C1*temp[j+2] - C0*temp[j+3]);
    }
    data[i] = (float) (C0*temp[nn-2] + C1*temp[nn-1] + C2*temp[0] + C3*temp[1]);
    data[i+ns2] = (float) (C3*temp[nn-2] - C2*temp[nn-1] + C1*temp[0] - C0*temp[1]);
    
    if( ns2*2 != nn ) data[nn-1] = 0;		// longueur impaire . suppression du dernier
    // attention, il faut laisser le dernier inchange pour faire la transformee inverse
  }

  if (isign<0) {
    // Calcul des valeurs
    data[0] = (float) (C2*temp[ns2-1] + C1*temp[nn-1] + C0*temp[0] + C3*temp[ns2]);
    data[1] = (float) (C3*temp[ns2-1] - C0*temp[nn-1] + C1*temp[0] - C2*temp[ns2]);
    for(i=0, j=2; i<ns2-1; i++) {
      data[j++] = (float) (C2*temp[i] + C1*temp[i+ns2] + C0*temp[i+1] + C3*temp[i+ns2+1]);
      data[j++] = (float) (C3*temp[i] - C0*temp[i+ns2] + C1*temp[i+1] - C2*temp[i+ns2+1]);
    }
  }
}


vector<float> Wavelet::daub4NS(vector<float> &FloatImage, int DimX, int DimY, int isign ) {
  int i,j,l,c;
  int Dim = max(DimX, DimY);
  int Dim2 = DimX * DimY;
  vector<float> temp(Dim);
					
  //Transformée de Haar non standard
  //int niv_x = (int) floor( log((double) DimX) / log((double) 2));
  //int niv_y = (int) floor( log((double) DimY) / log((double) 2));
		
  //if( niv_x != niv_y ) System.out.println("WT_2D : Transfo non standard ("+DimX+","+DimY+") . ("+niv_x+","+niv_y+") niveaux");
  int dx, dy;
  
  if (isign>0) {
    for(i=0; i<Dim2; i++) FloatImage[i]=FloatImage[i]/Dim;
    
    dx = DimX; dy = DimY;
    while(dx>1 && dy>1) {
      for(l=0; l<dy; l++) {
        //System.out.println("f["+l+",300] = "+FloatImage[l*Dim + 300]+" ; d="+d);
        for(c=0; c<dx; c++) temp[c]=FloatImage[l*DimX + c];
        daub4_step( temp, dx, isign);
        for(c=0; c<dx; c++) FloatImage[l*DimX + c]=temp[c];
        //System.out.println("l="+l+" ; d="+d);
      }
      
      for(c=0; c<dx; c++) {
        for(l=0; l<dy; l++) temp[l]=FloatImage[l*DimX + c];
        daub4_step( temp, dy, isign);
        for(l=0; l<dy; l++) FloatImage[l*DimX + c]=temp[l];
        //System.out.println("c="+c+" ; d="+d);
      }
      dx/=2; dy/=2;
    }
    
    // pour les images dont les dimensions correspondent a des nombres de niveaux differents
    // on continue la transformee en 1 dimension autant qu'il le faut
    while(dx>1) {
      l=0;
      for(c=0; c<dx; c++) temp[c]=FloatImage[l*DimX + c];
      daub4_step( temp, dx, isign);
      for(c=0; c<dx; c++) FloatImage[l*DimX + c]=temp[c];
      dx/=2;
    }
    while(dy>1) {
      c=0;
      for(l=0; l<dy; l++) temp[l]=FloatImage[l*DimX + c];
      daub4_step( temp, dy, isign);
      for(l=0; l<dy; l++) FloatImage[l*DimX + c]=temp[l];
      dy/=2;
    }
  }
  
  if (isign<0) {
    
    int limX[32];
    int limY[32];
    int    kx  = (int) floor(log((double) DimX) / log((double) 2));
    int    ky  = (int) floor(log((double) DimY) / log((double) 2));
    for ( j=DimX, i=kx; i>0; i--, j/=2 ) limX[i]=j;
    limX[0]=1;
    for ( j=DimY, i=ky; i>0; i--, j/=2 ) limY[i]=j;
    limY[0]=1;
  	
    int dkx=0,dky=0;
    // permet de gerer la reruction d'une image avec dimensions X et Y de niveaux differents
    if( kx>ky) dkx = kx-ky;
    else dky = ky-kx;
    
    int ix,iy;
    // si l'image a des dimensions correspondant a des nombres de niveaux differents
    // on commence par traiter en 1 dimension autant qu'il le faut
    if( dkx>0 )
      for( ix=0; ix<dkx; ix++ ) {
        l=0;
        for(c=0; c<limX[ix+1]; c++) temp[c]=FloatImage[l*DimX + c];
        daub4_step( temp, limX[ix+1], isign );
        for(c=0; c<limX[ix+1]; c++)	FloatImage[l*DimX + c]=temp[c];
      }
    if (dky > 0)
      for (iy = 0; iy < dky; iy++) {
        c = 0;
        for(l=0; l<limY[iy+1]; l++)	temp[l]=FloatImage[l*DimX + c];
        daub4_step( temp, limY[iy+1], isign);
        for(l=0; l<limY[iy+1]; l++) FloatImage[l*DimX + c]=temp[l];
      }
    
    for(ix=dkx,iy=dky;ix<kx && iy<ky;ix++,iy++) {
      
      for(l=0; l<limY[iy+1]; l++) {
        for(c=0; c<limX[ix+1]; c++) temp[c]=FloatImage[l*DimX + c];
        daub4_step( temp, limX[ix+1], isign );
        for(c=0; c<limX[ix+1]; c++)	FloatImage[l*DimX + c]=temp[c];
      }
      for(c=0; c<limX[ix+1]; c++) {
        for(l=0; l<limY[iy+1]; l++)	temp[l]=FloatImage[l*DimX + c];
        daub4_step( temp, limY[iy+1], isign);
        for(l=0; l<limY[iy+1]; l++) FloatImage[l*DimX + c]=temp[l];
      }
      
    }
    for(i=0; i<Dim2; i++)
      FloatImage[i]=FloatImage[i]*Dim;
  }
  
  return FloatImage;
}													       


vector<float> Wavelet::daub4Std(vector<float>& FloatImage, int DimX, int DimY, int isign) {
  int l,c;
  int Dim = (DimX>DimY)?DimX:DimY; //max(DimX, DimY);
  vector<float> temp(Dim);
  
  // Daub4_1D des lignes
  for(l=0; l<DimY; l++) {
    for(c=0; c<DimX; c++) temp[c]=FloatImage[l*DimX + c];
    Daub4_1D(temp, DimX, isign);
    for(c=0; c<DimX; c++) FloatImage[l*DimX + c]=temp[c];
  }
  // Daub4_1D des colonnes
  for(c=0; c<DimX; c++) {
    for(l=0; l<DimY; l++) temp[l]=FloatImage[l*DimX + c];
    Daub4_1D(temp, DimY, isign);
    for(l=0; l<DimY; l++) FloatImage[l*DimX + c]=temp[l];
  }
  
  return FloatImage;
}

void Wavelet::Daub4_1D(vector<float>& data, int nn, int isign) {	
  int n;
  if (isign > 0) {
    n = nn;
    while (n > 1) {
      daub4_step( data, n, isign );
      n/=2;
    }
  }

  if (isign<0) {
    n = 2;
    while(n<=nn) {
      daub4_step( data, n, isign );
      n*=2;
    }
  }
  
  return;
}

#include <vector>
#include <cmath>
#include <math.h>
#include "vectorfeature.hpp"
#include "imagelib.hpp"
#include "globalfeatureextraction.hpp"
#include <cmath>
using namespace std;
void calc_fractal_dim(const ImageFeature &img,double & fractal_dim) {
  double Block = 0.0;
  double b, e, xDurch,yDurch, Nenner,Zaehler, rNenner,rZaehler, yAchse, ymean, correlation;
  int m, n, M, N, pixelwert, max, min, a, f, g, i, j, k, l, p, s, t; 
  int 	lambda[11];
  double X[11], Y[11], yneu[11];
  
  for (i=0; i<11; ++i) {
    X[i] = 0;
    Y[i] = 0;
    yneu[i] = 0;
  }
	
  // Grauwertskalen-Einteilung mit Blockdicke = Lambda	
  lambda[0] = 3; lambda[1] = 4; lambda[2] = 5; lambda[3] = 6; lambda[4] = 7;
  lambda[5] = 8; lambda[6] = 9; lambda[7] =10; lambda[8] =12; lambda[9] =16; lambda[10]=256;
  
  xDurch = yDurch = Nenner = Zaehler = yAchse = ymean = rNenner = rZaehler = 0;

  N = img.xsize();
  M = img.ysize();
  
  // Berechnung
  for (p=0; p<10; ++p)		//  10 Durchläufe fur 10 verschiedene Lambda
	{
      Block = 0;			// Anzahl der Blöcke zu Lambda
      
      if ( M < 170 || N < 170 )	// Anpassung für kleine Bilder
        lambda[p] = lambda[p] + 5;	// Blockgröße gößer als eins
      
      if ( M < 64 || N < 64 )	        // Anpassung für kleine Bilder
        lambda[p] = lambda[p] + 8;	// Blockgröße gößer als eins
      
      a = 256 / lambda[p] ;    	// Anzahl der Quader 
      k = N/a ;              		// Breite der Blöcke
      if(k==0)              // verhindert Endlosschleife (n = n + k)  
        k=1;
      s = (N-a*k)/2;         		// Startspalte nach Zentrierung
      if ((N-a*k)%2 > 0 )
		s = s+1;
      
      
      l = M/a ;              		// Hoehe der Blöcke
      if(l==0)              // verhindert Endlosschleife (m = m + l)  
        l=1;
      
      t = (M-a*l)/2;         		// Startzeile nach Zentrierung
      if ((M-a*l)%2 > 0 )
		t = t+1;
      
      if(t<0) t=0;          // verhindert NAN-Fehler
      if(s<0) s=0;
      
      for (m = t; m < M-l+1; m = m + l)      // "Start"-Ecke der Blöcke
        {
          for (n = s; n < N-k+1; n = n + k)
            {
              max = 0;
              min = 256;
              f = g = 0;
              b = e = 0;
              
              for (j = 0; j < l; ++j)   // innerhalb der Blöcke
                {
                  for (i = 0; i < k; ++i)
                    {
                      // setzte pixelwert auf Bildpixel (n,m)
                      if(!(n+i>=0 && n+i<N && m+j>=0 && m+j<M))   // berechne pixelwert nur, wenn er im Bild liegt
                        continue;
                      
                      pixelwert = int(img(n+i,m+j,0)); //->Pointer(n+i,m+j));
                      
                      if (pixelwert > 255)
                        {
                          DBG(15) <<"pixelwert="<< pixelwert << " > 255" << endl;
                        }
                      
                      if (max < pixelwert)    // lokales Maximum bestimmen
					    max = pixelwert; 			
                      
                      if (min > pixelwert)    // lokales Minimum bestimmen
					    min = pixelwert;
                    }
                }
              
              b = (max+1) / lambda[p];        //   Graustufen / Dicke der Blöcke
              f =  (int) rint(b);		//  obere Blockgrenze
              
              if ( ((max+1) - (lambda[p]*f)) > 0 )
                f = f+1;
              
              e = (min+1) / lambda[p];
              g = (int) rint(e);		//   untere Blockgrenze
              
              if ( ((min+1) - (lambda[p]*g)) > 0 )
                g = g+1;

              Block = Block + f - g + 1;  	 	// Anzahl der Blöcke zu Lambda
            }
        }
      
      Y[p] = log10(Block);            // log der Block-Summen
      X[p] = log10(double(lambda[p]));	// log von lambda		
      
	}
  
  Y[10] = log10(double(1.0));    		// feste Randbedingungen
  X[10] = log10(double(lambda[10]));
  
  
  //  Berechnung der Ausgleichsgeraden zu X[i] und Y[i]
  
  
  for (i=0; i<11; ++i)	{
    xDurch = xDurch + X[i];
    yDurch = yDurch + Y[i];
  }
  
  xDurch = xDurch / 11;
  yDurch = yDurch / 11;
  
  for (i=0; i<11; ++i)	{
    Nenner = Nenner + (xDurch - X[i]) * (yDurch - Y[i]);
    Zaehler = Zaehler + (xDurch - X[i]) * (xDurch - X[i]);
  }
  
  fractal_dim = (-1) * (Nenner / Zaehler) ;	// Wert von fractal dim
  //  Berechnung der Correlation
  
  yAchse = fractal_dim * xDurch + yDurch;

  for (i=0; i<11; ++i)	{
    yneu[i] = (-1)*fractal_dim * X[i]+ yAchse;
    ymean = ymean + yneu[i];
  }
  
  ymean = ymean / 11;
  
  for (i=0; i<11; ++i)	{
    rZaehler = rZaehler + (yneu[i] - ymean) * (yneu[i] - ymean);
    rNenner = rNenner + (Y[i] - ymean) * (Y[i] - ymean);
  }
  
  correlation = 0;
  correlation = sqrt(rZaehler / rNenner);   	// Abweichung der Ausgleichsgeraden
  
  DBG(15) << "correlation=" <<correlation << endl;
  DBG(15) << "fractal_dim=" <<fractal_dim << endl;
}


void calc_histogram(const ImageFeature& img, double &Coar, double &Ent){
  
  int	m, n, M, N, pixelwert, mean, i; 
  double	G, disp, summe, help; 
  double	h[256];
  
  for ( i=0 ; i<256; ++i) {
    h[i] = 0.0;
  }
  
  N = img.xsize();//->DimensionList[0]->cValues;	// Breite
  M = img.ysize();//->DimensionList[1]->cValues;	// Hoehe

  //   Berechnung
  mean = 0;
  summe = disp = 0.0;
  Coar = Ent = 0;
  n = N;
  m = M;
  G = M * N;		// Anzahl der Pixelpaare
  for (m = 0; m < M; ++m) {
    for (n = 0; n < N; ++n) {
      // setzte pixelwert auf Bildpixel (n,m)
      pixelwert = int(img(n,m,0));
      if (pixelwert > 255) {
        DBG(20) << "pixelwert="<<pixelwert<< " > 255" << endl;
      }
      h[pixelwert]++;				// Counter für Histogramm
      summe = summe + pixelwert + 1;		// Addition aller Grauwerte
    }
  }
  mean = (int) rint (summe / G);			//  Berechnung des mittleren Grauwertes
  for ( i = 0; i < 256; ++i)
    h[i] = h[i] / G;		// Verteilung
  
  for ( i=0; i<256; ++i ) {
    disp = disp + h[i] * ((i+1)-mean) * ((i+1)-mean);   // Dispersion
  }
  DBG(20) << "disp=" << disp << endl;
  
  Coar = 1.0 - (1.0 / (1.0 + disp));		//Coarseness bzw. Körnung 
  
  for ( i=0; i<256; ++i ) {
    help=h[i];
    if ( help == 0 )
      help = 1;
    
    Ent =  Ent + (-1) * help * log(help) / log(2.0);	// Entropie
  }
  
  DBG(15) << "Coar=" << Coar << " Ent =" <<  Ent << " mean=" << mean;
}

void calc_circ_moran(const ImageFeature &img,
                     double & Circ1,  double & Circ2,
                     double & Circ3,  double & Circ4,
                     double & Circ1m, double & Circ2m,
                     double & Circ3m, double & Circ4m){
  int m, n, M, N, pixelwert, Mitte, mean, z, i, j, k,l;
  double Nenner, Nennerm, P, sum, E3, E4;
  int d[4];
  double ZW[4], Summe[8];
  
  N = img.xsize();//image->DimensionList[0]->cValues;	// Breite
  M = img.ysize();//image->DimensionList[1]->cValues;	// Hoehe
  //   Berechnung
  n = N;
  m = M;
  Nenner = Nennerm = 0;
  mean = 0;
  sum = P = 0;
  Circ1 = Circ2 = Circ3 = Circ4 = 0;
  Circ1m = Circ2m = Circ3m = Circ4m = 0;
  
  if ( M < 9 || N < 9 )	// Anpassung für kleine Bilder
	{
      DBG(15) << "Eingabebild "<< M<< "x" << N << " kleiner als Template min(9x9)" << endl;
      Circ1=Circ2=Circ3=Circ4=Circ1m=Circ2m=Circ3m=Circ4m=0;    // Feature-Extraktion wird trotzdem weitergeführt
      return;
	}
  
  for (i = 0; i < M; ++i) {  //  Berechnung des mittleren Grauwertes
    for (j = 0; j < N; ++j) {
      // setzte pixelwert auf Bildpixel (n,m)
      pixelwert = int(img(j,i,0));//(int) *((unsigned char*) image->Pointer(j,i));
      
      if (pixelwert > 255) {
        DBG(15) << "pixelwert=" << pixelwert <<" > 255" << endl;
      }
      sum = sum + pixelwert + 1;	//  Summe über alle 256 Grauwerte
    }
  }	
  
  P = M * N;	
  mean = (int) rint (sum / P);		//  Durchschnitt der Grauwerte
  for (i=0; i<8; ++i)
    Summe[i] = 0;
  
  for (m = 4; m < M-5; ++m) {
    for (n = 4; n < N-5; ++n) {
      // setzte pixelwert auf Bildpixel (n,m)
      pixelwert = int(img(n,m,0));//(int) *((unsigned char*) image->Pointer(n,m));
      if (pixelwert > 255) {
        DBG(20) << "pixelwert=" << pixelwert <<" > 255" << endl;
      }
      Nenner = Nenner + (pixelwert+1) * (pixelwert+1);		//  für alle Circ 
      
      Nennerm = Nennerm + (pixelwert+1-mean)*(pixelwert+1-mean);	//  für alle Circm 
      
      Mitte = pixelwert+1;		//  zentrales Pixel für alle Bänder
      
      for (i=0; i<4; ++i)
        ZW[i] = 0;
      
      for (j=0; j<4; ++j) {	// 4 Bänder mit Abstand 1,2,3,4
        d[j] = j+1;
        z = d[j];
        
        for (k =(-1)*z; k < z+1; k=k+2*z) {
          for (l =(-1)*z; l < z+1; ++l) {
            ZW[j] = ZW[j] +int(img(n+k,m+l,0))+1;
          }
        }
              
        for (k =(-1)*z+1; k < z; ++k) {
          for (l =(-1)*z; l < z+1; l=l+2*z) {
            ZW[j] = ZW[j] +int(img(n+k,m+l,0))+1;
          }				
        }
      }
      Summe[0] = Summe[0] + Mitte * ZW[0] ;			//  für Circ1
      Summe[1] = Summe[1] + Mitte * ZW[1] ;			//  für Circ2
      E3 = 0;			//  überzählige Ecken für Circ3 und Circ3m
      
      E3=
        int(img(n-3,m-3,0))+
        int(img(n-3,m+3,0))+
        int(img(n+3,m-3,0))+
        int(img(n+3,m+3,0))+ 4;
      
      Summe[2] = Summe[2] + Mitte * (ZW[2] - E3);		//  für Circ3
      
      E4 = 0;			//  überzählige Ecken für Circ4 und Circ4m
      E4=
        int(img(n-4,m-4,0))+
        int(img(n-4,m-3,0))+
        int(img(n-3,m-4,0))+
        int(img(n-4,m+4,0))+
        int(img(n-3,m+4,0))+
        int(img(n-4,m+3,0))+
        int(img(n+4,m-4,0))+
        int(img(n+4,m-3,0))+
        int(img(n+3,m-4,0))+
        int(img(n+4,m+4,0))+
        int(img(n+4,m+3,0))+
        int(img(n+3,m+4,0))+ 12;
      
      Summe[3] = Summe[3] + Mitte * (ZW[3] + E3 - E4);	//  für Circ4
      
      
      Summe[4] = Summe[4] + (Mitte-mean) * (ZW[0] - 8*mean);		//  für Circ1m
      Summe[5] = Summe[5] + (Mitte-mean) * (ZW[1] - 16*mean);		//  für Circ2m
      Summe[6] = Summe[6] + (Mitte-mean) * (ZW[2] - E3 - 20*mean);	//  für Circ3m
      Summe[7] = Summe[7] + (Mitte-mean) * (ZW[3] + E3 - E4 - 24*mean);//  für Circ4m
    }
  }	
  
  if (Nenner != 0)	{
    Circ1  = Summe[0] / ( 8 * Nenner);
    Circ2  = Summe[1] / (16 * Nenner);
    Circ3  = Summe[2] / (20 * Nenner);
    Circ4  = Summe[3] / (24 * Nenner);
  }

  if (Nennerm != 0) {
    Circ1m = Summe[4] / ( 8 * Nennerm);
    Circ2m = Summe[5] / (16 * Nennerm);
    Circ3m = Summe[6] / (20 * Nennerm);
    Circ4m = Summe[7] / (24 * Nennerm);
  }
  
  DBG(15) << "Circ1 =" << Circ1 << ",Circ2 ="<< Circ2 << ",Circ3 =" <<Circ3<< ",Circ4 ="<<Circ4<< endl;
  DBG(15) << "Circ1m =" << Circ1m << ",Circ2m ="<< Circ2m << ",Circ3m =" <<Circ3m<< ",Circ4m ="<<Circ4m<< endl;
  
}


void calc_sgld_statistics(const ImageFeature &img,
                          vector<double> &mean,
                          vector<double> &contrast,
                          vector<double> &ang_sec_mom,
                          vector<double> &entropy){

  mean        =vector<double>(8,0.0);
  contrast    =vector<double>(8,0.0);
  ang_sec_mom =vector<double>(8,0.0);
  entropy     =vector<double>(8,0.0);
  
  int	i,j,k, d, m, n, M, N, G, Diff, pixelwert;
  int	dist[2], P[4];
  
  double	h0[256], h1[256], h2[256], h3[256];
  
  N = img.xsize();	// Breite
  M = img.ysize();	// Hoehe

  // init distance array
  dist[0] = 1;
  dist[1] = 2;
	
  //  Berechnung der Histogramme und der stat. Funktionen
  for (j = 0; j < 2; ++j) {
    d = dist[j];
    for (k=0;k<256; ++k) {
      h0[k] = 0.0;
      h1[k] = 0.0;
      h2[k] = 0.0;
      h3[k] = 0.0;
    }
    
    //  Anzahl der Pixelpaare, abhängig von Abstand und Winkel
    G = N * M;
	
    for (i=0; i<4; ++i) {
      P[i] = 0; 
    }
    
    P[0] = (N-d)* M;
    P[1] = (N-d)*(M-d);
    P[2] =  N   *(M-d);	
    P[3] = (N-d)*(M-d);	
		
    //printf("G=%d,P[0]=%d,P[1]=%d,P[2]=%d,P[3]=%d\n",G,P[0],P[1],P[2],P[3]);
    
	
    for (m = 0; m < M; ++m) {		//  Differenz zum Winkel 0° (0°=0)
      for (n = 0; n < N-d; ++n) {
        // setzte pixelwert auf Bildpixel (n,m)
        pixelwert = int(img(n,m,0));
		
        if (pixelwert > 255) {
          DBG(20) << "pixelwert=" << pixelwert <<" > 255" << endl;
        }
				
        Diff = pixelwert - int(img(n+d,m,0));
        if (Diff < 0)
          Diff = (-1)*Diff;
        h0[Diff]++;
      }
    }	

    for (m = d; m < M ; ++m) {	//  Differenz zum Winkel 45°  (45°=1)
      for (n = 0; n < N-d; ++n) {
        // setzte pixelwert auf Bildpixel (n,m)
        pixelwert = int(img(n,m,0));
        if (pixelwert > 255) {
          DBG(20) << "pixelwert=" << pixelwert <<" > 255" << endl;
        }
        Diff = pixelwert - int(img(n+d,m-d,0));
        if (Diff < 0)
          Diff = (-1)*Diff;
        h1[Diff]++;
      }
    }	
    
    for (m = d; m < M; ++m) {		//  Differenz zum Winkel 90°  (90°=2)
      for (n = 0; n < N; ++n) {
        // setzte pixelwert auf Bildpixel (n,m)
        pixelwert = int(img(n,m,0));
        if (pixelwert > 255) {
          DBG(20) << "pixelwert=" << pixelwert <<" > 255" << endl;
        }
		
        Diff = pixelwert - int(img(n,m-d,0));
        if (Diff < 0)
          Diff = (-1)*Diff;
        h2[Diff]++;
      }
    }	

    for (m = d; m < M; ++m) {		//  Differenz zum Winkel 135°  (135°=3)
      for (n = d; n < N; ++n) {
        // setzte pixelwert auf Bildpixel (n,m)
        pixelwert = int(img(n,m,0));
        if (pixelwert > 255) {
          DBG(20) << "pixelwert=" << pixelwert <<" > 255" << endl;
        }
        Diff = pixelwert - int(img(n-d,m-d,0));
        if (Diff < 0)
          Diff = (-1)*Diff;
        
        h3[Diff]++;
      }
    }	
		
    for (i=0; i<256; ++i) {
      h0[i] = h0[i] / P[0];
      h1[i] = h1[i] / P[1];
      h2[i] = h2[i] / P[2];
      h3[i] = h3[i] / P[3];
    }
    //  Berechnung der Funktionswerte für den Abstand d
    mean[j*4+0] = mean[j*4+1] = mean[j*4+2]= mean[j*4+3] = 0.0;
	
    for (i=0; i<256; ++i) {
      mean[j*4+0] += ((i*h0[i])/256);
      mean[j*4+1] += ((i*h1[i])/256);
      mean[j*4+2] += ((i*h2[i])/256);
      mean[j*4+3] += ((i*h3[i])/256);			
    }
		
    DBG(20) << "mean("<<d<<",  0°)=" << mean[j*4+0]<< endl;
    DBG(20) << "mean("<<d<<", 45°)=" << mean[j*4+1]<< endl;
    DBG(20) << "mean("<<d<<", 90°)=" << mean[j*4+2]<< endl;
    DBG(20) << "mean("<<d<<",135°)=" << mean[j*4+3]<< endl;
	
    contrast[j*4+0] = 0.0;
    contrast[j*4+1] = 0.0;
    contrast[j*4+2] = 0.0;
    contrast[j*4+3] = 0.0;
    
    for (i=0; i<256; ++i) {
      contrast[j*4+0] += ((i*i)*h0[i]);
      contrast[j*4+1] += ((i*i)*h1[i]);
      contrast[j*4+2] += ((i*i)*h2[i]);
      contrast[j*4+3] += ((i*i)*h3[i]);						
    }
		
    DBG(20) << "contrast("<<d<<",  0°)=" << contrast[j*4+0]<< endl;
    DBG(20) << "contrast("<<d<<", 45°)=" << contrast[j*4+1]<< endl;
    DBG(20) << "contrast("<<d<<", 90°)=" << contrast[j*4+2]<< endl;
    DBG(20) << "contrast("<<d<<",135°)=" << contrast[j*4+3]<< endl;
		
    ang_sec_mom[j*4+0] = 0.0;
    ang_sec_mom[j*4+1] = 0.0;
    ang_sec_mom[j*4+2] = 0.0;
    ang_sec_mom[j*4+3] = 0.0;
    
    for (i=0; i<256; ++i) {
      ang_sec_mom[j*4+0] += (h0[i]*h0[i]);			
      ang_sec_mom[j*4+1] += (h1[i]*h1[i]);
      ang_sec_mom[j*4+2] += (h2[i]*h2[i]);
      ang_sec_mom[j*4+3] += (h3[i]*h3[i]);
    }
    
    DBG(20) << "ang_sec_mom("<<d<<",  0°)=" << ang_sec_mom[j*4+0]<< endl;
    DBG(20) << "ang_sec_mom("<<d<<", 45°)=" << ang_sec_mom[j*4+1]<< endl;
    DBG(20) << "ang_sec_mom("<<d<<", 90°)=" << ang_sec_mom[j*4+2]<< endl;
    DBG(20) << "ang_sec_mom("<<d<<",135°)=" << ang_sec_mom[j*4+3]<< endl;
    
    entropy[j*4+0] = 0.0;
    entropy[j*4+1] = 0.0;
    entropy[j*4+2] = 0.0;
    entropy[j*4+3] = 0.0;
	
    for (i=0; i<256; ++i) {
      if (h0[i]==0)
        h0[i] = 1;
      entropy[j*4+0] += ((-1) * (h0[i] * log(h0[i])));
      
      if (h1[i]==0)
        h1[i] = 1;
      entropy[j*4+1] += ((-1) * (h1[i] * log(h1[i])));
      
      if (h2[i]==0)
        h2[i] = 1;
      entropy[j*4+2] += ((-1) * (h2[i] * log(h2[i])));
      
      if (h3[i]==0)
        h3[i] = 1;
      entropy[j*4+3] += ((-1) * (h3[i] * log(h3[i])));		
    }
    
    DBG(20) << "entropy("<<d<<",  0°)=" << entropy[j*4+0]<< endl;
    DBG(20) << "entropy("<<d<<", 45°)=" << entropy[j*4+1]<< endl;
    DBG(20) << "entropy("<<d<<", 90°)=" << entropy[j*4+2]<< endl;
    DBG(20) << "entropy("<<d<<",135°)=" << entropy[j*4+3]<< endl;
  } /* for j: index in distance-array */ 
}


VectorFeature getGlobalTextureFeature(const ImageFeature &inImg) {


  vector<double> resultVec(43);
  
  ImageFeature img=makeGray(inImg);
  multiply(img,255);
  
  calc_fractal_dim(img,resultVec[0]);
  calc_histogram(img,resultVec[1],resultVec[2]);
  calc_circ_moran(img, resultVec[3],resultVec[4], resultVec[5],resultVec[6], resultVec[7],resultVec[8], resultVec[9],resultVec[10]);
  
  vector<double> meanVec, contrastVec, angular_second_momentVec, entropyVec;
  calc_sgld_statistics(img, meanVec, contrastVec, angular_second_momentVec, entropyVec);
  
  // now copy gld statistics to result
  int cnt=11;
  for(int i=0;i<8;++i) {
    resultVec[cnt]=meanVec[i];
    resultVec[cnt+8]=contrastVec[i];
    resultVec[cnt+16]=angular_second_momentVec[i];
    resultVec[cnt+24]=entropyVec[i];
    ++cnt;
  }

  // create the resulting vectorfeature
  VectorFeature result(resultVec);
  return result;
}

//#include "hungarian.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#define DBGI(level,instruction) if(DEBUG_LEVEL>=level) instruction

  
void solveAssignmentProblemDoubleRect(double **Array, int **Result, int m, int n) {
  // IMPORTANT! The values of Array[size1][size2] are changed in this routine
  //
  // adopted from KNUTH 1994 (The Stanford GraphBase, pp. 74ff)
  // <numbers> refer to the paragraphs in Knuth' text 

  // local variables <14>
  int i;
  int k; // current row of interest
  int l; // current column of interest
  int j; // another interesting column
  int*col_mate;  // column matching a row or -1
  int*row_mate;  // row matching a column or -1
  int*parent_row;  // parent in the forest / ancestor of column's mate or -1
  int*unchosen_row;  // node in the forest
  int*slack_row;  // row where minimum slack[] occurs
  int t; // total number of nodes in the forest
  int q; // total number of explored nodes in the forest
  int unmatched; 
  double s;  // current matrix element of interest
  double*row_dec;  // subtracted from row (\sigma_k)
  double*col_inc;  // added to column (\tau_l)
  double*slack;  // minimum uncovered entry in column
  double del;
  double*Array_k;

  
//   void printstate(){
//     printf("\n");
//     printf("    ");
//     for(j=0;j<n;++j) {
//       if(parent_row[j]>=0) printf("  +  ");
//       else printf("     ");
//     }
//     printf("\n");
//     for(i=0;i<m;++i) {
//       if(col_mate[i]>=0 && parent_row[col_mate[i]]<0) printf("  + ");
//       else printf("    ");
//       for(j=0;j<n;++j) {
//         printf("%3d",(int)(Array[i][j]-row_dec[i]+col_inc[j]));
//         if (col_mate[i]==j) printf("* ");
//         else if (parent_row[j]==i) printf("' ");
//         else printf("  ");
//       }
//       printf("\n");
//     }
//     printf("\n");
//   }
 
  if(m>n){
    fprintf(stderr,"cannot have m>n in %s, line %d",__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  // allocate intermediate data structures <15>
  col_mate=(int*)calloc(sizeof(int),m);
  row_mate=(int*)calloc(sizeof(int),n);
  parent_row=(int*)calloc(sizeof(int),n);
  unchosen_row=(int*)calloc(sizeof(int),m);
  slack_row=(int*)calloc(sizeof(int),n);

  row_dec=(double*)calloc(sizeof(double),m);
  col_inc=(double*)calloc(sizeof(double),n);
  slack=(double*)calloc(sizeof(double),n);
      
  // initialize result matrix to zero                                                                                      
  for(i=0;i<m;++i) {
    for(j=0;j<n;++j) {
      Result[i][j]=0; 
    } 
  }
  
  // initialize t, row_mate, parent_row, col_inc, slack <16>
  t=0; // forest starts empty
  for(l=0;l<n;++l){
    row_mate[l]=-1;
    parent_row[l]=-1;
    col_inc[l]=0.0;
    slack[l]=DBL_MAX;
  }

  // only needed for visualiztion: init col_mate
  for(k=0;k<m;++k){
    col_mate[k]=-1;
  }

  //  DBGI(25,printf("input matrix\n"));
  //DBGI(25,printstate());

  // subtract row minimum <Keysers>
  DBGI(25,printf ("subtract row minimum from each row\n"));
  for(k=0;k<m;++k){
    s=Array[k][0];
    for(l=1;l<n;++l) {
      if(Array[k][l]<s) {s=Array[k][l];}
    }
    for(l=0;l<n;++l) {
      Array[k][l]-=s;
    }
  }

  //  DBGI(25,printstate());

  if(m==n){ // otherwise heuristic does not work
    // subtract column minimum from each column <12>
    for(l=0;l<n;++l){
      s=Array[0][l];
      for(k=1;k<n;++k) {
	if(Array[k][l]<s) {s=Array[k][l];}
      }
      for(k=0;k<n;++k) {
	Array[k][l]-=s;
      }
    }
  }
  // end <12>
  
  // choose starting assignment <16>
  DBGI(25,printf("choose starting assignment\n"));
  for(k=0;k<m;++k){
    s=Array[k][0];
    Array_k=Array[k];
    for(l=1;l<n;++l)
      if(Array_k[l]<s) s=Array_k[l];
    row_dec[k]=s;  // row_dec[k] = row minimum
    for(l=0;l<n;++l) // at minimum set row_mate and col_mate if column has no mate yet
      if((s==Array_k[l])&&(row_mate[l]<0)){
        col_mate[k]=l;  
        row_mate[l]=k;
        DBGI(25,printf("matching col %d==row %d\n",l,k));
        goto row_done;
      }
    col_mate[k]=-1; // if column already has a mate, row is unchosen
    DBGI(25,printf("node %d: unmatched row %d\n",t,k));
    unchosen_row[t++]=k;
  row_done:;
  }

  //  DBGI(25,printstate());

  // Hungarian algorithm <18>
  // at most m stages with O(mn) operations -> O(m^2 n) runtime
  if(t==0) goto done;
  unmatched=t;
  while(1){
    DBGI(25,printf("we have matched %d of %d rows\n",m-t,m));
    q=0;
    while(1){
      while(q<t) {
	// explore node q of the forest <19>
        //        DBGI(25,printstate());
	k=unchosen_row[q]; // k iterates over unchosen rows, indexed by q
	DBGI(25,printf("explore node %d of the forest (row %d)\n",q,k));
	s=row_dec[k];      
	for(l=0;l<n;++l)
	  if(slack[l]>0.0){
	    del=Array[k][l]-s+col_inc[l]; // this is the "real" array value (-dec+inc)
	    // ??? if (del<0.0) del=0.0;
	    if(del<slack[l]){
	      if(del<=0.0){ // we found a new zero at [k][l]  //del==0 -> DBL_EPS??
		// changed == to <= since it can only be smaller than zero due to numerical "problems"
		if(row_mate[l]<0) goto breakthru; // matching can be increased
		slack[l]=0.0; // this column will now be chosen
		parent_row[l]=k;
		DBGI(25,printf("node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k));
		unchosen_row[t++]=row_mate[l];
	      }
	      else{
		slack[l]=del;
		slack_row[l]=k;
	      }
	    }
	  } // end <19>
        q++;
      }
      // introduce new zero into the matrix by modifying row_dec and col_inc <21>
      // we have explored the entire forest; none of the unchosen rows
      // has led to a breakthrough
      // an unchosen column with smallest slack will allow further progress
      //      DBGI(25,printstate());
      //      DBGI(25,printf("we tested all candidate rows, so now all zeroes are covered\n"));
      //DBGI(25,printf("now introduce new zero into the matrix\n"));
      s=DBL_MAX;
      for(l=0;l<n;++l)
        if(slack[l]>0.0 && slack[l]<s)
          s=slack[l];  // find minimum non-zero slack
      for(q=0;q<t;++q)
        row_dec[unchosen_row[q]]+=s; // and decrease all unchosen rows
      for(l=0;l<n;++l) {
        if(slack[l]>0.0){ // column l is not chosen   
          slack[l]-=s;
          if(slack[l]<=0.0) // slack[l]==0 -> DBL_EPS??  
	    // changed == to <= since it can only be smaller than zero due to numerical "problems"
	    // look at new zero <22>
            {
              k=slack_row[l];
              DBGI(25,printf("decreasing uncovered elements by %d\n  (and increasing all doubly covered elements)\n  produces zero at [%d,%d]\n",(int)s,k,l));
	      // the next few lines are only for the "visualization"
	      for(j=l+1;j<n;++j)
		if(slack[j]<=0.0) col_inc[j]+=s;
              //	      DBGI(25,printstate());
	      for(j=l+1;j<n;++j)
            if(slack[j]<=0.0) col_inc[j]-=s; 
          if(row_mate[l]<0){
            for(j=l+1;j<n;++j)
              if(slack[j]<=0.0) col_inc[j]+=s; // slack[l]==0 -> DBL_EPS??
                goto breakthru; // matching can be increased 
          }
	      else { // not a breakthrough, but the forest continues to grow
            parent_row[l]=k;
            DBGI(25,printf("node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k));
                unchosen_row[t++]=row_mate[l];
		//printf("not a breakthrough, but the forest continues to grow\n");
              }
            }
        }
	else {
	  col_inc[l]+=s;
	}
      }
      // end <21>
    }
  breakthru:  
    // update matching by pairing row k with column l <20>
//     DBGI(25,printf("breakthrough: column %d is not matched\n",l));
//     DBGI(25,printstate());
//     DBGI(25,printf("update matching by pairing row %d with column %d\n",k,l));
    while(1){
      j=col_mate[k];
      col_mate[k]=l;
      row_mate[l]=k;
      DBGI(25,printf("(re)matching col %d==row %d\n",l,k));
      if(j<0)break;
      k=parent_row[j];
      l=j;
    }
    // end <20>
    //    DBGI(25,printstate());
    if(--unmatched==0) goto done;
    // prepare for new stage by reinitializing the forest <17>
    DBGI(25,printf("\nget ready for a new stage\n"));
    t=0;
    for(l=0;l<n;l++){
      parent_row[l]=-1;
      slack[l]=DBL_MAX;
    }
    for(k=0;k<m;++k)
      if(col_mate[k]<0){
        DBGI(25,printf("node %d: unmatched row %d\n",t,k));
        unchosen_row[t++]=k;
      }
    // end <17>
  }
 done:
  // end <18>
  
  DBGI(25,printf("done\n"));
  for (i=0;i<m;++i){
    Result[i][col_mate[i]]=1;
  }

  free(col_mate);
  free(row_mate);
  free(parent_row);
  free(unchosen_row);
  free(row_dec);
  free(col_inc);
  free(slack);
  free(slack_row);
}


void solveMinWeightEdgeCover(double **Array, int **Result, int size1, int size2)
{ 
  // this routine does NOT change the Array; 
  // it is duplicated here, because the
  // original values are needed

  // the reduction is due to a description by R. Zens following the paper
  // @misc{ keijsper96efficient,
  //  author = "J. Keijsper and R. Pendavingh",
  //  title = "An efficient algorithm for minimum-weight bibranching",
  //  text = "J. Keijsper and R. Pendavingh, An efficient algorithm for minimum-weight
  //    bibranching, Report 96-12, Faculteit WINS, Universiteit van Amsterdam, 1996.",
  //  year = "1996",
  //  url = "citeseer.nj.nec.com/keijsper96efficient.html" }

  double **DuplArray;
  double *min1,*min2;
  double min= DBL_MAX;
  double *Array_y;
  double value;
  int *argmin1,*argmin2;
  int *cover1,*cover2;
  int x,y;

  min1=(double*)malloc(size1*sizeof(double));
  min2=(double*)malloc(size2*sizeof(double));
  argmin1=(int*)malloc(size1*sizeof(int));
  argmin2=(int*)malloc(size2*sizeof(int));
  cover1=(int*)calloc(sizeof(int),size1);
  cover2=(int*)calloc(sizeof(int),size2);

  // reduction step 1: determine minimal edge weights for each node
  // at the same time: reduction step 4: determine edges that have minimal weights for each node
  for (y=0;y<size1;++y){
    min1[y]=Array[y][0];
    argmin1[y]=0;
  }
  for (x=0;x<size2;++x){
    min2[x]=Array[0][x];
    argmin2[x]=0;
  }
  for (y=1;y<size1;++y){
    Array_y=Array[y];
    for (x=1;x<size2;++x){
      value=Array_y[x];
      if(value<min1[y]){
	min1[y]=value;
	argmin1[y]=x;
      }
      if(value<min2[x]){
        min2[x]=value;
        argmin2[x]=y;
      }
    }
  }

  // reduction step 2: adjust edge weights such that we can reduce to min weight matching
  for (y=0;y<size1;++y){
    Array_y=Array[y];
    for (x=0;x<size2;++x){
      Array_y[x]-=min1[y]+min2[x];
      if(Array_y[x]<min) min=Array_y[x];
    }
  }

  // reduction step 3: use Hungarian to determine min weight matching
  if(size1<=size2){
    // duplicate array
    DuplArray=(double**)malloc(size1*sizeof(double*));
    for(y=0;y<size1;++y){
      Array_y=Array[y];
      DuplArray[y]=(double*)malloc(size2*sizeof(double));
      for (x=0;x<size2;++x){
	DuplArray[y][x]=(Array_y[x]-min);
      }
    }
    // call Hungarian
    solveAssignmentProblemDoubleRect(DuplArray, Result, size1, size2);
    // free DuplArray
    for(y=0;y<size1;++y){ 
      free(DuplArray[y]);
    }
    free(DuplArray);
  }
  else{ // this case is not yet tested
    // transpose array if sizes are wrong way
    // duplicate array
    DuplArray=(double**)malloc(size2*sizeof(double*));
    for(y=0;y<size2;++y){
      DuplArray[y]=(double*)malloc(size1*sizeof(double));
      for (x=0;x<size1;++x){
	DuplArray[y][x]=(Array[x][y]-min);
      }
    }
    // call Hungarian
    solveAssignmentProblemDoubleRect(DuplArray, Result, size2, size1);
    // free DuplArray
    for(y=0;y<size2;++y){ 
      free(DuplArray[y]);
    }
    free(DuplArray);
  }

  // the algorithm used here matches all nodes; 
  // the reduction assumes that nodes with a weight>0 are not chosen
  // therefore: delete all edges that have a weight>0
  // (the original reduction is to "max weight matching" and assumes no edges<0 are chosen)
  // also: compute how many times each node is covered (step 5)
  for (y=0;y<size1;++y){
    Array_y=Array[y];
    for (x=0;x<size2;++x){
      if (Result[y][x]) {
	++cover1[y];
	++cover2[x];
	if (Array_y[x]>0) { 
	  Result[y][x]=0; // delete edge
	  Result[y][argmin1[y]]=1; // and add alternative targets
	  Result[argmin2[x]][x]=1; 
	  ++cover1[argmin2[x]];
	  ++cover2[argmin1[y]];
	}
      }
    }
  }

  // reduction step 5: look for vertices not yet covered and cover them
  for (y=0;y<size1;++y) {
    if(!cover1[y]){
      Result[y][argmin1[y]]=1;
      if(!(cover2[argmin1[y]]++)) {
        printf("Something went wrong in the Hungarian...\n");
      } 
    }
  }
  for (x=0;x<size2;++x) {
    if(!cover2[x]){
      Result[argmin2[x]][x]=1;
      if(!(cover1[argmin2[x]]++)) {
        printf("Something went wrong in the Hungarian...\n");
      } 
    }
  }
  
  // restore edge weights (undo step 2)
  for (y=0;y<size1;++y){
    Array_y=Array[y];
    for (x=0;x<size2;++x){
      Array_y[x]+=min1[y]+min2[x];
    }
  }


  free(min1);
  free(min2);
  free(argmin1);
  free(argmin2);
  free(cover1);
  free(cover2);

}

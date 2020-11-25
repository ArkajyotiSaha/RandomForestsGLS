#include <string>
#include <limits>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}


double dist2(double &a1, double &a2, double &b1, double &b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}

void getNNIndx(int i, int m, int &iNNIndx, int &iNN){

  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  }
}

void mkNNIndx(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){

  int i, j, iNNIndx, iNN;
  double d;

  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);

  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }

#ifdef _OPENMP
#pragma omp parallel for private(j, iNNIndx, iNN, d)
#endif
  for(i = 0; i < n; i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;
    if(i != 0){
      for(j = 0; j < i; j++){
	d = dist2(coords[i], coords[n+i], coords[j], coords[n+j]);
	if(d < nnDist[iNNIndx+iNN-1]){
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;
	  rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	}
      }
    }
  }

}

std::string getCorName(int i){

  if(i == 0){
    return "exponential";
  }else if(i == 1){
    return "spherical";
  }else if(i == 2){
    return "matern";
  }else if(i == 3){
    return "gaussian";
  }else{
    error("c++ error: cov.model is not correctly specified");
  }

}

double spCor(double &D, double &phi, double &nu, int &covModel, double *bk){

  //0 exponential
  //1 spherical
  //2 matern
  //3 gaussian

  if(covModel == 0){//exponential

    return exp(-phi*D);

  }else if(covModel == 1){//spherical

    if(D > 0 && D <= 1.0/phi){
      return 1.0 - 1.5*phi*D + 0.5*pow(phi*D,3);
    }else if(D >= 1.0/phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == 2){//matern

    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

    if(D*phi > 0.0){
      return pow(D*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k_ex(D*phi, nu, 1.0, bk);//thread safe bessel
    }else{
      return 1.0;
    }
  }else if(covModel == 3){//gaussian

    return exp(-1.0*(pow(phi*D,2)));

  }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

//Description: computes the quadratic term.
double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU){

  double a, b, q = 0;
  int i, j;

#ifdef _OPENMP
#pragma omp parallel for private(a, b, j) reduction(+:q)
#endif
  for(i = 0; i < n; i++){
    a = 0;
    b = 0;
    for(j = 0; j < nnIndxLU[n+i]; j++){
      a += B[nnIndxLU[i]+j]*u[nnIndx[nnIndxLU[i]+j]];
      b += B[nnIndxLU[i]+j]*v[nnIndx[nnIndxLU[i]+j]];
    }
    q += (u[i] - a)*(v[i] - b)/F[i];
  }

  return(q);
}

//trees
Node *miniInsert(Node *Tree, double *coords, int index, int d,int n){

  int P = 2;

  if(Tree==NULL){
    return new Node(index);
  }

  if(coords[index]<=coords[Tree->index]&&d==0){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }

  if(coords[index]>coords[Tree->index]&&d==0){
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }

  if(coords[index+n]<=coords[Tree->index+n]&&d==1){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }

  if(coords[index+n]>coords[Tree->index+n]&&d==1){
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }

  return Tree;
}

void get_nn(Node *Tree, int index, int d, double *coords, int n, double *nnDist, int *nnIndx, int iNNIndx, int iNN, int check){

  int P = 2;

  if(Tree==NULL){
    return;
  }

  double disttemp= dist2(coords[index],coords[index+n],coords[Tree->index],coords[Tree->index+n]);

  if(index!=Tree->index && disttemp<nnDist[iNNIndx+iNN-1]){
    nnDist[iNNIndx+iNN-1]=disttemp;
    nnIndx[iNNIndx+iNN-1]=Tree->index;
    //fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
    rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
  }

  Node *temp1=Tree->left;
  Node *temp2=Tree->right;

  if(d==0){

    if(coords[index]>coords[Tree->index]){
      std::swap(temp1,temp2);
    }

    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);

    if(fabs(coords[Tree->index]-coords[index])>nnDist[iNNIndx+iNN-1]){
      return;
    }

    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);
  }

  if(d==1){

    if(coords[index+n]>coords[Tree->index+n]){
      std::swap(temp1,temp2);
    }

    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);

    if(fabs(coords[Tree->index+n]-coords[index+n])>nnDist[iNNIndx+iNN-1]){
      return;
    }

    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);
  }

}


void mkNNIndxTree0(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){

  int i, iNNIndx, iNN;
  double d;
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  int BUCKETSIZE = 10;


  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }

  Node *Tree=NULL;
  int time_through=-1;

  for(i=0;i<n;i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;
    if(time_through==-1){
      time_through=i;
    }

    if(i!=0){
      for(int j = time_through; j < i; j++){
	getNNIndx(i, m, iNNIndx, iNN);
	d = dist2(coords[i], coords[i+n], coords[j], coords[n+j]);
	if(d < nnDist[iNNIndx+iNN-1]){
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;

	  //fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	  rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	}
      }


      if(i%BUCKETSIZE==0){

#ifdef _OPENMP
#pragma omp parallel for private(iNNIndx, iNN)
#endif
	for(int j=time_through;j<time_through+BUCKETSIZE;j++){

	  getNNIndx(j, m, iNNIndx, iNN);
	  get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	}


	for(int j=time_through;j<time_through+BUCKETSIZE;j++){
	  Tree=miniInsert(Tree,coords,j,0, n);
	}

	time_through=-1;
      }
      if(i==n-1){

#ifdef _OPENMP
#pragma omp parallel for private(iNNIndx, iNN)
#endif
	for(int j=time_through;j<n;j++){
	  getNNIndx(j, m, iNNIndx, iNN);
	  get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	}

      }
    }
    if(i==0){
      Tree=miniInsert(Tree,coords,i,0,n);
      time_through=-1;
    }
  }

  delete Tree;
}

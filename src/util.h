#include <string>

void zeros(double *a, int n);

double dist2(double &a1, double &a2, double &b1, double &b2);

//Description: given a location's index i and number of neighbors m this function provides the index to i and number of neighbors in nnIndx
void getNNIndx(int i, int m, int &iNNIndx, int &iNN);

//Description: creates the nearest neighbor index given pre-ordered location coordinates.
//Input:
//n = number of locations
//m = number of nearest neighbors
//coords = ordered coordinates for the n locations
//Output:
//nnIndx = set of nearest neighbors for all n locations (on return)
//nnDist = euclidean distance corresponding to nnIndx (on return)
//nnIndxLU = nx2 look-up matrix with row values correspond to each location's index in nnIndx and number of neighbors (columns 1 and 2, respectively)
//Note: nnIndx and nnDist must be of length (1+m)/2*m+(n-m-1)*m on input. nnIndxLU must also be allocated on input.
void mkNNIndx(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);


std::string getCorName(int i);

double spCor(double &D, double &phi, double &nu, int &covModel, double *bk);

double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU);

//trees
struct Node{
	int index; // which point I am
	Node *left;
	Node *right;
	Node (int i) { index = i; left = right = NULL; }
};

Node *miniInsert(Node *Tree, double *coords, int index, int d,int n);

void get_nn(Node *Tree, int index, int d, double *coords, int n, double *nnDist, int *nnIndx, int iNNIndx, int iNN, int check);

void mkNNIndxTree0(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);

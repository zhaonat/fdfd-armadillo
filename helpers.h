#include <armadillo>
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

sp_cx_mat createDws(string w, int s, double* dL,int*  N);
mat bwdmean_w(mat array, string w)

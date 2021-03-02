#include <armadillo>
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

sp_cx_mat createDws(string w, int s, double* dL,int*  N);
mat bwdmean_w(mat array, string w);
cx_mat create_sfactor(double * wrange,
                   string s,
                   double omega,
                   double eps0,
                   double mu0,
                   int Nw,
                   int Nw_pml);

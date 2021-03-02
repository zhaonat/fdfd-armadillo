#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include <helpers.h>
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

int main(){

    double L0 = 1e-6; //baseline units
    int N [2] = {200,200};
    double dL [2] = {0.01, 0.01};
    int n = N[0];  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)
    double wvlen = 1.0; //microns
    double eps0 = 8.85*1e-12*L0;
    double mu0 = 4*M_PI*1e-7*L0;
    double c0 = 1/sqrt(eps0*mu0);
    double omega = 2*M_PI*c0/(wvlen);  // angular frequency in rad/sec

    //specify the dielectric function
    cx_mat eps_r = ones<cx_mat>(n,n);
    // specify source;
    cx_mat Mz = zeros<cx_mat>(n,n);
    Mz[n/2, n/2] = 1; // simple point source
    cx_mat b = 1i*omega*vectorise(Mz);

    //make a diagonal matrix using sp_cx_mat...no diag function, whihc sucks
    mat ind_cur = linspace(0,m-1,m);

    mat AB = join_rows(ind_cur, ind_cur).t();
    umat locations =  conv_to<umat>::from(AB);
    cx_mat cvals = conv_to<cx_mat>::from(vectorise(eps0*eps_r));
    sp_cx_mat Teps(locations, cvals);
    sp_cx_mat Teps_inv(locations, 1/cvals);
    sp_cx_mat Dyb = createDws("y",-1,dL,N);
    sp_cx_mat Dxf = createDws("x",1,dL,N);
    sp_cx_mat Dxb = createDws("x",-1,dL,N);
    sp_cx_mat Dyf = createDws("y",1,dL,N);

    sp_cx_mat A = Dxf*(1/mu0)*Dxb + Dyf*(1/mu0)*Dyb + omega*omega*Teps;

    // different types of solves... superlu
    cx_mat X = spsolve(A,b);

    // reshape to a field
    //cx_mat Hz = reshape(X,n,n);

    cout << X << endl;

    // save this field into a data files
    //Hz.save("test_hz", hdf5_binary);
    X.save("test_hz.csv",csv_ascii);
    X.save("test_hz.h5",hdf5_binary);

    //cout<<Teps<<endl;
    //cout << diff(mat(N)) << endl;

    return 0;


}

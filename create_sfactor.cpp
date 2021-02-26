//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include <complex>
#include <cmath>

using namespace std;
using namespace arma;

//(wrange, s, omega, eps0, mu0, Nw, Nw_pml, lnR, m)
int create_sfactor(float omega, float eps0, float mu0, int Nw, int Nw_pml, string s){
  return 2;
}

complex<double> pmlfunc(double L, double d_pml, double sigma_max, double omega, double eps0){
  // evaluate the pml scaling funciton of the complex conductivity
  //
  double m = 3.5;
  complex<double> ansr = 1. - 1i*sigma_max/(omega*eps0)*pow((L/d_pml),m);
  return ansr;
}

int main(){
    /// it appears complex does not support float, only double
    // parameters for pmlfunc
    double L = 1;
    double d_pmli = 1;
    double omega = 1;
    double eps0 = 1;

    // =============================================
    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image
    string s = "f";
    double lnR = -12.0;
    double m = 3.5;
    double eta0 = 377.0;

    double wrange [2] = {-1.0, 1.0};
    int Nw = 20;
    int Nw_pml = 5;

    // double armadillo support scalar/vector??


    double sigma_max = -(m+1)*lnR/(2*eta0)/d_pmli;


    // play with complex values
    complex<double> z1 = 1i * m*lnR;
    cout << pow(2,2)<< endl;

    //would I prefer this being a vector?
    mat w_array = linspace(wrange[0], wrange[1], Nw+1); // outputs a column vector?
    //using vector or vec from armadillo?
    // d_pml = abs(wrange - loc_pml);
    // sigma_max = -(m+1)*lnR/(2*eta0)/(d_mpl);

    // specifies where the pml begins on each side
    double loc_pml [2] = {w_array[Nw_pml], w_array[Nw-Nw_pml]};

    //d_pml = abs(wrange - loc_pml); % pml thickness

    // we may need to do this manually;
    double d_pml [2] = {abs(wrange[0]+loc_pml[0]),abs(wrange[1]-loc_pml[1])};
    // mat test = 1./w_array;
    // cout << test << endl;


    if(s == "b"){
      mat ws = w_array.rows(0,Nw-1);
    }else{
      mat ws = w_array.rows(1,Nw-1);
    }

    // extract all indices in w_array that needs to be evaluated;
    vector<double> w_left;
    vector<double> w_right;
    for(int i = 0; i< Nw; i++){
      if(w_array[i] < d_pml[0]){
        w_left.push_back(i);
      }
      if(w_array[i] > d_pml[i]){
        w_right.push_back(i);
      }
    }

    // Construct the final sfactor array
    mat sfactor_array = ones(1, Nw);
    cx_mat y = conv_to< cx_mat >::from(sfactor_array);
    // assign subset of sfactor array to the values denoted in.

    //so some casting does indeed happen
    y.cols(1,3) = {1,2,3};
    cout << y << endl;

    // let's evaluate the pml function on w_left and w_right;
    for(int i = 0; i < w_right.size(); i++){
      cout << w_right[i] << endl;
      complex<double> tval = pmlfunc(L, d_pmli, sigma_max, omega, eps0);

    }

    complex<double> tval = pmlfunc(L, d_pmli, sigma_max, omega, eps0);
    cout << tval << endl;
    //(double L, double d_pml, double sigma_max, double omega, double eps0)

    //assign left half of the w_array

    // assign right half of the array

    return 0;


}

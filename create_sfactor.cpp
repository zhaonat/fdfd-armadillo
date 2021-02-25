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
    //
    // float lnR = -12.0;
    // float m = 3.5;
    // float eta0 = 377.0;
    //
    // w_array = linspace(wrange[0], wrange[1], Nw+1);
    // loc_pml = [];
    // //using vector or vec from armadillo?
    // d_pml = abs(wrange - loc_pml);
    // sigma_max = -(m+1)*lnR/(2*eta0)/(d_mpl);
    // if(s == 'b'){
    //   mat ws = w_array[1:-1];
    // }else{
    //   mat ws = w_array[1:-1]+w_array[2:end];
    // }

//    %% Output Parameter
//    % sfactor_array: 1D array with Nw elements containing PML s-factors for Dws

//    w_array = linspace(wrange(1), wrange(2), Nw+1);
//
//    loc_pml = [w_array(1 + Nw_pml), w_array(end - Nw_pml)]; % specifies where the pml begins on each side
//    d_pml = abs(wrange - loc_pml); % pml thickness
//
//    %% what happens when we have a 0 in the denominator
//    sigma_max = -(m+1)*lnR/(2*eta0) ./ d_pml; %usually the pml is the same thickness on both sides
//
//    %% forward or backward...idk what this is requiring
//    if s == 'b'
//        ws = w_array(1:end-1);
//    else  % s == 'f'
//        assert(s == 'f');
//        ws = (w_array(1:end-1) + w_array(2:end)) / 2;
//    end
//
//    ind_pml = {ws < loc_pml(1), ws > loc_pml(2)};  % {} means cell array
//
//    sfactor_array = ones(1, Nw);
//    for n = 1:2
//        sfactor = @(L) 1 - 1i * sigma_max(n)/(omega*eps0) * (L/d_pml(n)).^m;
//        sfactor_array(ind_pml{n}) = sfactor(abs(loc_pml(n) - ws(ind_pml{n})));
//    end
  return 2;
}

complex<double> pml(double L, double d_pml, double sigma_max, double omega, double eps0){
  complex<double> ansr = 1. - 1i*sigma_max/(omega*eps0)*(L/d_pml);
  return ansr;
}

int main(){
    /// it appears complex does not support float, only double

    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image
    string s = "f";
    double lnR = -12.0;
    double m = 3.5;
    double eta0 = 377.0;
    double wrange [2] = {-1.0, 1.0};
    int Nw = 20;

    // play with complex values
    complex<double> z1 = 1i * m*lnR;
    cout << z1 << endl;

    //would I prefer this being a vector?
    mat w_array = linspace(wrange[0], wrange[1], Nw+1); // outputs a column vector?
    //using vector or vec from armadillo?
    // d_pml = abs(wrange - loc_pml);
    // sigma_max = -(m+1)*lnR/(2*eta0)/(d_mpl);

    if(s == "b"){
      mat ws = w_array.rows(0,Nw-1);
    }else{
      mat ws = w_array.rows(1,Nw-1);
    }

    // Construct the final sfactor array
    mat sfactor_array = ones(1, Nw);
    cx_mat y = conv_to< cx_mat >::from(sfactor_array);
    // assign subset of sfactor array to the values denoted in.

    //so some casting does indeed happen
    y.cols(1,3) = {1,2,3};
    cout << y << endl;

    //assign left half of the w_array

    // assign right half of the array

    return 0;


}

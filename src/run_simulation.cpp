#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include "helpers.h"
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

// should accept user inputs externally...
// once we generate an executable, we want to feed stuff to the executable
// otherwise everything is hardcoded in the .cpp file which gets compiled
// and has to be recompiled
// run executable for input ouput
// ./a.out < input.txt > output.txt
// compilation:
// g++ run_simulation.cpp create_sfactor.cpp derivatives.cpp -std=c++14 -O2 -larmadillo


int main(int argc, char *argv[]){

    if( argc == 2 ) {
        printf("The argument supplied is %s\n", argv[1]);
     }
     else if( argc > 2 ) {
        printf("Too many arguments supplied.\n");
     }
     else {
        printf("One argument expected.\n");
     }

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

    //specify the dielectric function // should be given as an input
    cx_mat eps_r = ones<cx_mat>(n,n);

    // map field_solutions = solveTM(double L0, double dL, double wvlen,
    //                 cx_mat eps_r, cx_mat Mz, vector<double> Npml);
    //
    // Hz = field_solutions["Hz"];
    // Ex = field_solutions["Hz"];
    // Ey = field_solutions["Hz"];


    return 0;


}

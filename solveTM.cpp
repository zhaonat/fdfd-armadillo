#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>

cx_mat solveTM(double L0, double wvlen, vector<double> xrange, vector<double> yrange,
                cx_mat eps_r, cx_mat Mz, vector<double> Npml){

    double L0 = 1e-6; //baseline units
    double eps0 = 8.85*1e-12*L0;
    double mu0 = 4*M_PI*1e-7*L0;
    double c0 = 1/sqrt(eps0*mu0);

    // it's up to us to figure out dL, Lx, Ly, and N
    double dL [2] = {0.01, 0.01};
    double omega = 2*M_PI*c0/(wvlen);  // angular frequency in rad/sec

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

    return X;

}

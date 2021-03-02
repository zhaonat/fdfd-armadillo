#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include "helpers.h"

map solveTM(double L0, double dL, double wvlen,
                cx_mat eps_r, cx_mat Mz, vector<double> Npml){


    double eps0 = 8.85*1e-12*L0;
    double mu0 = 4*M_PI*1e-7*L0;
    double c0 = 1/sqrt(eps0*mu0);
    // it's up to us to figure out dL, Lx, Ly, and N
    double omega = 2*M_PI*c0/(wvlen);  // angular frequency in rad/sec

    // get size of Mz
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

    cx_mat s_vector_x_f = create_sfactor(xrange,"f",omega,eps0, mu0,N[0],npml[0]);
    cx_mat s_vector_x_b = create_sfactor(xrange,"b",omega,eps0, mu0,N[0],npml[0]);
    cx_mat s_vector_y_f = create_sfactor(yrange,"f",omega,eps0, mu0,N[1],npml[1]);
    cx_mat s_vector_y_b = create_sfactor(yrange,"b",omega,eps0, mu0,N[1],npml[1]);

    // convert these to matrices 2d, not yet flattened
    cx_mat Sx_f_2D = repmat(s_vector_x_f, N[0],1);
    cx_mat Sy_f_2D = repmat(trans(s_vector_y_f), 1,N[1]);
    cx_mat Sx_b_2D = repmat(s_vector_x_b,  N[0],1);
    cx_mat Sy_b_2D = repmat(trans(s_vector_y_b), 1,N[1]);
    cout << Sy_f_2D.n_rows << " "<<Sy_f_2D.n_cols << endl;

    //convert these to the final flattened diagonal matrices
    sp_cx_mat Sxf = spdiags(m, vectorise(Sx_f_2D));
    sp_cx_mat Sxb = spdiags(m, vectorise(Sx_f_2D));
    sp_cx_mat Syf = spdiags(m, vectorise(Sx_f_2D));
    sp_cx_mat Syb = spdiags(m, vectorise(Sx_f_2D));


    sp_cx_mat A = Dxf*(Teps_inv)*Dxb + Dyf*(Teps_inv)*Dyb + omega*omega*Teps;

    // different types of solves... superlu
    cx_mat X = spsolve(A,b); //Hz

    // let's determine the other fields;
    hx = -1/(1i*omega)*(1/mu0)*Dyf)*X;
    hy = ((1/mu0)*Dxf)*X*(1/(1i*omega));
    // Hy = reshape(hy,N);
    // Hx = reshape(hx,N);

    //store fields as a dictionary?
    map<string, cx_mat> field_solutions;
    field_solutions['Ez'] = X;
    field_solutions['Hy'] = hx;
    field_solutions['Hx'] = hy;


    return field_solutions;

}

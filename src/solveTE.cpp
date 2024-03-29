//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include "derivatives.h"
#include "helpers.h"
#include <chrono>

// should we have .h files or what?
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

//#define ARMA_DONT_USE_WRAPPER


//what should inputs be
// int solveTE(double L0, double wvlen, vector<double> xrange, vector<double> yrange,
//                 cx_mat eps_r, cx_mat Mz, vector<double> Npml){
//
//     double eps0 = 8.854e-12 * L0;  // vacuum permittivity in farad/L0
//     double mu0 = pi * 4e-7 * L0;  // vacuum permeability in henry/L0
//     double c0 = 1/sqrt(eps0*mu0);  // speed of light in vacuum in L0/sec
//
//     N = size(eps_r);  // [Nx Ny]
//
//     L = [diff(xrange) diff(yrange)];  // [Lx Ly]
//     dL = L./N;  // [dx dy]
//     M = prod(N);
//
//     omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec
//
//     //%% Deal with the s_factor
//     //%% Deal with the s_factor
//     //% Create the sfactor in each direction and for 'f' and 'b'
//     s_vector_x_f = create_sfactor(xrange, 'f', omega, eps0, mu0, N(1), Npml(1));
//     s_vector_x_b = create_sfactor(xrange, 'b', omega, eps0, mu0, N(1), Npml(1));
//     s_vector_y_f = create_sfactor(yrange, 'f', omega, eps0, mu0, N(2), Npml(2));
//     s_vector_y_b = create_sfactor(yrange, 'b', omega, eps0, mu0, N(2), Npml(2));
//
//
//     // Fill the 2D space with layers of appropriate s-factors
//     Sx_f_2D = zeros(N[0],N[1]);
//     Sx_b_2D = zeros(N[0],N[1]);
//     Sy_f_2D = zeros(N[0],N[1]);
//     Sy_b_2D = zeros(N[0],N[1]);
//
//     for j = 1:N(2)
//         Sx_f_2D(:, j) = s_vector_x_f .^-1;
//         Sx_b_2D(:, j) = s_vector_x_b .^-1;
//     end
//
//     for i = 1:N(1)
//         Sy_f_2D(i, :) = s_vector_y_f .^-1;
//         Sy_b_2D(i, :) = s_vector_y_b .^-1;
//     end
//
//
//     // Reshape the 2D s-factors into a 1D s-array
//     Sx_f_vec = reshape(Sx_f_2D, M, 1);
//     Sx_b_vec = reshape(Sx_b_2D, M, 1);
//     Sy_f_vec = reshape(Sy_f_2D, M, 1);
//     Sy_b_vec = reshape(Sy_b_2D, M, 1);
//
//     // Construct the 1D total s-array into a diagonal matrix
//     // NO spdiags in armadillo...looksl like we have to use sp_cx_mat;
//     //Sxf = spdiags(Sx_f_vec, 0, M, M);
//     Sxf = spdiags(M, Sx_f_vec);
//     Sxb = spdiags(M, Sx_b_vec);
//     Syf = spdiags(M, Sy_f_vec);
//     Syb = spdiags(M, Sy_b_vec);
//
//     // Set up the permittivity and permeability in the domain.
//     // eps_x = bwdmean_w(eps0*eps_r, 'x');
//     // eps_y = bwdmean_w(eps0*eps_r, 'y');
//
//     // Setup the Teps_x, Teps_y, and Tmu_z matrices
//     eps_flatten = vectorise(eps_r);
//
//     // ============================================
//     T_eps_x = spdiags(vectorise(eps_r), 0, M, M);
//     T_eps_y = spdiags(vectorise(eps_r(:)), 0, M, M);
//
//     // Construct derivate matrices
//     sp_cx_mat Dyb = createDws("y",-1,dL,N);
//     sp_cx_mat Dxf = createDws("x",1,dL,N);
//     sp_cx_mat Dxb = createDws("x",-1,dL,N);
//     sp_cx_mat Dyf = createDws("y",1,dL,N);
//
//     // Dyb = Syb*createDws("y", "b", dL, N);
//     // Dxb = Sxb*createDws("x", "b", dL, N);
//     // Dxf = Sxf*createDws("x", "f", dL, N);
//     // Dyf = Syf*createDws("y", "f", dL, N);
//
//     // Reshape Mz into a vector
//     mz = reshape(Mz, M, 1);
//
//     //matrix multipication
//     // Construct A matrix and b vector
//     // % A = Sx_f *Dxf* T_eps_y^-1 *Sx_b *Dxb + Sy_f *Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z;
//     // % A = Sx_f*Dxf* T_eps_y^-1 *Sx_b*Dxb + Sy_f*Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z;
//
//     I = speye(M);
//     // assume we can invert diagonal matrices with no problem;
//     A_mode = Dxf*(T_eps_x.i())*Dxb + Dyf*(T_eps_y.i())*Dyb;
//     A = A_mode + omega^2*mu0*I;
//     % % A = Dxf* T_eps_y^-1 *Dxb + Dyf* T_eps_x^-1* Dyb + omega^2*T_mu_z;
//     b = 1i * omega * mz;
//
//
//     %% Solve the equation.
//     if all(b==0)
//         hz = zeros(size(b));
//     else
//         hz = spsolve(A,b);
//     end
//     Hz = reshape(hz, N);
//
//
//     ex = 1/(1i*omega) * T_eps_y.i() * Dyb * hz;
//     ey = 1/(1i*omega) * T_eps_x.i() * (-Dxb * hz);
//
//     Ex = reshape(ex, N);
//     Ey = reshape(ey, N);
//
//
// }

sp_cx_mat spdiags(int m, cx_mat cvals){
  // m: size of the matrix mxm
  // cvals: mx1 cx matrix of the vlaues to put on the diagonal
  // only supports main diagonal, which is enough for fdfd
  mat ind_cur = linspace(0,m-1,m);
  mat AB = join_rows(ind_cur, ind_cur).t();
  umat locations =  conv_to<umat>::from(AB);
  sp_cx_mat Teps(locations, cvals);
  return Teps;
}

int main(){

    double L0 = 1e-6; //baseline units
    int N [2] = {200,200};
    int npml [2] = {20,20};
    double dL [2] = {0.01, 0.01};
    double xrange [2] = {-1,1};
    double yrange [2] = {-1,1};

    int n = N[0];  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)
    double wvlen = 0.5; //microns
    double eps0 = 8.85*1e-12*L0;
    double mu0 = 4*M_PI*1e-7*L0;
    double c0 = 1/sqrt(eps0*mu0);
    double omega = 2*M_PI*c0/(wvlen);  // angular frequency in rad/sec

    //specify the dielectric function
    cx_mat eps_r = ones<cx_mat>(n,n);
    // specify source;
    cx_mat Mz = zeros<cx_mat>(n,n);
    int ind_src [2] = {n/2,n/2};
    cout << n/2 << endl;
    // for(int i = 0; i<n; i++){
    //   Mz(i,100) = 1; // simple point source
    // }
    Mz(n/2,n/2) = 1;
    cx_mat b = 1i*omega*vectorise(Mz,0);

    //make a diagonal matrix using sp_cx_mat...no diag function, whihc sucks
    mat ind_cur = linspace(0,m-1,m);

    mat AB = join_rows(ind_cur, ind_cur).t();
    umat locations =  conv_to<umat>::from(AB);
    cx_mat cvals = conv_to<cx_mat>::from(vectorise(eps0*eps_r,1));
    sp_cx_mat Teps(locations, cvals);
    sp_cx_mat Teps_inv(locations, 1/cvals);
    sp_cx_mat Dyb = createDws("y",-1,dL,N);
    sp_cx_mat Dxf = createDws("x",1,dL,N);
    sp_cx_mat Dxb = createDws("x",-1,dL,N);
    sp_cx_mat Dyf = createDws("y",1,dL,N);

    sp_cx_mat A = Dxf*(1/mu0)*Dxb + Dyf*(1/mu0)*Dyb + omega*omega*Teps;

    // different types of solves... superlu
    // is there anyway to time this line?
    cx_mat X = spsolve(A,b);

    // reshape to a field
    //cx_mat Hz = reshape(X,n,n);

    // cout << X << endl;

    // save this field into a data files
    //Hz.save("test_hz", hdf5_binary);
    X.save("test_hz.csv",csv_ascii);
    X.save("test_hz.h5",hdf5_binary);

    //cout<<Teps<<endl;
    //cout << diff(mat(N)) << endl;

    // add a pml
    cx_mat s_vector_x_f = create_sfactor(xrange,"f",omega,eps0, mu0,N[0],npml[0]);
    cx_mat s_vector_x_b = create_sfactor(xrange,"b",omega,eps0, mu0,N[0],npml[0]);
    cx_mat s_vector_y_f = create_sfactor(yrange,"f",omega,eps0, mu0,N[1],npml[1]);
    cx_mat s_vector_y_b = create_sfactor(yrange,"b",omega,eps0, mu0,N[1],npml[1]);
    s_vector_x_f.save("sfactorxf.csv",csv_ascii);
    s_vector_x_b.save("sfactorxb.csv",csv_ascii);
    s_vector_y_f.save("sfactoryf.csv",csv_ascii);
    s_vector_y_b.save("sfactoryb.csv",csv_ascii);

    // convert these to matrices 2d, not yet flattened
    //
    cx_mat Sx_f_2D = repmat(1/s_vector_x_f, N[0],1);
    cx_mat Sy_f_2D = repmat(1/strans(s_vector_y_f), 1,N[1]); //trans does hermitian
    cx_mat Sx_b_2D = repmat(1/s_vector_x_b,  N[0],1);
    cx_mat Sy_b_2D = repmat(1/strans(s_vector_y_b), 1,N[1]);
    // how do we check all these matrices were correctly made?

    //save these matrices
    Sx_f_2D.save("Sx_f_2D.csv",csv_ascii);
    Sy_f_2D.save("Sy_f_2D.csv",csv_ascii);
    Sx_b_2D.save("Sx_b_2D.csv",csv_ascii);
    Sy_b_2D.save("Sy_b_2D.csv",csv_ascii);


    cout << Sx_f_2D << endl;

    cout << Sy_f_2D.n_rows << " "<<Sy_f_2D.n_cols << endl;

    //convert these to the final flattened diagonal matrices
    //     Sxf = spdiags(M, Sx_f_vec);
    //     Sxb = spdiags(M, Sx_b_vec);
    //     Syf = spdiags(M, Sy_f_vec);
    //     Syb = spdiags(M, Sy_b_vec);
    //concatenate Sx_f_2D row-wise?
    sp_cx_mat Sxf = spdiags(m, vectorise(Sx_f_2D,1));
    sp_cx_mat Sxb = spdiags(m, vectorise(Sx_b_2D,1));
    sp_cx_mat Syf = spdiags(m, vectorise(Sy_f_2D,1));
    sp_cx_mat Syb = spdiags(m, vectorise(Sy_b_2D,1));

    // extract spdiags
    sp_cx_mat sxfdiag = Sxf.diag();
    sp_cx_mat syfdiag = Syf.diag();
    // convert back to dense
    cx_mat sxfd = strans(conv_to<cx_mat>::from(sxfdiag));
    cx_mat syfd = strans(conv_to<cx_mat>::from(syfdiag));
    cout << sxfd << endl;
    sxfd.save("sxf_diagonal.csv", csv_ascii);
    syfd.save("syf_diagonal.csv", csv_ascii);



    sp_cx_mat Dybs = Syb*Dyb;
    sp_cx_mat Dxfs = Sxf*Dxf;
    sp_cx_mat Dxbs = Sxb*Dxb;
    sp_cx_mat Dyfs = Syf*Dyf;

    // matrix with PML
    // A =  PEC_mask_y*PEC_mask_x*(Dxb*(Tmy^-1)*Dxf+ ...
    //     Dyb*(Tmx^-1)*Dyf)*PEC_mask_x*PEC_mask_y+ omega^2*Tepz;
    sp_cx_mat As = Dxbs*(1/mu0)*Dxfs + Dybs*(1/mu0)*Dyfs + omega*omega*Teps;

    cx_mat Xs = spsolve(As,b);
    // can't even save complex sparse matrices
    //As.save("A_matrix.csv", csv_ascii);

    // cout << Xs << endl;
    // cout << X<<endl;

    // save this field into a data files
    //Hz.save("test_hz", hdf5_binary);
    Xs.save("pml_test_hz.csv",csv_ascii);

    return 0;


}

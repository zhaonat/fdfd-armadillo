//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
# include "derivatives.cpp"

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
    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)
    double wvlen = 1.0; //microns

    //specify the dielectric function
    cx_mat eps_r = ones<cx_mat>(n,n);
    // specify source;
    cx_mat Mz = zeros<cx_mat>(n,n);
    cx_mat b = 1i*vectorise(Mz);

    //make a diagonal matrix using sp_cx_mat...no diag function, whihc sucks
    mat ind_cur = linspace(0,m-1,m);

    mat AB = join_rows(ind_cur, ind_cur).t();
    umat locations =  conv_to<umat>::from(AB);
    cx_mat cvals = conv_to<cx_mat>::from(vectorise(eps_r));
    sp_cx_mat Teps(locations, cvals);
    sp_cx_mat Teps_inv(locations, 1/cvals);
    sp_cx_mat Dyb = createDws("y",-1,dL,N);
    sp_cx_mat Dxf = createDws("x",1,dL,N);
    sp_cx_mat Dxb = createDws("x",-1,dL,N);
    sp_cx_mat Dyf = createDws("y",1,dL,N);

    sp_cx_mat A = Dxf*(Teps_inv)*Dxb + Dyf*(Teps_inv)*Dyb + omega^2*Teps;

    spsolve(A,b);

    //cout<<Teps<<endl;
    //cout << diff(mat(N)) << endl;

    return 0;


}

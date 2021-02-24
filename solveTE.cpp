//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>

// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

//#define ARMA_DONT_USE_WRAPPER


// what should inputs be
int solveTE(float L0, float wvlen, vector<float> xrange, vector<float> yrange, mat eps_r, mat Mz, vector<float> Npml){
    
    float eps0 = 8.854e-12 * L0;  // vacuum permittivity in farad/L0
    float mu0 = pi * 4e-7 * L0;  // vacuum permeability in henry/L0
    float c0 = 1/sqrt(eps0*mu0);  // speed of light in vacuum in L0/sec

    N = size(eps_r);  % [Nx Ny]
    L = [diff(xrange) diff(yrange)];  % [Lx Ly]
    dL = L./N;  % [dx dy]

    M = prod(N);

    omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

    %% Deal with the s_factor
    %% Deal with the s_factor
    % Create the sfactor in each direction and for 'f' and 'b'
    s_vector_x_f = create_sfactor(xrange, 'f', omega, eps0, mu0, N(1), Npml(1));
    s_vector_x_b = create_sfactor(xrange, 'b', omega, eps0, mu0, N(1), Npml(1));
    s_vector_y_f = create_sfactor(yrange, 'f', omega, eps0, mu0, N(2), Npml(2));
    s_vector_y_b = create_sfactor(yrange, 'b', omega, eps0, mu0, N(2), Npml(2));


    % Fill the 2D space with layers of appropriate s-factors
    Sx_f_2D = zeros(N);
    Sx_b_2D = zeros(N);
    Sy_f_2D = zeros(N);
    Sy_b_2D = zeros(N);

    for j = 1:N(2)
        Sx_f_2D(:, j) = s_vector_x_f .^-1;
        Sx_b_2D(:, j) = s_vector_x_b .^-1;
    end

    for i = 1:N(1)
        Sy_f_2D(i, :) = s_vector_y_f .^-1;
        Sy_b_2D(i, :) = s_vector_y_b .^-1;
    end

    % surf(abs(Sy_f_2D)); pause

    % Reshape the 2D s-factors into a 1D s-array
    Sx_f_vec = reshape(Sx_f_2D, M, 1);
    Sx_b_vec = reshape(Sx_b_2D, M, 1);
    Sy_f_vec = reshape(Sy_f_2D, M, 1);
    Sy_b_vec = reshape(Sy_b_2D, M, 1);

    % Construct the 1D total s-array into a diagonal matrix
    Sxf = spdiags(Sx_f_vec, 0, M, M);
    Sxb = spdiags(Sx_b_vec, 0, M, M);
    Syf = spdiags(Sy_f_vec, 0, M, M);
    Syb = spdiags(Sy_b_vec, 0, M, M);

    %% Set up the permittivity and permeability in the domain.
    eps_x = bwdmean_w(eps0*eps_r, 'x');
    eps_y = bwdmean_w(eps0*eps_r, 'y');

    % Setup the Teps_x, Teps_y, and Tmu_z matrices
    T_eps_x = spdiags(eps_x(:), 0, M, M);
    T_eps_y = spdiags(eps_y(:), 0, M, M);


    %% Construct derivate matrices
    Dyb = Syb*createDws('y', 'b', dL, N);
    Dxb = Sxb*createDws('x', 'b', dL, N);
    Dxf = Sxf*createDws('x', 'f', dL, N);
    Dyf = Syf*createDws('y', 'f', dL, N);

    %% Reshape Mz into a vector
    mz = reshape(Mz, M, 1);


    %% Construct A matrix and b vector
    % A = Sx_f *Dxf* T_eps_y^-1 *Sx_b *Dxb + Sy_f *Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z;
    % A = Sx_f*Dxf* T_eps_y^-1 *Sx_b*Dxb + Sy_f*Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z;

    %% =======================================================================%%
    %tic;[L,U,P,Q] = lu(A);toc; for umf pack
    %% =======================================================================%%

    I = speye(M);
    A_mode = Dxf*(T_eps_x^-1)*Dxb + Dyf*(T_eps_y^-1)*Dyb;
    A = A_mode + omega^2*mu0*I;
    % % A = Dxf* T_eps_y^-1 *Dxb + Dyf* T_eps_x^-1* Dyb + omega^2*T_mu_z;
    b = 1i * omega * mz;


    %% Solve the equation.
    if all(b==0)
        hz = zeros(size(b));
    else
        hz = A\b;
    end
    Hz = reshape(hz, N);


    ex = 1/(1i*omega) * T_eps_y^-1 * Dyb * hz;
    ey = 1/(1i*omega) * T_eps_x^-1 * (-Dxb * hz);

    Ex = reshape(ex, N);
    Ey = reshape(ey, N);

    
}

int main(){

    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)

    //
    int sign = 1;
    float dw = 1;

    return 0;


}

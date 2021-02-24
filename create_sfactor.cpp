//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>

int create_sfactor(){

    float lnR = -12.0;
    float m = 3.5;
    float eta0 = 377.0;

    w_array = linspace(wrange[0], wrange[1], Nw+1);
    loc_pml = [];
    //using vector or vec from armadillo?
    d_pml = abs(wrange - loc_pml);
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
}

int main(){

    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image

    float lnR = -12.0;
    float m = 3.5;
    float eta0 = 377.0;

    w_array = linspace(wrange[0], wrange[1], Nw+1);
    loc_pml = [];
    //using vector or vec from armadillo?
    d_pml = abs(wrange - loc_pml);

    return 0;


}

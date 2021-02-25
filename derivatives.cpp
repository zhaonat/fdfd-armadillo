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
mat circshift(mat A, int r, int c){
    return shift(arma::shift(A,r,0),c,1);
}

mat circshift2(mat A, int r, int c){
    return shift(arma::shift(A,c,0),r,1);
}
// what should inputs be
sp_cx_mat createDws(string w, int s, double* dL,int*  N){
   // w: either "x", "y", "z"
   // s: integer sign -1 or 1, -1 for backward, 1 for forward
   // dL is two dimensional;
   string xyz ("xyz");
   int idx;
   idx = xyz.find(w);
   int sign = s;

   int shift_inds [2];
   //process w
   if(w == "x"){
     shift_inds[0] = sign*1;
     shift_inds[1] = 0;
   }else{
     shift_inds[0] = 0;
     shift_inds[1] = sign*1;
   }

   float dw = dL[idx];
   int m = N[0]*N[1];
   mat ind_cur = linspace(0,m-1,m);
   mat ind_adj = linspace(0,m-1,m);

   ind_adj = reshape(ind_cur, N[0],N[1]);
   ind_adj = circshift(ind_adj, shift_inds[0],shift_inds[1]);
   //can I now do something like: roll or circshift

   cout<<ind_adj <<endl;

   ind_adj = reshape(ind_adj, m,1);

   mat off_diag = (sign/dw)*ones(m,1);
   mat on_diag  = -(sign/dw)*ones(m,1);

   mat x = join_cols(ind_cur, ind_cur);
   mat y = join_cols(ind_cur, ind_adj);
   mat AB = join_rows(x,y).t();
   umat locations =  conv_to<umat>::from(AB);

   mat vals = join_cols(off_diag, on_diag);
   cx_mat cvals = conv_to<cx_mat>::from(vals);

   sp_cx_mat dws(locations, cvals);

   return dws;
}

int main(){

    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int m = N[0]*N[1];  // number of unknows (=number of pixels)

    //cout << X << endl;
    sp_cx_mat Dyb = createDws("y",-1,dL,N);
    sp_cx_mat Dxf = createDws("x",1,dL,N);
    sp_cx_mat Dxb = createDws("x",-1,dL,N);
    sp_cx_mat Dyf = createDws("y",1,dL,N);

    cout << Dxf*Dxb << endl;

    return 0;

  }

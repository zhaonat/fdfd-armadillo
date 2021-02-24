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
// what should inputs be
sp_mat dws(int w,int s,vector<float> dL,vector<float> N){
   // w:
   // s: integer sign -1 or 1
   // dL is two dimensional;
   float dw = dL[w];
   int sign =s;
   int m = N[0]*N[1];
   mat ind_cur = linspace(1,m,m);
   mat ind_adj = linspace(1,m,m);
   ind_adj = reshape(ind_cur, N[0],N[1]);
   ind_adj = circshift(ind_adj, 1,0);
   //can I now do something like: roll or circshift
   ind_adj = reshape(ind_adj, m,1);

   mat off_diag = (sign/dw)*ones(m,1);
   mat on_diag  = -(sign/dw)*ones(m,1);

   mat x = join_cols(ind_cur, ind_cur);
   mat y = join_cols(ind_cur, ind_adj);
   mat AB = join_rows(x,y).t();
   umat locations =  conv_to<umat>::from(AB);
   mat vals = join_cols(off_diag, on_diag);

   sp_mat dws(locations, vals);

   return dws;
}

int main(){
//    mat A(4, 5, fill::randu);
//    mat B(4, 5, fill::randu);
//
//    cout << A*B.t() << endl;
    int N [2] = {20,20};
    double dL [2] = {0.1, 0.1};
    int n = 20;  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)

    //Mat<int> ind_adj(ind_cur);
    //ind_adj.reshape(N[0],N[1]);
    //ind_adj.t(); //transpose to row order
    //shift(ind_adj, +1);
    //cout << shift(ind_adj,1) <<" " << size(ind_adj)<< endl;
    //declare sparse matrix

    //
    int sign = 1;
    float dw = 1;
    // outputs a vector type
    mat ind_cur = linspace(1,m,m);
    mat ind_adj = linspace(1,m,m);
    //vec ind_adj = linspace<vec>(1,m,m);

    // ind_cur = reshape(ind_cur, n,n);
    ind_adj = reshape(ind_cur, n,n);
    ind_adj = circshift(ind_adj, 1,0);
    //can I now do something like: roll or circshift
    ind_adj = reshape(ind_adj, m,1);

    mat off_diag = (sign/dw)*ones(m,1);
    mat on_diag  = -(sign/dw)*ones(m,1);

    // what's a umat?
    // umat locations = { { 1, 7, 9 },
    //                { 2, 8, 9 } };
    //concatenate ind_cur and ind_adj;
    // sdf

    mat x = join_cols(ind_cur, ind_cur);
    mat y = join_cols(ind_cur, ind_adj);
    mat AB = join_rows(x,y).t();
    umat locations =  conv_to<umat>::from(AB);

    cout <<"shape: " << AB.n_rows <<" "<< AB.n_cols<< endl;

    //Dws = sparse([ind_cur;ind_cur], [ind_adj;ind_cur], [off_diag;on_diag]);
    mat vals = join_cols(off_diag, on_diag);
    // shape;
    cout << "vals shape: " << vals.n_rows << endl;

    //cout << << endl;

    //cout << locations << endl;

    sp_mat dws(locations, vals);
    cout << dws << endl;

    //write sp_mat to saveable format
    string filename = "test_dws";
    dws.save( filename, arma_binary);



    //cout << X << endl;


    sp_mat eye = speye<sp_mat>(5,5);


    return 0;

//    // we want to create matrices as we did in the matlab file
//    // create range from 1:M
//    //mat ind_adj = Mat<int>();
//    // convert c++ array to matrix;

//    // for(int i=0; i < ind_cur.size(); i++){
//    //   std::cout << ind_cur.at(i) << ' ';}
//    // ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
//
//    // How to reshape into matrix; ind_adj = reshape(ind_adj, N);
//
//    // ind_adj = circshift(ind_adj, -sign * ('xyz' == w));
//    // ind_adj = ind_adj(:);

}

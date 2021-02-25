//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include <complex>
#include <cmath>

//function avg_array = bwdmean_w(center_array, w)
// %% Input Parameters
// % center_array: 2D array of values defined at cell centers
// % w: 'x' or 'y', direction in which average is taken
//
// %% Out Parameter
// % avg_array: 2D array of averaged values
//
// center_shifted = circshift(center_array, 1*('xyz'==w)); %doe sthis generalize easily into 3 Dimensions, CHECK!
// avg_array = (center_shifted + center_array) / 2;
//
// end
// use armadillo, it directly maps matlab syntax
using namespace std;
using namespace arma;

// mat circshift(mat A, int r, int c){
//     return shift(shift(A,r,0),c,1);
// }

mat bwdmean_w(mat array, string w){
  string xyz ("xyz");
  int shift_param;
  shift_param = xyz.find(w);
  mat center_shifted = shift(array, shift_param); //doe sthis generalize easily into 3 Dimensions, CHECK!
  mat avg_array = (center_shifted + array) / 2;
  return avg_array;
}

int main(){
  //center_array = 2;
  string w = "x";
  int n = 100;
  vector<double> center_array(n,1.);
  string xyz ("xyz");
  int shift_param;
  shift_param = xyz.find(w);
  cout <<"shift: " <<shift_param << endl;
  mat ctr = mat(center_array);
  mat center_shifted = shift(ctr, shift_param); //doe sthis generalize easily into 3 Dimensions, CHECK!
  mat avg_array = (center_shifted + ctr) / 2;
  // cout << avg_array << endl;
  bwdmean_w(ctr, w);
  return 0;
}

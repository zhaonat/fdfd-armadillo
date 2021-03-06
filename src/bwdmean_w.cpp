//#include <armadillo>
#include <armadillo>
#include <vector>
#include <iostream>
#include <numeric> //std::iota
#include <algorithm>
#include <complex>
#include <cmath>

using namespace std;
using namespace arma;

mat bwdmean_w(mat array, string w){
  string xyz ("xyz");
  int shift_param;
  shift_param = xyz.find(w);
  mat center_shifted = shift(array, shift_param); //doe sthis generalize easily into 3 Dimensions, CHECK!
  mat avg_array = (center_shifted + array) / 2;
  return avg_array;
}

// int main(){
//   //center_array = 2;
//   string w = "x";
//   int n = 100;
//   vector<double> center_array(n,1.);
//   string xyz ("xyz");
//   int shift_param;
//   shift_param = xyz.find(w);
//   cout <<"shift: " <<shift_param << endl;
//   mat ctr = mat(center_array);
//   mat center_shifted = shift(ctr, shift_param); //doe sthis generalize easily into 3 Dimensions, CHECK!
//   mat avg_array = (center_shifted + ctr) / 2;
//   // cout << avg_array << endl;
//   bwdmean_w(ctr, w);
//   return 0;
// }

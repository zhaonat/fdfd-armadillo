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

//const double eps0;

//(wrange, s, omega, eps0, mu0, Nw, Nw_pml, lnR, m)

complex<double> pmlfunc(double L, double d_pml, double sigma_max, double omega, double eps0){
  // evaluate the pml scaling funciton of the complex conductivity
  //
  double m = 3.5;
  complex<double> ansr = 1. - 1i*sigma_max/(omega*eps0)*pow((L/d_pml),m);
  return ansr;
}

cx_mat create_sfactor(double * wrange,
                   string s,
                   double omega,
                   double eps0,
                   double mu0,
                   int Nw,
                   int Nw_pml){
   // special constants
   double lnR = -12.0;
   double m = 3.5;
   double eta0 = sqrt(mu0/eps0);


   mat w_array = linspace(wrange[0], wrange[1], Nw+1); // outputs a column vector?
   // specifies where the pml begins on each side in units of L0 or microns
   double loc_pml [2] = {w_array[Nw_pml], w_array[Nw-Nw_pml]};
   // we may need to do this manually;
   double d_pml [2] = {abs(wrange[0]-loc_pml[0]),abs(wrange[1]-loc_pml[1])};
   // double armadillo support scalar/vector??
   double sigma_max [2]= {-(m+1)*lnR/(2*eta0)/d_pml[0], -(m+1)*lnR/(2*eta0)/d_pml[1]};

   // we use ws to evaluate the pml coordinates
   mat ws; //forward or backward
   if(s == "b"){
     ws = w_array.rows(0,Nw-1);
   }else{
     ws = (w_array.rows(0,Nw-1)+w_array.rows(1,Nw))/2;
   }

   // extract all indices in w_array that needs to be evaluated;
   vector<double> w_left; //this is the corresponding thing in ind_pml
   vector<double> w_right;

   for(int i = 0; i< Nw; i++){
     if(ws[i] < loc_pml[0]){
       w_left.push_back(ws[i]);
     }
     if(ws[i] >= loc_pml[1]){
       w_right.push_back(ws[i]);
     }
   }
   // Construct the final sfactor array
   cx_mat sfactor_array = ones<cx_mat>(1, Nw);

   for(int i = 0; i < w_right.size(); i++){
     //cout << w_right[i] << endl;
     complex<double> tval = pmlfunc(abs(loc_pml[1]-w_right[i]), d_pml[1], sigma_max[1], omega, eps0);
     sfactor_array[Nw-Nw_pml+i] = tval;
   }

   // left side
   for(int i = 0; i < w_right.size(); i++){
     complex<double> tval = pmlfunc(abs(loc_pml[0]-w_left[i]), d_pml[0], sigma_max[0], omega, eps0);
     sfactor_array[i] = tval;
   }
  return sfactor_array;
}

// int main(){
//     /// it appears complex does not support float, only double
//     // parameters for pmlfunc
//     double L0 = 1e-6;
//
//     double eps0 = 8.85*1e-12*L0;
//     double mu0 = 4*M_PI*1e-7*L0;
//     double c0 = 1/sqrt(mu0*eps0);
//
//     double wvlen = 3.0;
//     double omega = 2*M_PI*c0/wvlen;
//
//     // =============================================
//     string s = "b";
//     double lnR = -12.0;
//     double m = 3.5;
//
//     double eta0 = sqrt(mu0/eps0);
//
//     double wrange [2] = {-1.0, 1.0}; // specified units of L0 or microns typically
//     int Nw = 100;
//     int Nw_pml = 10;
//
//     //would I prefer this being a vector? // it's Nw+1...because we will later
//     //cut it down
//     mat w_array = linspace(wrange[0], wrange[1], Nw+1); // outputs a column vector?
//     //using vector or vec from armadillo?
//     // d_pml = abs(wrange - loc_pml);
//     // sigma_max = -(m+1)*lnR/(2*eta0)/(d_mpl);
//
//     // specifies where the pml begins on each side in units of L0 or microns
//     double loc_pml [2] = {w_array[Nw_pml], w_array[Nw-Nw_pml]};
//
//     //d_pml = abs(wrange - loc_pml); % pml thickness
//
//     // we may need to do this manually;
//     double d_pml [2] = {abs(wrange[0]-loc_pml[0]),abs(wrange[1]-loc_pml[1])};
//     // mat test = 1./w_array;
//     cout << d_pml[0] <<" "<< d_pml[1] << endl;
//     // double armadillo support scalar/vector??
//     double sigma_max [2]= {-(m+1)*lnR/(2*eta0)/d_pml[0], -(m+1)*lnR/(2*eta0)/d_pml[1]};
//
//
//     // we use ws to evaluate the pml coordinates
//     mat ws;
//     if(s == "b"){
//       ws = w_array.rows(0,Nw-1);
//     }else{
//       ws = (w_array.rows(0,Nw-1)+w_array.rows(1,Nw))/2;
//     }
//
//     cout << ws.n_rows << ", "<<Nw <<endl;
//     cout << ws << endl;
//     // extract all indices in w_array that needs to be evaluated;
//     vector<double> w_left; //this is the corresponding thing in ind_pml
//     vector<double> w_right;
//
//     for(int i = 0; i< Nw; i++){
//       if(ws[i] < loc_pml[0]){
//         w_left.push_back(ws[i]);
//       }
//       if(ws[i] >= loc_pml[1]){
//         w_right.push_back(ws[i]);
//       }
//     }
//     cout <<"right: "<< mat(w_right) << endl;
//     //reverse(w_left.begin(),myvector.end())
//     cout <<w_left.size() << " " << w_right.size() << endl;
//     // Construct the final sfactor array
//     cx_mat sfactor_array = ones<cx_mat>(1, Nw);
//     // assign subset of sfactor array to the values denoted in.
//
//     //so some casting does indeed happen
//     // sfactor_array.cols(1,3) = {1,2,3};
//     //cout << sfactor_array << endl;
//
//     // let's evaluate the pml function on w_left and w_right;
//     // right side
//     for(int i = 0; i < w_right.size(); i++){
//       //cout << w_right[i] << endl;
//       complex<double> tval = pmlfunc(abs(loc_pml[1]-w_right[i]), d_pml[1], sigma_max[1], omega, eps0);
//       sfactor_array[Nw-Nw_pml+i] = tval;
//     }
//
//     // left side
//     for(int i = 0; i < w_right.size(); i++){
//       complex<double> tval = pmlfunc(abs(loc_pml[0]-w_left[i]), d_pml[0], sigma_max[0], omega, eps0);
//       cout << omega<<" "<< eps0<<" " << d_pml[0]<<" " << w_left[i]<<" " << endl;
//       cout <<tval << endl;
//       sfactor_array[i] = tval;
//     }
//     cout << sfactor_array << endl;
//
//     //test save of the sfactor;
//     sfactor_array.save("sfactor.csv",csv_ascii);
//
//     double xrange [2] = {-1.0, 1.0};
//     double yrange [2] = {-1.0, 1.0};
//
//     int N [2] = {100,100};
//     int npml [2] = {10,10};
//     //test creation of sfactor matrices;
//     cx_mat s_vector_x_f = create_sfactor(xrange,"f",omega,eps0, mu0,N[0],npml[0]);
//     cx_mat s_vector_x_b = create_sfactor(xrange,"b",omega,eps0, mu0,N[0],npml[0]);
//     cx_mat s_vector_y_f = create_sfactor(yrange,"f",omega,eps0, mu0,N[1],npml[1]);
//     cx_mat s_vector_y_b = create_sfactor(yrange,"b",omega,eps0, mu0,N[1],npml[1]);
//
//     s_vector_x_f.save("sfactorxf.csv",csv_ascii);
//     s_vector_x_b.save("sfactorxb.csv",csv_ascii);
//     s_vector_y_f.save("sfactoryf.csv",csv_ascii);
//     s_vector_y_b.save("sfactoryb.csv",csv_ascii);
//
//     // put s_vector into matrices
//
//
//     return 0;
//
//
// }

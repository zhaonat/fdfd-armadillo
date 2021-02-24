#include <iostream>
#include <Eigen/Dense>
// #include <Eigen/Sparse>
#include "include/add.h"
using namespace std;

using Eigen::MatrixXd;

int main() {
  cout << "Hello World" << endl; // call the function
  cout << add(2,2) <<endl;
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  cout << m << endl;
  return 0;
}

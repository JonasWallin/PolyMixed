#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

extern "C" {
  void F77_NAME(dchud)(
      double* R,
      int* ldr,
      int* p,
      double* x,
      double* z,
      int* ldz,
      int* nz,
      double* y,
      double* rho,
      double* c,
      double* s);
}



//' WoodburySpecial
//' Computes the inverse of D + UU^T where
//' D is a diagonal matrix
//' U is a low dimensional matrix
//'
//' @param D is a vector of the diagonal of D
//' @param U is the low rank matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd WoodburySpecial(Eigen::VectorXd D, Eigen::MatrixXd U){

  Eigen::VectorXd Dinv = D.cwiseInverse();
  Eigen::MatrixXd Uh = Dinv.asDiagonal() * U;
  Eigen::MatrixXd UhTu = Uh.transpose() * U;
  UhTu.diagonal().array() += 1;
  Eigen::MatrixXd innerMatrix = UhTu.inverse();
  Eigen::MatrixXd res = -(Uh * innerMatrix) * Uh.transpose();
  res.diagonal() += Dinv;
  return(res);
}


//' Update a Cholesky factor R of Q=R'R
//' to R of Q+XX'= R'R
//'
//' @param R upper triangel matrix
//' @param x is a matrix
//'
// [[Rcpp::export]]
Eigen::MatrixXd CholUpdate(Eigen::MatrixXd R, Eigen::MatrixXd U){

  int p = R.rows();
  int m = U.cols();
  Eigen::VectorXd s,c;
  c.setZero(p);
  s.setZero(p);
  int ldz = 0;
  int nz = 0;
  for(int i =0; i < m; i++){
    dchud_(&R(0,0), &p ,&p, &U(0,i), 0,
           &ldz,
           &nz,
           0,
           0,
           &c(0),
           &s(0));
  }
         return(R);
}


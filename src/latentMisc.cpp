/*
	misc functions for latent class



*/


#include <RcppEigen.h>
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]


using Eigen::MatrixXd;
using Eigen::VectorXd;


	// [[Rcpp::export]]
Eigen::VectorXd getTValues(Eigen::VectorXi& index,
						   Eigen::MatrixXd& Xu,
						   Eigen::MatrixXd& uX,
					       Eigen::MatrixXd& Sigma_U,
					       Eigen::VectorXd& Y){

  int p_i = index.size();
	int p   = uX.cols();
  int np  = Xu.cols();
	Eigen::VectorXd tValues(p);
	tValues.setZero();
	Eigen::MatrixXd Xu_i(Xu.rows(), Xu.cols() + 1);
	Xu_i.block(0,0, Xu.rows(), Xu.cols()) =  Xu;
	for(int j = 0; j < p_i; j++){
	  int i = index(j) - 1;
		Xu_i.col(Xu.cols())      = uX.col(i);
		Eigen::MatrixXd  SinvX ;
	  SinvX   = Sigma_U.llt().solve(Xu_i);
	  Eigen::MatrixXd  XtSinvX = Xu_i.transpose() * SinvX;
	  Eigen::MatrixXd CovBeta  = XtSinvX.inverse();
	  Eigen::VectorXd beta_est = CovBeta * ( Xu_i.transpose() * Y);
	  double cb = CovBeta(np, np);
	  tValues(i)               = beta_est(np)/sqrt(cb);
	}
    return(tValues);
}

// [[Rcpp::export]]
Eigen::VectorXd getTValuesDiagS(Eigen::VectorXi& index,
                           Eigen::MatrixXd& Xu,
                           Eigen::MatrixXd& uX,
                           Eigen::VectorXd& Sigma_U,
                           Eigen::VectorXd& Y){

  int p_i = index.size();
  int p   = uX.cols();
  int np  = Xu.cols();
  Eigen::VectorXd tValues(p);
  tValues.setZero();
  Eigen::MatrixXd Xu_i(Xu.rows(), Xu.cols() + 1);
  Xu_i.block(0,0, Xu.rows(), Xu.cols()) =  Xu;
  for(int j = 0; j < p_i; j++){
    int i = index(j) - 1;
    Xu_i.col(Xu.cols())      = uX.col(i);
    Eigen::MatrixXd  SinvX ;
    SinvX   =  Sigma_U.asDiagonal().inverse() * Xu_i;
    Eigen::MatrixXd  XtSinvX = Xu_i.transpose() * SinvX;
    Eigen::MatrixXd CovBeta  = XtSinvX.inverse();
    Eigen::VectorXd beta_est = CovBeta * ( Xu_i.transpose() * Y);
    double cb = CovBeta(np, np);
    tValues(i)               = beta_est(np)/sqrt(cb);
  }
  return(tValues);
}

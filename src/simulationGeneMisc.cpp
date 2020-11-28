#include <Rcpp.h>

#include <Rcpp/sugar/functions/sample.h>
using namespace Rcpp;
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


/*
 *
 *  Code for random splitting crossing of Genes
 *  X1    - ( N x n) gene 1 for indivuals
 *  X2    - ( N x n) gene 2 for indivuals
 *  cross - (N  x1) number of crosses for each indivual
 *
 */
// [[Rcpp::export]]
void crossSelect(Eigen::Map<Eigen::MatrixXd> X1,
                 Eigen::Map<Eigen::MatrixXd> X2,
                 IntegerVector cross) {

  int N = X1.rows();
  int n = X1.cols();

  Eigen::Map<Eigen::MatrixXd> X2t = X2;
  IntegerVector index_n = seq(1, n-2);
  for(int i =0; i < N; i++){
    IntegerVector ind;
    ind = sample(index_n,cross[i]);
    std::sort(ind.begin(), ind.end());
    ind.push_back(n);
    ind.push_front(0);
    int test =  R::runif(0,1)>0.5;
    for(int ii=1; ii < ind.length(); ii++){
      if(test){
        X2.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = X1.block(i,ind(ii-1),1,ind(ii)-ind(ii-1));
        X1.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = X2t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1));
        test = 0;
      }else{
        test = 1;
      }
    }
  }
}



/*
 *
 *  Code for random splitting crossing of Genes
 *  X1    - ( N x n) gene 1 for indivuals
 *  X2    - ( N x n) gene 2 for indivuals
 *  Z1    - ( N x n) ancesteror gene 1 for indivuals
 *  Z2    - ( N x n) ancesteror gene 2 for indivuals
 *  n_gen - (1  x1) number of genrations
 *
 */
// [[Rcpp::export]]
void mixing_population(Eigen::Map<Eigen::MatrixXd> X1,
                 Eigen::Map<Eigen::MatrixXd> X2,
                 Eigen::Map<Eigen::MatrixXi> Z1,
                 Eigen::Map<Eigen::MatrixXi> Z2,
                 int n_gen,
                 double lambda_cross) {

  int N = X1.rows();
  int n = X1.cols();
  IntegerVector index_n = seq(1, n-2);
  Eigen::Map<Eigen::MatrixXd> X1t = X1;
  Eigen::Map<Eigen::MatrixXd> X2t = X2;
  Eigen::Map<Eigen::MatrixXi> Z1t = Z1;
  Eigen::Map<Eigen::MatrixXi> Z2t = Z2;
  //IntegerVector index_N = seq(0, N-1);
  //IntegerVector ind_N = sample(index_N,N, true);
  NumericVector Cross = Rcpp::rpois( N, lambda_cross );


  Eigen::RowVectorXd XM(n);
  Eigen::RowVectorXd XF(n);
  Eigen::RowVectorXi ZM(n);
  Eigen::RowVectorXi ZF(n);
  for(int g = 0; g < n_gen; g++){
    for(int i = 0; i < N; i++){

      int i_1 = floor(N*R::runif(0,1));
      int i_2 = floor(N*R::runif(0,1));

      // select mother
      if(R::runif(0,1)>0.5){
        XM = X1.row(i_1);
        ZM = Z1.row(i_1);
      }else{
        XM = X2.row(i_1);
        ZM = Z2.row(i_1);
      }
      // select father
      if(R::runif(0,1)>0.5){
        XF = X1.row(i_2);
        ZF = Z1.row(i_2);
      }else{
        XF = X2.row(i_2);
        ZF = Z2.row(i_2);
      }

      IntegerVector ind;
      ind = sample(index_n, (int)std::min(Cross[i], (double)index_n.length()));
      std::sort(ind.begin(), ind.end());
      ind.push_back(n);
      ind.push_front(0);
      int test =  R::runif(0,1)>0.5;
      for(int ii=1; ii < ind.length(); ii++){
        if(test){
          X1t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = XF.segment(ind(ii-1), ind(ii)-ind(ii-1));
          X2t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = XM.segment(ind(ii-1), ind(ii)-ind(ii-1));
          Z1t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = ZF.segment(ind(ii-1), ind(ii)-ind(ii-1));
          Z2t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = ZM.segment(ind(ii-1), ind(ii)-ind(ii-1));
          test = 0;
        }else{
          X1t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = XM.segment(ind(ii-1), ind(ii)-ind(ii-1));
          X2t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = XF.segment(ind(ii-1), ind(ii)-ind(ii-1));
          Z1t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = ZM.segment(ind(ii-1), ind(ii)-ind(ii-1));
          Z2t.block(i,ind(ii-1),1,ind(ii)-ind(ii-1)) = ZF.segment(ind(ii-1), ind(ii)-ind(ii-1));
          test = 1;
        }
      }

    }
    X1 = X1t;
    X2 = X2t;
    Z1 = Z1t;
    Z2 = Z2t;
  }
}




#ifndef ACEDLNMCPPADCLASS_HPP
#define ACEDLNMCPPADCLASS_HPP


// include CppAD
#include <cppad/cppad.hpp>

// autodiff include
// from https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
// #include "autodiff/common/meta.hpp"
// #include "autodiff/common/numbertraits.hpp"
// #include "autodiff/common/binomialcoefficient.hpp"
// #include "autodiff/common/vectortraits.hpp"
// #include "autodiff/forward/dual/dual.hpp"
// #include "autodiff/common/eigen.hpp"
// #include "autodiff/common/classtraits.hpp"
// #include "autodiff/forward/utils/derivative.hpp"
// #include "autodiff/forward/real/real.hpp"
// #include "autodiff/forward/utils/taylorseries.hpp"
// #include "autodiff/forward/real.hpp"
// #include "autodiff/forward/real/eigen.hpp"
// #include "autodiff/forward/utils/gradient.hpp"
// using namespace autodiff;



#include "defheader.h"


#include <cmath>

#include <iostream>
using namespace std;



// Eigen::VectorXd convertToDouble(const Vec input) {
//   int n = input.size();
//   Eigen::VectorXd output(n);

//   for (int i = 0; i < n; ++i) {
//       output(i) = CppAD::Value(input(i));
//   }
//   return output;
// }


// Eigen::MatrixXd convertToDouble(const Mat input) {
//   int nrow = input.rows();
//   int ncol = input.cols();
//   Eigen::MatrixXd output(nrow, ncol);

//   for (int i = 0; i < nrow; ++i) {
//       for (int j = 0; j < ncol; ++j) {
//           output(i, j) = CppAD::Value(input(i, j));
//       }
//   }

//   return output;
// }


// ************** PART 1: Define some functions **********************


// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindex(Scalar x,const Vec t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

Scalar weight(Scalar x, const Vec& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


Scalar Bspline(Scalar x, int j, const Vec& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  Scalar w1 = weight(x,t,j,p-1);
  Scalar w2 = weight(x,t,j+1,p-1);
  Scalar b1 = Bspline(x,j,t,p-1);
  Scalar b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

Scalar Bspline1st(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline(x,j+1,t,p-1);
  Scalar bb2 = Bspline(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Scalar Bspline2nd(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline1st(x,j+1,t,p-1);
  Scalar bb2 = Bspline1st(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Vec Bsplinevec(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec1st(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon1st(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec2nd(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}


Vec BsplinevecCon2nd(Scalar x, const Vec& t,int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}



// Lanczos approximation
// Source code from https://github.com/brianmartens/BetaFunction/blob/master/BetaFunction/bmath.h
Scalar lanczos_lgamma(Scalar z) {
    const Scalar LG_g = 7.0;
    const int LG_N = 9;

    const Scalar ln_sqrt_2_pi = 0.91893853320467274178;

    Vec lct(LG_N+1);
    lct << 0.9999999999998099322768470047347,
    676.520368121885098567009190444019,
   -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
     12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7;

    Scalar sum;
    Scalar base;

    // To avoid if condition z < 0.5, we calculate gamma(z+1) which is equal to z*gamma(z).
    // WAS:
    // z = z - 1.0;
    // if (z < 0.5) {
    //   // Use Euler's reflection formula:
    //   // Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
    //   out = log(g_pi / sin(g_pi * z)) - lanczos_lgamma(1.0 - z);
    //   return out;
    // }
    // gamma(z) ...

    // New: indeed gamma(z+1)
    base = z + LG_g + 0.5;  // Base of the Lanczos exponential
    sum = 0;
    // We start with the terms that have the smallest coefficients and largest
    // denominator.
    for(int i=LG_N; i>=1; i--) {
      sum += lct[i] / (z + ((double) i));
    }
    sum += lct[0];
    Scalar gammazplus1 = ((ln_sqrt_2_pi + log(sum)) - base) + log(base)*(z+0.5);
    Scalar out = gammazplus1 - log(z);

    return out;
}

// 1st derivative of log gamma function
// https://math.stackexchange.com/questions/481253/differentiate-log-gamma-function
// VERY SLOW!!!
// Scalar lgamma1st (Scalar z) {
//   // Euler's constant https://en.wikipedia.org/wiki/Euler%27s_constant#
//   const Scalar Euler = 0.57721566490153286060651209008240243104215933593992;
//   Scalar out = -1.0 * Euler;
//   int K = 1e6;
//   for (int i = 1; i < K; i++) {
//     out += 1.0/i - 1.0/(i+z-1.0);
//   }
//   std::cout << "z" << CppAD::Value(z) << std::endl;
//   std::cout << "lgamma1st" << CppAD::Value(out) << std::endl;
//   return out;
// }



// https://github.com/tminka/lightspeed/blob/master/digamma.m
Scalar lgamma1st (Scalar x) {
  const Scalar pi = 3.141592653589793238462643383279;
  const Scalar large = 9.5;
  const Scalar d1 = -0.5772156649015328606065121;
  const Scalar d2 = pi*pi/6.0;
  const Scalar small = 1e-6;
  const Scalar s3 = 1.0/12.0;
  const Scalar s4 = 1.0/120.0;
  const Scalar s5 = 1.0/252.0;
  const Scalar s6 = 1.0/240.0;
  const Scalar s7 = 1.0/132.0;
  const Scalar s8 = 691.0/32760.0;
  const Scalar s9 = 1.0/12.0;
  const Scalar s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  Scalar xplus10 = x + 10.0;
  Scalar y = 0.0;
  Scalar r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

Scalar lgamma2nd (Scalar x) {
  const Scalar pi = 3.141592653589793238462643383279;
  const Scalar c = pi*pi/6;
  const Scalar c1 = -2.404113806319188570799476;
  const Scalar b2 =  1.0/6.0;
  const Scalar b4 = -1.0/30.0;
  const Scalar b6 =  1.0/42.0;
  const Scalar b8 = -1.0/30.0;
  const Scalar b10 = 5.0/66.0;

  // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

  Scalar z = 1./(x*x);
  Scalar y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

  // std::cout << "x" << CppAD::Value(x) << std::endl;
  // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
  return y;
}

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Modelcppad {
  // The DLNM model

private:
  // DATA
  const Vec y; // Response
  const Mat Sw_large; // penalty matrix for w(l)
  const Mat Sf; // penalty matrix for f(E)
  const std::vector<Mat> B_inner_list;
  const Vec knots_f; // knots for f(E) B-spline
  const Mat Dw;  // \int w(l)^2 dl = 1

  const Mat Xfix; // fixed effects
  const Mat Xrand; // random effects
  const Vec r; // rank of each smooth

  const Mat Zf;

  const Vec Xoffset; // offset
  
public:
  const Mat K;
  const Vec a;
  const Mat Bindex; 
  
  int n;
  int kw;
  int kwopt;
  int kE;
  int kbetaR;
  int kbetaF;
  int mE;
  int p; // number of smooth terms in Xrand

  // PARAMETERS
  Vec alpha_f;
  Vec phi;
  Vec con_index_par; // unconstrained index parameter
  Vec phiKa;
  Scalar log_theta;
  Scalar log_smoothing_f;
  Scalar log_smoothing_w;

  Vec betaF; // parameters for fixed effects
  Vec betaR; // parameters for random effects
  Vec logsmoothing; // log smoothing parameters for random effects

  // Components generated
  Mat B_inner; // B_inner = sum_m gamma_m B_inner_list.at(m)
  Scalar theta;
  Scalar smoothing_f;
  Scalar smoothing_w;
  Vec smoothing;
  Vec phi_long;
  Scalar alpha_w_C_denominator;
  Vec alpha_w_C;
  Vec alpha_w_C_pen;
  Vec con_index_par_long;
  Scalar index_par_denominator;
  Vec index_par; // index parameter = con_index_par_long/index_par_denominator
  Mat Bf_matrix;
  Vec E;
  Vec eta;
  Vec eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Vec mu; // log(mu) = eta + eta_remaining + Xoffset
  Mat BtBindex;

  // Components for derivatives

  Mat dlogmu_dw_mat;
  Mat dlogmu_dindex_par_mat;
  Vec dlogdensity_dmu_vec;
  Mat dmu_dw_mat;
  Mat dmu_dindex_par_mat;
  Mat dw_dphi_mat;
  Mat dindex_par_dcon_index_par_mat;
  
  Vec Bf;

  Vec d2logdensity_dmudmu_vec;
  std::vector<Mat> d2w_dphidphi_list;
  std::vector<Mat> d2index_par_dcon_index_pardcon_index_par_list;
  Mat he_alpha_w_mat;
  Mat he_index_par_mat;
  Mat he_alpha_f_alpha_w_mat;
  Mat he_alpha_f_index_par_mat;


  // gradient and hessian for updating alpha_f, betaR and betaF
  Mat he_inner_mat;



  // full hessian
  Mat he_alpha_f_mat;
  Mat he_betaR_mat;
  Mat he_betaF_mat;
  Mat he_phi_mat;
  Mat he_con_index_par_mat;
  Mat he_alpha_f_phi_mat;
  Mat he_alpha_f_con_index_par_mat;
  Mat he_alpha_f_betaF_mat;
  Mat he_alpha_f_betaR_mat;
  Mat he_phi_betaF_mat;
  Mat he_phi_betaR_mat;
  Mat he_betaR_betaF_mat;
  Mat he_con_index_par_betaF_mat;
  Mat he_con_index_par_betaR_mat;
  Mat he_phi_index_par_mat;
  Mat he_alpha_w_index_par_mat;
  Mat he_phi_con_index_par_mat;

  Mat he_s_u_mat;
  
  // for AIC
  Vec he_alpha_f_log_theta_vec;
  Vec he_phi_log_theta_vec;
  Vec he_betaR_log_theta_vec;
  Vec he_con_index_par_log_theta_vec;
  Vec he_betaF_log_theta_vec;
  Vec d2logdensity_dmudtheta_vec;
  Scalar he_log_theta_scalar;
  Mat I_mat;
  Mat I_alpha_f_mat;
  Mat I_betaR_mat;
  Mat I_alpha_w_mat;
  Mat I_phi_mat;



  Mat K_alpha_f_mat;
  Mat K_betaR_mat;
  Mat K_betaF_mat;
  Mat K_alpha_w_mat;
  Mat K_phi_mat;
  Mat K_index_par_mat;
  Mat K_con_index_par_mat;
  Vec K_log_theta_vec;

  Mat Kleft;
  Mat Khat; // Khat = Kleft * Kleft.transpose()


  // AD tape for LAML
  // CppAD::ADFun<double> gr;
  bool ifhastape = false;

  // Constructor
  Modelcppad(const Vec& y_,
            const std::vector<Mat>& B_inner_list_,
            const Vec& knots_f_,
            const Mat& Sw_large_,
            const Mat& Sf_,
            const Mat& Dw_,
            const Mat& Xrand_,
            const Mat& Xfix_,
            const Mat& Zf_,
            const Vec& Xoffset_,
            const Vec& r_,
            const Mat& K_,
            const Vec& a_,
            const Mat& Bindex_,
            Vec& alpha_f_,
            Vec& phi_,
            Vec& con_index_par_,
            Scalar log_theta_,
            Scalar log_smoothing_f_,
            Scalar log_smoothing_w_,
            Vec& betaR_,
            Vec& betaF_,
            Vec& logsmoothing_) :
    y(y_), B_inner_list(B_inner_list_), knots_f(knots_f_), Sw_large(Sw_large_), Sf(Sf_), Dw(Dw_), Xrand(Xrand_), Xfix(Xfix_), Zf(Zf_), Xoffset(Xoffset_), r(r_),
    K(K_), a(a_), Bindex(Bindex_),
    alpha_f(alpha_f_), phi(phi_), con_index_par(con_index_par_), log_theta(log_theta_), log_smoothing_f(log_smoothing_f_), log_smoothing_w(log_smoothing_w_), betaR(betaR_), betaF(betaF_), logsmoothing(logsmoothing_) {

      n = y.size(); // sample size
      kw = a.size() + 1;
      kwopt = K.cols(); // conL = TRUE: kwopt = kw-2; conL = FALSE: kwopt = kw-1
      kE = alpha_f.size();
      kbetaR = betaR.size();
      kbetaF = betaF.size();
      mE = B_inner_list.size(); // number of exposures, which is the length of
      p = r.size();

      B_inner.resize(n, kw);
      B_inner.setZero();

      theta = exp(log_theta);
      smoothing_f = exp(log_smoothing_f);
      smoothing_w = exp(log_smoothing_w);

      smoothing.resize(p);
      for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));

      phiKa = K * phi + a;

      phi_long.resize(kw); // phi_long = c(1, phiKa)
      phi_long(0) = 1.0;
      for (int j = 0; j < (kw - 1); j++) {
        phi_long(j + 1) = phiKa(j);
      }
      alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));
      alpha_w_C = phi_long / alpha_w_C_denominator;
      alpha_w_C_pen = phi / alpha_w_C_denominator;

      
      BtBindex = Bindex.transpose() * Bindex;

      con_index_par_long.resize(mE); // con_index_par_long = c(1,con_index_par)
      con_index_par_long(0) = 1.0;
      for (int i = 0; i < (mE - 1); i++) {
        con_index_par_long(i + 1) = con_index_par(i);
      }
      index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
      index_par = Bindex * con_index_par_long/index_par_denominator;

      for (int i = 0; i < mE; i++) {
        B_inner += B_inner_list.at(i) * index_par(i);
      }

      E = B_inner * alpha_w_C;
      

      Bf_matrix.resize(n, kE);
      eta.resize(n);
      eta_remaining.resize(n);
      mu.resize(n);
      
      for (int i = 0; i < n; i++) {
        Bf = BsplinevecCon(E(i), knots_f, 4, Zf);
        Bf_matrix.row(i) = Bf;
        eta(i) = Bf.dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }

      // Initialize the derivative components and NegativeLogLikelihood
      dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
      d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi
      

      dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
      d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();


      he_s_u_mat.resize(kwopt+kE+mE-1+kbetaR+kbetaF, kwopt+kE+mE-1+kbetaR+kbetaF);

      I_mat.resize(kwopt+kE+mE-1+kbetaR+kbetaF + 1, kwopt+kE+mE-1+kbetaR+kbetaF + 1);

      derivative_coef();
      derivative_he();
    }

  // Functions to set parameters
  void setAlphaF(const Vec alpha_f_) {
    alpha_f = alpha_f_;
  }

  void setPhi(const Vec phi_) {
    phi = phi_;
  }
  void setCon_Index_Par(const Vec con_index_par_){
    con_index_par = con_index_par_;
  }
  void setBetaF(const Vec betaF_) {
    betaF = betaF_;
  }
  void setBetaR(const Vec betaR_) {
    betaR = betaR_;
  }
  void setLogTheta(const Scalar log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Scalar log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    smoothing_f = exp(log_smoothing_f);
  }

  void setLogSmoothingW(const Scalar log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    smoothing_w = exp(log_smoothing_w);
  }
  void setLogsmoothing(const Vec logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }


  void derivative_coef() {
    phiKa = K * phi + a;
    for (int j = 0; j < (kw - 1); j++) {
      phi_long(j + 1) = phiKa(j);
    }
    alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));
    alpha_w_C = phi_long / alpha_w_C_denominator;
    alpha_w_C_pen = phi / alpha_w_C_denominator;
    

    for (int i = 0; i < (mE - 1); i++) {
      con_index_par_long(i + 1) = con_index_par(i);
    }
    index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
    index_par = Bindex * con_index_par_long/index_par_denominator;
    B_inner.setZero();
    for (int i = 0; i < mE; i++) {
      B_inner += B_inner_list.at(i) * index_par(i);
    }
    
    
    E = B_inner * alpha_w_C;


    for (int i = 0; i < n; i++) {
      Bf = BsplinevecCon(E(i), knots_f, 4, Zf);
      Bf_matrix.row(i) = Bf;
      eta(i) = Bf.dot(alpha_f);
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }

    dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
    d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi

    dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
    d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();


    dlogmu_dw_mat = dlogmu_dw();
    dlogmu_dindex_par_mat = dlogmu_dindex_par();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_dw_mat = dmu_dw();
    dmu_dindex_par_mat = dmu_dindex_par();
    
  }

  // update full gradient and hessian of alpha_f, phi and betaR and betaF
  void derivative_he () {
    

    // obtain hessian
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
   
    he_alpha_w_mat = he_alpha_w();
    he_index_par_mat = he_index_par();

    he_alpha_f_alpha_w_mat = he_alpha_f_alpha_w();
    he_alpha_f_index_par_mat = he_alpha_f_index_par();

    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_phi_mat = he_phi();
    he_con_index_par_mat = he_con_index_par();

    he_alpha_f_phi_mat = he_alpha_f_phi();
    he_alpha_f_con_index_par_mat = he_alpha_f_con_index_par();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_phi_betaF_mat = he_phi_betaF();
    he_phi_betaR_mat = he_phi_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();
    he_con_index_par_betaF_mat = he_con_index_par_betaF();
    he_alpha_w_index_par_mat = he_alpha_w_index_par();
    he_phi_index_par_mat = he_phi_index_par();
    he_phi_con_index_par_mat = he_phi_con_index_par();
    he_con_index_par_betaR_mat = he_con_index_par_betaR();


    he_s_u_mat.setZero();

    
    he_s_u_mat.block(0, 0, kE, kE)  = he_alpha_f_mat;
    he_s_u_mat.block(0, kE, kE, kwopt) = he_alpha_f_phi_mat;
    he_s_u_mat.block(kE, 0, kwopt, kE) = he_alpha_f_phi_mat.transpose();
    he_s_u_mat.block(kE, kE, kwopt, kwopt) = he_phi_mat;


    he_s_u_mat.block(0, kE+kwopt, kE, mE-1) = he_alpha_f_con_index_par_mat;
    he_s_u_mat.block(kE+kwopt, 0, mE-1, kE) = he_alpha_f_con_index_par_mat.transpose();
    he_s_u_mat.block(kE, kE+kwopt, kwopt, mE-1) = he_phi_con_index_par_mat;
    he_s_u_mat.block(kE+kwopt, kE, mE-1, kwopt) = he_phi_con_index_par_mat.transpose();
    he_s_u_mat.block(kE+kwopt, kE+kwopt, mE-1, mE-1) = he_con_index_par_mat;


    he_s_u_mat.block(kE+kwopt+mE-1, kE+kwopt+mE-1, kbetaR, kbetaR) = he_betaR_mat;
    he_s_u_mat.block(kE+kwopt+kbetaR+mE-1, kE+kwopt+kbetaR+mE-1, kbetaF, kbetaF) = he_betaF_mat;

    he_s_u_mat.block(0,kE+kwopt+mE-1,kE,kbetaR) = he_alpha_f_betaR_mat;
    he_s_u_mat.block(kE,kE+kwopt+mE-1,kwopt,kbetaR) = he_phi_betaR_mat;
    he_s_u_mat.block(kE+kwopt,kE+kwopt+mE-1,mE-1,kbetaR) = he_con_index_par_betaR_mat;
    he_s_u_mat.block(0,kE+kwopt+mE-1+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_s_u_mat.block(kE,kE+kwopt+mE-1+kbetaR,kwopt,kbetaF) = he_phi_betaF_mat;
    he_s_u_mat.block(kE+kwopt,kE+kwopt+mE-1+kbetaR,mE-1,kbetaF) = he_con_index_par_betaF_mat;
    he_s_u_mat.block(kE+kwopt+mE-1, kE+kwopt+mE-1+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_mat;

    he_s_u_mat.block(kE+kwopt+mE-1,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1,kE,kbetaR,kwopt) = he_phi_betaR_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1,kE+kwopt,kbetaR,mE-1) = he_con_index_par_betaR_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1+kbetaR,kE,kbetaF,kwopt) = he_phi_betaF_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1+kbetaR,kE+kwopt,kbetaF,mE-1) = he_con_index_par_betaF_mat.transpose();
    he_s_u_mat.block(kE+kwopt+mE-1+kbetaR, kE+kwopt+mE-1, kbetaF, kbetaR) = he_betaR_betaF_mat.transpose();

  }

  void prepare_AIC () {
    // I_alpha_f_mat = I_alpha_f();
    // I_betaR_mat = I_betaR();
    I_alpha_w_mat = I_alpha_w();
    I_phi_mat = I_phi();

    Mat Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    Ones.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      double ri = CppAD::Value(r(i));
      int ki = static_cast<int>(ri);
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }

    I_mat.setZero();
    I_mat.block(0,0,kwopt+kE+mE-1+kbetaR+kbetaF, kwopt+kE+mE-1+kbetaR+kbetaF) = he_s_u_mat;
    I_mat.block(0, 0, kE, kE) = he_s_u_mat.block(0, 0, kE, kE) - smoothing_f*Sf;
    I_mat.block(kE, kE, kwopt, kwopt) = I_phi_mat;
    I_mat.block(kE+kwopt+mE-1, kE+kwopt+mE-1, kbetaR, kbetaR) = he_s_u_mat.block(kE+kwopt+mE-1, kE+kwopt+mE-1, kbetaR, kbetaR) - Ones;


    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_phi_log_theta_vec = he_phi_log_theta();
    he_betaR_log_theta_vec = he_betaR_log_theta();
    he_con_index_par_log_theta_vec = he_con_index_par_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();
    he_log_theta_scalar = he_log_theta();

    I_mat.row(kwopt+kE+mE-1+kbetaR+kbetaF) << he_alpha_f_log_theta_vec.transpose(), he_phi_log_theta_vec.transpose(),he_con_index_par_log_theta_vec.transpose(), he_betaR_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose(), he_log_theta_scalar;
    I_mat.col(kwopt+kE+mE-1+kbetaR+kbetaF) = I_mat.row(kwopt+kE+mE-1+kbetaR+kbetaF).transpose();

  }



  void prepare_AIC_proposed () {
    
    // Matrix K
    K_alpha_f_mat = K_alpha_f();
    K_betaR_mat = K_betaR();
    K_betaF_mat = K_betaF();
    K_alpha_w_mat = K_alpha_w();
    K_phi_mat = K_phi();
    K_index_par_mat = K_index_par();
    K_con_index_par_mat = K_con_index_par();
    K_log_theta_vec = K_log_theta();

    Kleft.resize(kwopt+kE+mE-1+kbetaR+kbetaF+1, n);
    Kleft.setZero();

    Khat.resize(kwopt+kE+mE-1+kbetaR+kbetaF+1, kwopt+kE+mE-1+kbetaR+kbetaF+1);

    Kleft.block(0, 0, kE, n) = K_alpha_f_mat;
    Kleft.block(kE, 0, kwopt, n) = K_phi_mat;
    Kleft.block(kE+kwopt,0,mE-1,n) = K_con_index_par_mat;
    Kleft.block(kE+kwopt+mE-1,0,kbetaR,n) = K_betaR_mat;
    Kleft.block(kE+kwopt+mE-1+kbetaR,0,kbetaF,n) = K_betaF_mat;
    Kleft.row(kwopt+kE+mE-1+kbetaR+kbetaF) = K_log_theta_vec.transpose();

    Khat = Kleft * Kleft.transpose();
  }



  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Vec dlogdensity_dmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Vec d2logdensity_dmudmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
 


  
  // 2. mean model
  // d log(mu) / d alpha_w
  Mat dlogmu_dw () {
    Mat out(n, kw);
    Vec Bf1st;
    for (int i = 0; i < n; i++) {
      Bf1st = BsplinevecCon1st(E(i), knots_f, 4, Zf);
      out.row(i) = B_inner.row(i) * (Bf1st.dot(alpha_f));
    }
    return out;
  }
  // d mu / d alpha_w
  Mat dmu_dw () {
    Mat out(n, kw);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dw_mat.row(i) * mu(i);
    }
    return out;
  }

  // d log mu / d index_par
  Mat dlogmu_dindex_par () {
    Mat out(n, mE);
    Vec Bf1st;
    Vec tmp;

    for (int i = 0; i < n; i++) {
      Bf1st = BsplinevecCon1st(E(i), knots_f, 4, Zf);
      tmp = alpha_w_C*Bf1st.dot(alpha_f);
      for (int j = 0; j < mE; j++) {
        out(i,j) = B_inner_list.at(j).row(i).dot(tmp);
      }
    }
    return out;
  }
  // d mu / d index_par
  Mat dmu_dindex_par () {
    Mat out(n, mE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dindex_par_mat.row(i) * mu(i);
    }
    return out;
  }

  // 3. Re-parameterization
  // d alpha_w / d phi
  Mat dw_dphi () {
    // deriv_g <- diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw) - phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    // deriv_g <- deriv_g[1:(kw), 2:kw]

    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));

    // diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw)
    Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(kw);
    D.diagonal().setConstant(1.0/alpha_w_C_denominator);
    Mat Ddense = D.toDenseMatrix();

    // phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    Mat deriv_g2 = (phi_long * (Dw * phi_long).transpose()) * pow(alpha_w_C_denominator, -3);

    Mat deriv_g = Ddense - deriv_g2;
    // Remove the first column
    Mat out = deriv_g.block(0, 1, kw, kw - 1) * K;
    return out;
  }


  std::vector<Mat> d2w_dphidphi () {
    std::vector<Mat> out;
    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long)) = tmp^(1/2);
    Scalar tmp1 = pow(alpha_w_C_denominator, -3); // pow(tmp, -1.5)
    Scalar tmp2 = pow(alpha_w_C_denominator, -5); // pow(tmp, -2.5)
    Vec Dwphi = Dw * phi_long;
    Mat outlarge(kw, kw);
    for (int s = 0; s < kw; s++) {
      if (s == 0) {
        outlarge = -1.0 * tmp1 * Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
      } else {
        Mat m1(kw, kw);
        m1.setZero();
        m1.row(s) = Dwphi.transpose()*tmp1;
        m1.col(s) = m1.col(s) + Dwphi*tmp1;
        Mat m2 = -1.0 * tmp1* Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
        outlarge = -1.0 * m1 + m2; // or m1 - m2 or m2 - m1
      }
      out.push_back(K.transpose() * outlarge.block(1, 1, kw-1, kw-1) * K);
    }
    return out;
  }


  // d index_par / d con_index_par
  // Mat dindex_par_dcon_index_par () {
  //   // index_par_denominator = sqrt(con_index_par_long.dot(con_index_par_long));
  //   // index_par = con_index_par_long/index_par_denominator;

  //   Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(mE);
  //   D.diagonal().setConstant(1.0/index_par_denominator);
  //   Mat Ddense = D.toDenseMatrix();

  //   // identity matrix
  //   Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> Dp(mE);
  //   Dp.diagonal().setConstant(1.0);
  //   Mat Dpdense = Dp.toDenseMatrix();


  //   Mat deriv_g2 = (con_index_par_long * (Dpdense * con_index_par_long).transpose()) * pow(index_par_denominator, -3);

  //   Mat deriv_g = Ddense - deriv_g2;
  //   // Remove the first column
  //   Mat out = deriv_g.block(0, 1, mE, mE - 1);
  //   return out;
  // }
  Mat dindex_par_dcon_index_par () {

    Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(mE);
    D.diagonal().setConstant(1.0/index_par_denominator);
    Mat Ddense = D.toDenseMatrix();

    Mat deriv_g2 = (con_index_par_long * (BtBindex * con_index_par_long).transpose()) * pow(index_par_denominator, -3);

    Mat deriv_g = Ddense - deriv_g2;
    // Remove the first column
    Mat out = Bindex * (deriv_g.block(0, 1, mE, mE - 1));
    return out;
  }



  // std::vector<Mat> d2index_par_dcon_index_pardcon_index_par () {
  //   std::vector<Mat> out;

  //   // identity matrix
  //   Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> Dp(mE);
  //   Dp.diagonal().setConstant(1.0);
  //   Mat Dpdense = Dp.toDenseMatrix();

  //   Scalar tmp1 = pow(index_par_denominator, -3); // pow(tmp, -1.5)
  //   Scalar tmp2 = pow(index_par_denominator, -5); // pow(tmp, -2.5)
  //   // Vec Dwphi = Dw * phi_long;
  //   for (int s = 0; s < mE; s++) {
  //     Mat outlarge(mE, mE);
  //     if (s == 0) {
  //       outlarge = -1.0 * tmp1 * Dpdense + 3.0 * tmp2 * con_index_par_long * con_index_par_long.transpose();
  //     } else {
  //       Mat m1(mE, mE);
  //       m1.setZero();
  //       m1.row(s) = con_index_par_long.transpose()*tmp1;
  //       m1.col(s) = m1.col(s) + con_index_par_long*tmp1;
  //       Mat m2 = -1.0 * tmp1* Dpdense + 3.0 * tmp2 * con_index_par_long * con_index_par_long.transpose();
  //       outlarge = -1.0 * m1 + m2; // or m1 - m2 or m2 - m1
  //     }
  //     out.push_back(outlarge.block(1, 1, mE-1, mE-1));
  //   }
  //   return out;
  // }
  
  std::vector<Mat> d2index_par_dcon_index_pardcon_index_par () {
    std::vector<Mat> out;
    std::vector<Mat> outreal;
    
    Scalar tmp1 = pow(index_par_denominator, -3); // pow(tmp, -1.5)
    Scalar tmp2 = pow(index_par_denominator, -5); // pow(tmp, -2.5)
    Vec BtBindexcon_index_par_long = BtBindex * con_index_par_long;
    Mat outlarge(mE, mE);

    
    for (int s = 0; s < mE; s++) {
      if (s == 0) {
        outlarge = -1.0 * tmp1 * BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
      } else {
        Mat m1(mE, mE);
        m1.setZero();
        m1.row(s) = BtBindexcon_index_par_long.transpose()*tmp1;
        m1.col(s) = m1.col(s) + BtBindexcon_index_par_long*tmp1;
        Mat m2 = -1.0 * tmp1* BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
        outlarge = -1.0 * m1 + m2;
      }

      out.push_back(outlarge.block(1, 1, mE-1, mE-1));
    }

    Mat outrealtmp(mE-1, mE-1);
    for (int s = 0; s < mE; s++) {
      outrealtmp.setZero();
      for (int i = 0; i < mE; i++) {
        outrealtmp += Bindex(s, i) * out.at(i);
      }
      outreal.push_back(outrealtmp);
    }

    return outreal;
  }

  // components for K

  Mat K_alpha_f () {
    // dmu_df_mat: n * kE
    // dlogdensity_dmu_vec: n * 1
    // out: kE * n
    Mat out(kE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Bf_matrix.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }


  Mat K_betaR () {
    // dmu_dbetaR_mat: n * kbetaR
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaR * n
    Mat out(kbetaR, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Xrand.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }


  Mat K_betaF () {
    // dmu_dbetaF_mat: n * kbetaF
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaF * n
    Mat out(kbetaF, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Xfix.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

   Mat K_alpha_w () {
    // dmu_dw_mat: n * kw
    // dlogdensity_dmu_vec: n * 1
    // out: kw * n
    Mat out(kw, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Mat K_phi () {
    // K_alpha_w_mat: kw * n
    // dw_dphi_mat: kw * kwopt
    // out: kwopt * n
    return dw_dphi_mat.transpose() * K_alpha_w_mat;
  }


  Mat K_index_par () {
    // dmu_dindex_par_mat: n * mE
    // dlogdensity_dmu_vec: n * 1
    // out: mE * n
    Mat out(mE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dindex_par_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }


  Mat K_con_index_par () {
    // K_index_par_mat: mE * n
    // dindex_par_dcon_index_par_mat: mE * (mE-1)
    // out: (mE-1) * n
    return dindex_par_dcon_index_par_mat.transpose() * K_index_par_mat;
  }

  Vec K_log_theta () {
    // out: n * 1. DO NOT FORGET TO TRANSPOSE WHEN USING
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = -1.0 * theta * (log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta));
    }
    return out;
  }


  // *** Hessian ***
  Mat he_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Bf_matrix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * Bf_matrix.row(i));
    }
    return - out1 - out2 + smoothing_f*Sf;
  }

  // Mat I_alpha_f () {
  //   Mat out1(kE, kE);
  //   Mat out2(kE, kE);
  //   out1.setZero();
  //   out2.setZero();
  //   for (int i = 0; i < n; i++) {
  //     out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Bf_matrix.row(i)) * mu(i) * mu(i);
  //     out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * Bf_matrix.row(i));
  //   }
  //   return - out1 - out2;
  // }



  Mat he_betaR () {
    Mat out1(kbetaR, kbetaR);
    Mat out2(kbetaR, kbetaR);
    Mat Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xrand.row(i)).transpose() * (Xrand.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Xrand.row(i).transpose() * Xrand.row(i));
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      double ri = CppAD::Value(r(i));
      int ki = static_cast<int>(ri);
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }

  // Mat I_betaR () {
  //   Mat out1(kbetaR, kbetaR);
  //   Mat out2(kbetaR, kbetaR);
  //   out1.setZero();
  //   out2.setZero();
  //   for (int i = 0; i < n; i++) {
  //     out1 += d2logdensity_dmudmu_vec(i) * (Xrand.row(i)).transpose() * (Xrand.row(i)) * mu(i) * mu(i);
  //     out2 += dlogdensity_dmu_vec(i) * (mu(i) * Xrand.row(i).transpose() * Xrand.row(i));
  //   }
  //   return - out1 - out2;
  // }

  Mat he_betaF () {
    Mat out1(kbetaF, kbetaF);
    Mat out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xfix.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Xfix.row(i).transpose() * Xfix.row(i));
    }
    return - out1 - out2;
  }

  Mat he_alpha_w () {
    Mat out1(kw, kw);
    Mat out2(kw, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * (BsplinevecCon2nd(E(i), knots_f, 4, Zf).dot(alpha_f)) * B_inner.row(i).transpose() * B_inner.row(i));
    }
    return - out1 - out2 + smoothing_w*Sw_large;
  }

  Mat I_alpha_w () {
    Mat out1(kw, kw);
    Mat out2(kw, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * (BsplinevecCon2nd(E(i), knots_f, 4, Zf).dot(alpha_f)) * B_inner.row(i).transpose() * B_inner.row(i));
    }
    return - out1 - out2;
  }
  

  Mat he_phi () {
    Mat out1 = dw_dphi_mat.transpose() * he_alpha_w_mat * dw_dphi_mat;
    Mat out2(kwopt, kwopt);
    out2.setZero();

    Vec gr_pen_w = smoothing_w * Sw_large * alpha_w_C;

    Vec gr_alpha_w_vec = - dmu_dw_mat.transpose() * dlogdensity_dmu_vec + gr_pen_w;

    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }

  Mat I_phi () {
    Mat out1 = dw_dphi_mat.transpose() * I_alpha_w_mat * dw_dphi_mat;
    Mat out2(kwopt, kwopt);
    out2.setZero();


    Vec gr_alpha_w_vec = - dmu_dw_mat.transpose() * dlogdensity_dmu_vec;

    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }
  

  Mat he_index_par () {
    
    Mat out1(mE, mE);
    Mat out2(mE, mE);
    Mat out_tmp(mE, mE);
    out1.setZero();
    out2.setZero();
    
    Vec Bf2nd;
    Scalar tmp;
    for (int i = 0; i < n; i++) {
      Bf2nd = BsplinevecCon2nd(E(i), knots_f, 4, Zf);
      tmp = Bf2nd.dot(alpha_f);
      out_tmp.setZero();
      for (int ii = 0; ii < mE; ii++) {
        for (int jj = 0; jj < mE; jj++) {
          out_tmp(ii,jj) = tmp * (B_inner_list.at(ii).row(i).dot(alpha_w_C)) * (B_inner_list.at(jj).row(i).dot(alpha_w_C));
        }
      }

      out1 += d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dindex_par_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i) * out_tmp);
    }

    return - out1 - out2;
  }


  Mat he_con_index_par () {
    Mat out1 = dindex_par_dcon_index_par_mat.transpose() * he_index_par_mat * dindex_par_dcon_index_par_mat;
    Mat out2(mE-1, mE-1);
    out2.setZero();
    Vec tmp = - dmu_dindex_par_mat.transpose() * dlogdensity_dmu_vec;
    for (int s = 0; s < mE; s++) {
      out2 = out2 + tmp(s) * d2index_par_dcon_index_pardcon_index_par_list.at(s);
    }
    return out1 + out2;
  }



  Mat he_alpha_f_alpha_w () {
    Mat out1(kE, kw);
    Mat out2(kE, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i) * mu(i)).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i)*BsplinevecCon1st(E(i), knots_f, 4, Zf) * B_inner.row(i));
    }
    return - out1 - out2;
  }

  Mat he_alpha_f_phi () {
    Mat out = he_alpha_f_alpha_w_mat * dw_dphi_mat;
    return out;
  }

  Mat he_alpha_f_index_par () {
    Mat out1(kE, mE);
    Mat out2(kE, mE);
    Mat out_mat(kE, mE);
    out1.setZero();
    out2.setZero();
    Vec Bf1st;
    for (int i = 0; i < n; i++) {
      out_mat.setZero();
      Bf1st = BsplinevecCon1st(E(i), knots_f, 4, Zf);
      for (int jj = 0; jj < mE; jj++) {
        out_mat.col(jj) = Bf1st * (B_inner_list.at(jj).row(i).dot(alpha_w_C));
      }
      out1 += d2logdensity_dmudmu_vec(i) * Bf_matrix.row(i).transpose() * dmu_dindex_par_mat.row(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i)*out_mat);
    }
    return - out1 - out2;
  }


  Mat he_alpha_f_con_index_par () {
    Mat out = he_alpha_f_index_par_mat * dindex_par_dcon_index_par_mat;
    return out;
  }

  Mat he_alpha_f_betaF () {
    Mat out1(kE, kbetaF);
    Mat out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Bf_matrix.row(i).transpose() * (Xfix.row(i) * mu(i));
    }
    return - out1 - out2;
  }
  Mat he_alpha_f_betaR () {
    Mat out1(kE, kbetaR);
    Mat out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Xrand.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Bf_matrix.row(i).transpose() * (Xrand.row(i) * mu(i));
    }
    return - out1 - out2;
  }
  Mat he_phi_betaR () {
    Mat out1(kwopt, kbetaR);
    Mat out2(kwopt, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
    }
    return - out1 - out2;
  }
  Mat he_phi_betaF () {
    Mat out1(kwopt, kbetaF);
    Mat out2(kwopt, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
    }
    return - out1 - out2;
  }
  Mat he_betaR_betaF () {
    Mat out1(kbetaR, kbetaF);
    Mat out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xrand.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Xrand.row(i).transpose() * (Xfix.row(i) * mu(i));
    }
    return - out1 - out2;
  }


  Mat he_con_index_par_betaF () {
    Mat out1(mE-1, kbetaF);
    Mat out2(mE-1, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
    }
    return - out1 - out2;
  }




  Mat he_alpha_w_index_par () {
    Mat out1(kw, mE);
    Mat out2(kw, mE);
    Mat out_mat(kw,mE);
    Vec Bf1st;
    Vec Bf2nd;
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out_mat.setZero();
      Bf1st = BsplinevecCon1st(E(i), knots_f, 4, Zf);
      Bf2nd = BsplinevecCon2nd(E(i), knots_f, 4, Zf);
      for (int jj = 0; jj < mE; jj++) {
        out_mat.col(jj) = B_inner_list.at(jj).row(i).transpose() * (Bf1st.dot(alpha_f)) + B_inner.row(i).transpose() * (Bf2nd.dot(alpha_f)) * (B_inner_list.at(jj).row(i).dot(alpha_w_C));
      }
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i)*out_mat);
    }
    return - out1 - out2;
  }


  Mat he_phi_index_par () {
    Mat out = dw_dphi_mat.transpose() * he_alpha_w_index_par_mat;
    return out;
  }

  Mat he_phi_con_index_par () {
    Mat out = he_phi_index_par_mat * dindex_par_dcon_index_par_mat;
    return out;
  }


  Mat he_con_index_par_betaR () {
    Mat out1(mE-1, kbetaR);
    Mat out2(mE-1, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
    }
    return - out1 - out2;
  }

  // for AIC

  Vec d2logdensity_dmudtheta () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }

  Vec he_alpha_f_log_theta () {
    Mat dmu_df_mat(n, kE);
    for (int i = 0; i < n; i++) {
      dmu_df_mat.row(i) = Bf_matrix.row(i) * mu(i);
    }

    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }

  Vec he_phi_log_theta () {
    return -1.0*theta * dw_dphi_mat.transpose() * ( dmu_dw_mat.transpose() * d2logdensity_dmudtheta_vec );
  }

  Vec he_betaR_log_theta () {
    Mat dmu_dbetaR_mat(n, kbetaR);
    for (int i = 0; i < n; i++) {
      dmu_dbetaR_mat.row(i) = Xrand.row(i) * mu(i);
    }

    return -1.0*dmu_dbetaR_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }

  Vec he_con_index_par_log_theta () {
    return -1.0*theta * dindex_par_dcon_index_par_mat.transpose() * ( dmu_dindex_par_mat.transpose() * d2logdensity_dmudtheta_vec );
  }


  Vec he_betaF_log_theta () {
    Mat dmu_dbetaF_mat(n, kbetaF);
    for (int i = 0; i < n; i++) {
      dmu_dbetaF_mat.row(i) = Xfix.row(i) * mu(i);
    }
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }


  Scalar he_log_theta () {

    Scalar dlogdensity_dtheta_scalar = 0.0;
    for (int i = 0; i < n; i++) {
      dlogdensity_dtheta_scalar += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    Scalar d2logdensity_dthetadtheta_scalar = 0.0;
    for (int i = 0; i < n; i++) {
      d2logdensity_dthetadtheta_scalar += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    }

    return -1.0 * theta * dlogdensity_dtheta_scalar - theta * theta * d2logdensity_dthetadtheta_scalar;
  }




  // *********** LAML ***********
  Scalar logdetH05() {
    Scalar logdetH05 = 0.0;
    Eigen::FullPivLU<Mat> lu(he_s_u_mat);
    Mat LU = lu.matrixLU();
    // Scalar c = lu.permutationP().determinant(); // -1 or 1
    Scalar lii;
    for (int i = 0; i < (kwopt+mE-1+kE+kbetaR+kbetaF); i++) {
      lii = LU(i,i);
      // std::cout << "lii : " << CppAD::Value(lii) << std::endl;
      // std::cout << "c : " << CppAD::Value(c) << std::endl;
      // if (lii < 0.0) c *= -1;

      logdetH05 += log(CppAD::abs(lii));
    }
    // logdetH05 += log(c);



    return logdetH05/2.0;

  }
};







#endif
#RCPP based functions
###### Rcpp Functions for Normal kernel ####
{
  Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
  Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
  Sys.setenv("PKG_LIBS"="-fopenmp")
  
  
  sourceCpp(code = '
            
            // [[Rcpp::depends(RcppArmadillo)]]
            #include <RcppArmadillo.h>
            #include <omp.h>
            #include <cmath>
            // [[Rcpp::plugins(openmp)]]
            
            using namespace Rcpp;
            using namespace arma;
            
            // [[Rcpp::export]]
            arma::mat pnormpar_arma(arma::mat X){
            arma::mat p=normcdf(X);
            return p;
            }
            
            // [[Rcpp::export]]
            NumericVector dnormpar(NumericVector x){
            double c = 1/sqrt(2*M_PI);
            int n = x.size();
            NumericVector ret(n);
            omp_set_num_threads(10);
            #pragma omp parallel for
            for(int i=0; i<n; ++i)
            ret[i] = exp(-x[i]*x[i]/2)*c;
            return ret;
            }
            
            // [[Rcpp::export]]
            NumericVector pnormpar(NumericVector x){
            NumericVector ret = pnorm(x);
            return ret;
            }
            
            // [[Rcpp::export]]
            NumericMatrix dnormpar_symneg_mat(NumericMatrix x){
            double c = 1/sqrt(2*M_PI);
            int nr = x.rows();
            int nc= x.cols();
            NumericMatrix ret(nr,nc);
            #pragma omp parallel for if(nc> 50)
            for(int i=0; i<nr; ++i){
            for(int j=0; j<=i; ++j){
            
            if ((x(i,j)<(-5)) | (x(i,j)>5)) {
            ret(i,j) = 0;
            } else {
            ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
            }
            ret(j,i)=ret(i,j);
            }
            }
            return ret;
            }
            
            
            // [[Rcpp::export]]
            arma::mat dnormpar_mat_arma(arma::mat x){
            double c = 1/sqrt(2*M_PI);
            arma::mat ret = exp(-x*x/2)*c;
            return ret;
            }
            
            // [[Rcpp::export]]
            NumericMatrix dnormpar_mat(NumericMatrix x){
            double c = 1/sqrt(2*M_PI);
            int nr = x.rows();
            int nc= x.cols();
            NumericMatrix ret(nr,nc);
            #pragma omp parallel for if(nc> 50)
            for(int i=0; i<nr; ++i){
            for(int j=0; j<nc; ++j){
            
            ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
            }
            }
            return ret;
            }
            
            // [[Rcpp::export]]
            NumericMatrix pnormpar_mat(NumericMatrix x){
            int nr= x.rows();
            int nc= x.cols();
            NumericMatrix p(nr,nc);
            #pragma omp parallel for if(nc> 50)
            for(int i=0; i<nr; ++i){
            for(int j=0; j<nc; ++j){
            
            if (x(i,j)<(-5)) {
            p(i,j) = 0;
            } else if (x(i,j)>5) {
            p(i,j) = 1;
            } else {
            //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
            p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
            }
            }
            }
            return p;
            }
            
            // [[Rcpp::export]]
            arma::mat pnormpar_mat_arma(arma::mat x){
            arma::mat ret= 0.5*(1+erf(x/pow(2,0.5)));
            return ret;
            }
            ')
}


#lhood0102
{
  
  # sourceCpp(code ='
  #           // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  #           #include <cmath>
  #           using namespace Rcpp;
  #           using namespace arma;
  
  likelihood.include0102 <- ' 
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  
  // [[Rcpp::export]]
  Rcpp::NumericVector lhood_normker0102(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec d = Rcpp::as<arma::vec>(e["delta"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec V0 = Rcpp::as<arma::vec>(e["V0"]);
  int n=X.n_rows;
  int nf=X0.n_rows;
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(V) -Xbeta;
  arma::vec Rf = log(V0) -X0beta;
  arma::vec In= arma::ones(n,1);
  arma::mat r= In*trans(R);
  arma::mat RtR=(trans(r)-r);
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  arma::mat RtR2= RtR.rows(find(d==1));
  arma::mat pnormRRn=pnormpar_mat(-RtR2/a[1]);
  
  double s1s2= sum(log(V0))/n;
  double s3= sum(log(sum(dnormRRn,0)/(n*a[0])))/n;
  arma::vec s_4=pnormRRn*g;
  double s4=sum((log(s_4/n))/n);
  double s=s1s2-s3+s4;
  Rcpp::NumericVector out(1);
  out[0] = s;
  return out;
  }'
  #)
}

#grad0102
{
  
  # sourceCpp(code = '
  #           // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  # 
  #           using namespace Rcpp;
  #           using namespace arma;
  
  gradient.include0102 <-'  
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat dnormpar_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 100)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  }
  }
  return ret;
  };
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  NumericVector grad_normker0102(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec d = Rcpp::as<arma::vec>(e["delta"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec V0 = Rcpp::as<arma::vec>(e["V0"]);
  int n=X.n_rows;
  int nf=X0.n_rows;
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(V) -Xbeta;
  arma::vec Rf = log(V0) -X0beta;
  arma::vec In= arma::ones(n,1);
  arma::mat r= In*trans(R);
  arma::mat RtR=(trans(r)-r);
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  
  arma::mat RtR2= -RtR.rows(find(d==1));
  arma::mat pnormRRn= pnormpar_mat(RtR2/a[1]);
  
  arma::vec denom= pnormRRn*g;
  arma::mat dnorm_Rf=dnormRRn%(RtRf/a[0]);
  arma::vec rowsums_dnormRRn = sum(dnormRRn,1);
  arma::mat dnormpar_RtR2=dnormpar_mat(RtR2/a[1]);
  int nb=beta.n_elem;
  NumericVector grad(nb);
  for (int k=0; k<nb; ++k){
  arma::mat Xk=arma::ones(nf,1)*trans(X0.col(k));
  arma::mat XktXk=(trans(Xk)-Xk)/a[0];
  arma::mat Xk2=arma::ones(n,1)*trans(X.col(k));
  arma::mat XktXk2=(trans(Xk2)-Xk2)/a[1];
  double rrr= -(accu((sum(dnorm_Rf%XktXk,1))/rowsums_dnormRRn))/n;
  double rrr2= (accu(((((dnormpar_RtR2))%(XktXk2.rows(find(d==1))))*g)/(denom)))/n;
  grad[k] =  (rrr+rrr2);
  } 
  return grad;
  }
  '
  #)
}

#body+inline0102
{
  likelihood.body0102 <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&lhood_normker0102)));
  '
  gradient.body0102 <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&grad_normker0102)));
  '
  
  settings <- getPlugin("RcppArmadillo")
  settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
  settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
  
  l_s_m_01_02_f.CPP <- cxxfunction(signature(), body=likelihood.body0102,
                                   inc=likelihood.include0102, plugin="RcppArmadillo",settings = settings)
  grad_beta01_02_f.CPP <- cxxfunction(signature(), body=gradient.body0102,
                                      inc=gradient.include0102, plugin="RcppArmadillo",settings = settings)
}

#lhood0102pert
{
  
  # sourceCpp(code ='
  #           // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  #           #include <cmath>
  #           using namespace Rcpp;
  #           using namespace arma;
  
  likelihood.include0102pert <- ' 
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  Rcpp::NumericVector lhood_normker0102pert(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::vec G = Rcpp::as<arma::vec>(e["G"]);
  arma::vec G0 = Rcpp::as<arma::vec>(e["G0"]);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec d = Rcpp::as<arma::vec>(e["delta"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec V0 = Rcpp::as<arma::vec>(e["V0"]);
  int n=X.n_rows;
  int nf=X0.n_rows;
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(V) -Xbeta;
  arma::vec Rf = log(V0) -X0beta;
  arma::vec In= arma::ones(n,1);
  arma::mat r= In*trans(R);
  arma::mat RtR=(trans(r)-r);
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  arma::mat RtR2= RtR.rows(find(d==1));
  arma::mat pnormRRn=pnormpar_mat(-RtR2/a[1]);
  
  double s1s2   = as_scalar(trans(G0)*log(V0)/n);
  arma::vec c = (dnormRRn*G0)/(n*a[0]);
  double s3     = as_scalar(trans(log(c))*G0/n);
  arma::vec c1 = pnormRRn*(g%G)/n;
  double s_4    = as_scalar(trans(log(c1))*G0/n);
  double s      = s1s2-s3+s_4;
  Rcpp::NumericVector out(1);
  out[0] = s;
  return out;
  
  }'
  #)
}

#grad0102pert
{
  
  # sourceCpp(code = '
  #           // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  # 
  #           using namespace Rcpp;
  #           using namespace arma;
  
  gradient.include0102pert <-'  
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat dnormpar_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 100)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  }
  }
  return ret;
  };
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  NumericVector grad_normker0102pert(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::vec G = Rcpp::as<arma::vec>(e["G"]);
  arma::vec G0 = Rcpp::as<arma::vec>(e["G0"]);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec d = Rcpp::as<arma::vec>(e["delta"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec V0 = Rcpp::as<arma::vec>(e["V0"]);
  int n=X.n_rows;
  int nf=X0.n_rows;
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(V) -Xbeta;
  arma::vec Rf = log(V0) -X0beta;
  arma::vec In= arma::ones(n,1);
  arma::mat r= In*trans(R);
  arma::mat RtR=(trans(r)-r);
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  
  arma::mat RtR2= -RtR.rows(find(d==1));
  arma::mat pnormRRn= pnormpar_mat(RtR2/a[1]);
  
  arma::vec denom= pnormRRn*(g%G);
  arma::mat dnorm_Rf= dnormRRn%(RtRf/a[0]);
  arma::vec rowsums_dnormRRn = dnormRRn*(G0);
  arma::mat dnormpar_RtR2=dnormpar_mat(RtR2/a[1]);
  int nb=beta.n_elem;
  NumericVector grad(nb);
  for (int k=0; k<nb; ++k){
  arma::mat Xk=arma::ones(nf,1)*trans(X0.col(k));
  arma::mat XktXk=(trans(Xk)-Xk)/a[0];
  arma::mat Xk2=arma::ones(n,1)*trans(X.col(k));
  arma::mat XktXk2=(trans(Xk2)-Xk2)/a[1];
  double rrr=  as_scalar(-((trans(G0)*(dnorm_Rf%XktXk))/trans(rowsums_dnormRRn)/n)*(G0));
  double rrr2= as_scalar(trans(G0)*(((((dnormpar_RtR2))%(XktXk2.rows(find(d==1))))*(g%G))/denom/n));
  grad[k] =  (rrr+rrr2);
  } 
  return grad;
  }
  '
  #)
}

#body+inline0102pert
{
  likelihood.body0102pert <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&lhood_normker0102pert)));
  '
  gradient.body0102pert <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&grad_normker0102pert)));
  '
  
  settings <- getPlugin("RcppArmadillo")
  settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
  settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
  
  l_s_m_01_02pert_f.CPP <- cxxfunction(signature(), body=likelihood.body0102pert,
                                       inc=likelihood.include0102pert, plugin="RcppArmadillo",settings = settings)
  grad_beta01_02pert_f.CPP <- cxxfunction(signature(), body=gradient.body0102pert,
                                          inc=gradient.include0102pert, plugin="RcppArmadillo",settings = settings)
}

#lhood12
{
  
  #sourceCpp(code ='
  # // [[Rcpp::depends(RcppArmadillo)]]
  # #include <RcppArmadillo.h>
  # 
  # using namespace Rcpp;
  # using namespace arma;
  likelihood.include12 <- '  
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  
  // [[Rcpp::export]]
  Rcpp::NumericVector lhood_normker12(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec W = Rcpp::as<arma::vec>(e["W"]);
  arma::vec W0 = Rcpp::as<arma::vec>(e["W0"]);
  arma::vec d3 = Rcpp::as<arma::vec>(e["delta3"]);
  int n = Rcpp::as<int>(e["n"]);
  int ne=X.n_rows;
  
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(W) -Xbeta;
  arma::vec Rv = log(V) -Xbeta;
  arma::vec Rf = log(W0) -X0beta;
  arma::vec In= arma::ones(ne,1);
  arma::mat r= In*trans(R);
  arma::mat rv= In*trans(Rv);
  arma::mat RtR=(trans(r)-r);
  arma::mat RtRv=(trans(rv)-r);
  
  int nf=X0.n_rows;
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  arma::mat RtR2W= RtR.cols(find(d3==1));
  arma::mat RtR2V= RtRv.cols(find(d3==1));
  arma::mat pnormRRnW=pnormpar_mat(-RtR2W/a[1]);
  arma::mat pnormRRnV=pnormpar_mat(-RtR2V/a[1]);
  
  double s1s2= sum(log(W0))/n;
  double s3= sum(log(sum(dnormRRn,0)/(n*a[0])))/n;
  arma::vec s_4= -trans(pnormRRnW-pnormRRnV)*g;
  double s4=sum((log(s_4/n))/n);
  double s=s1s2-s3+s4;
  Rcpp::NumericVector out(1);
  out[0] = s;
  return out;
  }'
  #)
}


#grad12
{
  
  #sourceCpp(code = '
  # // [[Rcpp::depends(RcppArmadillo)]]
  # #include <RcppArmadillo.h>
  # 
  # using namespace Rcpp;
  # using namespace arma;
  
  gradient.include12 <-'
  
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat dnormpar_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 100)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  }
  }
  return ret;
  };
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  NumericVector grad_normker12(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec W = Rcpp::as<arma::vec>(e["W"]);
  arma::vec W0 = Rcpp::as<arma::vec>(e["W0"]);
  arma::vec d3 = Rcpp::as<arma::vec>(e["delta3"]);
  int n = Rcpp::as<int>(e["n"]);
  int ne=X.n_rows;
  
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(W) -Xbeta;
  arma::vec Rv = log(V) -Xbeta;
  arma::vec Rf = log(W0) -X0beta;
  arma::vec In= arma::ones(ne,1);
  arma::mat r= In*trans(R);
  arma::mat rv= In*trans(Rv);
  arma::mat RtR=(trans(r)-r);
  arma::mat RtRv=(trans(rv)-r);
  
  int nf=X0.n_rows;
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  
  arma::mat RtR2W= RtR.cols(find(d3==1));
  arma::mat RtR2V= RtRv.cols(find(d3==1));
  arma::mat pnormRRnW=pnormpar_mat(-RtR2W/a[1]);
  arma::mat pnormRRnV=pnormpar_mat(-RtR2V/a[1]);
  
  arma::mat dnorm_Rf=dnormRRn%(RtRf/a[0]);
  arma::vec rowsums_dnormRRn = sum(dnormRRn,1);
  arma::mat denom= -trans(pnormRRnW-pnormRRnV)*g;
  arma::mat dnormpar_RtR2W= dnormpar_mat(RtR2W/a[1]);
  arma::mat dnormpar_RtR2V= dnormpar_mat(RtR2V/a[1]);
  arma::mat dnormpar_RtR2= dnormpar_mat(RtR2W/a[1])-dnormpar_mat(RtR2V/a[1]);
  int nb=beta.n_elem;
  NumericVector grad(nb);
  for (int k=0; k<nb; ++k){
  arma::mat Xk=arma::ones(nf,1)*trans(X0.col(k));
  arma::mat XktXk=(trans(Xk)-Xk)/a[0];
  arma::mat Xk2=arma::ones(ne,1)*trans(X.col(k));
  arma::mat XktXk2=(trans(Xk2)-Xk2)/a[1];
  double rrr= -(accu((sum(dnorm_Rf%XktXk,1))/rowsums_dnormRRn))/n;
  //arma::mat tt = (trans(g)*(-dnormpar_RtR2%(XktXk2.cols(find(d3==1)))))/trans(denom);
  double rrr2= (accu((trans(g)*(-dnormpar_RtR2%(XktXk2.cols(find(d3==1)))))/trans(denom)))/n;
  grad[k] =  (rrr+rrr2);
  } 
  return grad;
  }
  '
  #)
}

#body+inline12
{
  likelihood.body12 <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&lhood_normker12)));
  '
  gradient.body12 <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&grad_normker12)));
  '
  settings <- getPlugin("RcppArmadillo")
  settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
  settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
  
  
  l_s_m_12_f.CPP <- cxxfunction(signature(), body=likelihood.body12,
                                inc=likelihood.include12, plugin="RcppArmadillo",settings = settings)
  grad_beta12_f.CPP <- cxxfunction(signature(), body=gradient.body12,
                                   inc=gradient.include12, plugin="RcppArmadillo",settings = settings)
}

#lhood12pert
{
  
  # sourceCpp(code ='
  #           // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  #           
  #           using namespace Rcpp;
  #           using namespace arma;
  likelihood.include12pert <- '  
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  Rcpp::NumericVector lhood_normker12pert(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::vec G = Rcpp::as<arma::vec>(e["G"]);
  arma::vec G0 = Rcpp::as<arma::vec>(e["G0"]);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec W = Rcpp::as<arma::vec>(e["W"]);
  arma::vec W0 = Rcpp::as<arma::vec>(e["W0"]);
  arma::vec d3 = Rcpp::as<arma::vec>(e["delta3"]);
  int n = Rcpp::as<int>(e["n"]);
  int ne=X.n_rows;
  
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(W) -Xbeta;
  arma::vec Rv = log(V) -Xbeta;
  arma::vec Rf = log(W0) -X0beta;
  arma::vec In= arma::ones(ne,1);
  arma::mat r= In*trans(R);
  arma::mat rv= In*trans(Rv);
  arma::mat RtR=(trans(r)-r);
  arma::mat RtRv=(trans(rv)-r);
  
  int nf=X0.n_rows;
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  arma::mat RtR2W= RtR.cols(find(d3==1));
  arma::mat RtR2V= RtRv.cols(find(d3==1));
  arma::mat pnormRRnW=pnormpar_mat(-RtR2W/a[1]);
  arma::mat pnormRRnV=pnormpar_mat(-RtR2V/a[1]);
  
  double s1s2           = as_scalar(trans(G0)*log(W0)/n);
  arma::vec c           = (dnormRRn*G0)/(n*a[0]);
  double s3             = as_scalar((trans(log(c))*G0)/n);
  arma::vec c1          = trans(pnormRRnV-pnormRRnW)*(g%G)/n;
  double s4             = as_scalar(trans(log(c1))*G0/n);
  double s              = s1s2-s3+s4;
  Rcpp::NumericVector out(1);
  out[0] = s;
  return out;
  }'
  #)
}

#grad12pert
{
  
  
  #sourceCpp(code = '
  # // [[Rcpp::depends(RcppArmadillo)]]
  # #include <RcppArmadillo.h>
  # 
  # using namespace Rcpp;
  # using namespace arma;
  
  gradient.include12pert <-'
  arma::mat dnormpar_symneg_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<=i; ++j){
  
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  ret(j,i)=ret(i,j);
  }
  }
  return ret;
  }
  
  
  arma::mat dnormpar_mat(arma::mat x){
  double c = 1/sqrt(2*M_PI);
  int nr = x.n_rows;
  int nc= x.n_cols;
  arma::mat ret(nr,nc);
  #pragma omp parallel for if(nc> 100)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  //if (x(i,j)<(-5) | x(i,j)>5) {
  //ret(i,j) = 0;
  //} else {
  ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  //}
  }
  }
  return ret;
  };
  
  arma::mat pnormpar_mat(arma::mat x){
  int nr= x.n_rows;
  int nc= x.n_cols;
  arma::mat p(nr,nc);
  #pragma omp parallel for if(nc> 50)
  for(int i=0; i<nr; ++i){
  for(int j=0; j<nc; ++j){
  
  //if (x(i,j)<(-5)) {
  //p(i,j) = 0;
  //} else if (x(i,j)>5) {
  //p(i,j) = 1;
  //} else {
  //p(i,j) = R::pnorm(x(i,j),0.0, 1.0, 1, 0);
  p(i,j)=0.5*(1+erf(x(i,j)/pow(2,0.5)));
  //}
  }
  }
  return p;
  }
  
  // [[Rcpp::export]]
  NumericVector grad_normker12pert(SEXP xs, SEXP env){
  arma::vec beta = Rcpp::as<arma::vec>(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);
  arma::vec G = Rcpp::as<arma::vec>(e["G"]);
  arma::vec G0 = Rcpp::as<arma::vec>(e["G0"]);
  arma::mat X = Rcpp::as<arma::mat>(e["X"]);
  arma::mat X0 = Rcpp::as<arma::mat>(e["X0"]);
  arma::vec g = Rcpp::as<arma::vec>(e["gamma_m"]);
  arma::vec a = Rcpp::as<arma::vec>(e["a"]);
  arma::vec V = Rcpp::as<arma::vec>(e["V"]);
  arma::vec W = Rcpp::as<arma::vec>(e["W"]);
  arma::vec W0 = Rcpp::as<arma::vec>(e["W0"]);
  arma::vec d3 = Rcpp::as<arma::vec>(e["delta3"]);
  int n = Rcpp::as<int>(e["n"]);
  int ne=X.n_rows;
  
  arma::vec Xbeta = X * beta;
  arma::vec X0beta = X0 * beta;
  arma::vec R = log(W) -Xbeta;
  arma::vec Rv = log(V) -Xbeta;
  arma::vec Rf = log(W0) -X0beta;
  arma::vec In= arma::ones(ne,1);
  arma::mat r= In*trans(R);
  arma::mat rv= In*trans(Rv);
  arma::mat RtR=(trans(r)-r);
  arma::mat RtRv=(trans(rv)-r);
  
  int nf=X0.n_rows;
  arma::vec In0= arma::ones(nf,1);
  arma::mat rf= In0*trans(Rf);
  arma::mat RtRf=(trans(rf)-rf);
  
  arma::mat dnormRRn=dnormpar_symneg_mat(RtRf/a[0]);
  
  arma::mat RtR2W= RtR.cols(find(d3==1));
  arma::mat RtR2V= RtRv.cols(find(d3==1));
  arma::mat pnormRRnW=pnormpar_mat(-RtR2W/a[1]);
  arma::mat pnormRRnV=pnormpar_mat(-RtR2V/a[1]);
  
  arma::mat dnorm_Rf=dnormRRn%(RtRf/a[0]);
  arma::vec rowsums_dnormRRn = dnormRRn*G0;
  arma::mat denom= trans(pnormRRnV-pnormRRnW)*(g%G);
  arma::mat dnormpar_RtR2W= dnormpar_mat(RtR2W/a[1]);
  arma::mat dnormpar_RtR2V= dnormpar_mat(RtR2V/a[1]);
  arma::mat dnormpar_RtR2=dnormpar_RtR2W-dnormpar_RtR2V;
  int nb=beta.n_elem;
  NumericVector grad(nb);
  for (int k=0; k<nb; ++k){
  arma::mat Xk=arma::ones(nf,1)*trans(X0.col(k));
  arma::mat XktXk=(trans(Xk)-Xk)/a[0];
  arma::mat Xk2=arma::ones(ne,1)*trans(X.col(k));
  arma::mat XktXk2=(trans(Xk2)-Xk2)/a[1];
  double rrr      = as_scalar(-((trans(G0)*(dnorm_Rf%XktXk))/trans(rowsums_dnormRRn))*G0/n);
  double rrr2     = as_scalar(trans(G0)*((trans((-dnormpar_RtR2%(XktXk2.cols(find(d3==1)))))*(g%G))/(denom))/n);
  grad[k]         =  rrr+rrr2;
  } 
  return grad;
  }'
  #)
}

#body+inline12
{
  likelihood.body12pert <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&lhood_normker12pert)));
  '
  gradient.body12pert <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&grad_normker12pert)));
  '
  settings <- getPlugin("RcppArmadillo")
  settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
  settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
  
  
  l_s_m_12pert_f.CPP <- cxxfunction(signature(), body=likelihood.body12pert,
                                    inc=likelihood.include12pert, plugin="RcppArmadillo",settings = settings)
  grad_beta12pert_f.CPP <- cxxfunction(signature(), body=gradient.body12pert,
                                       inc=gradient.include12pert, plugin="RcppArmadillo",settings = settings)
}

{
  # sourceCpp(code = '
  # // [[Rcpp::depends(RcppArmadillo)]]
  #           // [[Rcpp::depends(RcppEigen)]]
  #           // [[Rcpp::depends(RcppNumerical)]]
  #           #include <RcppArmadillo.h>
  #           #include <RcppNumerical.h>
  #           
  #           using namespace Numer;
  #           using namespace Rcpp;
  #           using namespace arma;
  #           
  #           using Rcpp::NumericVector;
  #           using Rcpp::NumericMatrix;
  #           using arma::vec;
  #           typedef Eigen::Map<Eigen::MatrixXd> MapMat;
  #           typedef Eigen::Map<Eigen::VectorXd> MapVec;
  #           
  #           arma:: vec matmult_mv(NumericMatrix M, NumericVector c) {
  #           arma:: vec mult= as<arma::mat>(M)*as<arma::vec>(c) ;
  #           return mult ;
  #           };
  #           
  #           double matmult_vv(NumericVector m, NumericMatrix c) {
  #           double mult= arma::as_scalar(as<arma::vec>(m).t()*as<arma::mat>(c)) ;
  #           return mult ;
  #           }
  # 
  #           NumericMatrix dnormpar_mat(NumericMatrix x){
  #           double c = 1/sqrt(2*M_PI);
  #           int nr = x.rows();
  #           int nc= x.cols();
  #           NumericMatrix ret(nr,nc);
  #           #pragma omp parallel for if(nc> 50)
  #           for(int i=0; i<nr; ++i){
  #           for(int j=0; j<nc; ++j){
  #           
  #           ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  #           }
  #           }
  #           return ret;
  #           } 
  # 
  #           NumericMatrix pnorm_ltbl_cpp(NumericMatrix x,NumericVector tbl){
  #           int nr = x.rows();
  #           int nc= x.cols();
  #           NumericMatrix p(nr,nc);
  #           #pragma omp parallel for if(nc> 50)
  #           for(int i=0; i<nr; ++i){
  #           for(int j=0; j<nc; ++j){
  # 
  #           if (x(i,j)<(-7)) {
  #           p(i,j) = 0;
  #           } else if (x(i,j)>7) {
  #           p(i,j) = 1;
  #           } else {
  #           p(i,j) = tbl[round((x(i,j)-(-7))*pow(10,4))];
  # 
  #           }
  #           }
  #           }
  #           return p;
  #           }
  # 
  #           class h01: public Numer::Func
  #           {
  #           private:
  #           const Rcpp::NumericVector tbl;
  #           const Rcpp::NumericVector b;
  #           const Rcpp::NumericVector g;
  #           const Rcpp::NumericVector a;
  #           const Rcpp::NumericMatrix X;
  #           const Rcpp::NumericVector V;
  #           const Rcpp::NumericVector D;
  #           
  #           public:
  #           
  #           h01(const Rcpp::NumericVector tbl_,const Rcpp::NumericVector b_,const Rcpp::NumericVector g_, const Rcpp::NumericVector a_, const Rcpp::NumericMatrix x_, const Rcpp::NumericVector v_, const Rcpp::NumericVector d_) : tbl(tbl_),b(b_), g(g_), a(a_), X(x_), V(v_), D(d_) {}
  #           
  #           double operator()(const double& t) const
  #           {
  #           int n=X.rows();
  #           arma:: vec Xb=matmult_mv(X,b);
  #           arma:: vec R= (log(as<arma::vec>(V))-Xb-log(t));
  #           NumericMatrix Rr= as<NumericMatrix>(wrap(R));
  #           NumericMatrix dnormRRn= dnormpar_mat(Rr/a[0]);
  #           double nomi= (matmult_vv(D,dnormRRn))/(n*a[0]*t);
  #           NumericMatrix pnormRRn= pnorm_ltbl_cpp(Rr/a[1],tbl);
  #           double denomi= (matmult_vv(g,pnormRRn))/n;
  #           return nomi/denomi;
  #           }
  #           };
  #           
  #           // [[Rcpp::export]]
  #           double integrate_test01(const double &lower,const double &upper,const Rcpp::NumericVector &tbl,const Rcpp::NumericVector &b, const Rcpp::NumericVector &g,const Rcpp::NumericVector &a, Rcpp::NumericMatrix &x, Rcpp::NumericVector &v, Rcpp::NumericVector &d)
  #           {
  #           h01 f(tbl,b,g,a,x,v,d);
  #           double err_est;
  #           int err_code;
  #           const double res = integrate(f, lower, upper, err_est, err_code);
  #           return res;
  #           //return Rcpp::List::create(
  #           //Rcpp::Named("approximate") = res,
  #           //Rcpp::Named("error_estimate") = err_est,
  #           //Rcpp::Named("error_code") = err_code
  #           //);
  #           }
  # 
  #           // [[Rcpp::export]]
  #           
  #           NumericVector H01(NumericVector& t,const NumericVector& table,const NumericVector& beta, const NumericVector& gamma,const NumericVector& a_n, const NumericMatrix& X, const NumericVector& V, const NumericVector& delta){
  #           int n = t.size();
  #           NumericVector tbl = table;
  #           NumericVector b = beta;
  #           NumericVector g = gamma;
  #           NumericVector a = a_n;
  #           NumericMatrix x = X;
  #           NumericVector v = V;
  #           NumericVector d = delta;
  #           double lower;
  #           double upper;
  #           NumericVector ret(n);
  #           for(int i=0; i<(n-1); ++i){
  #           ret[i] = integrate_test01(lower=t[i],upper=t[i+1],tbl,b,g,a,x,v,d);
  #           }
  #           return ret;
  #           }
  #           ')
}

{
  # sourceCpp(code='
  # // [[Rcpp::depends(RcppArmadillo)]]
  #           #include <RcppArmadillo.h>
  # 
  #           using namespace Rcpp;
  #           using namespace arma;
  #           
  # arma:: vec matmult_mv(NumericMatrix M, NumericVector c) {
  #           arma:: vec mult= as<arma::mat>(M)*as<arma::vec>(c) ;
  #           return mult ;
  #           };
  #           
  #           double matmult_vv(NumericVector m, NumericMatrix c) {
  #           double mult= arma::as_scalar(as<arma::vec>(m).t()*as<arma::mat>(c)) ;
  #           return mult ;
  #           }
  #           
  #           NumericMatrix dnormpar_mat(NumericMatrix x){
  #           double c = 1/sqrt(2*M_PI);
  #           int nr = x.rows();
  #           int nc= x.cols();
  #           NumericMatrix ret(nr,nc);
  #           #pragma omp parallel for if(nc> 50)
  #           for(int i=0; i<nr; ++i){
  #           for(int j=0; j<nc; ++j){
  #           
  #           ret(i,j) = exp(-x(i,j)*x(i,j)/2)*c;
  #           }
  #           }
  #           return ret;
  #           } 
  #           
  #           NumericMatrix pnorm_ltbl_cpp(NumericMatrix x,NumericVector tbl){
  #           int nr = x.rows();
  #           int nc= x.cols();
  #           NumericMatrix p(nr,nc);
  #           //#pragma omp parallel for if(nc> 50)
  #           for(int i=0; i<nr; ++i){
  #           for(int j=0; j<nc; ++j){
  #           
  #           if (x(i,j)<(-5)) {
  #           p(i,j) = 0;
  #           } else if (x(i,j)>5) {
  #           p(i,j) = 1;
  #           } else {
  #           p(i,j) = tbl[round((x(i,j)-(-5))*pow(10,4))];
  #           
  #           }
  #           }
  #           }
  #           return p;
  #           }
  #           
  #           // [[Rcpp::export]]
  #           double h01(double t,NumericVector tbl, NumericMatrix X, NumericVector b,NumericVector a, NumericVector g, NumericVector V, NumericVector D){
  #           int n=X.rows();
  #           arma:: vec Xb=matmult_mv(X,b);
  #           arma:: vec R= (log(as<arma::vec>(V))-Xb-log(t));
  #           NumericMatrix Rr= as<NumericMatrix>(wrap(R));
  #           NumericMatrix dnormRRn= dnormpar_mat(Rr/a[0]);
  #           double nomi= (matmult_vv(D,dnormRRn))/(n*a[0]*t);
  #           NumericMatrix pnormRRn= pnorm_ltbl_cpp(Rr/a[1],tbl);
  #           double denomi= (matmult_vv(g,pnormRRn))/n;
  #           return nomi/denomi;
  #           }
  #           ')
}

{
  # sourceCpp(code = '#include <Rcpp.h>
  #           //#include <omp.h>
  #           using namespace Rcpp;
  #            // [[Rcpp::export]]
  #            SEXP pnorm_ltbl_cpp(SEXP xx,SEXP table){
  #               NumericMatrix x = as<Rcpp::NumericMatrix>(xx);
  #               NumericVector tbl = as<Rcpp::NumericVector>(table);
  #           int nr = x.rows();
  #           int nc= x.cols();
  #           NumericMatrix p(nr,nc);
  #           #pragma omp parallel for if(nc> 50)
  #           for(int i=0; i<nr; ++i){
  #           for(int j=0; j<nc; ++j){
  # 
  #           if (x(i,j)<(-5)) {
  #           p(i,j) = 0;
  #           } else if (x(i,j)>5) {
  #           p(i,j) = 1;
  #           } else {
  #           p(i,j) = tbl[round((x(i,j)-(-5))*pow(10,4))];
  # 
  #           }
  #           }
  #           }
  #           return p;
  #           }
  #           ')
}

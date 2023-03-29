// Estimating parameters in a Tweedie distribution.
// Modified to include keep vector


#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep, y);
  PARAMETER(mu);
  PARAMETER(ln_phi);
  PARAMETER(logit_p);
  
  Type phi = exp(ln_phi);
  Type p = invlogit(logit_p) + 1;
  Type ans = 0;
  vector<Type> cdf(y.size());
  cdf.fill(0);

  //Compound Poisson-Gamma parameter
  Type lambda = pow(mu,(2-p))/((2-p)*phi);
  
  for(int i=0; i<y.size(); i++){
    ans -= keep(i) * dtweedie(y(i), mu, phi, p, true);

    if(y(i) == 0){
      cdf(i) = squeeze( ppois(y(i),lambda) );
    }
    if(y(i) > 0){
      cdf(i) = squeeze( ppois(Type(0),lambda) );
    }
    
    ans -= keep.cdf_lower(i) * log( cdf(i) );
    ans -= keep.cdf_upper(i) * log( 1.0 - cdf(i) );
    
  }

  REPORT(cdf);
  ADREPORT(phi);
  ADREPORT(p);
  
  return ans;
}

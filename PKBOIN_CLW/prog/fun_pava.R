require(Rcpp)
require(RcppArmadillo)
require(inline)

pava_cpp <- 
"
Rcpp::NumericVector ya(y);
Rcpp::NumericVector wa(w);
int n = ya.size(); 
Rcpp::NumericVector r (n); 
for(int i = 0; i < n; i++)
 r(i) = 1; 
 
int stable = 0; 
double www = 0; 
double ttt = 0; 
Rcpp::NumericVector result (n); 
while(stable == 0)
{
  stable = 1; 
  for(int i = 1; i< n; i++)
  {
    if(ya(i-1) > ya(i))
    {
      stable = 0; 
      www = wa(i-1) + wa(i); 
      ttt = (wa(i-1)*ya(i-1) + wa(i)*ya(i))/www; 
      ya(i) = ttt; 
      wa(i) = www; 
      ya.erase(i-1); 
      wa.erase(i-1); 
      r(i) += r(i-1); 
      r.erase(i-1); 
      n = n-1; 
    }
  }
}

int l = r.size(); 
int k = 1; 
for(int i = 1; i<= l; i++)
{
  for(int j = 1; j <= r(i-1); j++)
  {
    result(k-1) = ya(i-1); 
    k += 1; 
  }
}

return Rcpp::wrap(result); 
"

pava_rcpp <- cxxfunction(signature(y = "numeric", w = "numeric"), 
                         pava_cpp, plugin = "RcppArmadillo", verbose = F)
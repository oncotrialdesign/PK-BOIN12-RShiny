# function for selecting the OBD
fun_TITE_PKBOIN12_OBD <- function(doseDT, pT, u11, u00, zeta1)
{
  used <- which(doseDT$n > 0)
  phat <- rep(pT, nrow(doseDT))
  phat[used] <- doseDT$x[used]/doseDT$n[used]

  ptilde <- pava_rcpp(phat, w = 1/((doseDT$x + 0.05) * (doseDT$n - doseDT$x + 0.05)/((doseDT$n + 0.1)^2 * (doseDT$n + 0.1 + 1)))) + 0.001*(1:nrow(doseDT))
  MTD <- which.min(abs(ptilde - pT))

  ubar <- (u11 * doseDT$y + u00 * (doseDT$n - doseDT$x))/100
  ubar <- (ubar + 1)/(doseDT$n + 2)   

  ubar[doseDT$keep == 0] <- -100
  ubar[doseDT$n == 0] <- -100


  r_hat <- doseDT$r_d[used]
  n_used <- doseDT$n[used]
  x_used <- doseDT$x[used]
  rtilde <- pava_rcpp(r_hat, w = 1/((x_used + 0.05) * (n_used - x_used + 0.05)/((n_used + 0.1)^2 * (n_used + 0.1 + 1)))) + 0.001*(length(used))
  if(sum(rtilde < 1.25 * zeta1) > 0)
  {
    PKmin <-  which.min(abs(rtilde - 1.25 * zeta1)[1:sum(rtilde < 1.25 * zeta1)])
  } else
  {
    PKmin <-  1
  }

  if(PKmin > MTD)
  {
    OBD <- MTD
  } else
  {
    OBD <- which.max(ubar[PKmin:MTD]) + PKmin - 1
  }

  return(OBD)

}

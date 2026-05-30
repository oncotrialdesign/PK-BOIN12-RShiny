# decision table for the dosing decision
fun_TITE_PKBOIN12dec <- function(pT, qE, lambda_e, lambda_d, csize, cN, cutoff_tox=0.95, cutoff_eff = 0.90)
{

  result <- data.frame(n = (1:cN) * csize,
                       DU_T = rep(NA, cN),
                       DU_E = rep(NA, cN),
                       D = rep(NA, cN),
                       E = rep(NA, cN))
  result$D = ceiling(result$n * lambda_d)
  result$E = floor(result$n * lambda_e)

  for(i in 1:cN)
  {
    if(1 - pbeta(pT, result$n[i] + 1, 1) >= cutoff_tox)
      result$DU_T[i] = min(which(sapply(0:result$n[i], FUN = function(x) 1 - pbeta(pT, x + 1, result$n[i] + 1 - x)) >= cutoff_tox)) - 1
    if(pbeta(qE, 1, result$n[i] + 1) >= cutoff_eff)  # otherwise is NA
      result$DU_E[i] = max(which(sapply(0:result$n[i], FUN = function(x) pbeta(qE, x + 1, result$n[i] + 1 - x)) >= cutoff_eff)) - 1
  }
  return(result)
}


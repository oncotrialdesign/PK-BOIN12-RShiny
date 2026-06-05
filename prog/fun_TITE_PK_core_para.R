# core function for PKBOIN-12 and TITE-PKBOIN-12
fun_TITE_PK_core_para <-
  function(index, rM, pM, qM, OBDlist, decisionM, pT, qE, pqcorr, csize, cN, design, mindelta,
           size, lambda_d, lambda_e, zeta1, CV, g_P,
           alphaTO, betaTO, alphaEO, betaEO, current,
           doselimit, ub, u11, u00,
           accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
           use_susp, accrual_random, considerPK  # TITE parameters
           )
{
  rlist <- rM[index,]
  plist <- pM[index,]
  qlist <- qM[index,]
  trueOBD <- OBDlist[[index]]  

  Result <- fun_TITE_PKBOIN12(index, rlist, plist, qlist, trueOBD, pT, qE, pqcorr, lambda_e, lambda_d, zeta1, CV, g_P,
                              csize, cN, decisionM, ub, u11, u00, current, doselimit,
                              accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
                              use_susp, accrual_random, considerPK)

  return(Result)
}

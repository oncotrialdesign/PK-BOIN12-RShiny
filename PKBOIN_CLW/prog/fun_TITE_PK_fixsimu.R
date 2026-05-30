# function for main simulation
fun_TITE_PK_fixsimu <-
  function(dN, rV, pV, qV, pT, qE, pqcorr = 0, psi0PK, CV, g_P, csize, cN, design, utility = FALSE,
           alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,
           alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5, current = 1,
           doselimit = Inf,
           u11 = 100, u00 = 0, cutoff_tox = 0.95, cutoff_eff = 0.9,
           repsize = 10000, n_cores = 10,                                 # replication parameters
           accrual = 10, susp = 0.5, tox_win = 30, eff_win = 60, tox_dist = "Uniform",   #TITE parameters
           eff_dist = "Uniform", tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = TRUE, accrual_random = FALSE,
           considerPK = TRUE)
{

  if(utility == FALSE)
  {
    trueOBD <- findOBD(pV, qV, pT, qE)
  } else
  {
    trueOBD <- findOBD_RDS(pV, qV, pT, qE, u11, u00)
  }

  if(design %in% c("PKBOIN-12", "TITE-PKBOIN-12"))
  {
    phi1 = 0.6 * pT
    phi2 = 1.4 * pT

    lambda_e = log((1-phi1)/(1-pT))/log(pT * (1-phi1)/phi1/(1-pT))  # lower boundary
    lambda_d = log((1-pT)/(1-phi2))/log(phi2 * (1-pT)/pT/(1-phi2))  # upper boundary

    psi1PK = 0.6 * psi0PK
    zeta1 = (psi0PK + psi1PK)/2

    u = 100 * qE * (1 - pT) + u00 * (1 - pT) * (1 - qE) + u11 * pT * qE
    ub = (u + (100 - u)/2)/100

    decisionM <- fun_TITE_PKBOIN12dec(pT, qE, lambda_e, lambda_d, csize, cN, cutoff_tox, cutoff_eff)   # the same as BOIN12

  }


  rM <- kronecker(matrix(rV, nrow = 1), rep(1, repsize))
  pM <- kronecker(matrix(pV, nrow = 1), rep(1, repsize))
  qM <- kronecker(matrix(qV, nrow = 1), rep(1, repsize))

  OBDlist <- lapply(1:repsize, FUN = function(x) trueOBD)

  # replicate the simulation study
  ResultDF <- mclapply(1:repsize, FUN = function(x)
    fun_TITE_PK_core_para(x, rM, pM, qM, OBDlist, decisionM, pT, qE, pqcorr, csize, cN, design, mindelta,
                  size, lambda_d, lambda_e, zeta1, CV, g_P,
                  alphaTO, betaTO, alphaEO, betaEO, current,
                  doselimit, ub, u11, u00,
                  accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
                  use_susp, accrual_random, considerPK
                  ),
    mc.cores = n_cores) %>%
    do.call(rbind, .)

  return(data.frame(n = csize * cN,
                    design = design,
                    trueOBD = trueOBD[1],
                    OBD = ResultDF$OBD,
                    rN = ResultDF$rN,
                    select_OBD = ResultDF$select_OBD,
                    num_at_OBD = ResultDF$num_at_OBD,
                    risk_allocate = ResultDF$risk_allocate,
                    num_overdose_OBD = ResultDF$num_overdose_OBD,
                    num_overdose_nOBD = ResultDF$num_overdose_nOBD,
                    duration = ResultDF$duration                            # TITE
                    ) %>%
           cbind(., ResultDF[, (ncol(ResultDF) - dN + 1):ncol(ResultDF)])
  )
}

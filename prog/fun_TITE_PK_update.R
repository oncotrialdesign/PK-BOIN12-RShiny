# function for determining the dosing determination time for the next cohort and the statistics for the decision
fun_TITE_PK_update <- function(cid, CV, g_P, patDT, doseDT, current, time_current, tox_win, eff_win, csize,
                            tox_dist, eff_dist, accrual, susp = 0.5,
                            u11 = 60, u10 = 0, u01 = 100, u00 = 40,
                            tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = TRUE, accrual_random = FALSE,
                            pqcorr = 0)
{
  dN <- nrow(doseDT)

### Step 1: update the cid'th cohort information in patDT #####
  ###### Accrual (1.1) ######################
  #   1.1 enrollment time for each patient
  if(accrual_random == TRUE)  # can also consider exponential distribution!
  {
    patDT$enroll[patDT$cid == cid] <- time_current + cumsum(runif(csize, 0, 2 * accrual)) + 1    # add 1 here
  } else
  {
    patDT$enroll[patDT$cid == cid] <- time_current + (0:(csize - 1)) * accrual + 1    # add 1 here
  }
  #   1.2 toxicity and efficacy assessment end time
  patDT$toxend[patDT$cid == cid] <- patDT$enroll[patDT$cid == cid] + tox_win
  patDT$effend[patDT$cid == cid] <- patDT$enroll[patDT$cid == cid] + eff_win
  #   1.2.2 generate PK individually
  patDT$PK[patDT$cid == cid] <- truncnorm::rtruncnorm(csize, a = 0, b = Inf, mean = doseDT$PK[current], sd = doseDT$PK[current] * CV)
  #   1.2.3 generate true toxicity and efficacy
  patDT$tox[patDT$cid == cid] <- pmin(pmax(doseDT$tox[current] + doseDT$tox[current] * g_P * (patDT$PK[patDT$cid == cid] - doseDT$PK[current])/doseDT$PK[current], 0), 1)
  patDT$eff[patDT$cid == cid] <- pmin(pmax(doseDT$eff[current] + doseDT$eff[current] * g_P * (patDT$PK[patDT$cid == cid] - doseDT$PK[current])/doseDT$PK[current], 0), 1)

  #   1.3 toxicity and efficacy result at the end (0 or 1)
  if(pqcorr == 0)
  {
    patDT$toxobs[patDT$cid == cid] <- rbinom(n = csize, size = 1, prob = patDT$tox[patDT$cid == cid])
    patDT$effobs[patDT$cid == cid] <- rbinom(n = csize, size = 1, prob = patDT$eff[patDT$cid == cid])
  } else
  {
    joint_prob <- pqcorr * sqrt(patDT$tox[patDT$cid == cid] * (1 - patDT$tox[patDT$cid == cid]) *
                                  patDT$eff[patDT$cid == cid] * (1 - patDT$eff[patDT$cid == cid])) +
      patDT$tox[patDT$cid == cid] * patDT$eff[patDT$cid == cid]
    rand <- runif(csize, min = 0, max = 1)
    patDT$toxobs[patDT$cid == cid] <- as.numeric(rand <= patDT$tox[patDT$cid == cid])
    patDT$effobs[patDT$cid == cid] <- as.numeric((rand <= joint_prob) | (rand > 1 - patDT$eff[patDT$cid == cid] + joint_prob))
  }


  ######## DLT and Response Time (1.4)  Very Important #######
  #   1.4 generate the toxicity and efficacy response time
  new_tox = sum(patDT$cid == cid & patDT$toxobs == 1)   # number of DLT in the current cohort
  new_eff = sum(patDT$cid == cid & patDT$effobs == 1)   # number of response in the current cohort
  if(tox_dist == "Uniform") {
    patDT$toxdat[patDT$cid == cid & patDT$toxobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$toxobs == 1] + runif(new_tox, max = tox_win)
  } else if(tox_dist == "UnifCateg") {
    tox_set <- 1:(tox_win/5-1)
    patDT$toxdat[patDT$cid == cid & patDT$toxobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$toxobs == 1] + sample(tox_set, size = new_tox, replace = T) * 5
  }

  if(eff_dist == "Uniform") {
    patDT$effdat[patDT$cid == cid & patDT$effobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$effobs == 1] + runif(new_eff, max = eff_win)
  } else if(eff_dist == "UnifCateg") {
    eff_set <- 1:(eff_win/5)
    patDT$effdat[patDT$cid == cid & patDT$effobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$effobs == 1] + sample(eff_set, size = new_eff, replace = T) * 5
  }


  #   1.5 update the dose level of the current cohort
  patDT$d[patDT$cid == cid] <- current

  #   1.6 update tox_confirm, eff_confirm
  patDT$tox_confirm[patDT$cid == cid] <- pmin(patDT$toxend[patDT$cid == cid], patDT$toxdat[patDT$cid == cid])
  patDT$eff_confirm[patDT$cid == cid] <- pmin(patDT$effend[patDT$cid == cid], patDT$effdat[patDT$cid == cid])


### Step 2: derive time_next ####

  tmp <- patDT %>% filter(d == current)
  n_d <- nrow(tmp)                 # number of patients assigned to dose level d (including current cohort)

  if(use_susp == TRUE)   # consider suspension rule
  {
    min_n <- min(floor(n_d * susp + 1), n_d)        # susp default value = 0.5 in (0, 1]; if susp = 1, it is the same as non-TITE designs
    tox_next <- sort(tmp$tox_confirm)[min_n]
    eff_next <- sort(tmp$eff_confirm)[min_n]
    if(accrual_random == TRUE)
    {
      time_next <- max(tox_next, eff_next, max(patDT$enroll, na.rm = T) + runif(1, 0, 2 * accrual) + 1)
    } else
    {
      time_next <- max(tox_next, eff_next, time_current + csize * accrual + 1)
    }
    # the earliest time fitting the accrual suspension rule and also the first patient in the next cohort is ready for enrollment
  } else
  {
    tox_next <- max(tmp$tox_confirm)
    eff_next <- max(tmp$eff_confirm)
    if(accrual_random == TRUE)
    {
      time_next <- max(tox_next, eff_next, max(patDT$enroll, na.rm = T)) + runif(1, 0, 2 * accrual) + 1
    } else
    {
      time_next <- max(tox_next, eff_next, time_current + csize * accrual + 1)
    }
  }

### Step 3: update doseDT ####

  patDT_current <- patDT %>% filter(!is.na(d)) %>%
    mutate(delta_t = -1, delta_e = -1, t_i = time_next - enroll)

  # 3.1  delta_t, delta_e (for each patient)
  patDT_current$delta_t[patDT_current$toxdat <= time_next] <- 1     # DLT observed
  patDT_current$delta_t[patDT_current$toxend <= time_next & patDT_current$toxdat > time_next] <- 0 # no DLT (DLT assessment end before time_next, and no DLT observed)
    # if delta_t = -1, the toxicity assessment has not finished at time_next

  patDT_current$delta_e[patDT_current$effdat <= time_next] <- 1     # response observed
  patDT_current$delta_e[patDT_current$effend <= time_next & patDT_current$effdat > time_next] <- 0 # no response
    # if delta_e = -1, the efficacy assessment has not finished at time_next

  ####  weight adjusting accounts for the partial information (3.2) very important #######
  # 3.2  w_t, w_e (for each patient)
  patDT_current <- patDT_current %>%
    mutate(w_t = ifelse(t_i < tox_win & delta_t < 0, t_i/tox_win, 0),     # t_i/A_t
           w_e = ifelse(t_i < eff_win & delta_e < 0, t_i/eff_win, 0))     # t_i/A_e



  # 3.3 ESS_t, ESS_e, pi_t_hat, pi_e_hat 
  for(i in 1:dN)
  {
    if(sum(patDT_current$d == i) > 0)   # at least one cohort has been assigned to dose i
    {
      # update n, x, y
      doseDT$n[i] <- sum(patDT_current$d == i)
      delta_t_i = patDT_current$delta_t[patDT_current$d == i]
      delta_e_i = patDT_current$delta_e[patDT_current$d == i]

      doseDT$x[i] <- sum(delta_t_i == 1)    # number of observed DLT
      doseDT$y[i] <- sum(delta_e_i == 1)    # number of observed responses

      # update the effective sample size
      w_t_i = patDT_current$w_t[patDT_current$d == i]
      w_e_i = patDT_current$w_e[patDT_current$d == i]

      doseDT$ESS_t[i] <- sum(delta_t_i == 1) + sum(delta_t_i == 0) + sum(w_t_i)
      doseDT$pi_t_hat[i] <- doseDT$x[i]/doseDT$ESS_t[i]     # update unless n[i] > 0

      doseDT$ESS_e[i] <- sum(delta_e_i == 1) + sum(delta_e_i == 0) + sum(w_e_i)
      doseDT$pi_e_hat[i] <- doseDT$y[i]/doseDT$ESS_e[i]     # update unless n[i] > 0
    }
  }

  # 3.4 tox_exp, eff_exp (for each patient)
  patDT_current <- patDT_current %>% mutate(tox_exp = delta_t, eff_exp = delta_e)

  patDT_current$tox_exp[patDT_current$delta_t == -1] <-
    (doseDT$pi_t_hat[patDT_current$d] * (1 - patDT_current$w_t) / (1 - doseDT$pi_t_hat[patDT_current$d] * patDT_current$w_t))[patDT_current$delta_t == -1]

  patDT_current$eff_exp[patDT_current$delta_e == -1] <-
    (doseDT$pi_e_hat[patDT_current$d] * (1 - patDT_current$w_e) / (1 - doseDT$pi_e_hat[patDT_current$d] * patDT_current$w_e))[patDT_current$delta_e == -1]

  #### quasi-event x_d for TITE-BOIN12 ####
  # 3.5 x_d (for each dose level)
  for(i in 1:dN) {
    if(sum(patDT_current$d == i) > 0)   # at least one cohort has been assigned to dose i
    {
      tox_exp_i <- patDT_current$tox_exp[patDT_current$d == i]
      eff_exp_i <- patDT_current$eff_exp[patDT_current$d == i]
      PK_i <- patDT_current$PK[patDT_current$d == i]

      doseDT$x_d[i] <-
        (u11 * sum(tox_exp_i * eff_exp_i) + u10 * sum(tox_exp_i * (1 - eff_exp_i)) +
        u01 * sum((1 - tox_exp_i) * eff_exp_i) + u00 * sum((1 - tox_exp_i) * (1 - eff_exp_i)))/100  # u10: DLT + no response
      doseDT$r_d[i] <- mean(PK_i)
      doseDT$r_sd[i] <- sd(PK_i)
    }
  }

### Step 4 return the outputs #####
  return(list(patDT = patDT, doseDT = doseDT, time_next = time_next ))

}




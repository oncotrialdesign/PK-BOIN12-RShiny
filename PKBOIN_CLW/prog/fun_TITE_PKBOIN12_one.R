# function for making the dosing decision
TITE_PKBOIN12_one <- function(doseDT, current, pT, qE, lambda_e, lambda_d, zeta1, csize, decisionM, ub, u11, u00)
{
  n <- doseDT$n[current]
  phat <- doseDT$pi_t_hat[current]        
  qhat <- doseDT$pi_e_hat[current]
  rhat <- doseDT$r_d[current]        # observed PK value
  cn <- round(n/csize)   # number of cohort assigned to the dose level
  dN <- nrow(doseDT)
  
  above <- ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # NA or number
  below <- ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1])   # NA or number
  
  PKmin <-  suppressWarnings(min(which(doseDT$r_d > zeta1)))   # Inf means no dose has PK value > zeta1
  
  lower_exist <- FALSE
  if(!is.infinite(PKmin) & !is.na(below))
  {
    if(PKmin < below)
    {
      lower_exist <- TRUE      # TRUE means the admissible set need to include doses {PKmin, ... d-2}
    }
  }
  
  if(phat * n >= decisionM$DU_T[cn])    # DU_T  
  {
    doseDT$keep[current:dN] = 0
    newdose <- below
    return(list(doseDT = doseDT, newdose = newdose))    # skip the following part
  }
  
  if(is.na(decisionM$DU_E[cn]) == FALSE)   # futility
  {
    if(qhat * n <= decisionM$DU_E[cn])    
      doseDT$keep[current] = 0
  }
  
  if(is.na(above) == FALSE)
  {
    if(n >= 9 & phat < lambda_d & doseDT$n[above] == 0)  # dose exploration
    {
      newdose <- above
      return(list(doseDT = doseDT, newdose = newdose))    # skip the following part
    }
  }
  
  if(current == dN & n >= 6)    # concentration
  {
    if(pnorm(1.25 * zeta1, mean = doseDT$r_d[current], sd = doseDT$r_sd[current]/sqrt(n)) > 0.95)  # 1.25 *
    {
      doseDT$keep <- 0
      newdose <- NA
      return(list(doseDT = doseDT, newdose = newdose))
    }
  }
  
  if(n >= 6)    # concentration
  {
    if(pnorm(1.25 * zeta1, mean = doseDT$r_d[current], sd = doseDT$r_sd[current]/sqrt(n)) > 0.95)
    {
      if(sum(doseDT$keep) > 0 & which.max(doseDT$keep) < current & current >= 2)
      {
        doseDT$keep[which.max(doseDT$keep)] <- 0
      }
    }
  }

  if(sum(doseDT$keep == 0) == dN)   # no dose left
  {
    newdose <- NA
    return(list(doseDT = doseDT, newdose = newdose))
  }
  
  # posterior prob, Inf means eliminated
  if(doseDT$keep[current] == 0)   # eliminate current dose
  {
    post_current <- -Inf
  } else {
    xd_current <- doseDT$x_d[current]        
    post_current <- 1 - pbeta(ub, 1 + xd_current, n + 1 - xd_current)
  }
  
  if(is.na(below))
  {
    post_below <- -Inf
  } else {
    xd_below <- doseDT$x_d[below]             
    post_below <- 1 - pbeta(ub, 1 + xd_below, doseDT$n[below] + 1 - xd_below)
  }
  
  if(is.na(above))
  {
    post_above <- -Inf
  } else{
    xd_above <-  doseDT$x_d[above]           
    post_above <- 1 - pbeta(ub, 1 + xd_above, doseDT$n[above] + 1 - xd_above)
  }
  
  lower_list <- NULL
  if(lower_exist == TRUE)
  {
    lower_list <- (PKmin:(below - 1))[doseDT$keep[PKmin:(below - 1)] == 1]
    if(length(lower_list) > 0)
    {
      xd_lower <- doseDT$x_d[lower_list]
      post_lower <- 1 - pbeta(ub, 1 + xd_lower, doseDT$n[lower_list] + 1 - xd_lower)
    }
  }
  
  if(lambda_d <= phat)  # de-escalate
  {
    if(!is.na(below))
    {
      if(rhat > zeta1 & length(lower_list) > 0)  # {PKmin, ... , d-1}
      {
        remain <- c(lower_list, below)
        post <- c(post_lower, post_below)
        newdose <- remain[which.max(post + (1:length(post)) * 1e-6)]
      } else {                                   # d-1
        newdose <- below
      }
      return(list(doseDT = doseDT, newdose = newdose))
    } else {
      if(is.na(current))                         # NA
      {
        newdose <- NA
        return(list(doseDT = doseDT, newdose = newdose))
      } else                                     # current
      {
        newdose <- current    # if current exist, then we use current
        return(list(doseDT = doseDT, newdose = newdose))
      }
    }
  } else   # not de-escalate
  {
    post <- c(post_below, post_current)
    remain <- c(below, current)
    if((n < 6 & !is.na(above)) | (n >= 6 & phat <= lambda_e & !is.na(above)))
    {
      post <- c(post, post_above)
      remain <- c(remain, above)
    }
    if(rhat > zeta1 & length(lower_list) > 0)  # include {PKmin, ... , d-2}
    {
      remain <- c(lower_list, remain)
      post <- c(post_lower, post)
    }
    
    newdose <- remain[which.max(post + (1:length(post)) * 1e-6)]
    if(max(post) == -Inf)
    {
      return(list(doseDT = doseDT, newdose = NA))
    } else
    {
      return(list(doseDT = doseDT, newdose = newdose))
    }
  }
  
}

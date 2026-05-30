# function for identifying the true OBD with RDS
findOBD_RDS <- function(tox, eff, pT, qE, u11, u00)
{
  MTD <- sum(tox <= pT)
  if(MTD == 0)
  {
    return(-1)
  } else if(max(eff[1:MTD]) < qE)
  {
    return(-1)
  } else
  {
    RDS <- 100 * eff * (1 - tox) + u00 * (1-tox) * (1-eff) + u11 * tox * eff
    OBDlist <- which.max(RDS[1:MTD])
    return(OBDlist)
  }
}

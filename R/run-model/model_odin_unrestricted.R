## GONORRHOEA TRANSMISSION MODEL
## STRATIFIED BY SEX, RISK, AGE, GESTATIONAL STATUS FOR HETEROSEXUAL POPULATION

# Difference from model_odin:
# Instead of testing proportion of those accessing syndromic management or screening, 
# intervention set up to test entire target population

# Main changes made:
# line 449 when defining delta_init, which has been set to 1 for each target pop
# line 537 when defining tau_hc, set as 1 for each target group (tested once per eligible person per year)

n_sex  <- 2  # male = 1, female = 2
n_risk <- 3 # low risk = 1, med risk = 2, high risk = 3
n_age  <- 2  # young (15-24 years) = 1, old (25+ years) = 2
n_gest <- 2 # non-pregnant = 1, pregnant = 2

## CORE EQUATIONS --------------------------------------------------------------

# Transmission between compartments

# Never sexually active, not susceptible
deriv(U[,,,]) <- ent[i,j,k,l]  -
  psi[i,j,k,l] * U[i,j,k,l] - 
  # aging
  omega[k] * U[i,j,k,l] + (if (k==2) omega[k-1] * U[i,j,k-1,l] else 0) -
  mu[i,k] * U[i,j,k,l]
# no pregnancy for non-sexually active

# Ever sexually active, susceptible
deriv(S[,,,]) <- psi[i,j,k,l] * U[i,j,k,l] +
  pi[i] * X[i,j,k,l]  + pi[i] * Y[i,j,k,l]  + pi[i] * Z[i,j,k,l] +
  pi[i] * Ry[i,j,k,l] + 
  pi[i] * Rh_symp[i,j,k,l] + pi[i] * Rh_asymp[i,j,k,l] +
  sigma[i] * Ty[i,j,k,l] + 
  sigma[i] * Th_symp[i,j,k,l] + sigma[i] * Th_asymp[i,j,k,l] + 
  sigma[i] * My[i,j,k,l] -
  lambda[i,j,k,l] * S[i,j,k,l] -
  # aging
  omega[k] * S[i,j,k,l] + (if (k==2) omega[k-1] * S[i,j,k-1,l] else 0) +
  # pregnancy (f[i,j,k,l] = 0 when i==1 and l==2)
  (if (l==1) - f[i,j,k,l] * S[i,j,k,l] + alpha * S[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * S[i,j,k,l-1] - alpha * S[i,j,k,l] else 0) -
  # death
  mu[i,k] * S[i,j,k,l]

# Asymptomatic NG 
deriv(X[,,,]) <- lambda[i,j,k,l] * (1 - phi[i]) * S[i,j,k,l] -
  tau_hc[i,j,k,l] * X[i,j,k,l] - 
  pi[i] * X[i,j,k,l] -
  omega[k] * X[i,j,k,l] + (if (k==2) omega[k-1] * X[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * X[i,j,k,l] + alpha * X[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * X[i,j,k,l-1] - alpha * X[i,j,k,l] else 0) -
  mu[i,k] * X[i,j,k,l]

# Symptomatic NG that will be treated (STI clinic)
deriv(Y[,,,]) <- lambda[i,j,k,l] * phi[i] * gamma[i,j,k,l] * S[i,j,k,l] -
  tau_m[i,j,k,l] * Y[i,j,k,l] - 
  tau_hc[i,j,k,l] * Y[i,j,k,l] -
  pi[i] * Y[i,j,k,l] -
  omega[k] * Y[i,j,k,l] + (if (k==2) omega[k-1] * Y[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Y[i,j,k,l] + alpha * Y[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Y[i,j,k,l-1] - alpha * Y[i,j,k,l] else 0) -
  mu[i,k] * Y[i,j,k,l]

# Symptomatic NG that will not be treated (STI clinic)
deriv(Z[,,,]) <- lambda[i,j,k,l] * phi[i] * (1 - gamma[i,j,k,l]) * S[i,j,k,l] -
  tau_hc[i,j,k,l] * Z[i,j,k,l] -
  pi[i] * Z[i,j,k,l] -
  omega[k] * Z[i,j,k,l] + (if (k==2) omega[k-1] * Z[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Z[i,j,k,l] + alpha * Z[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Z[i,j,k,l-1] - alpha * Z[i,j,k,l] else 0) -
  mu[i,k] * Z[i,j,k,l]

# Symptomatic NG: syndromic diagnosis and correctly treated (STI clinic)
deriv(My[,,,]) <- tau_m[i,j,k,l] * zeta[i] * (1 - delta[i,j,k,l]) * Y[i,j,k,l] - 
  sigma[i] * My[i,j,k,l] -
  omega[k] * My[i,j,k,l] + (if (k==2) omega[k-1] * My[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * My[i,j,k,l] + alpha * My[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * My[i,j,k,l-1] - alpha * My[i,j,k,l] else 0) -
  mu[i,k] * My[i,j,k,l]

# Symptomatic NG: POCT diagnosis and correctly treated (STI clinic)
deriv(Ty[,,,]) <- tau_m[i,j,k,l] * eta[i] * delta[i,j,k,l] * Y[i,j,k,l] - 
  sigma[i] * Ty[i,j,k,l] -
  omega[k] * Ty[i,j,k,l] + (if (k==2) omega[k-1] * Ty[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Ty[i,j,k,l] + alpha * Ty[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Ty[i,j,k,l-1] - alpha * Ty[i,j,k,l] else 0) -
  mu[i,k] * Ty[i,j,k,l]

# Symptomatic NG: not correctly treated at STI clinic
deriv(Ry[,,,]) <- tau_m[i,j,k,l] * (1 - zeta[i] * (1 - delta[i,j,k,l]) -
                                      eta[i] * delta[i,j,k,l]) * Y[i,j,k,l] - 
  pi[i] * Ry[i,j,k,l] -
  omega[k] * Ry[i,j,k,l] + (if (k==2) omega[k-1] * Ry[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Ry[i,j,k,l] + alpha * Ry[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Ry[i,j,k,l-1] - alpha * Ry[i,j,k,l] else 0) -
  mu[i,k] * Ry[i,j,k,l]

# Any NG: POCT diagnosis and correctly treated (other HC setting)
deriv(Th_symp[,,,]) <- tau_hc[i,j,k,l] * eta[i] * Y[i,j,k,l] +
  tau_hc[i,j,k,l] * eta[i] * Z[i,j,k,l] - 
  sigma[i] * Th_symp[i,j,k,l] -
  omega[k] * Th_symp[i,j,k,l] + (if (k==2) omega[k-1] * Th_symp[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Th_symp[i,j,k,l] + alpha * Th_symp[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Th_symp[i,j,k,l-1] - alpha * Th_symp[i,j,k,l] else 0) -
  mu[i,k] * Th_symp[i,j,k,l]

deriv(Th_asymp[,,,]) <- tau_hc[i,j,k,l] * eta[i] * X[i,j,k,l] -
  sigma[i] * Th_asymp[i,j,k,l] -
  omega[k] * Th_asymp[i,j,k,l] + (if (k==2) omega[k-1] * Th_asymp[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Th_asymp[i,j,k,l] + alpha * Th_asymp[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Th_asymp[i,j,k,l-1] - alpha * Th_asymp[i,j,k,l] else 0) -
  mu[i,k] * Th_asymp[i,j,k,l]

# Any NG: POCT diagnosis and not correctly treated (other HC setting)
deriv(Rh_symp[,,,]) <-  tau_hc[i,j,k,l] * (1 - eta[i]) * Y[i,j,k,l] +
  tau_hc[i,j,k,l] * (1 - eta[i]) * Z[i,j,k,l] - 
  pi[i] * Rh_symp[i,j,k,l] -
  omega[k] * Rh_symp[i,j,k,l] + (if (k==2) omega[k-1] * Rh_symp[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Rh_symp[i,j,k,l] + alpha * Rh_symp[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Rh_symp[i,j,k,l-1] - alpha * Rh_symp[i,j,k,l] else 0) -
  mu[i,k] * Rh_symp[i,j,k,l]

deriv(Rh_asymp[,,,]) <- tau_hc[i,j,k,l] * (1 - eta[i]) * X[i,j,k,l] -
  pi[i] * Rh_asymp[i,j,k,l] -
  omega[k] * Rh_asymp[i,j,k,l] + (if (k==2) omega[k-1] * Rh_asymp[i,j,k-1,l] else 0) +
  (if (l==1) - f[i,j,k,l] * Rh_asymp[i,j,k,l] + alpha * Rh_asymp[i,j,k,l+1] else 0) +
  (if (l==2) f[i,j,k,l-1] * Rh_asymp[i,j,k,l-1] - alpha * Rh_asymp[i,j,k,l] else 0) -
  mu[i,k] * Rh_asymp[i,j,k,l]

# Other symptomatic RTI that will be treated (STI clinic) (NO LONGER IMPLEMENTED)
# deriv(O[,,,]) <- 
#   #rate_other[i,j,k,l] * NS[i,j,k,l] - 
#   rate_other[i,j,k,l] * S[i,j,k,l] -
#   pi[i] * O[i,j,k,l] -
#   tau_m[i,j,k,l] * O[i,j,k,l] -
#   omega[k] * O[i,j,k,l] + (if (k==2) omega[k-1] * O[i,j,k-1,l] else 0) +
#   (if (l==1) - f[i,j,k,l] * O[i,j,k,l] + alpha * O[i,j,k,l+1] else 0) +
#   (if (l==2) f[i,j,k,l-1] * O[i,j,k,l-1] - alpha * O[i,j,k,l] else 0) -
#   mu[i,k] * O[i,j,k,l]

# Other symptomatic RTI: syndromic diagnosis and correctly treated (STI clinic) (NO LONGER IMPLEMENTED)
# deriv(Mo[,,,]) <- tau_m[i,j,k,l] * sens_o[i] * O[i,j,k,l] - 
#   sigma[i] * Mo[i,j,k,l] -
#   omega[k] * Mo[i,j,k,l] + (if (k==2) omega[k-1] * Mo[i,j,k-1,l] else 0) +
#   (if (l==1) - f[i,j,k,l] * Mo[i,j,k,l] + alpha * Mo[i,j,k,l+1] else 0) +
#   (if (l==2) f[i,j,k,l-1] * Mo[i,j,k,l-1] - alpha * Mo[i,j,k,l] else 0) -
#   mu[i,k] * Mo[i,j,k,l]

# Other symptomatic RTI: syndromic diagnosis and not correctly treated (STI clinic) (NO LONGER IMPLEMENTED)
# deriv(Ro[,,,]) <- tau_m[i,j,k,l] * (1 - sens_o[i]) * O[i,j,k,l] - 
#   pi[i] * Ro[i,j,k,l] -
#   omega[k] * Ro[i,j,k,l] + (if (k==2) omega[k-1] * Ro[i,j,k-1,l] else 0) +
#   (if (l==1) - f[i,j,k,l] * Ro[i,j,k,l] + alpha * Ro[i,j,k,l+1] else 0) +
#   (if (l==2) f[i,j,k,l-1] * Ro[i,j,k,l-1] - alpha * Ro[i,j,k,l] else 0) -
#   mu[i,k] * Ro[i,j,k,l]

## < Cumulative numbers ----

# Number individuals moving into each compartment per year
## NG
# new infections: S -> I
n_X[,,,] <- lambda[i,j,k,l] * (1 - phi[i]) * S[i,j,k,l] 
n_Y[,,,] <- lambda[i,j,k,l] * phi[i] * gamma[i,j,k,l] * S[i,j,k,l]    
n_Z[,,,] <- lambda[i,j,k,l] * phi[i] * (1 - gamma[i,j,k,l]) * S[i,j,k,l]    
n_I[,,,] <- n_X[i,j,k,l] + n_Y[i,j,k,l] + n_Z[i,j,k,l]
# receive POCT + correct treatment: I -> T
n_Th_symp[,,,] <- tau_hc[i,j,k,l] * eta[i] * (Y[i,j,k,l] + Z[i,j,k,l])
n_Th_asymp[,,,] <- tau_hc[i,j,k,l] * eta[i] * X[i,j,k,l]
n_Ty[,,,] <- tau_m[i,j,k,l] * eta[i] * delta[i,j,k,l] * Y[i,j,k,l]
# receive syndromic diagnosis + correct treatment (I -> M)
n_My[,,,] <- tau_m[i,j,k,l] * zeta[i] * (1 - delta[i,j,k,l]) * Y[i,j,k,l]
# receive incorrect diagnosis + not treated (I -> R)
n_Rh_symp[,,,] <- tau_hc[i,j,k,l] * (1 - eta[i]) * (Y[i,j,k,l] + Z[i,j,k,l])
n_Rh_asymp[,,,] <- tau_hc[i,j,k,l] * (1 - eta[i]) * X[i,j,k,l]
n_Ry[,,,] <- tau_m[i,j,k,l] * (1 - zeta[i] * (1 - delta[i,j,k,l]) - eta[i] * delta[i,j,k,l]) * Y[i,j,k,l]


## Other symptomatic RTI receiving treatment (NO LONGER IMPLEMENTED)
# new infections 
# n_O[,,,] <- rate_other[i,j,k,l] * NS[i,j,k,l] 
# n_O[,,,] <- rate_other[i,j,k,l] * S[i,j,k,l] # assumes no co-infection for ng and other symp infections
# receiving syndromic diagnosis + correct treatment
# n_Mo[,,,] <- tau_m[i,j,k,l] * sens_o[i] * O[i,j,k,l]
# receiving incorrect diagnosis and not treated
# n_Ro[,,,] <- tau_m[i,j,k,l] * (1 - sens_o[i]) * O[i,j,k,l]

# Cumulative number of individuals moving into each compartment
# deriv(cum_X[,,,])  <- n_X[i,j,k,l]
# deriv(cum_Y[,,,])  <- n_Y[i,j,k,l]
# deriv(cum_Z[,,,])  <- n_Z[i,j,k,l]
# deriv(cum_Th_symp[,,,]) <- n_Th_symp[i,j,k,l]
# deriv(cum_Th_asymp[,,,]) <- n_Th_asymp[i,j,k,l]
# deriv(cum_Ty[,,,]) <- n_Ty[i,j,k,l]
# deriv(cum_My[,,,]) <- n_My[i,j,k,l]
# deriv(cum_Rh_symp[,,,]) <- n_Rh_symp[i,j,k,l]
# deriv(cum_Rh_asymp[,,,]) <- n_Rh_asymp[i,j,k,l]
# deriv(cum_Ry[,,,]) <- n_Ry[i,j,k,l]

# deriv(cum_O[,,,])  <- n_O[i,j,k,l]
# deriv(cum_Mo[,,,]) <- n_Mo[i,j,k,l]
# deriv(cum_Ro[,,,]) <- n_Ro[i,j,k,l]

# Person-years spent in each compartment
deriv(py_X[,,,])  <- X[i,j,k,l]
deriv(py_Y[,,,])  <- Y[i,j,k,l]
deriv(py_Z[,,,])  <- Z[i,j,k,l]
deriv(py_Th_symp[,,,]) <- Th_symp[i,j,k,l]
deriv(py_Th_asymp[,,,]) <- Th_asymp[i,j,k,l]
deriv(py_Ty[,,,]) <- Ty[i,j,k,l]
deriv(py_My[,,,]) <- My[i,j,k,l]
deriv(py_Rh_symp[,,,]) <- Rh_symp[i,j,k,l]
deriv(py_Rh_asymp[,,,]) <- Rh_asymp[i,j,k,l]
deriv(py_Ry[,,,]) <- Ry[i,j,k,l]

## ADDITIONAL EQUATIONS --------------------------------------------------------

## < Population size ----

# total infections 
I[,,,] <- X[i,j,k,l] + Y[i,j,k,l] + Z[i,j,k,l] +
  Ty[i,j,k,l] + Ry[i,j,k,l] + My[i,j,k,l] +
  Th_symp[i,j,k,l] + Th_asymp[i,j,k,l] + 
  Rh_symp[i,j,k,l] + Rh_asymp[i,j,k,l]

# sexually active population size
NS[,,,] <- S[i,j,k,l] + I[i,j,k,l]

# total population size
N[,,,] <- U[i,j,k,l] + S[i,j,k,l] + I[i,j,k,l]

## < Demographic parameters ----

## Entrants to never sexually active population (U) by sex, risk, age, and gest
# enter as age = 1 (younger) and gest = 1 (non-pregnant)

# number entrants by sex
# calculated from number turning age 15 years
ent_m <- interpolate(ent_t, entm_y, "linear") 
ent_f <- interpolate(ent_t, entf_y, "linear") 

# calculate entrants per risk group, according to proportion in each risk group
ent[,,,] <- 
  # men, young, non-preg
  if (i==1 && k==1 && l==1) ent_m * upsilon[i,j] else 
    # women, young, non-preg
    if (i==2 && k==1 && l==1) ent_f * upsilon[i,j] else
      # older age group and pregnant 
      0

## Pregnancy rate by sex and age

# pregnancy rate calculated for full population
f_y <- interpolate(f_t, fy_y, "linear")
f_o <- interpolate(f_t, fo_y, "linear")

# convert pregnancy rate among full pop to pregnancy rate among sexually active pop
# i.e. originally calculated for full population but need to account for having 
# only sexually active individuals get pregnant
# f_input = pregnant/N
# f_model = pregnant/S = pregnant/(q.N)
# therefore: f_input.N = f_model.q.N --> f_model = f_input / q
f[,,,] <- # women, young
  if (i==2 && k==1) f_y / q[i,j,k,l] else 
    # women, old
    # q = 1 for older pop (all sexually active)
    if (i==2 && k==2) f_o / q[i,j,k,l] else 
      # men
      0

## Rate of becoming sexually active

# if q = 1 (100% of pop is sexually active)
# move from U -> S at fast rate to maintain proportion of 100% (NS/N)

# if q < 1 (varies by risk group)
# move from U -> S at rate that approximates proportion q 
# manually calibrated

psi[,,,] <- if(q[i,j,k,l] == 1) 100 else -0.3*log(1 - q[i,j,k,l])

## < Sexual partnership probability ----

# rho[i,j,k,l,i5,i6,i7,i8] has 8 dimensions, defined as:
# i: my sex      (1st dim)
# j: my risk     (2nd dim)
# k: my age      (3rd dim)
# l: my gest     (etc... )
# i5: their sex
# i6: their risk
# i7: their age 
# i8: their gest

# Conditions for which sexual partnerships are possible between srag and s'r'a'g'
rho_set[,,,,,,,] <- 
  # no same sex partnerships
  if(i == i5) 0 else 
    # no partnerships between FSW and low or med risk men
    # if i = 1, j = 1 or 2, i5 = 2, i6 = 3      -> rho = 0
    if ((i == 1 && (j == 1 || j == 2) && i5 == 2 && i6 == 3) ||
        (i == 2 && j == 3 && i5 == 1 && (i6 == 1 || i6 == 2))) 0 else 
          # no partnerships with pregnant men
          if ((i == 1 && l == 2) || (i5 == 1 && i8 == 2)) 0 else
            # otherwise, partnerships are allowed
            1

# Number of sexual partnerships offered by s'r'a'g' (p2)
# rho_sp refers to partner, but needs 4 dimensions (first 4 dim: i,j,k,l) 
rho_sp[,,,] <- c[i,j,k,l] * NS[i,j,k,l] 

# Number of sexual partnerships offered to srag (p1) by s'r'a'g' (p2)
# take conditions in rho_set into account
# allocate total sp offered by p2 to any interaction with p2 
rho_sub[,,,,,,,] <- rho_set[i,j,k,l,i5,i6,i7,i8] * rho_sp[i5,i6,i7,i8]

# Probability of sexual partnerships between srag and s'r'a'g'
# numerator: partnerships offered to srag (p1) by s'r'a'g' (p2) according to conditions
# denominator: total partnerships offered to srag by s' according to conditions
# p2 offers partnerships to multiple p1 groups
# denom accounts for this by summing across total
rho[,,,,,,,] <-
  # if statement avoids dividing by 0
  if (rho_sub[i,j,k,l,i5,i6,i7,i8] == 0) 0 else
    rho_sub[i,j,k,l,i5,i6,i7,i8] / sum(rho_sub[i,j,k,l,i5,,,])

## < Partnership balancing ----

B[,,,,,,,] <-
  # if statement avoids diving by 0
  if (rho[i,j,k,l,i5,i6,i7,i8] == 0 || rho[i5,i6,i7,i8,i,j,k,l] == 0) 0 else
    (c[i,j,k,l] * rho[i,j,k,l,i5,i6,i7,i8] * NS[i,j,k,l]) /
  (c[i5,i6,i7,i8] * rho[i5,i6,i7,i8,i,j,k,l] * NS[i5,i6,i7,i8])

# theta_c: extent to which compromise balancing for women
# If theta_c = 1: c_ad = c for p1, c_ad changes for p2
# If theta_c = 0: c_ad = c for p2, c_ad changes for p1

bal_exp[] <- theta_c - 1

# Adjusted total partners for p1
c_ad[,,,,,,,] <-
  # if statement avoids dividing by 0 [0^(-0.5) = Inf]
  if (B[i,j,k,l,i5,i6,i7,i8] == 0) 0 else
    c[i,j,k,l] * B[i,j,k,l,i5,i6,i7,i8] ^ bal_exp[i]

# Adjusted total partners for p1, stratified across p2 
c_ad_strat[,,,,,,,] <- c_ad[i,j,k,l,i5,i6,i7,i8] * rho[i,j,k,l,i5,i6,i7,i8]

# Number of sex partners per year
sp[,,,,,,,] <- c_ad[i,j,k,l,i5,i6,i7,i8] * rho[i,j,k,l,i5,i6,i7,i8] * NS[i,j,k,l]

## < Transmission probability ----

# Number sex acts per partnership per year
# Distribute total sex acts across adjusted number of partners 
# e.g. low risk: 52 sex acts/year and 1.2 partners/year = 43.3 sex acts/partner
# e.g. med risk: 57.2 sex acts/year and 3.7 partners/year = 15.5 sex acts/partner
# Sex acts per partner = n_weekly * 52 / sum(c_ad_strat)

# distribute sex acts evenly across partners
n[,,,,,,,] <- 
  # if c_ad_strat = 0, then no partnership and no sex acts
  if (c_ad_strat[i,j,k,l,i5,i6,i7,i8] == 0) 0 else
    # if low risk
    if (j == 1) n_l * 52 / sum(c_ad_strat[i,j,k,l,i5,,,]) else
      # if med risk
      if (j == 2) n_m * 52 / sum(c_ad_strat[i,j,k,l,i5,,,]) else  
        # if FSW and client, 1 sex act per partnership
        if (j == 3 && i6 == 3) 1 else
          # if high risk men, distribute total sex acts across low and med risk partners 
          # FSW sex acts are then additional to the low and med risk sex acts
          if (i == 1 && j == 3 && (i6 == 1 || i6 == 2)) 
            n_h * 52 / (sum(c_ad_strat[i,j,k,l,i5,1,,]) + sum(c_ad_strat[i,j,k,l,i5,2,,])) else
              0

# Consistency of condom use per partnership 
# average of what is reported by each partner type, weighted in favour of women's reporting

chi_avg[1,,,,,,,] <- rho_set[i,j,k,l,i5,i6,i7,i8] * 
  (chi[i,j,k,l] * (1 - theta_chi) + chi[i5,i6,i7,i8] * theta_chi)

chi_avg[2,,,,,,,] <- rho_set[i,j,k,l,i5,i6,i7,i8] * 
  (chi[i,j,k,l] * theta_chi + chi[i5,i6,i7,i8] * (1 - theta_chi))

# Per partnership transmission probability
kappa_p[,,,,,,,] <- 1 -
  ((1 - (1 - e) * kappa[i]) ^ (n[i,j,k,l,i5,i6,i7,i8] * chi_avg[i,j,k,l,i5,i6,i7,i8])) *
  ((1 - kappa[i]) ^ (n[i,j,k,l,i5,i6,i7,i8] * (1 - chi_avg[i,j,k,l,i5,i6,i7,i8])))

## < Force of infection ----

lambda_sub[,,,,,,,] <-
  # if statement avoids dividing by 0 (NS = 0 for pregnant males)
  if (i5 == 1 && i8 == 2) 0 else
    kappa_p[i,j,k,l,i5,i6,i7,i8] * c_ad[i,j,k,l,i5,i6,i7,i8] * rho[i,j,k,l,i5,i6,i7,i8] * I[i5,i6,i7,i8] / NS[i5,i6,i7,i8]

lambda[,,,] <- sum(lambda_sub[i,j,k,l,,,,])

## < Diagnosis and treatment ----

# Values for test_pop ----  THESE ARE FIXED VALUES 
# 0: no one tested
# 1: agyw all (not possible for symptomatic group)
# 2: agyw sexually active
# 3: pregnant women (any risk and age)
# 4: fsw (any age)
# 5: male (any risk and age)
# 6: male high risk (any age)

### <> 1. POCT: STI service ----

# Tests available per year for those with discharge (any cause) who access STI services 
test_ty <- interpolate(test_t, test_y, "constant")

# Eligible population: number to develop symptoms and seek STI care per year
# Symptoms due to NG and symptoms due to other infections
# Assume care seeking proportion is same for NG and other STIs

# Rate of developing symptoms due to other infections:
# Aetiologic proportion = rate_ng / (rate_ng + rate_othersti)
# Assume aetiologic proportion is constant at equilibrium
# Once testing is implemented, rate_ng decreases and aetiologic prop changes
# At equilibrium (prior to testing): rate_othersti = rate_ng * (1 - ap) / ap
# Once testing is implemented:       rate_othersti = constant

# Old version uses compartment for O (number of individuals presenting with symptoms due to non-NG causes)
# rate_ng[,,,] <- lambda[i,j,k,l] * phi[i] * gamma[i,j,k,l]
# 
# rate_other[,,,] <- 
#   if(i == 1 && l == 2) 
#     0 else 
#   if(test_ty == 0) 
#     lambda[i,j,k,l] * phi[i] * gamma[i,j,k,l] * (1 - ap[i]) / ap[i] else rate_other[i,j,k,l]

# Annual number exiting Y and exiting O 
# n_symp[,,,] <- Y[i,j,k,l] * tau_m[i,j,k,l] + O[i,j,k,l] * tau_m[i,j,k,l]

# New version uses aetiologic proportion without creating a new compartment

# Number of symptomatic cases presenting to care, caused by other RTIs and not NG 
# Note rate_other has different definition here to version above
rate_other[,,,] <- if(test_ty == 0) 
  Y[i,j,k,l] * tau_m[i,j,k,l] * (1 - ap[i]) / ap[i] else rate_other[i,j,k,l]

# Total number symptomatic cases presenting to care
n_symp[,,,] <-  Y[i,j,k,l] * tau_m[i,j,k,l] + rate_other[i,j,k,l]

# Proportion of eligible population tested 
# delta is probability of being tested when accessing STI services
# delta is equivalent for all population strata

# ANY CHANGES TO TEST_POP DEFINITION SHOULD BE CONSIDERED FOR SCENARIO 2
delta_init[,,,] <-
  # at t = 0, n_sym = 0 so to avoid dividing by 0
  if(n_symp[i,j,k,l] == 0) 0 else
    # agyw all
    if (test_pop == 1 && i == 2 && k == 1 && l == 1) 0 else
      # agyw sexually active
      if (test_pop == 2 && i == 2 && k == 1 && l == 1 && test_ty > 0) 1 else
        # pregnant women
        if (test_pop == 3 && i == 2 && l == 2 && test_ty > 0) 1 else
          # fsw
          if (test_pop == 4 && i == 2 && j == 3 && test_ty > 0) 1 else
            # male (any risk and age)
            if (test_pop == 5 && i == 1 && l == 1 && test_ty > 0) 1 else
              # male high risk (any age)
              if (test_pop == 6 && i == 1 && j == 3 && l == 1 && test_ty > 0) 1 else
                # no one tested
                0

# If tests > number eligible (delta > 1), set delta = 1
delta[,,,] <- if (delta_init[i,j,k,l] > 1) 1 else delta_init[i,j,k,l]

# Number eligible for test in scenario (i.e. denominator that actually receive test)
n_elig_y[,,,] <-
  if(Y[i,j,k,l] == 0) 0 else
    # agyw all
    if (test_pop == 1 && i == 2 && k == 1 && l == 1) 0 else
      # agyw sexually active
      if (test_pop == 2 && i == 2 && k == 1 && l == 1) Y[i,j,k,l] * tau_m[i,j,k,l] else
        # pregnant women
        if (test_pop == 3 && i == 2 && l == 2) Y[i,j,k,l] * tau_m[i,j,k,l] else
          # fsw
          if (test_pop == 4 && i == 2 && j == 3) Y[i,j,k,l] * tau_m[i,j,k,l] else
            # male (any risk and age)
            if (test_pop == 5 && i == 1 && l == 1) Y[i,j,k,l] * tau_m[i,j,k,l] else
              # male high risk (any age)
              if (test_pop == 6 && i == 1 && j == 3 && l == 1) Y[i,j,k,l] * tau_m[i,j,k,l] else
                # no one tested
                0

n_elig_o[,,,] <-
  # agyw all
  if (test_pop == 1 && i == 2 && k == 1 && l == 1) 0 else
    # agyw sexually active
    if (test_pop == 2 && i == 2 && k == 1 && l == 1) rate_other[i,j,k,l] else
      # pregnant women
      if (test_pop == 3 && i == 2 && l == 2) rate_other[i,j,k,l] else
        # fsw
        if (test_pop == 4 && i == 2 && j == 3) rate_other[i,j,k,l] else
          # male (any risk and age)
          if (test_pop == 5 && i == 1 && l == 1) rate_other[i,j,k,l] else
            # male high risk (any age)
            if (test_pop == 6 && i == 1 && j == 3 && l == 1) rate_other[i,j,k,l] else
              # no one tested
              0

# receive POCT
# n_test_y[,,,] <- delta[i,j,k,l] * tau_m[i,j,k,l] * Y[i,j,k,l]
# n_test_o[,,,] <- delta[i,j,k,l] * rate_other[i,j,k,l]
n_test_y[,,,] <- delta[i,j,k,l] * n_elig_y[i,j,k,l]
n_test_o[,,,] <- delta[i,j,k,l] * n_elig_o[i,j,k,l]

### <> 2. POCT: Other health service ----

# Tests available per year for individuals who access "other" health services
# Other = non-STI related services, which varies per population group
# e.g. ANC for pregnant women, PHC for men, TO BE FURTHER DEFINED
test_th  <- interpolate(test_t, test_h, "constant")

# Rate of testing among population accessing other health services
# tau is rate at which access health services and get tested
# tau = test_th /  (hc_rate * N)
# hc_rate is rate at which population accesses health services
# hc_rate * N is number accessing health services per year
# hc_rate = -log(1 - hc_probability)

# Eligible population
N_hc[,,,] <-
  if(test_pop == 2)  
    # use sexually active population for agyw sex active
    -log(1 - hc[i,j,k,l]) * (S[i,j,k,l] + X[i,j,k,l] + Y[i,j,k,l] + Z[i,j,k,l]) else
      # use full population for other test pops
      -log(1 - hc[i,j,k,l]) * (U[i,j,k,l] + S[i,j,k,l] + X[i,j,k,l] + Y[i,j,k,l] + Z[i,j,k,l])


# ANY CHANGES TO TEST_POP DEFINITION SHOULD BE CONSIDERED FOR SCENARIO 1
tau_hc[,,,] <-
  # no one tested
  if (test_pop == 0) 0 else
    # agyw all
    if (test_pop == 1 && i == 2 && k == 1 && l == 1 && test_th > 0) 1 else
      # agyw sexually active
      if (test_pop == 2 && i == 2 && k == 1 && l == 1 && test_th > 0) 1 else
        # pregnant women
        if (test_pop == 3 && i == 2 && l == 2 && test_th > 0) 1 else
          # fsw
          if (test_pop == 4 && i == 2 && j == 3 && test_th > 0) 1 else
            # male (any risk and age)
            if (test_pop == 5 && i == 1 && l == 1 && test_th > 0) 1 else
              # male high risk (any age)
              if (test_pop == 6 && i == 1 && j == 3 && l == 1 && test_th > 0) 1 else
                # no one tested
                0


# Number eligible for test in scenario (i.e. denominator that actually receive test)
n_elig_h[,,,] <- 
  # no one tested
  if (test_pop == 0) 0 else
    # agyw all
    if (test_pop == 1 && i == 2 && k == 1 && l == 1) N_hc[i,j,k,l] else
      # agyw sexually active
      if (test_pop == 2 && i == 2 && k == 1 && l == 1) N_hc[i,j,k,l] else
        # pregnant women
        if (test_pop == 3 && i == 2 && l == 2) N_hc[i,j,k,l] else
          # fsw
          if (test_pop == 4 && i == 2 && j == 3) N_hc[i,j,k,l] else
            # male (any risk and age)
            if (test_pop == 5 && i == 1 && l == 1) N_hc[i,j,k,l] else
              # male high risk (any age)
              if (test_pop == 6 && i == 1 && j == 3 && l == 1) N_hc[i,j,k,l] else
                # no one tested
                0

# receive POCT
n_test_h[,,,] <- tau_hc[i,j,k,l] * n_elig_h[i,j,k,l]

### <> 3. Over-treatment ----

# Number of false positive results among those not infected with NG

fp_sm[,,,] <- rate_other[i,j,k,l] * (1 - delta[i,j,k,l]) * (1 - spec_sm[i])

fp_poct[,,,] <- rate_other[i,j,k,l] * delta[i,j,k,l] * (1 - spec_poct[i]) +
  (U[i,j,k,l] + S[i,j,k,l]) * tau_hc[i,j,k,l] * (1 - spec_poct[i])


## COMPARTMENTS ----------------------------------------------------------------

## < Assign initial states -----
initial(U[,,,])   <- U0[i,j,k,l]
initial(S[,,,])   <- S0[i,j,k,l]
initial(X[,,,])   <- X0[i,j,k,l]
initial(Y[,,,])   <- Y0[i,j,k,l]
initial(Z[,,,])   <- Z0[i,j,k,l]
initial(My[,,,])  <- 0
initial(Ty[,,,])  <- 0
initial(Th_symp[,,,])  <- 0
initial(Th_asymp[,,,])  <- 0
initial(Ry[,,,])  <- 0
initial(Rh_symp[,,,])  <- 0
initial(Rh_asymp[,,,])  <- 0
# initial(O[,,,])   <- 0
# initial(Mo[,,,])  <- 0
# initial(Ro[,,,])  <- 0

# initial(cum_X[,,,])  <- 0
# initial(cum_Y[,,,])  <- 0
# initial(cum_Z[,,,])  <- 0
# initial(cum_My[,,,]) <- 0
# initial(cum_Ty[,,,]) <- 0
# initial(cum_Th_symp[,,,]) <- 0
# initial(cum_Th_asymp[,,,]) <- 0
# initial(cum_Ry[,,,]) <- 0
# initial(cum_Rh_symp[,,,]) <- 0
# initial(cum_Rh_asymp[,,,]) <- 0

# initial(cum_O[,,,])  <- 0
# initial(cum_Mo[,,,]) <- 0
# initial(cum_Ro[,,,]) <- 0

initial(py_X[,,,])  <- 0
initial(py_Y[,,,])  <- 0
initial(py_Z[,,,])  <- 0
initial(py_My[,,,]) <- 0
initial(py_Ty[,,,]) <- 0
initial(py_Th_symp[,,,]) <- 0
initial(py_Th_asymp[,,,]) <- 0
initial(py_Ry[,,,]) <- 0
initial(py_Rh_symp[,,,]) <- 0
initial(py_Rh_asymp[,,,]) <- 0

## < Define initial states ----
U0[,,,]  <- user()
S0[,,,]  <- user()
X0[,,,]  <- user()
Y0[,,,]  <- user()
Z0[,,,]  <- user()

## < Set dimensions ----
dim(U)    <- c(n_sex, n_risk, n_age, n_gest)
dim(S)    <- c(n_sex, n_risk, n_age, n_gest)
dim(X)    <- c(n_sex, n_risk, n_age, n_gest)
dim(Y)    <- c(n_sex, n_risk, n_age, n_gest)
dim(Z)    <- c(n_sex, n_risk, n_age, n_gest)
dim(My)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Ty)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Th_symp)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Th_asymp)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Ry)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Rh_symp)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Rh_asymp)   <- c(n_sex, n_risk, n_age, n_gest)
# dim(O)    <- c(n_sex, n_risk, n_age, n_gest)
# dim(Mo)   <- c(n_sex, n_risk, n_age, n_gest)
# dim(Ro)   <- c(n_sex, n_risk, n_age, n_gest)

dim(U0)   <- c(n_sex, n_risk, n_age, n_gest)
dim(S0)   <- c(n_sex, n_risk, n_age, n_gest)
dim(X0)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Y0)   <- c(n_sex, n_risk, n_age, n_gest)
dim(Z0)   <- c(n_sex, n_risk, n_age, n_gest)

dim(n_X)  <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Y)  <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Z)  <- c(n_sex, n_risk, n_age, n_gest)
dim(n_I)  <- c(n_sex, n_risk, n_age, n_gest)
dim(n_My)  <- c(n_sex, n_risk, n_age, n_gest)
dim(n_test_h) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_test_y) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_test_o) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_elig_h) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_elig_y) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_elig_o) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Ty) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Th_symp) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Th_asymp) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Ry) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Rh_symp) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_Rh_asymp) <- c(n_sex, n_risk, n_age, n_gest)
# dim(n_O)  <- c(n_sex, n_risk, n_age, n_gest)
# dim(n_Mo) <- c(n_sex, n_risk, n_age, n_gest)
# dim(n_Ro) <- c(n_sex, n_risk, n_age, n_gest)

# dim(cum_X)  <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Y)  <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Z)  <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_My) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Ty) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Th_symp) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Th_asymp) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Ry) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Rh_symp) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Rh_asymp) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_O)  <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Mo) <- c(n_sex, n_risk, n_age, n_gest)
# dim(cum_Ro) <- c(n_sex, n_risk, n_age, n_gest)

dim(py_X)  <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Y)  <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Z)  <- c(n_sex, n_risk, n_age, n_gest)
dim(py_My) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Ty) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Th_symp) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Th_asymp) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Ry) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Rh_symp) <- c(n_sex, n_risk, n_age, n_gest)
dim(py_Rh_asymp) <- c(n_sex, n_risk, n_age, n_gest)

## PARAMETERS -------------------------------------------------------

## < Define input parameters ----

# demographic
ent_t[]         <- user()
entm_y[]        <- user()
entf_y[]        <- user()
f_t[]           <- user()
fy_y[]          <- user()
fo_y[]          <- user()
q[,,,]          <- user() 
upsilon[,]      <- user()
omega[]         <- user()
alpha           <- user()
mu[,]           <- user()

# infection
c[,,,]          <- user() 
n_l             <- user()
n_m             <- user()
n_h             <- user()
kappa[]         <- user()
theta_c         <- user()
theta_chi       <- user()
chi[,,,]        <- user()
e               <- user()
phi[]           <- user()
gamma[,,,]      <- user()
pi[]            <- user()
sigma[]         <- user()

# testing and treatment
test_t[]        <- user()
test_h[]        <- user()
test_y[]        <- user()
test_pop        <- user()
ap[]            <- user()
zeta[]          <- user()
eta[]           <- user()
# sens_o[]        <- user()    
spec_sm[]       <- user()
spec_poct[]     <- user()
tau_m[,,,]      <- user()
hc[,,,]         <- user()

## < Set dimensions ----

# compartments
dim(NS)         <- c(n_sex, n_risk, n_age, n_gest)
dim(N)          <- c(n_sex, n_risk, n_age, n_gest)
dim(I)          <- c(n_sex, n_risk, n_age, n_gest)

# demographic
dim(ent_t)      <- user()
dim(entm_y)     <- user()
dim(entf_y)     <- user()
dim(f_t)        <- user()
dim(fy_y)       <- user()
dim(fo_y)       <- user()
dim(ent)        <- c(n_sex, n_risk, n_age, n_gest)
dim(f)          <- c(n_sex, n_risk, n_age, n_gest)
dim(q)          <- c(n_sex, n_risk, n_age, n_gest) 
dim(upsilon)    <- c(n_sex, n_risk)
dim(psi)        <- c(n_sex, n_risk, n_age, n_gest)
dim(omega)      <- c(n_age)
dim(mu)         <- c(n_sex, n_age)

# sexual partnership probability
dim(rho_set)    <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(rho_sp)     <- c(n_sex, n_risk, n_age, n_gest)
dim(rho_sub)    <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(rho)        <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(B)          <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(bal_exp)    <- c(n_sex)
dim(c_ad)       <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(c_ad_strat) <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(sp)         <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(n)          <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(kappa_p)    <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(chi_avg)    <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(lambda_sub) <- c(n_sex, n_risk, n_age, n_gest, n_sex, n_risk, n_age, n_gest)
dim(lambda)     <- c(n_sex, n_risk, n_age, n_gest)

# infection
dim(c)          <- c(n_sex, n_risk, n_age, n_gest) 
dim(kappa)      <- c(n_sex) 
dim(chi)        <- c(n_sex, n_risk, n_age, n_gest)
dim(phi)        <- c(n_sex)
dim(gamma)      <- c(n_sex, n_risk, n_age, n_gest)
dim(pi)         <- c(n_sex)
dim(sigma)      <- c(n_sex)

# testing and treatment
dim(test_t)     <- user()
dim(test_h)     <- user()
dim(test_y)     <- user()
dim(ap)         <- c(n_sex)
# dim(rate_ng)    <- c(n_sex, n_risk, n_age, n_gest)
dim(rate_other) <- c(n_sex, n_risk, n_age, n_gest)
dim(n_symp)     <- c(n_sex, n_risk, n_age, n_gest)
dim(delta_init) <- c(n_sex, n_risk, n_age, n_gest)
dim(delta)      <- c(n_sex, n_risk, n_age, n_gest)
dim(zeta)       <- c(n_sex)
dim(eta)        <- c(n_sex)
# dim(sens_o)     <- c(n_sex)
dim(spec_sm)    <- c(n_sex)
dim(spec_poct)  <- c(n_sex)
dim(hc)         <- c(n_sex, n_risk, n_age, n_gest)
dim(N_hc)       <- c(n_sex, n_risk, n_age, n_gest)
dim(tau_hc)     <- c(n_sex, n_risk, n_age, n_gest)
dim(tau_m)      <- c(n_sex, n_risk, n_age, n_gest)
dim(fp_sm)      <- c(n_sex, n_risk, n_age, n_gest)
dim(fp_poct)    <- c(n_sex, n_risk, n_age, n_gest)

## OUTPUT ----------------------------------------------------------------------

output(N)       <- N
output(NS)      <- NS
output(I)       <- I

#output(psi) <- psi
#output(ent) <- ent
#output(upsilon) <- upsilon
#output(q) <- q
#output(f) <- f

#output(rho)     <- rho
#output(c)       <- c
#output(c_ad)    <- c_ad
#output(c_ad_strat) <- c_ad_strat
#output(sp)      <- sp

#output(chi) <- chi
#output(chi_avg) <- chi_avg
#output(n) <- n
#output(kappa_p) <- kappa_p
#output(lambda) <- lambda
#output(lambda_sub) <- lambda_sub

# output(rate_ng) <- rate_ng
#output(rate_other) <- rate_other
output(n_symp) <- n_symp
#output(delta_init) <- delta_init
output(delta) <- delta
#output(gamma) <- gamma
#output(phi) <- phi
output(tau_m) <- tau_m
output(N_hc)  <- N_hc
output(tau_hc) <- tau_hc
#output(fp_sm) <- fp_sm
#output(fp_poct) <- fp_poct

output(n_X) <- n_X
output(n_Y) <- n_Y
output(n_Z) <- n_Z
output(n_I) <- n_I
output(n_My) <- n_My
output(n_Ty) <- n_Ty
output(n_Th_symp) <- n_Th_symp
output(n_Th_asymp) <- n_Th_asymp
output(n_test_h) <- n_test_h
output(n_test_y) <- n_test_y
output(n_test_o) <- n_test_o
output(n_elig_h) <- n_elig_h
output(n_elig_y) <- n_elig_y
output(n_elig_o) <- n_elig_o
output(n_Ry) <- n_Ry
output(n_Rh_symp) <- n_Rh_symp
output(n_Rh_asymp) <- n_Rh_asymp
# output(n_O) <- n_O
# output(n_Mo) <- n_Mo
# output(n_Ro) <- n_Ro


sir <- odin({
  
  # Equations for transitions between compartments by age group
  
  deriv(Sm) <- (mu_m*Nm) - (b*p_bm*(Ib_total/Nb)+mu_m)*Sm #assuming birth and death rate are the same for simplicity
  
  deriv(Em) <- (b*p_bm*(Ib_total/Nb)*Sm)-((theta_m+mu_m)*Em)
  
  deriv(Im) <- (theta_m*Em)-(mu_m*Im)
  
  deriv(Sb[]) <- gamma*Nb - (mu_b+lambda_b[i])*Sb[i]
  
  deriv(Eb[]) <- lambda_b[i]*Sb[i] - (mu_b+theta_b)*Eb[i]
  
  deriv(Ib[]) <- theta_b*Eb[i]-(mu_b+sigma)*Ib[i]
  
  deriv(Rb[]) <- (sigma*Ib[i]) - (mu_b*Rb[i])
  
  Ib_total <- sum(Ib[])
  
  Nb <- sum(Sb[]) + sum(Eb[]) + sum(Ib[]) + sum(Rb[])
  
  # Calculate force of infection (lambda)
  
  dim(beta) <- c(2,2)
  beta[1,1] <- 0.3
  beta[1,2] <- 0.6
  beta[2,1] <- 0.6
  beta[2,2] <- 0.3
  lambda_b[] <- beta[i]*Im/Nm

  # AI generated
  #P1 <- Im / Nm
  #P2 <- 0.5 * Im / Nm   # or any real second mechanism
  #lambda_b[] <- beta[i,1] * P1 + beta[i,2] * P2
  #beta <- parameter()
  #lambda_b[] <- sum(beta[i,])*Im/Nm
  
  
  #initial conditions
  
  initial(Sb[1]) <- 500
  initial(Sb[2]) <- 500
  initial(Eb[]) <- 0
  initial(Ib[]) <- 0
  initial(Rb[]) <- 0
  initial(Sm) <- Nm - Im0
  initial(Em) <- 0
  initial(Im) <- Im0
  
  # User defined parameters - default in parentheses:
  Im0 <- parameter(5)
  Nm <- parameter(1000) # unknown parameter
  # Nb <- parameter(1000) # unknown parameter
  mu_m <- parameter(0.02) # random value, should be related to temperature
  b <- parameter(0.3) # random value
  p_bm <- parameter(0.9) # random value, should be related to temperature
  theta_m <- parameter(0.1) # random value, should be related to temperature
  gamma <- parameter(0.0015) # changed to same as death rate  (avian fertility rate)
  mu_b <- parameter(0.0015) # death rate for mature birds (table)
  #p_mb <- parameter(0.9) # table
  theta_b <- parameter(0.15) # table
  sigma <- parameter(0.1) # avian recovery rate, table 0.57
  
  # Dimensions of arrays
  dim(Sb) <- 2
  dim(Eb) <- 2
  dim(Ib) <- 2
  dim(Rb) <- 2
  dim(lambda_b) <- 2
})

sys <- dust_system_create(sir, list())
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

#plot

matplot(t, t(y$Ib),
        type = "l", lty = 1,
        xlab = "Time",
        ylab = "Infected birds")

legend("topright",
       legend = c("Group 1", "Group 2"),
       col = 1:2, lty = 1)

# distance between populations (km)

daily_temperatures <- bologna.temperature.series

daily_temp <- daily_temperatures$tavg
str(daily_temp)

BriereNorm <- function(temp, T_minB, T_maxB, mB) {
  
  briere <- (temp - T_minB) * (T_maxB - temp)^mB
  
  briere[briere < 0] <- 0
  briere[is.na(briere)] <- 0
  
  briere <- briere / max(briere, na.rm = TRUE)
  
  return(briere)
}


briere_value <- rep(BriereNorm(
  temp = daily_temp, 
  T_minB = 14, 
  T_maxB = 34.3, 
  mB = 0.7), #24 peak
  10
)

briere_time <- seq(0, length(briere_value)-1, by=1)

temp_seq <- seq(0.3, 31.15, by = 0.1)


plot(temp_seq,
     BriereNorm(temp_seq, 14, 34.3, 0.7), 
     type = "l",
     xlab = "Temperature (°C)",
     ylab = "Briere response")


# Scenario 13 (metres)
distance <- matrix(c(
  1000, 200, 20, 100, 120,
  200, 1000, 120, 100, 100,
  20, 120, 1000, 40, 20,
  100, 100, 40, 1000, 160,
  120, 100, 20, 160, 1000
), nrow = 5, byrow = TRUE)


sir <- odin({
  
  # Equations for transitions between compartments by group
  
  deriv(S[]) <- (gamma[i]*N[i]) - (lambda[i]*S[i]) - (mu*S[i]) - seed[i]
  
  deriv(E[]) <- (lambda[i]*S[i]) - (theta*E[i]) - (mu*E[i])
  
  deriv(I[]) <- (theta*E[i]) - (mu*I[i]) - (sigma*I[i]) + seed[i]
  
  deriv (R[]) <- (sigma*I[i]) - (mu*R[i])
  
  # beta
  
  beta0 <- parameter(0.5)
  betaT <- briere * beta0
  h0 <- parameter(0.05)
  r0 <- parameter(100)
  alpha <- parameter(2)
  beta[,] <- (h0 / (1 + (distance[i,j] / r0)^alpha))
  new_beta[,] <- if (i == j) betaT else beta[i, j]

  
  lambda[] <- (
    new_beta[i,1] * I[1]/N[i] +
    new_beta[i,2] * I[2]/N[i] +
    new_beta[i,3] * I[3]/N[i] +
    new_beta[i,4] * I[4]/N[i] +
    new_beta[i,5] * I[5]/N[i]
  )
  
  
  
  # initial conditions
  
  N[] <- S[i] + E[i] + I[i] + R[i]
  
  initial(S[]) <- 300 
  
  initial(E[]) <- 0
  
  initial(I[]) <- 0 
  
  initial(R[]) <- 0
  
  # seed in march every year
  year <- floor(time/365) + 1
  seed[] <- if(
    i == seed_pop[year] &&
    time%%365 >= 60 && 
    time%%365 < 61 
  ) I0 else 0
  
  # User defined parameters:
  
  briere <- interpolate(briere_time, briere_value, "linear")
  
  briere_time <- parameter(constant = TRUE)
  briere_value <- parameter(constant = TRUE)
  
  seed_pop <- parameter(constant = TRUE)
  
  I0 <- parameter(1) # initial infected
  
  mu <- parameter(0.000913242) # death rate
  
  theta <- parameter(0.5) # incubation period
  
  sigma <- parameter(0.25) # recovery rate
  
  distance <- parameter(constant = TRUE)

  
  
  # birth rate
  gamma[] <- if (time %% 365 >= 90 && time %% 365 <= 212) 0.0027548209 else 0 
  # 2 per bird may-july
  
  # Dimensions of arrays
  
  dim(S) <- 5
  
  dim(E) <- 5
  
  dim(I) <- 5
  
  dim(R) <- 5
  
  dim(N) <- 5
  
  
  dim(new_beta) <- c(5, 5)
  
  dim(beta) <- c(5, 5)
  
  dim(distance) <- c(5,5)
  
  dim(briere_time, briere_value) <- parameter(rank = 1)
  
  dim(lambda) <- 5
  
  dim(seed_pop) <- 10
  
  dim(seed) <- 5
  
  dim(gamma) <- 5
})

##############
# Simulation #
##############

sys <- dust_system_create(
  sir,
  list(
    beta = beta,
    briere_time = briere_time,
    briere_value = briere_value,
    seed_pop=seed_pop,
    distance=distance
    
  )
)

dust_system_set_state_initial(sys)

t <- seq(0, 3649) # 10 years

y <- dust_system_simulate(sys, t)

y <- dust_unpack_state(sys, y)

##################
# infected birds #
##################

# data frames

df <- as.data.frame(t(y$I))
df$t <- t
colnames(df)[1:5] <- c("P1", "P2", "P3", "P4", "P5")
df_long <- pivot_longer(df, cols = c("P1", "P2", "P3", "P4", "P5"),
                        names_to = "Population",
                        values_to = "I")


prop <- y$I / (y$S + y$E + y$I + y$R) # proportion infected

df <- as.data.frame(t(prop))
df$t <- t
colnames(df)[1:5] <- c("P1", "P2", "P3", "P4", "P5")

df_propI_long <- pivot_longer(df,
                              cols = c(P1, P2, P3, P4, P5),
                              names_to = "Population",
                              values_to = "PropI")

df_plotI <- df_propI_long[df_propI_long$t >= (365 * 4), ]  # plot last 6 years

## 6 years ##


month_breaks <- c(0,31,59,90,120,151,181,212,243,273,304,334)

plotI6 <- ggplot(data = subset(df_plotI), aes(x = t %% 365, y = PropI, color = Population)) +
  geom_line() +
  facet_wrap(~(floor(t/365)+1), ncol = 1) +
  scale_x_continuous(
    breaks = month_breaks,
    labels = c("Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion Infected"
  ) +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

plotI6

## 10 years ##

plotI10 <- ggplot(data = subset(df_propI_long), aes(x = t %% 365, y = PropI, color = Population)) +
  geom_line() +
  facet_wrap(~(floor(t/365)+1), ncol = 1) +
  scale_x_continuous(
    breaks = month_breaks,
    labels = c("Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec")
  ) +
  labs(
    x = NULL,
    y = "Proportion Infected"
  ) +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )



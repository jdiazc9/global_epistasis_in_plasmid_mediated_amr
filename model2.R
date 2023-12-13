rm(list = ls())
library(deSolve)
library(ggplot2)
library(tidyr)
library(scales)

# model formulation: dynamics of intracellular antibiotic
dAdt <- function(t, state, params) with(as.list(c(state, params)),
                                        {
                                          dA_in <- v0*(A_out - A_in) - v1*E*A_in/(K1 + A_in) - v2*P*A_in/(K2 + A_in)       # internal antibiotic dynamics
                                          return(list(c(dA_in)))
                                        })

params <- c(
  v0 = 1,         # passive rate of internalization
  A_out = 500,    # extracellular antibiotic concentration
  v1 = 1,         # active degradation rate
  K1 = 50,        # internal concentration for half-max rate of degradation
  E = 0,          # enzyme concentration
  v2 = 1,         # active pumping rate
  K2 = 50,        # internal concentration for half-max rate of pumping
  P = 0           # pump concentration
)

state <- c(
  A_in = 0        # initial intracellular concentration
)

# how does fitness depend on intracellular antibiotic?
f_vs_Ain <- function(A_in, Fmax = 1, Kf = 50) Fmax * Kf/(Kf + A_in)

# background strains (enzyme expression normally distributed with avg. 250 and sd 50)
df <- data.frame(strain_id = 1:50,
                 E_background = rnorm(50, 250, 50),
                 P_background = 0)
df$E_plasmid <- df$E_background # + rnorm(nrow(df), 100, 10) # effect of plasmid acquisition: enzyme expression remains the same
df$P_plasmid <- df$P_background + rnorm(nrow(df), 100, 5) # effect of plasmid acquisition: pump expression increases by 100 +- 20 units

df$fitness_background <- NA
df$fitness_plasmid <- NA

# integrate equations
times <- seq(0, 24, by = 0.1)
for (i in 1:nrow(df)) {
  
  params['E'] <- df$E_background[i]
  params['P'] <- df$P_background[i]
  out <- as.data.frame(ode(y = state, times = times, func = dAdt, parms = params))
  out <- gather(out, variable, value, 2:ncol(out), factor_key=TRUE)
  
  df$fitness_background[i] <- f_vs_Ain(out$value[out$time == times[length(times)]])
  
  params['E'] <- df$E_plasmid[i]
  params['P'] <- df$P_plasmid[i]
  out <- as.data.frame(ode(y = state, times = times, func = dAdt, parms = params))
  out <- gather(out, variable, value, 2:ncol(out), factor_key=TRUE)
  
  df$fitness_plasmid[i] <- f_vs_Ain(out$value[out$time == times[length(times)]])
  
}

df$fitness_effect <- df$fitness_plasmid - df$fitness_background

ggplot(df,
       aes(x = fitness_background, y = fitness_effect)) +
  geom_point(color = 'black',
             shape = 16,
             cex = 2.5,
             alpha = 0.75) +
  geom_smooth(method = 'lm',
              formula = y ~ x,
              se = F,
              fullrange = T,
              color = 'gray') +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste('Background fitness, ', italic(F)[B], ' (a.u.)', sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = expression(paste('Fitness effect, ', Delta*italic(F), ' (a.u.)', sep = ''))) +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = paste('./plots/model2_dF_vs_FB.pdf', sep = ''),
       width = 85,
       height = 85,
       units = 'mm',
       limitsize = F)











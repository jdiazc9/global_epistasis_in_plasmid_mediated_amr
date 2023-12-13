rm(list = ls())
library(deSolve)
library(ggplot2)
library(tidyr)
library(scales)

# model formulation: dynamics of intracellular antibiotic
dAdt <- function(t, state, params) with(as.list(c(state, params)),
                                        {
                                          dA_in <- v0*(A_out - A_in) - v1*E*A_in/(K1 + A_in)       # internal antibiotic dynamics
                                          return(list(c(dA_in)))
                                        })

params <- c(
  v0 = 1,         # passive rate of internalization
  A_out = 100,    # extracellular antibiotic concentration
  v1 = 1,         # active degradation/externalization/inactivation rate
  K1 = 50,        # internal concentration for half-max rate of degradation/externalization/inactivation
  E = 0           # enzyme concentration
)

state <- c(
  A_in = 0        # initial intracellular concentration
)

# integrate equations
times <- seq(0, 24, by = 0.1)

plot_this <- data.frame(E = numeric(0),
                        time = numeric(0),
                        variable = character(0),
                        value = numeric(0))

for (Ei in seq(0, 500, length.out = 50)) {
  
  params['E'] <- Ei
  out <- as.data.frame(ode(y = state, times = times, func = dAdt, parms = params))
  out <- gather(out, variable, value, 2:ncol(out), factor_key=TRUE)
  
  plot_this <- rbind(plot_this,
                     cbind(E = Ei, out))
  
}



ggplot(plot_this, aes(x = time, y = value, group = E, color = E)) +
  geom_line()

ggplot(plot_this[plot_this$time == times[length(times)], ],
       aes(x = E, y = value)) +
  geom_point(color = '#b72d27',
             shape = 16) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste('Enzyme concentration, [', italic(E), '] (a.u.)', sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = expression('Intracellular'~'antibiotic'~'concentration, ['*italic(A)[`in`]*'] (a.u.)')) +
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

ggsave(filename = paste('./plots/Ain_vs_E.pdf', sep = ''),
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# how does fitness depend on intracellular antibiotic?
f_vs_Ain <- function(A_in, Fmax = 1, Kf = 50) Fmax * Kf/(Kf + A_in)

plot_this$fitness <- f_vs_Ain(plot_this$value)
ggplot(plot_this[plot_this$time == times[length(times)], ],
       aes(x = value, y = fitness)) +
  geom_point(color = '#9e509b',
             shape = 16) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = expression('Intracellular'~'antibiotic'~'concentration, ['*italic(A)[`in`]*'] (a.u.)')) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = 'Fitness (a.u.)') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = paste('./plots/F_vs_Ain.pdf', sep = ''),
       width = 50,
       height = 50,
       units = 'mm',
       limitsize = F)

ggplot(plot_this[plot_this$time == times[length(times)], ],
       aes(x = E, y = fitness)) +
  geom_point(color = 'black',
             shape = 16) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste('Enzyme concentration, [', italic(E), '] (a.u.)', sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = 'Fitness (a.u.)') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = paste('./plots/F_vs_E.pdf', sep = ''),
       width = 60,
       height = 60,
       units = 'mm',
       limitsize = F)

# background strains (enzyme expression normally distributed with avg. 250 and sd 50)
df <- data.frame(strain_id = 1:50,
                 E_background = rnorm(50, 250, 50))
df$E_plasmid <- df$E_background + rnorm(nrow(df), 100, 10) # effect of plasmid acquisition: enzyme expression increases by 100 units

df$fitness_background <- NA
df$fitness_plasmid <- NA

for (i in 1:nrow(df)) {
  
  params['E'] <- df$E_background[i]
  out <- as.data.frame(ode(y = state, times = times, func = dAdt, parms = params))
  out <- gather(out, variable, value, 2:ncol(out), factor_key=TRUE)
  
  df$fitness_background[i] <- f_vs_Ain(out$value[out$time == times[length(times)]])
  
  params['E'] <- df$E_plasmid[i]
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

ggsave(filename = paste('./plots/dF_vs_FB.pdf', sep = ''),
       width = 85,
       height = 85,
       units = 'mm',
       limitsize = F)











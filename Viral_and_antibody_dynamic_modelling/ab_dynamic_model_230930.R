# script to run viral dynamic models
library(ggplot2)
library(nlmixr2)
library(gridExtra)
#
# Read data
ab.dat <- read.csv("panoramic_virology_molnupiravir_230927.csv", 
                   head = TRUE, stringsAsFactors = FALSE)
#### Add covariates for the final model: ####
# 1.  Vaccine effect for not fully vaccinated
ab.dat$vac_dich <- 0
ab.dat$vac_dich[ab.dat$n_vacccine <= 2] <- 1
# 2. Dichotomous molnupiravir effect, drug_eff gives as-treated analysis
ab.dat$mol_all_dich <- ab.dat$drug_eff
#### run 7 sex vacc on baseline, only mol on slope - final model #####
si.mod.full.ab <- function() {
  ini({
    tdelta.ab <-  -2.5  # Increase in s antibody
    tab0    <- 3.3      # Viral load at treatment initiation
    eta.delta.ab + eta.ab.ab0  ~ c(0.1,
                                   0.001, 0.1)
    add.err.ab <- 0.1     # residual variability
    beta_sex <- - 1
    beta_vac <- -0.3
    beta_mol <- - 0.1
  })
  model({
    delta.ab <- exp(tdelta.ab + eta.delta.ab +
                      beta_mol * mol_all_dich)  # individual value of delta
    a0    <- 10^(tab0 + eta.ab.ab0 +
                   beta_sex * sex_male  +
                   beta_vac * vac_dich)        # individual value of v0
    A_a(0) = a0
    
    d/dt(A_a) =  delta.ab * A_a 
    
    s_antibody = log10(A_a) 
    s_antibody ~ add(add.err.ab)       # define error model
    
    treat_group = treat_group
    sample_group = sample_group
  })
}
#si.fit.7.ab<- nlmixr2(si.mod.full.ab, ab.dat[ab.dat$dvid == "s_antibody" &
#                                                    ab.dat$cens == 0, ],
#                               table = tableControl(cwres = TRUE), #, npde = TRUE),  
#                               est = "foce")
#bootstrapFit(si.fit.7.ab)
#saveRDS(si.fit.7.ab, "si.fit.7.ab.rds")
# final model saved as rds
si.fit.7.ab<- readRDS("si.fit.7.ab.rds")
###### make parameter estimate table #######
pe <- si.fit.7.ab$parFixedDf
# take exp of all covariates to get fractional change
pe$`Back-transformed`[4:6] <- exp(pe$`Back-transformed`[4:6])
pe$`CI Lower`[4:6] <- exp(pe$`CI Lower`[4:6])
pe$`CI Upper`[4:6] <- exp(pe$`CI Upper`[4:6])
pe$`Back-transformed` <- round(pe$`Back-transformed`, 3)
pe$`CI Lower` <- round(pe$`CI Lower`, 3)
pe$`CI Upper` <- round(pe$`CI Upper`, 3)
pe$est <- paste0(pe$`Back-transformed`, " (", pe$`CI Lower`, ", ", pe$`CI Upper`, ")")
pe$est[3] <- pe$`Back-transformed`[3]
pe$`BSV(CV% or SD)`[1] <- signif(pe$`BSV(CV% or SD)`[1], 3)
pe$`BSV(CV% or SD)`[2] <- signif(sqrt(pe$`BSV(CV% or SD)`[2]) * 100, 3)
pe$`BSV(CV% or SD)`[is.na(pe$`BSV(CV% or SD)`)] <- "-"

ab.parest <- data.frame(matrix(ncol = 3, nrow = 6))
colnames(ab.parest) <- c("Parameter", "Estimate (95%CI)", "IIV (%CV)")
ab.parest$Parameter <- c("Delta (/d)", "A0 (log10 U/mL)", "Additive error",
                         "Beta_sex", "Beta_vac", "Beta_mol")
ab.parest$`Estimate (95%CI)` <- pe$est
ab.parest$`IIV (%CV)` <- pe$`BSV(CV% or SD)`
write.csv(ab.parest, "antibody_dynamic_model_parameter_table.csv", row.names = FALSE)

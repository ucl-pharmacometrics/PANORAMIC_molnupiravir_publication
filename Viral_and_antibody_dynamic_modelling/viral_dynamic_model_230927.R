# script to run viral dynamic models
library(ggplot2)
library(nlmixr2)
#
# Read data
vl.dat <- read.csv("panoramic_virology_molnupiravir_230927.csv", 
                   head = TRUE, stringsAsFactors = FALSE)
#### Add covariates for the final model: ####
# 1. Time varying treatment effect need on and off treatment as decline rate slows after treatment stops
vl.dat$mol_dich_on <- 0
vl.dat$mol_dich_off <- 0
vl.dat$mol_dich_on[vl.dat$treat_group == "Molnupiravir" & vl.dat$t_start_therapy < 5.9] <- 1
vl.dat$mol_dich_off[vl.dat$treat_group == "Molnupiravir" & vl.dat$t_start_therapy >= 5.9] <- 1
# multiply by drug_eff to get as-treated analysis
vl.dat$mol_dich_on <- vl.dat$mol_dich_on * vl.dat$drug_eff
# 2. Allometric effect for age
vl.dat$log_age_cov <- log(vl.dat$age_y / median(vl.dat$age_y))
# 3. Allometric effect for antibody
vl.dat$log_ab_cov <- log(vl.dat$bl_antibody / median(vl.dat$bl_antibody))
# 4. Allometric effect for time since symptom onset
vl.dat$t_symp_day <- vl.dat$t_symp_enrol + 1 
vl.dat$log_symp_cov <- log(vl.dat$t_symp_day  / mean(vl.dat$t_symp_day))
#### run 11 full on baseline, only mol on slope - final model #####
si.mod.mol.dich.full.ud  <- function() {
  ini({
    tdelta <-  -0.05    # Death rate of infected cells
    tv0    <- 7.16       # Viral load at treatment initiation
    eta.delta + eta.v0 ~ c(0.05,
                           -0.001, 0.5)
    add.err <- 0.9     # residual variability
    
    beta_mol_on <- 0.67
    beta_mol_off <- -1.3
    beta_sex <- 0.2
    beta_age <- 0.75
    beta_ab <- -1.0
    beta_sym <- -0.75   
    
  })
  model({
    delta <- exp(tdelta + eta.delta + 
                   beta_mol_on * mol_dich_on +
                   beta_mol_off * mol_dich_off)  # individual value of delta
    v0    <- 10^(tv0 + eta.v0 +
                   beta_sex * sex_male +
                   beta_age * log_age_cov +
                   beta_ab * log_ab_cov +
                   beta_sym * log_symp_cov)    # individual value of v0
    A_v(0) = v0
    
    d/dt(A_v) =  - delta * A_v 
    
    virus = log10(A_v) 
    virus ~ add(add.err)       # define error model
    
    treat_group = treat_group
    sample_group = sample_group
  })
}
#si.fit.11 <- nlmixr2(si.mod.mol.dich.full.ud, vl.dat[vl.dat$dvid == "virus", ],
#                               table = tableControl(cwres = TRUE, npde = TRUE),  
#                               est = "foce")
#bootstrapFit(si.fit.11)
#saveRDS(si.fit.11, "si.fit.11.rds")
# final model saved as rds
si.fit.11 <- readRDS("si.fit.11.rds")
###### make parameter estimate table #######
pe <- si.fit.11$parFixedDf
# take exp of beta mol on, off and sex to get fractional change
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

vl.parest <- data.frame(matrix(ncol = 3, nrow = 9))
colnames(vl.parest) <- c("Parameter", "Estimate (95%CI)", "IIV (%CV)")
vl.parest$Parameter <- c("Delta (/d)", "V0 (log10 cp/mL)", "Additive error",
                         "Beta_mol_on", "Beta_mol_off", "Beta_sex",
                         "Beta_age", "Beta_ab", "Beta_sym")
vl.parest$`Estimate (95%CI)` <- pe$est
vl.parest$`IIV (%CV)` <- pe$`BSV(CV% or SD)`
write.csv(vl.parest, "viral_dynamic_model_parameter_table.csv", row.names = FALSE)

# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw
draw_R_mass <- function(n,R_basic, intervention_efficacy,
                        vac_eff, p_detected, p_vacc
                        ) {
  #Four types of cases
  
  # 1: Naive case, not detected
  R_undetected_unvacc <-  R_basic
  
  # 2: Naive case, but detected
  R_detected_unvacc <- R_basic * (1 - intervention_efficacy)
  
  # 3: Vaccinated case, not detected
  R_undetected_vacc <- R_basic * (1 - vac_eff)
  
  # 4: Vaccinated, detected case
  R_detected_vacc <- R_basic * (1 - intervention_efficacy) * (1 - vac_eff)
  
  
  group_R = c(R_detected_vacc, R_detected_unvacc, 
              R_undetected_vacc,R_undetected_unvacc)
  
  group_p = c(p_detected*p_vacc, p_detected*(1-p_vacc),
              (1 - p_detected)*p_vacc, (1 - p_detected)*(1-p_vacc))
  sample(
    group_R,
    size = n,
    prob = group_p,
    replace = TRUE)
}

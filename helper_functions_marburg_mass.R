# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw
draw_R_mass <- function(n,R_detected_vacc,R_detected_unvacc,
                        R_undetected_vacc,R_undetected_unvacc,
                        p_detected,p_vacc
) {
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

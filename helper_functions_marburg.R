# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw
draw_R <- function(n,R_detected, R_undetected,p_detected
) {
  group_R = c(R_detected, R_undetected)
  group_p = c(p_detected, 1 - p_detected)
  sample(
    group_R,
    size = n,
    prob = group_p,
    replace = TRUE)
}

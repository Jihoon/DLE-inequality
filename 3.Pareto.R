k = seq(0.01, 2, 0.01)

base_gini = given_gini * avg_ene / (avg_ene - (1-k)*actual.min)

alpha = 0.5*(1/base_gini + 1)

base_mean = (avg_ene - (1-k)*actual.min) / k
base_mean2 = alpha * actual.min / (alpha-1)

sol = which((base_mean > base_mean2) ==FALSE)[1]

k_sol = k[sol]
alpha_sol = alpha[sol]
d_sol = (1-k_sol) * actual.min

x0 = dpareto(xval0, actual.min, shape=alpha_sol)   # Gini = 1/(2*alpha-1)= 0.54 Staus quo
# r0 = rpareto(1e7, actual.min, shape=alpha_sol)   # Gini = 1/(2*alpha-1)
# mean(r0)
# gini(r0)
# min(r0)
# z0 = k_sol*r0 + d_sol
# mean(z0)
# gini(z0)
# min(z0)
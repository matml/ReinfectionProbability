start_me <- function() {

	library(rPython)
	library(dplyr)
	python.load(file.path(dir_python, "Algorithm2_new_M_improved.py"))

	dir_project <<- path.expand("~/work/projects/ReinfectionProbability")
	dir_python <<- dir_project %>% file.path("python")
	dir_R <<- dir_project %>% file.path("R")

}

p_reinfection_SIR <- function(gamma = 0, gamma_A = gamma, beta = 0, beta_A_C = beta, beta_C_A = beta, alpha_I = 0, alpha_I_A = alpha_I, alpha_R = 0, alpha_R_A = alpha_R, sigma = 0, sigma_A = sigma, n_pop, n_reinfection, S_0, I_0, R_0) {

	python.call("p_reinfection_SIR", gam = gamma, gamA = gamma_A, bet = beta, betAc = beta_A_C, betcA = beta_C_A , alpI = alpha_I, alpIA = alpha_I_A, alpR = alpha_R, alpRA = alpha_R_A, sig = sigma, sigA = sigma_A, N = n_pop, M = n_reinfection, i0 = I_0, s0 = S_0, r0 = R_0)

}

p_reinfection_SIRS <- function(n_reinfection, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0) {

	gamma <- 1/D_infection
	beta <- R0/D_infection/n_pop
	alpha_R <- 1/D_immunity

	p_reinfection_SIR(n_reinfection = n_reinfection, gamma = gamma, beta = beta, alpha_R = alpha_R, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0)

}

p_reinfection_AoN <- function(n_reinfection, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0) {

	gamma <- prop_immunity/D_infection
	beta <- R0/D_infection/n_pop
	alpha_I <- (1-prop_immunity)/D_infection

	p_reinfection_SIR(n_reinfection = n_reinfection, gamma = gamma, beta = beta, alpha_I = alpha_I, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0)

}

p_reinfection_PPI <- function(n_reinfection, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0) {

	gamma <- 1/D_infection
	beta <- R0/D_infection/n_pop
	sigma <- (1 - partial_protection)

	p_reinfection_SIR(n_reinfection = n_reinfection, gamma = gamma, beta = beta, sigma = sigma, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0)

}

logLike_0123 <- function(data = c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0), p_reinfection, ..., echo = FALSE) {

	# prob n_reinfection >=1, >= 2, >=3
	p_1_or_more <- p_reinfection(n_reinfection = 1, ...)
	p_2_or_more <- p_reinfection(n_reinfection = 2, ...)
	p_3_or_more <- p_reinfection(n_reinfection = 3, ...)
	
	# prob exactly 0, 1, 2 reinfections
	p_0 <- 1 - p_1_or_more
	p_1 <- p_1_or_more - p_2_or_more
	p_2 <- p_2_or_more - p_3_or_more

	# multinomial likelihood observation
	logLike <- dmultinom(x = data, prob = c(p_0, p_1, p_2, p_3_or_more), log = TRUE)

	if(echo){
		E <- round(c(p_0, p_1, p_2, p_3_or_more)*sum(data))
		cat(paste(c("E[0] =", "E[1] =", "E[2] =", "E[3+] =", "logLike = "), c(E,round(logLike, 3)), collapse = "; "), "\t")
	}

	return(logLike)

}

logLike_SIRS <- function(x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and D_immunity
	R0 <- x[1]
	D_immunity <- x[2]

	logLike_0123(data, p_reinfection_SIRS, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0)

}


logLike_AoN <- function(x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and prop_immunity
	R0 <- x[1]
	prop_immunity <- x[2]

	logLike_0123(data, p_reinfection_AoN, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0)

}


logLike_PPI <- function(x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and partial_protection
	R0 <- x[1]
	partial_protection <- x[2]

	logLike_0123(data, p_reinfection_PPI, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0, echo = TRUE)

}



test <- function() {

	n_reinfection <- 2
	
	R0 <- 10
	D_infection <- 3
	n_pop <- 285
	S_0 <- n_pop - 2
	I_0 <- 1
	R_0 <- 0

	D_immunity <- 60
	prop_immunity <- 0.45
	partial_protection <- 0.95

	# p1 <- p_reinfection_SIRS(n_reinfection, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0)
	# p2 <- p_reinfection_AoN(n_reinfection, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0)
	# p3 <- p_reinfection_PPI(n_reinfection, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0)

	data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)
	logLike_0123(data, p_reinfection_SIRS, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0, echo = TRUE)
	logLike_0123(data, p_reinfection_AoN, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0, echo = TRUE)
	logLike_0123(data, p_reinfection_PPI, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0, echo = TRUE)

	df_PPI_grid <- expand.grid(R0 = 1:20, partial_protection = seq(0,1,0.1))
	df_PPI_logLike_grid <- df_PPI_grid %>% group_by(R0, partial_protection) %>% do(logLike = logLike_PPI(x = c(.$R0, .$partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))


	df_PPI_grid <- expand.grid(R0 = 1:20, partial_protection = seq(0,1,0.1))
	df_PPI_logLike_grid <- df_PPI_grid %>% group_by(R0, partial_protection) %>% do(logLike = logLike_PPI(x = c(.$R0, .$partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))

	# df_SIRS_grid <- expand.grid(R0 = 1:20, D_immunity = seq(5, 50, 5))
	# df_SIRS_logLike_grid <- df_SIRS_grid %>% group_by(R0, D_immunity) %>% do(logLike = logLike_SIRS(x = c(.$R0, .$D_immunity), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))

	df_AoN_grid <- expand.grid(R0 = 1:20, prop_immunity = seq(0,1,0.1))
	df_AoN_logLike_grid <- df_AoN_grid %>% group_by(R0, prop_immunity) %>% do(logLike = logLike_AoN(x = c(.$R0, .$prop_immunity), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))

	library(tidyr)
	df_AoN_logLike_grid2 <- df_AoN_logLike_grid %>% unnest(logLike) %>% mutate(model = "AoN") %>% rename(y = prop_immunity)
	df_PPI_logLike_grid2 <- df_PPI_logLike_grid %>% unnest(logLike) %>% mutate(model = "PPI") %>% rename(y = partial_protection)

	df_models_logLike <- df_AoN_logLike_grid2 %>% bind_rows(df_PPI_logLike_grid2)
	
	p <- ggplot(df_models_logLike, aes(x = R0, y = y)) + facet_wrap(~model, scales = "free_y")
	p <- p + geom_raster(aes(fill = log(-logLike)))
	p <- p + stat_density_2d(aes(fill = logLike), geom = "polygon")
	p
	# logLike_PPI(x = c(R0, partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0)
	# optim(par = c(R0, partial_protection), fn = logLike_PPI, data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0, method = "L-BFGS-B", lower = c(0,0), upper = c(20, 1))

}


main <- function() {

	start_me()

	# test()
}


main()
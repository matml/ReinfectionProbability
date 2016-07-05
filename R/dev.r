start_me <- function() {

	library(rPython)
	library(dplyr)
	python.load(file.path(dir_python, "Algorithm2_new_M_improved.py"))

	dir_project <<- path.expand("~/work/projects/ReinfectionProbability")
	dir_python <<- dir_project
	dir_R <<- dir_project %>% file.path("R")

}

p_reinfection_SIR <- function(gamma = 0, beta = 0, alpha_I = 0, alpha_R = 0, sigma = 0, n_pop, n_reinfection, S_0, I_0, R_0) {

	python.call("p_reinfection_SIR", gam = gamma, bet = beta, alpI = alpha_I, alpR = alpha_R, sig = sigma, N = n_pop, M = n_reinfection, i0 = I_0, s0 = S_0, r0 = R_0)

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

test <- function() {

	n_reinfection <- 5
	
	R0 <- 10
	D_infection <- 3
	n_pop <- 285
	S_0 <- n_pop - 2
	I_0 <- 1
	R_0 <- 0

	D_immunity <- 20
	prop_immunity <- 0.5
	partial_protection <- 0.9

	p1 <- p_reinfection_SIRS(n_reinfection, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0)
	p2 <- p_reinfection_AoN(n_reinfection, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0)
	p3 <- p_reinfection_PPI(n_reinfection, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0)

}


main <- function() {

	start_me()

	# test()
}


main()
init_packrat <- function(dir_project, force = FALSE) {

	library(packrat)

	setwd(dir = dir_project)

	if(file.exists("packrat") & !force){

		packrat::on()

	} else {

		packrat::init(options = list(vcs.ignore.src = TRUE))

	}
}

my_library <- function(x) {

	if(!require(x, character.only = TRUE)){
		install.packages(x)
	}

	library(x, character.only = TRUE)

}

start_me <- function() {


	dir_project <<- path.expand("~/work/projects/ReinfectionProbability")
	dir_python <<- file.path(dir_project, "python")
	dir_R <<- file.path(dir_project, "R")

	init_packrat(dir_R)

	my_library("rPython")
	# my_library("plyr")
	my_library("dplyr")
	my_library("tidyr")
	my_library("readr")
	my_library("ggplot2")
	# my_library("multidplyr")

	python.load(file.path(dir_python, "SIR_Reinfection.py"))
	python.load(file.path(dir_python, "SIS_Reinfection.py"))

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
	prob <- c(p_0, p_1, p_2, p_3_or_more)
	if(all(is.finite(prob)) & !all(prob == 0) & !any(prob < 0)){
		logLike <- dmultinom(x = data, prob = prob, log = TRUE)		
	} else {
		logLike <- as.numeric(NA)
	}

	if(echo){
		E <- round(c(p_0, p_1, p_2, p_3_or_more)*sum(data))
		cat(paste(c("E[0] =", "E[1] =", "E[2] =", "E[3+] =", "logLike = "), c(E,round(logLike, 3)), collapse = "; "), "\t")
	}

	ans <- data_frame(p_0, p_1, p_2, p_3_or_more, logLike)

	return(ans)

}

logLike_SIRS <- function(i = 1, x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and D_immunity
	R0 <- x[1]
	D_immunity <- x[2]

	x <- system.time(ans <- logLike_0123(data, p_reinfection_SIRS, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0))
	# ans <- logLike_0123(data, p_reinfection_SIRS, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0)

	message(paste("Iteration:", i, "R0 =", R0, "D_immunity =", D_immunity, "time:", as.numeric(x[["elapsed"]])))

	return(ans)

}


logLike_AoN <- function(i = 1, x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and prop_immunity
	R0 <- x[1]
	prop_immunity <- x[2]

	x <- system.time(ans <- logLike_0123(data, p_reinfection_AoN, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0))
	# ans <- logLike_0123(data, p_reinfection_AoN, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0)

	print(paste("Iteration:", i, "R0 =", R0, "prop_immunity =", prop_immunity, "time:", as.numeric(x[["elapsed"]])))

	return(ans)

}


logLike_PPI <- function(i = 1, x, data, D_infection, n_pop, S_0, I_0, R_0) {

	# x contains R0 and partial_protection
	R0 <- x[1]
	partial_protection <- x[2]

	x <- system.time(ans <- logLike_0123(data, p_reinfection_PPI, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0))
	# ans <- logLike_0123(data, p_reinfection_PPI, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0)

	print(paste("Iteration:", i, "R0 =", R0, "partial_protection =", partial_protection, "time:", as.numeric(x[["elapsed"]])))
	

	return(ans)

}


grid_logLike_SIRS <- function(R0, D_immunity, ...) {

	expand.grid(R0 = R0, D_immunity = D_immunity) %>% 
	group_by(R0, D_immunity) %>% 
	do(logLike_SIRS(i = 1, x = c(.$R0, .$D_immunity), ...)) %>% 
	ungroup

}

grid_logLike_AoN <- function(R0, prop_immunity, ...) {

	expand.grid(R0 = R0, prop_immunity = prop_immunity) %>% 
	group_by(R0, prop_immunity) %>% 
	do(logLike_AoN(i = 1, x = c(.$R0, .$prop_immunity), ...)) %>% 
	ungroup

}

grid_logLike_PPI <- function(R0, partial_protection, ...) {

	expand.grid(R0 = R0, partial_protection = partial_protection) %>% 
	group_by(R0, partial_protection) %>% 
	do(logLike_PPI(i = 1, x = c(.$R0, .$partial_protection), ...)) %>% 
	ungroup

}

grid_logLike_SIRS_parallel <- function(R0, D_immunity, i_job, n_job, ...) {

	expand.grid(R0 = R0, D_immunity = D_immunity) %>% 
	mutate(i = 1:n()) %>% 
	filter(i%in%seq(i_job, n(), n_job)) %>% 
	group_by(R0, D_immunity) %>% 
	do(logLike_SIRS(i = .$i, x = c(.$R0, .$D_immunity), ...)) %>% 
	ungroup

}

grid_logLike_AoN_parallel <- function(R0, prop_immunity, i_job, n_job, ...) {

	expand.grid(R0 = R0, prop_immunity = prop_immunity) %>% 
	mutate(i = 1:n()) %>% 
	filter(i%in%seq(i_job, n(), n_job)) %>% 
	group_by(R0, prop_immunity) %>% 
	do(logLike_AoN(i = .$i, x = c(.$R0, .$prop_immunity), ...)) %>% 
	ungroup

}


grid_logLike_PPI_parallel <- function(R0, partial_protection, i_job, n_job, ...) {

	expand.grid(R0 = R0, partial_protection = partial_protection) %>% 
	mutate(i = 1:n()) %>% 
	filter(i%in%seq(i_job, n(), n_job)) %>% 
	group_by(R0, partial_protection) %>% 
	do(logLike_PPI(i = .$i, x = c(.$R0, .$partial_protection), ...)) %>% 
	ungroup

}


grid_logLike_AoN_multidplyr <- function(R0, prop_immunity, data, D_infection, n_pop, S_0, I_0, R_0) {

	# n_reinfection <- 2
	
	# R0 <- 10
	# D_infection <- 3
	# n_pop <- 285
	# S_0 <- n_pop - 2
	# I_0 <- 1
	# R_0 <- 0

	# D_immunity <- 60
	# prop_immunity <- 0.45
	# partial_protection <- 0.95

	# data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)


	my_cluster <- create_cluster(cores = 2)
	set_default_cluster(my_cluster)
	
	# cluster_call(my_cluster, logLike_0123)
	# cluster_call(my_cluster, p_reinfection_AoN)
	# cluster_call(my_cluster, p_reinfection_SIR)
	# cluster_assign_value(my_cluster, "logLike_0123", logLike_0123)
	# cluster_assign_value(my_cluster, "D_infection", D_infection)
	# cluster_assign_value(my_cluster, "n_pop", n_pop)
	# cluster_assign_value(my_cluster, "S_0", S_0)
	# cluster_assign_value(my_cluster, "I_0", I_0)
	# cluster_assign_value(my_cluster, "R_0", R_0)
	# cluster_assign_value(my_cluster, "p_reinfection_SIR", p_reinfection_SIR)
	# cluster_assign_value(my_cluster, "start_me", start_me)
	# cluster_call(my_cluster, start_me)


	# cluster_call(my_cluster, logLike_AoN(x = c(2, 0.5), data, D_infection, n_pop, S_0, I_0, R_0))
	

	# cluster_ls(my_cluster)
	# cluster_assign_value(my_cluster, "data", data)

	# my_cluster %>% 
	# cluster_call(logLike_AoN) %>% 
	party_df <- expand.grid(R0 = R0, prop_immunity = prop_immunity) %>% 
	partition(R0, prop_immunity, cluster = my_cluster) 



	# %>% 
	cluster_library(party_df, "rPython")
	cluster_assign_value(party_df, "logLike_AoN", logLike_AoN)
	cluster_assign_value(party_df, "logLike_0123", logLike_0123)
	cluster_assign_value(party_df, "p_reinfection_AoN", p_reinfection_AoN)

	cluster_ls(party_df)

	party_df %>% 
	do(logLike_AoN(x = c(.$R0, .$prop_immunity), data, D_infection, n_pop, S_0, I_0, R_0)) 

	# %>% 
	# collect()

	# vignette(multidplyr)

}


test <- function() {

	n_reinfection <- 2
	
	R0 <- 47.5
	D_infection <- 3
	n_pop <- 284
	S_0 <- n_pop - 1
	I_0 <- 1
	R_0 <- 0

	D_immunity <- 42
	prop_immunity <- 0.45
	partial_protection <- 0.95

	data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)

	ans_SIR <- logLike_SIRS(i = 1, x = c(R0, D_immunity), data, D_infection, n_pop, S_0, I_0, R_0)
	print(ans_SIR)
	ans_AoN <- logLike_AoN(i = 1, x = c(R0, prop_immunity), data, D_infection, n_pop, S_0, I_0, R_0)
	print(ans_AoN)
	ans_PPI <- logLike_PPI(i = 1, x = c(R0, partial_protection), data, D_infection, n_pop, S_0, I_0, R_0)
	print(ans_PPI)



	# data_frame(
	# 	model = c("SIRS", "AoN", "PPI"), 
	# 	R0 = R0, 
	# 	D_infection = D_infection, 
	# 	D_immunity = c(D_immunity, NA, NA), 
	# 	prop_immunity = c(NA, prop_immunity, NA), 
	# 	partial_protection = c(NA, NA, partial_protection)
	# 	) %>% bind_cols(bind_rows(ans_SIR, ans_AoN, and_PPI)) %>% 
	# write_csv("../doc/check.csv")



	# p1 <- p_reinfection_SIRS(n_reinfection, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0)
	# system.time(p2 <- p_reinfection_AoN(n_reinfection, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0))
	# p3 <- p_reinfection_PPI(n_reinfection, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0)

	# data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)
	# logLike_0123(data, p_reinfection_SIRS, R0, D_infection, D_immunity, n_pop, S_0, I_0, R_0, echo = TRUE)
	# system.time(ll <- logLike_0123(data, p_reinfection_AoN, R0, D_infection, prop_immunity, n_pop, S_0, I_0, R_0, echo = TRUE))
	# logLike_0123(data, p_reinfection_PPI, R0, D_infection, partial_protection, n_pop, S_0, I_0, R_0, echo = TRUE)

	# df_PPI_grid <- expand.grid(R0 = 1:20, partial_protection = seq(0,1,0.1))
	# df_PPI_logLike_grid <- df_PPI_grid %>% group_by(R0, partial_protection) %>% do(logLike = logLike_PPI(x = c(.$R0, .$partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))


	# df_PPI_grid <- expand.grid(R0 = 1:20, partial_protection = seq(0,1,0.1))
	# df_PPI_logLike_grid <- df_PPI_grid %>% group_by(R0, partial_protection) %>% do(logLike = logLike_PPI(x = c(.$R0, .$partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))

	# df_SIRS_grid <- expand.grid(R0 = 1:20, D_immunity = seq(5, 50, 5))
	# df_SIRS_logLike_grid <- df_SIRS_grid %>% group_by(R0, D_immunity) %>% do(logLike = logLike_SIRS(x = c(.$R0, .$D_immunity), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))


	# grid_logLike_AoN_parallel(R0 = seq(1, 50, 0.5), prop_immunity = seq(0.01,1, 0.01), i_job = 2, n_job = 10,  data, D_infection, n_pop, S_0, I_0, R_0)

	# grid_logLike_SIRS_parallel(R0 = 3:4, D_immunity = 4, i_job = 1, n_job = 2,  data, D_infection, n_pop, S_0, I_0, R_0)

	# df_AoN_grid <- expand.grid(R0 = 1:20, prop_immunity = seq(0,1,0.1))
	# df_AoN_logLike_grid <- df_AoN_grid %>% group_by(R0, prop_immunity) %>% do(logLike = logLike_AoN(x = c(.$R0, .$prop_immunity), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0))

	
	# df_AoN_logLike_grid2 <- df_AoN_logLike_grid %>% unnest(logLike) %>% mutate(model = "AoN") %>% rename(y = prop_immunity)
	# df_PPI_logLike_grid2 <- df_PPI_logLike_grid %>% unnest(logLike) %>% mutate(model = "PPI") %>% rename(y = partial_protection)

	# df_models_logLike <- df_AoN_logLike_grid2 %>% bind_rows(df_PPI_logLike_grid2)
	
	# p <- ggplot(df_models_logLike, aes(x = R0, y = y)) + facet_wrap(~model, scales = "free_y")
	# p <- p + geom_raster(aes(fill = log(-logLike)), interpolate =TRUE)
	# # p <- p + stat_density_2d(aes(fill = logLike), geom = "polygon")
	# p
	# # logLike_PPI(x = c(R0, partial_protection), data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0)
	# # optim(par = c(R0, partial_protection), fn = logLike_PPI, data = data, D_infection = D_infection, n_pop = n_pop, S_0 = S_0, I_0 = I_0, R_0 = R_0, method = "L-BFGS-B", lower = c(0,0), upper = c(20, 1))

}


main_cluster <- function() {

	start_me()

	model <- Sys.getenv("model")
	analysis <- Sys.getenv("analysis")
	jobDir <- Sys.getenv("jobDir")
	i_job <- as.numeric(Sys.getenv("process")) + 1 + 3
	n_job <- as.numeric(Sys.getenv("replicate")) + 3

	message("Model:",sQuote(model), "\nanalysis:", analysis, "\ni_job:", i_job, "\nn_job:", n_job)

	D_infection <- 2
	n_pop <- 284
	S_0 <- n_pop - 1
	I_0 <- 1
	R_0 <- 0

	R0 <- seq(0.5, 50, 0.5)
	# R0 <- seq(9, 10, 1)
	# R0 <- 10
	prop_immunity <- seq(0.01,1, 0.01)
	# prop_immunity <- seq(0.5, 0.6, 0.1)
	# prop_immunity <- 0.8
	partial_protection <- prop_immunity
	D_immunity <- prop_immunity*100

	data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)
	
	if(model == "AoN"){

		ans <- grid_logLike_AoN_parallel(R0 = R0, prop_immunity = prop_immunity, i_job = i_job, n_job = n_job,  data, D_infection, n_pop, S_0, I_0, R_0)

	} else if (model == "PPI") {

		ans <- grid_logLike_PPI_parallel(R0 = R0, partial_protection = partial_protection, i_job = i_job, n_job = n_job,  data, D_infection, n_pop, S_0, I_0, R_0)


	} else if (model == "SIRS") {

		ans <- grid_logLike_SIRS_parallel(R0 = R0, D_immunity = D_immunity, i_job = i_job, n_job = n_job,  data, D_infection, n_pop, S_0, I_0, R_0)

	}

	saveRDS(ans, file.path(jobDir, sprintf("ans_%s_job_%s_of_%s.rds", analysis, i_job, n_job)))

}

main_local <- function() {

	model <- "SIRS"
	analysis <- "large"
	jobDir <- file.path("rds", model)
	i_job <- 1
	n_job <- 1

	message("Model:",sQuote(model), "\nanalysis:", analysis, "\ni_job:", i_job, "\nn_job:", n_job)

	D_infection <- 2
	n_pop <- 284
	S_0 <- n_pop - 1
	I_0 <- 1
	R_0 <- 0

	R0 <- seq(0.5, 50, 0.5)
	# R0 <- seq(9, 10, 1)
	# R0 <- 10
	prop_immunity <- seq(0.01,1, 0.01)
	# prop_immunity <- seq(0.5, 0.6, 0.1)
	# prop_immunity <- 0.8
	partial_protection <- prop_immunity
	D_immunity <- prop_immunity*100

	data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)
	
	ans <- grid_logLike_SIRS_parallel(R0 = R0, D_immunity = D_immunity, i_job = i_job, n_job = n_job,  data, D_infection, n_pop, S_0, I_0, R_0)

	saveRDS(ans, file.path(jobDir, sprintf("ans_%s_job_%s_of_%s.rds", analysis, i_job, n_job)))


}

combine_outputs <- function() {

	import_res <- function(model) {
		
		files <- list.files(file.path("rds", model))
		files <- files[grepl("ans", files)]

		df_ans <- NULL

		for(file in files){
			df_ans <- df_ans %>% bind_rows(readRDS(file.path("rds", model, file)))
		}

		x_var <- names(df_ans)[1]
		y_var <- names(df_ans)[2]

		names(df_ans)[1:2] <- c("x", "y")

		df_ans <- mutate(df_ans, x_var = x_var, y_var = y_var)

		return(df_ans)
	}



	df_res_raw <- data_frame(model = c("AoN", "PPI", "SIRS")) %>% group_by(model) %>% do(import_res(.$model)) %>% ungroup 

	# parameters without logLike (negative p_1)
	# put p_1 = 0 and renormalise
	# and compute log-like
	data <- c(n_0 = 11, n_1 = 181, n_2 = 92, n_3_or_more = 0)
	df_res_edit <- df_res_raw %>% filter(is.na(logLike)) %>% mutate(
		p_1 = 0, 
		sum_p = p_0 + p_1 + p_2 + p_3_or_more, 
		p_0 = p_0/sum_p,
		p_1 = p_1/sum_p,
		p_2 = p_2/sum_p,
		p_3_or_more = p_3_or_more/sum_p
		) %>% 
	select(-sum_p) %>% 
	group_by(x,y) %>% 
	mutate(logLike = dmultinom(x = data, prob = c(p_0, p_1, p_2, p_3_or_more), log = TRUE))

	df_res <- df_res_raw %>% filter(!is.na(logLike)) %>% bind_rows(df_res_edit) %>% gather(variable, value, c(p_0, p_1, p_2, p_3_or_more, logLike))


	df_res_raw %>% group_by(model) %>% filter(logLike == max(logLike, na.rm = TRUE))

	# TdC paper MLE
	df_MLE_TdC <- frame_data(
		~model, ~x, ~y,
		"AoN", 11.20, 0.53,
		"PPI", 6.92, 0.83
		)

	# optimal probabilities
	df_optimal_prob <- data_frame(variable = c("p_0", "p_1", "p_2", "p_3_or_more"), n = data) %>% mutate(p_optim  = n/sum(n))

	# logLike + confidence intervals
	df_plot_ll <- df_res %>% filter(variable == "logLike")
	df_plot_ll_max_CI <- df_plot_ll %>% group_by(model) %>% filter(value == max(value)) %>% ungroup %>% mutate(min_logLike_CI = value - qchisq(0.95, 2)/2)

	all_models <- c("AoN", "PPI", "SIRS")

	for(my_model in all_models) {

		df_plot_ll_model <- df_plot_ll %>% filter(model == my_model) #%>% filter(y > 0.9)
		df_plot_ll_max_CI_model <- df_plot_ll_max_CI %>% filter(model == my_model)
		df_MLE_TdC_model <- df_MLE_TdC %>% filter(model == my_model)

		p <- ggplot(df_plot_ll_model, aes(x = x, y = y)) #+ facet_grid(model~variable, scales = "free_y")
		p <- p + geom_tile(aes(fill = -log(-value))) 
		p <- p + geom_contour(aes(z = value), breaks = df_plot_ll_max_CI_model$min_logLike_CI, linetype = "dashed", color = "black")
		p <- p + geom_point(data = df_plot_ll_max_CI_model)
		p <- p + geom_point(data = df_MLE_TdC_model, aes(x = x, y = y), shape = 17)
		p <- p + scale_fill_distiller("Log-Likelihood", palette = "Reds", direction = 1, labels = function(x) {round(-exp(-x))})
		p <- p + theme_minimal() + xlab("R0") + ylab(switch(my_model, AoN = "Proportion with immunity", "PPI" = "Partial protection", SIRS = "Duration of immunity"))
		# print(p)
		ggsave(paste0("../Manuscript/fig/", my_model, "_loglike.pdf"), width = 8, height = 6)

		# probabilities
		df_plot_prob_model <- df_res %>% filter(variable != "logLike", model == my_model) %>% left_join(df_optimal_prob, "variable") %>% mutate(diff_prob = abs(value - p_optim))
		df_plot_ll_model_prob <-  data_frame(variable = df_plot_prob_model$variable %>% unique, df = list(df_plot_ll_model %>% select(-variable))) %>% unnest(df)
		df_plot_ll_max_CI_model_prob <- data_frame(variable = df_plot_prob_model$variable %>% unique, df = list(df_plot_ll_max_CI_model %>% select(-variable))) %>% unnest(df)

		p <- ggplot(df_plot_prob_model, aes(x = x, y = y)) + facet_wrap(~variable, nrow = 1)
		p <- p + geom_raster(aes(fill = diff_prob), interpolate = TRUE) 
		p <- p + geom_contour(data = df_plot_ll_model_prob, aes(z = value), breaks = df_plot_ll_max_CI_model$min_logLike_CI, linetype = "dashed", color = "black")
		p <- p + geom_point(data = df_plot_ll_max_CI_model_prob)
		p <- p + geom_point(data = df_MLE_TdC_model, aes(x = x, y = y), shape = 17)
		p <- p + scale_fill_distiller("Diff P", palette = "Reds", direction = -1)
		p <- p + theme_minimal() + xlab("R0") + ylab(switch(my_model, AoN = "Proportion with immunity", "PPI" = "Partial protection", SIRS = "Duration of immunity"))
		# print(p)
		ggsave(paste0("../Manuscript/fig/", my_model, "_prob.pdf"), width = 8, height = 3)

	}

}


check_PPI_traj <- function() {

	# my_library("devtools")
	# my_library("httr")
	# install_github("sbfnk/fitR")
	# library(fitR)
	data(SIR)


	PPI_state.names <- c("S", "I", "R", "I2", "R2", "I3", "R3", "incidence")	
	PPI_theta.names <- c("R0", "D_inf", "partial_protection")

	PPI_simulate <- function (theta, init.state, times) 
	{
		PPI_transitions <- list(
			c(S = -1, I = 1, incidence = 1), 
			c(I = -1, R = 1), 
			c(I2 = 1, R = -1, incidence = 1), 
			c(I2 = -1, R2 = 1), 
			c(R2 = -1, I3 = 1, incidence = 1), 
			c(I3 = -1, R3 = 1), 
			c(R3 = -1, I3 = 1, incidence = 1)
			)

		PPI_rateFunc <- function(x, parameters, t) {
			beta <- parameters[["R0"]]/parameters[["D_inf"]]
			nu <- 1/parameters[["D_inf"]]
			sigma <- parameters[["partial_protection"]]

			S <- x[["S"]]
			I <- x[["I"]]
			R <- x[["R"]]
			I2 <- x[["I2"]]
			R2 <- x[["R2"]]
			I3 <- x[["I3"]]
			R3 <- x[["R3"]]
			incidence <- x[["incidence"]]
			N <- S + I + R + I2 + R2 + I3 + R3
			
			infection_force <- beta * (I+I2+I3)/N

			return(c(infection_force * S, nu * I, (1-sigma)*infection_force*R, nu * I2, (1-sigma)*infection_force*R2, nu * I3, (1-sigma)*infection_force*R3))

		}

		df <- simulateModelStochastic(theta, init.state, times, 
			PPI_transitions, PPI_rateFunc)

		df <- df %>% mutate(incidence = c(0, diff(incidence)))

		return(df)
	}


	PPI_rPointObs <- function (model.point, theta) {
		obs.point <- rpois(n = 1, lambda = model.point[["incidence"]])
		return(c(obs = obs.point))
	}

	PPI_dPointObs <- function (data.point, model.point, theta, log = FALSE) {
		return(dpois(x = data.point[["obs"]], lambda = model.point[["incidence"]], 
			log = log))
	}

	PPI_dprior <- function (theta, log = FALSE) {
		log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 100, 
			log = TRUE)
		log.prior.D <- dunif(theta[["D_inf"]], min = 0, max = 30, 
			log = TRUE)
		log.prior.pp <- dunif(theta[["partial_protection"]], min = 0, max = 1, 
			log = TRUE)
		log.sum <- log.prior.R0 + log.prior.D + log.prior.pp
		return(ifelse(log, log.sum, exp(log.sum)))
	}


	PPI <- fitmodel(
		name = "PPI",
		state.names = PPI_state.names,
		theta.names = PPI_theta.names,
		simulate = PPI_simulate,
		dprior = PPI_dprior,
		rPointObs = PPI_rPointObs,
		dPointObs = PPI_dPointObs)


	data(FluTdC1971)
	df_data <- FluTdC1971 #%>% select(-time) %>% dplyr::rename(time = date) %>% mutate(time = as.Date(time)) %>% as_data_frame

	plotFit(PPI, theta = c("R0" = 23, "D_inf" = 2, "partial_protection" = 0.98), init.state = c(S = 283, I = 1, R = 0, I2 = 0, R2 = 0, I3 = 0, R3 = 0, incidence = 0), data = df_data, n.replicates = 100, all.vars = TRUE, summary = FALSE)


}

main <- function() {

	start_me()

	# test()
	# main_local()


}


main()
# main_cluster()



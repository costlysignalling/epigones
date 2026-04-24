library(rethinking)
source("00_functions.R")
source("01_process_data.R")
source("02_define_traditions.R")
traditions<-keep_first_pioneer(traditions) #Keep just the most important pioneer

# after loading your RDS files and computing D and d$namlong...
data <- prepare_stan_data(D, d$namlong, traditions)
dA<-data$dA

#there are also other calculated variables for nice visualizations of predictions later
str(data)

# now you can hand dA straight into cmdstan:
m   <- cmdstan_model("m_step_back7.stan")

# Create an init function that returns valid values. It facilitates beginning of the sampling.
init_fun <- function() {
  list(L_Rho_e = diag(1, dA$Nd),         # valid Cholesky factor
       sigma_e = rep(1, dA$Nd))}           # positive values for sigma_e

fit <- m$sample(data = dA, chains = 8, parallel_chains = 8, 
                iter_warmup = 1000, iter_sampling = 1000, init = init_fun)

save(fit, file = "fit_step_back7.Rdata")

post<-extract.samples(fit)
save(post, file = "post_step_back7.Rdata")

paste(names(post),collapse=",")
#reducing post for transfer to another computer (stripping of lam, mu, lsig, sig) is a good idea
post.r<-post[c("a","lsb","lnud","lnur","lnup","dsb","dnud","dnur","dnup","sb","nud","nur","nup")]
save(post.r, file = "post_reduced_step_back7.Rdata")

precis(fit)
precis(post.r)

tab<-precis(fit,depth=2)
write.table(tab, file = "tab_step_back7.txt", sep = "\t",col.names = NA)

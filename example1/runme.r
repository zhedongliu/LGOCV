#get the script directory
script_dir <- dirname(sys.frame(1)$ofile)
#load libraries
if(!require("rstan",quietly=TRUE)){
	install.packages("rstan")
	require("rstan",quietly=TRUE)
}

if(!require("INLA",quietly=TRUE)){
	install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
	require("INLA",quietly=TRUE)
}

#data generation and plot data
set.seed(2020)
beta = log(10)
group_size = 10
group_each = 10
n = group_size*group_each
group_mu = rnorm(group_size)
eta = beta + rep(group_mu,each = group_each)
id = rep(1:group_size,each= group_each)
paper_data = TRUE

#normal responses
y_normal = rnorm(n,eta,sd = .1)
if(paper_data){y_normal = readRDS(paste0(script_dir,"/paper_data/data_norm.RDS"))$y}

#binomial responses
N_trials = 20
y_binom = rbinom(n = n,size = N_trials,prob = 1/(1+exp(-eta)))
if(paper_data){y_binom = readRDS(paste0(script_dir,"/paper_data/data_binom.RDS"))$y}
#exponential responses
y_exp = rexp(n,exp(eta))
if(paper_data){y_exp = readRDS(paste0(script_dir,"/paper_data/data_exp.RDS"))$y}
#Run and store stan models
n_leave =  n - group_each
Niter = 1e5
Nwarmup = Niter/10
Nthin = 1
Nchains = 1
#This will take enormous amount of time. I will suggests to run it on a multi-core machine and parallelize it. 
rerun = TRUE
if(rerun){
	for(group_interest in 1:group_size){
  		y_leave = y_normal[-which(id == group_interest)]
  		id_leave = id[-which(id == group_interest)]
  		stan_model_dir = paste0(script_dir,"/rstan_model/group_normal.stan")
  		res_stan = stan(file = stan_model_dir,
                 		data = list(y = y_leave,id=id_leave,n=n_leave,m=group_size),
                  		iter = Niter,
                  		warmup = Nwarmup,
                  		thin = Nthin,
                  		chains = Nchains)
  		res_stan_list = extract(res_stan, permuted = TRUE)
  		file_name = paste0(script_dir,"/rstan_result/stan_res_normal",group_interest,".rds")
  		saveRDS(object = res_stan_list,file = file_name)
	}

	for(group_interest in 1:group_size){
  		y_leave = y_binom[-which(id == group_interest)]
  		id_leave = id[-which(id == group_interest)]
  		stan_model_dir = paste0(script_dir,"/rstan_model/group_binom.stan")
  		res_stan = stan(file = stan_model_dir,
                  	data = list(y = y_leave,id=id_leave,n=n_leave,m=group_size),
                  	iter = Niter,
                  	warmup = Nwarmup,
                  	thin = Nthin,
                  	chains = Nchains)
  		res_stan_list = extract(res_stan, permuted = TRUE)
  		file_name = paste0(script_dir,"/rstan_result/stan_res_binom",group_interest,".rds")
  		saveRDS(object = res_stan_list,file = file_name)
	}

	for(group_interest in 1:group_size){
  		y_leave = y_exp[-which(id == group_interest)]
  		id_leave = id[-which(id == group_interest)]
  		stan_model_dir = paste0(script_dir,"/rstan_model/group_exp.stan")
  		res_stan = stan(file = stan_model_dir,
                  		data = list(y = y_leave,id=id_leave,n=n_leave,m=group_size),
                  		iter = Niter,
                  		warmup = Nwarmup,
                  		thin = Nthin,
                  		chains = Nchains)
  		res_stan_list = extract(res_stan, permuted = TRUE)
  		file_name = paste0(script_dir,"/rstan_result/stan_res_exp",group_interest,".rds")
  		saveRDS(object = res_stan_list,file = file_name)
	}
}





#Fit INLA models, compare with Stan result
Design = matrix(0,100,2)
tau_g = seq(-3.5,1.5,length.out = 100)
Design[,1] = tau_g
Design[,2] = 1
formula = y~ 1 + f(id,model = "iid",hyper = list(prec = list(prior =  "normal",param = c(0, 1e-4))))


#####################
res_normal = inla(formula, 
           family= "gaussian", 
           data=list(id=id,y=y_normal),
           inla.mode = "experimental",
           control.family = list(hyper = list(prec = list(initial = log(100),fixed = TRUE))),
           control.compute = list(control.gcpo = list(enable = TRUE,num.level.sets = 1),config = T),
           control.fixed = list(prec.intercept = 1e-4),
           control.inla = list(int.strategy = "user",int.design = Design))

gcpo_stan_normal = c()
for(group_interest in 1:10){
   file_name = paste0(script_dir,"/rstan_result/stan_res_normal",group_interest,".rds")
   res_stan_list = readRDS(file = file_name)
   eta = res_stan_list$mu+res_stan_list$mu_group[,group_interest]
   #eta = res_stan_list$mu+res_stan_list$mu_group
   y_here = y_normal[which(id == group_interest)]
   for(i in 1:10){
     fy = dnorm(x = y_here[i],mean = eta,sd = .1)
     temp = mean(fy)
     gcpo_stan_normal = c(gcpo_stan_normal,temp)
   }
}




###################
res_binom = inla(formula, 
           family= "binomial",
           Ntrials = 20,
           data=list(id=id,y=y_binom),
           inla.mode = "experimental",
           control.compute = list(control.gcpo = list(enable = TRUE,num.level.sets = 1),config = T),
           control.fixed = list(prec.intercept = 1e-4),
           control.inla = list(int.strategy = "user",int.design = Design))

gcpo_stan_binom = c()
for(group_interest in 1:10){
   file_name = paste0(script_dir,"/rstan_result/stan_res_binom",group_interest,".rds")
   res_stan_list = readRDS(file = file_name)
   eta = res_stan_list$mu+res_stan_list$mu_group[,group_interest]
   #eta = res_stan_list$mu+res_stan_list$mu_group
   y_here = y_binom[which(id == group_interest)]
   for(i in 1:10){
     fy = dbinom(x = y_here[i],size = 20,prob = 1/(1+exp(-eta)))
     temp = mean(fy)
     gcpo_stan_binom = c(gcpo_stan_binom,temp)
   }
}




#####################
res_exp = inla(formula, 
           family= "exponential", 
           data=list(id=id,y=y_exp),
           inla.mode = "experimental",
           control.compute = list(control.gcpo = list(enable = TRUE,num.level.sets = 1),config = T),
           control.fixed = list(prec.intercept = 1e-4),
           control.inla = list(int.strategy = "user",int.design = Design))

gcpo_stan_exp = c()
 for(group_interest in 1:10){
   file_name = paste0(script_dir,"/rstan_result/stan_res_exp",group_interest,".rds")
   res_stan_list = readRDS(file = file_name)
   eta = res_stan_list$mu+res_stan_list$mu_group[,group_interest]
   #eta = res_stan_list$mu+res_stan_list$mu_group
   y_here = y_exp[which(id == group_interest)]
   for(i in 1:10){
     fy = dexp(x = y_here[i],rate = exp(eta))
     temp = mean(fy)
     gcpo_stan_exp = c(gcpo_stan_exp,temp)
   }
 }




pdf(file = paste0(script_dir,"/data_and_compare.pdf"),
    width = 10,
    height = 15)
par(mfrow = c(3,2))
plot(id,y_normal)
plot(gcpo_stan_normal,res_normal$gcpo$gcpo)
abline(a = 0,b = 1)
plot(id,y_binom)
plot(gcpo_stan_binom,res_binom$gcpo$gcpo)
abline(a = 0,b = 1)
plot(id,y_exp)
plot(gcpo_stan_exp,res_exp$gcpo$gcpo)
abline(a = 0,b = 1)
dev.off()















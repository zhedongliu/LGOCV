#get the script directory
script_dir <- dirname(sys.frame(1)$ofile)
#load libraries
if(!require("INLA",quietly=TRUE)){
	install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
	require("INLA",quietly=TRUE)
}

###simulate data
set.seed(2023)
n = 2000
mu = 2
rho = 0.9
var_add = 1
var_marginal = var_add/(1-rho^2)
prec_marginal = 1/var_marginal
eta = mu + arima.sim(model = list(ar = rho),n = n,sd = sqrt(var_add))
y = eta + rnorm(n,sd = 0.1)
id = 1:n
paper_data = TRUE
if(paper_data){y = readRDS(paste0(script_dir,"/paper_data.RDS"))}

###fit model
formula = y ~ 1 + f(id, model = "ar1",hyper = list(rho = list(initial = log((1+rho)/(1-rho)),fixed  = TRUE),prec = list(initial = log(prec_marginal),fixed = TRUE)))
res_observe = inla(formula = formula,
                   family = "gaussian",
                   inla.mode = "experimental",
                   control.family = list(hyper = list(prec=list(initial = log(100),fixed = T))),
                   data = list(id = id,y=y),
                   control.inla = list(int.strategy = "eb"),
                   control.compute = list(config = TRUE))

###compute CVs
lfocv = numeric(10)
days_start = 1500
res_refit = numeric(10)
for(days_ahead in 1:10){
    groups = list()
    for(i in days_start:n){
        groups[[i]] = (i-days_ahead+1):n
    }
    lfocv[days_ahead] = mean(log(inla.group.cv(result = res_observe,groups=groups)$cv[days_start:n]))
    print(days_ahead)
}

lgocv = numeric(10)
for(num.level.sets in 1:10){
    lgocv[num.level.sets] = mean(log(inla.group.cv(result = res_observe,num.level.sets = num.level.sets,strategy = "prior")$cv))
    print(num.level.sets)
}


###Plot
lfocv_to_ahead = splinefun(lfocv,1:10,method = "natural")
lgocv_to_m = splinefun(lgocv,1:10,method = "natural")


pdf(file = paste0(script_dir,"/plot.pdf"),
    width = 20,
    height = 10)
par(mfrow = c(1,2))
#################################################################
xx = seq(min(lfocv) -0.02,max(lfocv) + 0.02,length.out = 100)
plot(lfocv_to_ahead(xx),xx,type="l")
points(1:10,lfocv)
##########################################################################
xx = 1:10
yy = lfocv_to_ahead(lgocv)
plot(xx,yy,xlim=c(0,11))
dev.off()



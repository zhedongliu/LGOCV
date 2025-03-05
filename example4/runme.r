#get the script directory
script_dir <- dirname(sys.frame(1)$ofile)
#load libraries
if(!require("INLA",quietly=TRUE)){
	install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
	require("INLA",quietly=TRUE)
}
inla.setOption(num.threads ="1:8")

if(!require("sf",quietly=TRUE)){
	install.packages("sf")
	require("sf",quietly=TRUE)
}
if(!require("spdep",quietly=TRUE)){
	install.packages("spdep")
	require("spdep",quietly=TRUE)
}
if(!require("data.table",quietly=TRUE)){
	install.packages("data.table")
	require("data.table",quietly=TRUE)
}
if(!require("tsModel",quietly=TRUE)){
	install.packages("tsModel")
	require("tsModel",quietly=TRUE)
}
if(!require("dlnm",quietly=TRUE)){
	install.packages("dlnm")
	require("dlnm",quietly=TRUE)
}


#load data and preporcess
map <- read_sf(paste0(script_dir,"/data/shape_brazil.shp"))
nb.map <- poly2nb(as_Spatial(map$geometry))
g.file <- paste0(script_dir,"/data/map.graph")
nb2INLA(g.file, nb.map)
grid <- read.csv(paste0(script_dir,"/data/br_states_grid.csv"))
data <- fread(paste0(script_dir,"/data/data_2000_2019.csv"), header = T)
nlag = 6
lag_tmin <- Lag(data$tmin, group = data$micro_code, k = 0:nlag)
# Maximum temperature (Tmax)
lag_tmax <- Lag(data$tmax, group = data$micro_code, k = 0:nlag)
# Palmer drought severity index (PDSI)
lag_pdsi <- Lag(data$pdsi, group = data$micro_code, k = 0:nlag)
lag_tmin <- lag_tmin[data$year > 2000,]
lag_tmax <- lag_tmax[data$year > 2000,]
lag_pdsi <- lag_pdsi[data$year > 2000,]
data <- data[data$year > 2000,]
data$time <- data$time - 12
# total number of months
ntime <- length(unique(data$time))
# total number of years
nyear <- length(unique(data$year))
# total number of microregions
nmicro <- length(unique(data$micro_code))
# total number of states
nstate <- length(unique(data$state_code))
lagknot = equalknots(0:nlag, 2)
# Tmin
var <- lag_tmin
basis_tmin <- crossbasis(var,
                         argvar = list(fun = "ns", knots = equalknots(data$tmin, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))
# Tmax
var <- lag_tmax
basis_tmax <- crossbasis(var,
                         argvar = list(fun = "ns", knots = equalknots(data$tmax, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))
# PDSI
var <- lag_pdsi
basis_pdsi <- crossbasis(var,
                         argvar = list(fun="ns", knots = equalknots(data$pdsi, 2)),
                         arglag = list(fun="ns", knots = lagknot))

# set indicator to zero at point of interest (centring point)
# re-parameterise model to extract different predictions
urban_ind1 <- data$urban - quantile(data$urban, p = 0.75) # highly urbanised 
urban_ind2 <- data$urban - quantile(data$urban, p = 0.5) # intermediate
urban_ind3 <- data$urban - quantile(data$urban, p = 0.25) # more rural

water_ind1 <- data$water_shortage - quantile(data$water_shortage, p = 0.75) # high frequency shortages
water_ind2 <- data$water_shortage - quantile(data$water_shortage, p = 0.5) # intermediate
water_ind3 <- data$water_shortage - quantile(data$water_shortage, p = 0.25) # low frequency shortages

urban_basis1_pdsi <- basis_pdsi*urban_ind1
urban_basis2_pdsi <- basis_pdsi*urban_ind2
urban_basis3_pdsi <- basis_pdsi*urban_ind3

# multiply the PDSI cross-basis variables by the water shortage linear terms
water_basis1_pdsi <- basis_pdsi*water_ind1
water_basis2_pdsi <- basis_pdsi*water_ind2
water_basis3_pdsi <- basis_pdsi*water_ind3

# assign unique column names to cross-basis matrix for inla() model
# note: not necessary for glm(), gam() or glm.nb() models
colnames(basis_tmin) = paste0("basis_tmin.", colnames(basis_tmin))
colnames(basis_tmax) = paste0("basis_tmax.", colnames(basis_tmax))
colnames(basis_pdsi) = paste0("basis_pdsi.", colnames(basis_pdsi))

colnames(urban_basis1_pdsi) = paste0("urban_basis1_pdsi.", colnames(urban_basis1_pdsi))
colnames(urban_basis2_pdsi) = paste0("urban_basis2_pdsi.", colnames(urban_basis2_pdsi))
colnames(urban_basis3_pdsi) = paste0("urban_basis3_pdsi.", colnames(urban_basis3_pdsi))

colnames(water_basis1_pdsi) = paste0("water_basis1_pdsi.", colnames(water_basis1_pdsi))
colnames(water_basis2_pdsi) = paste0("water_basis2_pdsi.", colnames(water_basis2_pdsi))
colnames(water_basis3_pdsi) = paste0("water_basis3_pdsi.", colnames(water_basis3_pdsi))


data$micro_index <- rep(1:nmicro, ntime)

# create state index
# state length
k <- unique(data$state_code)

for (j in 1:nstate)
{
  data$state_index[data$state_code == k[j]] <- j 
}

# create year index
# set first year (in this case 2001) to 1
data$year_index <- data$year - 2000 

Y  <- data$dengue_cases # response variable
N  <- length(Y) # total number of data points
E  <- data$population/10^5 # model offset so that response is equivalent to an incidence rate per 100,000 people
T1 <- data$month # for random effect to account for annual cycle (seasonality)
T2 <- data$year_index # for random effect to account for inter-annual variability
S1 <- data$micro_index # for microregion spatial random effect
S2 <- data$state_index # for state interaction with month random effect
Vu <- data$urban # include level of urbanisation (% pop living in urban areas) variable along with linear urban interaction
Vw <- data$water_shortage # include frequency of water shortages along with linear water shortage interaction

# create dataframe for model testing
df <- data.frame(Y, E, T1, T2, S1, S2, Vu,Vw)
# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
# models
lab = c("basemodel", "model0.1", "model0.2", "model0.3", "model0.4", "model0.5","model1.1_urban", "model1.2_urban", "model1.3_urban","model1.1_water", "model1.2_water", "model1.3_water")
baseformula <- Y ~ 1 + f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE,scale.model = TRUE,  hyper = precision.prior) +  f(S1, model = "bym2", replicate = T2, graph = g.file, scale.model = TRUE, hyper = precision.prior)
formula0.1 <- update.formula(baseformula, ~. + basis_tmin)
formula0.2 <- update.formula(baseformula, ~. + basis_tmax)
formula0.3 <- update.formula(baseformula, ~. + basis_pdsi)
formula0.4 <- update.formula(baseformula, ~. + basis_tmin + basis_pdsi)
formula0.5 <- update.formula(baseformula, ~. + basis_tmax + basis_pdsi)
baseformula2 <- Y ~ 1 + f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE,scale.model = TRUE,  hyper = precision.prior) + f(S1, model = "bym2", replicate = T2, graph = g.file, scale.model = TRUE, hyper = precision.prior) + basis_tmin + basis_pdsi
formula1.1 <- update.formula(baseformula2, ~. + urban_basis1_pdsi + Vu)
formula1.2 <- update.formula(baseformula2, ~. + urban_basis2_pdsi + Vu)
formula1.3 <- update.formula(baseformula2, ~. + urban_basis3_pdsi + Vu)
formula1.4 <- update.formula(baseformula2, ~. + water_basis1_pdsi + Vw)
formula1.5 <- update.formula(baseformula2, ~. + water_basis2_pdsi + Vw)
formula1.6 <- update.formula(baseformula2, ~. + water_basis3_pdsi + Vw)
formulas <- list(baseformula, formula0.1, formula0.2, formula0.3, formula0.4, formula0.5,formula1.1, formula1.2, formula1.3, formula1.4, formula1.5, formula1.6)
lab <- c("basemodel", "model0.1", "model0.2", "model0.3", "model0.4", "model0.5","model1.1_urban", "model1.2_urban", "model1.3_urban","model1.1_water", "model1.2_water", "model1.3_water")
rerun = TRUE
if(rerun){
	models <- vector(mode = "list",length = 12)
	# Use the base model to generate groups
	print(paste0(lab[1]," started at ", Sys.time()))
        model <- inla(formula = formulas[[1]], 
  				data = df, 
  				family = "nbinomial", 
  				offset = log(E),
       				control.compute = list(dic = TRUE,cpo = TRUE,waic = TRUE, return.marginals = FALSE,control.gcpo=list(enable = T, num.level.sets = 3)),
       				inla.mode = "experimental",
       				control.fixed = list(prec.intercept = 1, prec = 1))
        save(model, file = paste0(script_dir,"/output/", lab[1],".RData"))
        print(paste0(lab[1]," finished at ", Sys.time()))
        models[[1]] <- model
        groups = lapply(1:N,FUN = function(i){models[[1]]$gcpo$groups[[i]]$idx})
	for(i in 2:12){
		 print(paste0(lab[i]," started at ", Sys.time()))
                 model <- inla(formula = formulas[[i]], 
  				data = df, 
  				family = "nbinomial", 
  				offset = log(E),
       				control.compute = list(dic = TRUE,cpo = TRUE,waic = TRUE, return.marginals = FALSE,control.gcpo=list(enable = T, groups = groups)),
       				inla.mode = "experimental",
       				control.fixed = list(prec.intercept = 1, prec = 1))
       		 save(model, file = paste0(script_dir,"/output/", lab[i],".RData"))
                 print(paste0(lab[i]," finished at ", Sys.time()))
                 models[[i]] <- model
	}
}else{
	models = lapply(1:12,FUN = function(i){load(file = paste0(script_dir,"/output/", lab[i],".RData"));return(model)})
}

table <- data.table(Model  = c("base", "tmin", "tmax", "pdsi", "tmin + pdsi", "tmax + pdsi","high urban", "intermediate urban", "low urban",
                                "high water shortage", "intermediate water shortage", "low water shortage"), 
                     DIC = NA,
                     WAIC = NA,
                     LOOCV = NA,
                     LGOCV = NA)




table$DIC = round(unlist(lapply(1:12,FUN = function(i){mean(models[[i]]$dic$dic)})),2)
table$WAIC = round(unlist(lapply(1:12,FUN = function(i){mean(models[[i]]$waic$waic)})),2)
table$LOOCV = round(unlist(lapply(1:12,FUN = function(i){-mean(log(models[[i]]$cpo$cpo))})),4)
table$LGOCV = round(unlist(lapply(1:12,FUN = function(i){-mean(log(models[[i]]$gcpo$gcpo))})),4)

print(paste("The minimum DIC,WAIC, LOOCV and LGOCV are", min(table$DIC),min(table$WAIC),min(table$LOOCV),min(table$LGOCV)))
table$DIC = table$DIC - min(table$DIC)
table$WAIC = table$WAIC - min(table$WAIC)
table$LOOCV = table$LOOCV - min(table$LOOCV)
table$LGOCV = table$LGOCV - min(table$LGOCV)
print(table)


groups = lapply(1:N,FUN = function(i){models[[1]]$gcpo$groups[[i]]$idx})
jan_idx = which(df$T1 == 1)
jan_hist = table(unlist(lapply(jan_idx,FUN = function(i){df$T1[groups[[i]]]})))
print(round(jan_hist/jan_hist[1],5))

april_idx = which(df$T1 == 4)
april_hist = table(unlist(lapply(april_idx,FUN = function(i){df$T1[groups[[i]]]})))
print(round(april_hist/april_hist[4],5))

july_idx = which(df$T1 == 7)
july_hist = table(unlist(lapply(july_idx,FUN = function(i){df$T1[groups[[i]]]})))
print(round(july_hist/july_hist[7],5))

november_idx = which(df$T1 == 11)
november_hist = table(unlist(lapply(november_idx,FUN = function(i){df$T1[groups[[i]]]})))
print(round(november_hist/november_hist[11],5))


#get the script directory
script_dir <- dirname(sys.frame(1)$ofile)
#load libraries
if(!require("INLA",quietly=TRUE)){
	install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
	require("INLA",quietly=TRUE)
}

#plot function and load data

germany.map <- function(data,
                        cutpoints=seq(min(data),max(data),length=256),
                        autoscale=FALSE,
                        legend=FALSE,
                        append=FALSE,
                        main = "")
{
  if (autoscale)
  {
    data = (data-min(data))/(max(data)-min(data)+1e-8)
  }
  #cutpoints = c(-1e9,cutpoints, 1e9)
  
  farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))
  
  xmin <- 1:length(germany)
  xmax <- 1:length(germany)
  ymin <- 1:length(germany)
  ymax <- 1:length(germany)
  
  for(i in 1:length(germany))
  {
    xmin[i] <- min(germany[[i]][,2],na.rm=T)
    xmax[i] <- max(germany[[i]][,2],na.rm=T)
    ymin[i] <- min(germany[[i]][,3],na.rm=T)
    ymax[i] <- max(germany[[i]][,3],na.rm=T)
  }
  
  breite <- c(min(xmin),max(xmax))
  hoehe <- c(min(ymin),max(ymax))
  
  if (TRUE) {
    ## correct the aspect-ratio,  this seems quite ok
    print(paste("aspect.ratio", diff(hoehe)/diff(breite)))
    fac = 0.75
    for(k in 1:length(germany))
      germany[[k]][, 2] = breite[1] + fac*(germany[[k]][, 2]-breite[1])
    ##breite[2] = breite[1] + fac*diff(breite)
  }
  
  breite[2] = breite[2] + 0.1 * diff(breite)
  
  if (!append) plot(breite,hoehe,type="n",axes=T, xlab=" ", ylab=" ",xaxt='n',yaxt='n',main= main,cex.main=4.3)
  
  for(k in length(germany):1)
  {
    polygon(germany[[k]][,2]+1250,germany[[k]][,3],col=farben[k])
  }
  
  if (legend)
  {
    nc = 20
    x1 = 5100 +1250
    x2 = x1 + 300
    y1 = 500
    dy = 350
    
    for(i in 1:nc)
    {
      polygon(c(x1, x1,x2, x2),c(y1+dy*(i-1),y1+dy*i,y1+dy*i,y1+dy*(i-1)),col=gray(i/(nc+1)))
    }
    for(i in 2:nc)
    {
      text(x2+dy/2, y1+dy*(i-1), as.character(round(cutpoints[(i/nc)*length(cutpoints)], 2)),
           cex=.7,col=rgb(0,0,0))
    }
  }
}
data("Germany")

source(paste0(script_dir,"/germany.map"))

SMR = Germany$Y/Germany$E

g = system.file("demodata/germany.graph", package="INLA")

Germany = cbind(Germany,region.struct=Germany$region)

n = 544

m = 10
#fit model
formula = Y ~ 1 + f(x, model="rw2") + f(region.struct,model="besag",graph=g) + f(region,model="iid") 

res = inla(formula,
            family="poisson",
            E = E,
            data=Germany,
            inla.mode = "experimental",
            control.fixed = list(prec.intercept = 1), 
            control.compute = list(control.gcpo = list(enable = TRUE,  num.level.sets = m, size.max=50 ,strategy = "prior",remove.fixed = F),config =T))
            

res_lgocv_prior = inla.group.cv(result = res, num.level.sets = m, strategy = "prior")
res_lgocv_prior_keep = inla.group.cv(result = res, num.level.sets = m, strategy = "prior", keep = "region.struct")
res_lgocv_post = inla.group.cv(result = res, num.level.sets = m)

#plot

id = 50
map1_1 = numeric(544) + 1
map1_2 = numeric(544) + 1
map1_3 = numeric(544) + 1
map1_1[res_lgocv_prior$groups[[id]]$idx] = .5
map1_2[res_lgocv_prior_keep$groups[[id]]$idx] = .5
map1_3[res_lgocv_post$groups[[id]]$idx] = .5
map1_1[id] = 0
map1_2[id] = 0
map1_3[id] = 0

id = 500
map2_1 = numeric(544) + 1
map2_2 = numeric(544) + 1
map2_3 = numeric(544) + 1
map2_1[res_lgocv_prior$groups[[id]]$idx] = .5
map2_2[res_lgocv_prior_keep$groups[[id]]$idx] = .5
map2_3[res_lgocv_post$groups[[id]]$idx] = .5
map2_1[id] = 0
map2_2[id] = 0
map2_3[id] = 0

pdf(file = paste0(script_dir,"/germany_map.pdf"),
    width = 40,
    height = 40)
par(mfrow = c(2,2))
germany.map(map1_2,main = "(a) Groups constructed by prior (spatial only) at location id:50")
germany.map(map2_2,main = "(b) Groups constructed by prior (spatial only) at location id:500")
germany.map(map1_3,main = "(c) Groups constructed by posterior at location id:50")
germany.map(map2_3,main = "(d) Groups constructed by posterior at location id:500")
dev.off()



































<h1> Load R packages </h1>


library(GWmodel)      ### GW models

library(sp)           ## Data management

library(spdep)        ## Spatial autocorrelation

library(RColorBrewer) ## Visualization

library(classInt)     ## Class intervals

library(raster)       ## spatial data

library(grid)         # plot

library(gridExtra)    # Multiple plot

library(ggplot2)      # Multiple plot

library(gtable)

<h1> Load Data </h1>

dataFolder<-"D:\\Dropbox\\Spatial Data Analysis and Processing in R\\Data_GWR\\"

county<-shapefile(paste0(dataFolder,"COUNTY_ATLANTIC.shp"))

state<-shapefile(paste0(dataFolder,"STATE_ATLANTIC.shp"))

mf<-read.csv(paste0(dataFolder,"data_atlantic_1998_2012.csv"), header=T)

<h1>Create a data frame</h1>

df=mf[c(1,4:9)]
head(df)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/c2b256cd-3cef-4a52-b23e-78435da10215)

<h1>Scale co-variates</h1>

df[, 3:7] = scale(df[, 3:7])

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/cccb5ad7-c6a9-4f90-a8ff-549f3ab4adda)


<h1>Merge data with county shape file</h1>

SPDF<-merge(county,df, by="FIPS")
names(SPDF)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/fd9148de-19f7-411e-9ba8-4eb6111a5339)

<h1>Bandwidth selection</h1>
DM<-gw.dist(dp.locat=coordinates(SPDF))

bw.gwr <- bw.ggwr(Rate ~ POV+SMOK+PM25+NO2+SO2,  
                 data = SPDF,
                 family = "poisson",
                 approach = "AICc",
                 kernel = "bisquare", 
                 adaptive = TRUE,
                 dMat = DM )
bw.gwr


<h1>Fit the model</h1>

bgwr.res <- ggwr.basic(Rate ~ POV+SMOK+PM25+NO2+SO2, 
                      data =SPDF,
                      family = "poisson",
                      bw = bw.gwr, 
                      kernel = "bisquare", 
                      adaptive = TRUE,
                      dMat = DM)
bgwr.res

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/20b81bb2-8c16-47e6-858a-e75cda8f202c)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/e7f6b705-1cef-4864-806d-ae7a89186d1b)




Extract GWPR results
### Create spatial data frame

county@data$y<-bgwr.res$SDF$y

county@data$yhat<-bgwr.res$SDF$yhat

county@data$residual<-bgwr.res$SDF$residual

rsd=sd(county@data$residual)

county@data$stdRes<-(county@data$residual)/sd(county@data$residual)

county@data$LLN=county@data$yhat-1.645*rsd

county@data$ULN=county@data$yhat+1.645*rsd


# Intercept
county@data$Intercept<-bgwr.res$SDF$Intercept

county@data$est_SMOK<-bgwr.res$SDF$SMOK

county@data$est_POV<-bgwr.res$SDF$POV

county@data$est_PM25<-bgwr.res$SDF$PM25

county@data$est_NO2<-bgwr.res$SDF$NO2

county@data$est_SO2<-bgwr.res$SDF$SO2

# T-values

county@data$t_Intercept<-bgwr.res$SDF$Intercept_TV

county@data$t_SMOK<-bgwr.res$SDF$SMOK_TV

county@data$t_POV<-bgwr.res$SDF$POV_TV

county@data$t_PM25<-bgwr.res$SDF$PM25_TV

county@data$t_NO2<-bgwr.res$SDF$NO2_TV

county@data$t_SO2<-bgwr.res$SDF$SO2_TV

# Calculate psudo-t values

county@data$p_SMOK<-2*pt(-abs(bgwr.res$SDF$SMOK_TV),df=3103)

county@data$p_POV<-2*pt(-abs(bgwr.res$SDF$POV_TV),df=3103)

county@data$p_PM25<-2*pt(-abs(bgwr.res$SDF$PM25_TV),df=3103)

county@data$p_NO2<-2*pt(-abs(bgwr.res$SDF$NO2_TV),df=3103)

county@data$p_SO2<-2*pt(-abs(bgwr.res$SDF$SO2_TV),df=3103)


county$sig_SMOK <-ifelse(county@data$est_SMOK > 0 &
                          county@data$p_SMOK <= 0.05 , 1, 0)
                          
county$sig_POV <-ifelse(county@data$est_POV > 0 &
                           county@data$p_POV <= 0.05 , 1, 0)
                           
county$sig_PM25 <-ifelse(county@data$est_PM25 > 0 &
                          county@data$p_PM25 <= 0.05 , 1, 0)
                          
county$sig_NO2 <-ifelse(county@data$est_NO2 > 0 &
                           county@data$p_NO2 <= 0.05 , 1, 0)
                           
county$sig_SO2 <-ifelse(county@data$est_SO2 > 0 &
                           county@data$p_SO2 <= 0.05 , 1, 0)

<h1>Plot GWRP Statistics</h1>

polys<- list("sp.lines", as(state, "SpatialLines"), col="grey", lwd=.8,lty=1)
  
col.palette<-colorRampPalette(c("blue",  "sky blue", "green","yellow", "red"),space="rgb",interpolate = "linear")

<h1>Plot Local Estimates</h1>

col.palette<-colorRampPalette(c("lightcyan","cyan","cyan1", "cyan2","cyan3","cyan4", "darkblue"),space="rgb",interpolate = "linear") 

est_smok<-spplot(county,"est_SMOK", main = "Smoking", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=col.palette(100))

est_pov<-spplot(county,"est_POV", main = "Poverty", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=col.palette(100))

est_pm25<-spplot(county,"est_PM25", main = "PM25", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=col.palette(100))

est_no2<-spplot(county,"est_NO2", main = "NO2", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=col.palette(100))

est_so2<-spplot(county,"est_SO2", main = "SO2", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=col.palette(100))
       
grid.arrange(est_smok, est_pov,est_pm25,est_no2, est_so2,ncol= 5, heights = c(30,6), top = textGrob("Local Estimates",gp=gpar(fontsize=25)))

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/206836b9-caa8-48e9-ad72-c978d24e5bf1)

<h1>Plot Local t-values</h1>
col.palette.t<-colorRampPalette(c("blue",  "sky blue", "green","yellow","pink", "red"),space="rgb",interpolate = "linear") 

t_smok<-spplot(county,"t_SMOK", main = "Smoking", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=rev(col.palette.t(100)))

t_pov<-spplot(county,"t_POV", main = "Poverty", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=rev(col.palette.t(100)))

t_pm25<-spplot(county,"t_PM25", main = "PM25", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=rev(col.palette.t(100)))

t_no2<-spplot(county,"t_NO2", main = "NO2", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=rev(col.palette.t(100)))

t_so2<-spplot(county,"t_SO2", main = "SO2", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=rev(col.palette.t(100)))
grid.arrange(t_smok, t_pov,t_pm25,t_no2, t_so2,ncol=5, heights = c(30,6), top = textGrob("Local t-values",gp=gpar(fontsize=25)))

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/9bbccde8-0aa9-49dd-bce5-c337ed38bb08)

<h1>Plot Std-Residuals</h1>

myPaletteRes <- colorRampPalette(c("lightseagreen","lightsteelblue1", "moccasin","hotpink", "red"))

std_res<-spplot(county,"stdRes", main = "GWRP Std. Residuals", 
       sp.layout=list(polys),
       col="transparent",
       col.regions=myPaletteRes(100))
       
#windows(width=4, height=3.5)

#tiff( file="FIG_GWRP_Std_Residuals.tif", 

#width=4, height=3.5,units = "in", pointsize = 12, res=1600,

#restoreConsole = T,bg="transparent")

print(std_res)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/ef7ef18a-991b-4a06-bd78-385caba6a4df)


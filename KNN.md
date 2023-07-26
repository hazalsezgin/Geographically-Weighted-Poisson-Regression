<h1>Spatial Regression Analysis, K-nearest neighbors</h1>
  
1) K-neighbors listesi oluşturma
   
Bu, her ilçe için  bir konum sağlar ve bu konumlara dayalı olarak mesafe değerleri hesaplayabiliriz. 

İlçe çokgenlerinden xy verilerini oluşturmak için basitçe geosphere paketindeki centroid() fonksiyonunu kullanabiliriz.

all.xy <-centroid(tx_sp)

colnames(all.xy) <- c("x","y")

2) Sonraki adım olarak, merkez noktalarımıza dayalı olarak k-en yakın değeri kullanarak bir komşu listesi oluşturmamız gerekiyor. 
Bu örnekte k = 1, k = 3 ve k = 5 değerlerini inceleyeceğiz.
Ardından, modelin komşuları içerecek bir yarıçap oluşturabilmesi için mesafe değerini hesaplamamız gerekiyor. 
Son olarak, komşuluk içindeki komşuların listesini oluşturmamız gerekiyor.

#Create neighbors
all.dist.k1 <- knn2nb(knearneigh(all.xy, k=1, longlat = TRUE))

all.dist.k3 <- knn2nb(knearneigh(all.xy, k=3, longlat = TRUE))

all.dist.k5 <- knn2nb(knearneigh(all.xy, k=5, longlat = TRUE))


#Determine max k distance value to neighbor

all.max.k1 <- max(unlist(nbdists(all.dist.k1, all.xy, longlat=TRUE)))

all.max.k3 <- max(unlist(nbdists(all.dist.k3, all.xy, longlat=TRUE)))

all.max.k5 <- max(unlist(nbdists(all.dist.k5, all.xy, longlat=TRUE)))

#Calculate neighbors based on distance

all.sp.dist.k1 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k1, longlat = TRUE)

all.sp.dist.k3 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k3, longlat = TRUE)

all.sp.dist.k5 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k5, longlat = TRUE)

#Create neighbor list

all.dist.neighb.k1 <- nb2listw(all.sp.dist.k1,style="W", zero.policy = TRUE)

all.dist.neighb.k3 <- nb2listw(all.sp.dist.k3,style="W", zero.policy = TRUE)

all.dist.neighb.k5 <- nb2listw(all.sp.dist.k5,style="W", zero.policy = TRUE)


Bu, bir (1), üç (3) ve beş (5) komşuya dayalı olarak yürütülecek mekansal gecikme (spatial lag) ve mekansal hata (spatial error) uzaklık modelleri için gereken mekansal matrisleri sağlar.
Gecikme ve hata modelleri, yukarıda tamamlanan modellere benzerdir, ancak mekansal matris için komşuluk temellendirmesi yerine mesafe hesaplamasını kullanır.
Ayrıca, yukarıda Moran Korelasyonu ve LaGrange Çarpanı Testlerini zaten tamamladığımız için, bu adımı bu çalışma için tekrarlamamıza gerek yoktur ve analize doğrudan geçebiliriz.

Distance Lag Model

Her k-mesafe değeri için bir mesafe gecikme modeli hesaplamak için aşağıdakileri kullanacağız:

all.dist.lag.k1 <- spatialreg::lagsarlm(equation, data = tx_pov, listw = all.dist.neighb.k1)

all.dist.lag.k3 <- spatialreg::lagsarlm(equation, data = tx_pov, listw = all.dist.neighb.k3)

all.dist.lag.k5 <- spatialreg::lagsarlm(equation, data = tx_pov, listw = all.dist.neighb.k5)

Bu örnek için, yalnızca K=1 gecikme modelinin özetini görüntüleyeceğiz:

summary(all.dist.lag.k1, Nagelkerke = TRUE)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/b9399eac-073a-491f-9247-6a1c3448639f)


![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/7fd2e7d4-5f36-4617-91e1-72cf69065fa1)


Distance Error Model

Her k-mesafe değeri için bir mesafe hata modeli hesaplamak için aşağıdakileri kullanacağız:

all.dist.err.k1 <- spatialreg::errorsarlm(equation, data = tx_pov, listw = all.dist.neighb.k1)

all.dist.err.k3 <- spatialreg::errorsarlm(equation, data = tx_pov, listw = all.dist.neighb.k3)

all.dist.err.k5 <- spatialreg::errorsarlm(equation, data = tx_pov, listw = all.dist.neighb.k5)

Bu örnek için, yalnızca K=1 hata modelinin özetini görüntüleyeceğiz:

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/79ab4df5-6066-40c2-9835-27c352cdd487)

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/df42619c-c345-41a1-bd32-2ec73df9cc0f)

 Comparing distance models, model selection

 ![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/5300bcf5-cd14-42e3-a828-71f8e95c44d6)

 Mapping the results

 Yoksulluk verilerini mekansal verilere bağlamak için, iki veri setini birleştirmek için ortak bir sütuna ihtiyacımız var. Bu örnek için FIPS kodlarını kullanacağız.4
 
 dist.err.data <- summary(all.dist.err.k1, correlation=TRUE, Nagelkerke = TRUE)

dist.err.output <- cbind.data.frame(tx_pov$FIPS,
                                    dist.err.data$fitted.values, 
                                    dist.err.data$residual, 
                                    tx_pov$child.pov.2016, 
                                    tx_pov$lnsinglemom, 
                                    tx_pov$lnuninsured, 
                                    tx_pov$lnlesshs, 
                                    tx_pov$lnincome_ratio,
                                    stringsAsFactors = FALSE)

#Renaming columns
colnames(dist.err.output) <- c("fips","fitted","resid","childpov",
                        "single_mom","uninsured","less_hs","income_ratio")

tx_fortify <- fortify(tx_sp)

tx_poly <- merge(x = tx_fortify, y = dist.err.output, 
                 by.x = "id", by.y = "fips", all = TRUE)

bivariate_data <- bi_class(tx_poly, x = childpov, y = single_mom, 
                           dim = 3, style = "quantile")

legend <- bi_legend(pal = "DkViolet",
                    dim = 3,
                    xlab = "Child Poverty",
                    ylab = "Single Mother\n Households",
                    size = 6)


world <- map_data("world")
states <- map_data("state")
southern_states <- subset(states, region %in% 
                            c("texas", "arkansas", "louisiana", "mississippi", 
                              "alabama", "georgia", "florida", "north carolina",
                              "south carolina", "tennessee", "oklahoma", 
                              "kentucky", "west virginia", "virginia", 
                              "maryland", "delaware", "district of columbia"))


                              


mom_pov_map <- ggplot() + 
  geom_polygon(data = world, aes(x=long,y=lat, group=group), fill = "gray95", color = "white") +
  geom_polygon(data = states, aes(x=long,y=lat, group=group), fill = "gray", color = "white") +
  geom_polygon(data = southern_states, aes(x=long,y=lat, group=group), fill = NA, size = 0.01, color = "white") +  
  geom_polygon(data = bivariate_data, aes(x=long, y=lat, group=group, fill = bi_class), color = "grey50", show.legend = FALSE) + 
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  coord_map("conic", lat0 = 30, xlim=c(-107,-92), ylim=c(25,37)) +
  theme_void() + theme(legend.title.align=0.5) +
  theme(panel.background = element_rect(fill = 'deepskyblue'),
        panel.grid.major = element_line(colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "Child Poverty", 
       title = "Bivariate Map of Child Poverty and Single Mother Households") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
mom_pov_map

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/60a48b9a-408a-4d16-a57c-890bc29f1acd)


final_map <- ggdraw() +
  draw_plot(mom_pov_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.60, 0.035, 0.25, 0.25)
final_map

![image](https://github.com/hazalsezgin/Geographically-Weighted-Poisson-Regression/assets/77546910/37c0abf3-edeb-4e7c-8430-d8a46ceb94fa)



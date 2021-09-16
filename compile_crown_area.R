######## Calculate crown areas from trees with different heights ########
####### corwns can be overlapping, contained within other crowns, etc. ##
###### AS April 2020, based on NEON data  ##############################

library(neonUtilities)
library(geoNEON)
library(sf)
library(rgeos)
library(viridis)
library(tidyverse)
library(gdata)

### Start with downloading data from NEON (see https://www.neonscience.org/neonDataStackR)
zipsByProduct(dpID="DP1.10098.001", site=c("ABBY","BART"),  
              startdate="2017-01", enddate="2018-11",
              package="basic", check.size=T)
stackByTable(filepath="./filesToStack10098/", savepath = "./",saveUnzippedFiles=T)

### Load data
maptag <- read.csv("./stackedFiles/vst_mappingandtagging.csv") ### location and plant ID
indiv <- read.csv("./stackedFiles/vst_apparentindividual.csv") 

maptag <- maptag %>% mutate(year=sapply(strsplit(as.character(maptag$date), "-"),"[",1))
maptag <- maptag %>%  mutate(plotyear= paste(maptag$plotID, maptag$year, sep="_"))

### Optional: Select desired site(s) and year(s)
patt <- "ABBY.*2017"
look <- maptag[grepl(patt, maptag$plotyear), ]

### OR: Run everything (can take a while depending on the no of site/years)
# look <- maptag

### Calculate tree coordinates
tree_coords <- getLocTOS(look, "vst_mappingandtagging")

### Merge with growth measurements
tree_coords_combi <- tree_coords %>%
  left_join(indiv[,c(9,12:27)],by = "individualID")%>%
  mutate_if(is.factor, as.character)


### Replace missing crown90Diameter with maxDiameter
tree_coords_combi[which(is.na(tree_coords_combi$ninetyCrownDiameter)),]$ninetyCrownDiameter <- tree_coords_combi[which(is.na(tree_coords_combi$ninetyCrownDiameter)),]$maxCrownDiameter
tree_coords_combi <- tree_coords_combi %>% mutate(individualID = paste(tree_coords_combi$taxonID,
                                  sapply(strsplit(as.character(tree_coords_combi$individualID), "\\.",),"[",5), sep = "_"))

### Reduce input (for easier overview later)
tree_combi2 <- tree_coords_combi %>% 
  select("plotID", "taxonID", "height", "adjNorthing","adjEasting", 
         "maxCrownDiameter", "ninetyCrownDiameter", "individualID")


### Spatial polygons for plots
plotsize <- read.csv("./stackedFiles/vst_perplotperyear.csv")
plotsize$year <- sapply(strsplit(as.character(plotsize$date), "-"),"[",1)
plotsize$plotyear <- paste(plotsize$plotID, plotsize$year, sep="_")

centroids <- plotsize %>% select(c(5,6,13,14,7,28,26,29,36,37)) %>% unique() %>%
  mutate_if(is.factor, as.character)

### Optional: Select site and year if needed 
centroids_sel <- centroids[grepl(patt, centroids$plotyear),]

### OR: run all
# centroids_sel <- centroids

### Spatial polygons for sites
radius <- sqrt(centroids_sel$totalSampledAreaTrees)/2 # radius (m) for plots
yPlus <- centroids_sel$northing+radius
xPlus <- centroids_sel$easting+radius
yMinus <- centroids_sel$northing-radius
xMinus <- centroids_sel$easting-radius

# Calculate polygon coordinates for each plot centroid. 
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus,yMinus,  # SE corner
             xMinus,yMinus, # SW corner
             xMinus,yPlus)  # NW corner again - close ploygon
# Extract the plot ID information
ID <- sort(as.character(centroids_sel$plotID))

polys <- SpatialPolygons(mapply(function(poly, id){
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
}, 
split(square, row(square)), ID))
polydf <- SpatialPolygonsDataFrame(Sr = polys, data = centroids_sel,match.ID = F)


### Start loop ###
### Calculate crown area per plot (ID) 

### Test: Start with one plot
k=1

# for(k in 1:length(ID)){  ### for loop
  keep(tree_combi2,ID,k,polydf,sure=T)
  # errfile <- file(paste0("./error_file_crowns_per_plot_",ID[k], ".txt")) ### useful for loop
  # tryCatch({
  seli <- tree_combi2[tree_combi2$plotID==ID[k]&!(is.na(tree_combi2$adjNorthing)),]
  seli <- seli[order(seli$height,decreasing = F),] ### order sel by height
  seli <- seli[!(is.na(seli$maxCrownDiameter)),]
  seli <- seli[!(is.na(seli$height)),]
  
  if(!(anyDuplicated(seli[,c("adjNorthing","adjEasting")])==0)){
  seli <- seli[-which(duplicated(seli[,c("adjNorthing","adjEasting")])),]
  }
  
  ### Spatial object
  sel <- seli
  sel$individualID <- as.character(sel$individualID)
  sel <- sel[order(sel$height,decreasing = T),]
  
  if(!(anyDuplicated(sel$individualID)==0)){
    sel <- sel[-which(duplicated(sel$individualID)),] 
  }
  
  coordinates(sel) <- c("adjEasting", "adjNorthing")### X before y
  
  elli_a <- ((sel@data$maxCrownDiameter)/2)*((sel@data$ninetyCrownDiameter)/2)*pi
  elli_r <- sqrt(elli_a/pi)
  
  buff <- gBuffer(sel,byid = T,width = elli_r,
                  id = sel@data$individualID, quadsegs=100) ### buffer by ellipse radius
  
  # plot(buff)
  buff@plotOrder <- order(buff@data$height)
  
  buffx <- buff
  
  remx <- gContainsProperly(buff,byid = T,returnDense = F)
  remx[sapply(remx, is.null)] <- NULL 
  
  #### Plot
  plot(buff, col=magma(nrow(sel),direction = -1, alpha=0.7,begin = 0.15))
  polygonsLabel(buff, buff$individualID,method = "centroid", cex=0.5)
  
  if(is_empty(remx)==F){
    ### check for crowns fully contained within other crowns
    try(for(i in 1:20){ ### run this a couple of times
      remx <- gContainsProperly(buff,byid = T,returnDense = F)
      remx[sapply(remx, is.null)] <- NULL 
      
      if(is_empty(remx)==F){
      rem_dat <- as.data.frame(unlist(remx))
      rem_dat$big <- row.names(rem_dat)
      names(rem_dat)[1] <- "small_ID"
      ssx <- substr(sapply(strsplit(rem_dat$big, "_"),"[",2),1,5)
      rem_dat$big <- paste(sapply(strsplit(rem_dat$big, "_"),"[",1), ssx, sep="_")
      
      for(i in 1:nrow(rem_dat)){
        rem_dat$small[i] <-  buff@data$individualID[rem_dat$small_ID[i]]
        rem_dat$height_small[i] <-buff@data$height[rem_dat$small_ID[i]]
        rem_dat$height_big[i] <- buff@data$height[which(buff$individualID==rem_dat$big[i])]  
      }
      
      rem_dat <- rem_dat[,c(2,5,3,4,1)] 
      rem <- rem_dat[rem_dat$height_big>rem_dat$height_small,] ### remove trees with smaller diameter that are less tall
      
      #### remove smaller plots
      buff <- subset(buff, !(buff$individualID %in% rem$small))
      }
    }, silent = T)
  }
  
  # plot(subset(buff,buff$individualID %in%rem$small), add=T, col=2)
  
  buffi <- buff
  
  plot(buffi,col=rev(rainbow(n=length(buffi),alpha=0.5)))
  polygonsLabel(buffi, buffi$individualID,method = "centroid", cex=0.5,doPlot = T)
  
  buffi_sf <- st_as_sf(buffi)
  plot(buffi_sf$geometry, col=rainbow(n=length(buff),alpha=0.5))
  
  ### Set precision (some won't run): 1000 == 1/1000 m (mm level)
  combi <- buffi_sf %>% st_set_precision(1000) %>% st_intersection 
 
 ### Or: Snap features to themselves
 # xx <- st_snap(buffi_sf, buffi_sf, tolerance = 5) 
 # combi <- xx %>% st_intersection
 
  plot(combi$geometry,col=rainbow(n=nrow(combi)))
  
  ### add origins for all tree segements
  unis <- buffi_sf %>% mutate(n.overlaps=1, origins=1:nrow(buffi_sf))  
  uu <- unis[!(unis$origins %in% combi$origins),]
  
  if(nrow(uu)>0){
    combi <- rbind(combi,uu)
  }
  
 ### Calculate area
  combi$area <- st_area(combi)
  if(nrow(uu)>0){
  combi[combi$origins %in% uu$origins,]$area <- 0
  }
  
  ### Add info about height as classes for tree selection
  combi2 <- as.data.frame(combi)
  combi2$group <- seq(1:nrow(combi2))

  for (i in 1: nrow(combi2)){
    combi2$groupx[i] <- paste(combi2[i,]$origins[[1]], collapse = "|")
  }

  combi2$heightmax <- numeric(nrow(combi2))
  for (i in 1:nrow(combi2)){
    if(combi2[i,]$n.overlaps==1){
      combi2[i,]$heightmax <- combi2[i,]$height
    } else{
      aa <- buffi[unlist(combi2[i,]$origins),]@data   
      # pp <- paste(paste0("^",unlist(strsplit(combi2[i,]$groupx,"\\|"))), collapse="|")
      #   selx <- combi2[grepl(pp, combi2$groupx) & combi$n.overlaps==1,] 
        combi2[i,]$heightmax <- max(aa$height, na.rm=T)
    }
  }
  
  combi2$indsel <- character(nrow(combi))
  combi2$IND <- character(nrow(combi))
  
  for(i in 1:nrow(combi2)){
    aa <- buffi[unlist(combi2[i,]$origins),]@data #### plotnames
    combi2[i,]$indsel <- paste(aa$individualID,collapse = ",")
    if(combi2[i,]$n.overlaps>1){
      pp <- paste(paste0("^",unlist(strsplit(combi2[i,]$groupx,"\\|")),"$"), collapse="|")
      patt <- paste(aa$individualID,collapse = "|")
      
      nnh <- combi2[grepl(pp, combi2$groupx),] ### for regular overlaps
      groupi <- nnh[nnh$heightmax==max(nnh$heightmax),"group"]
      
      # ppp <- gsub("\\^", "", pp) ### for overlaps within crowns
      # ppp <- gsub("\\$", "", ppp)
      # ppp <- paste0("^", ppp, "$")
      # nnh2 <- combi2[grepl(ppp, combi2$groupx),]
      # nnh <- bind_rows(nnh2,nnh)
      # groupi <- nnh[nnh$heightmax==max(nnh$heightmax),"group"]
      # aa$heightmax <- max(aa$height)  ### add if clause for crown contained in crowns
  
      if(length(groupi)>1){
        combi2$IND[i] <- patt
      }else{
        combi2[i,]$group <- groupi
        
        indi <- nnh[nnh$heightmax==max(nnh$heightmax),grep(pattern = "individualID",colnames(combi2))]
        combi2$IND[i] <- indi[!(is.na(indi))]
        
        # 
        # indi <- aa[aa$height==max(aa$height),grep(pattern = "individualID",colnames(combi2))]
        # combi2$IND[i] <- indi[!(is.na(indi))] 
        
      }
    }else{
      indi <-  as.character(combi2[i,grep(pattern = "individualID",colnames(combi2))])
      combi2$IND[i] <- indi[!(is.na(indi))]
    }
  }

### And another round for crowns countained in crowns  
outx <- combi2
bo <- st_as_sf(polydf[polydf$plotID==ID[k],])

# pdf(paste0("./crowns_",ID[k], ".pdf")) ### save if needed
plot(bo$geometry, border=4)
plot(outx$geometry, col=rainbow(n=length(outx),alpha=1), add=T)
plot(bo$geometry, border=4, add=T)
## polygonsLabel(buffi, buffi$individualID,method = "centroid", cex=0.5,doPlot = T)
# dev.off()

outxx <- st_as_sf(outx)

# plot(outxx$geometry)
st_agr(outxx) <- "constant"
st_agr(bo)  <- "constant"

outxxx <- st_intersection(outxx,bo)
plot(outxxx$geometry, col=rainbow(n=length(outx),alpha=0.5))
polygonsLabel(buffi, buffi$individualID,method = "centroid", cex=0.5,doPlot = T)

outxxx$area <- st_area(outxxx)

if(nrow(uu)>0){
  outxxx[outxxx$origins %in% uu$origins,]$area <- 0
}


### Add single tree info for trees containing trees
combi2 <- outxxx 

#### same height
if(exists("rem_dat")){
  sam <- rem_dat[rem_dat$height_big-rem_dat$height_small==0,]
  ss <- c(paste(sam$big, sam$small, sep=","),paste(sam$small, sam$big, sep=","))
  
  combi2[combi2$indsel %in% ss,]$IND <- gsub(",","\\|",combi2[combi2$indsel %in% ss,]$indsel)
  
  ### smaller tree is taller
  tall <- rem_dat[rem_dat$height_big-rem_dat$height_small<0,]
  tt <- c(paste(tall$big, tall$small, sep=","),paste(tall$small, tall$big, sep=","))
  
  xx <- combi2[combi2$indsel %in% tt,]
  
  for(i in 1:nrow(xx)){
    ll <- unlist(strsplit(xx$indsel[i],","))
    xx$IND[i] <- tall[tall$small%in%ll,]$small
  }
  
  # if(length(tt)>0 & any(combi2$indsel %in% tt)){
  #   combi2 <- combi2[-(which(combi2$indsel %in% tt)),]  
  #   combi2 <- rbind(combi2,xx) 
  # }
  if(length(tt)>0 & any(combi2$indsel %in% tt)){
    if(any(combi2$indsel %in%tt)){
      combi2 <- combi2[-(which(combi2$indsel %in% tt)),]
    }
  }
  
  if(nrow(xx)>0){
    for(u in 1:nrow(xx)){
      if(!(xx$IND[u] %in% combi2$IND))
        combi2 <- rbind(combi2,xx)
    }
  }
}

combi2 <- as.data.frame(as_tibble(combi2))
out <- combi2 %>% group_by(IND)%>% summarize(area_total=sum(area))
combi2$origins <- as.character(combi2$origins)
combi2 <- subset(combi2 , select=-c(geometry))

### Save crown segments if needed (usually not)
# write.csv(combi2, paste0("./" ,ID[k],"_crown_segments.csv"),row.names = F) 

### Overlaps with same with same height
cc <- out[grepl("\\|", out$IND),]$IND 

if(length(cc)>0){
  for (i in 1:length(cc)){ 
    namx <- unlist(strsplit(cc[i],"\\|"))
    # safe <- selx
    selx <- bind_rows(combi2[combi2$IND ==cc[i],],combi2[combi2$IND %in% namx,]) ### plots to work with
    # selx <- as_tibble(selx)
    oo <- selx[grepl("\\|", selx$IND),] ### split plots
    keep <- selx[!(grepl("\\|", selx$IND)),] ### decided plots to keep 
    
    ### three or more overlpas with same height: remove id's of trees that are not at max height from split
    nammax <- unique(selx[(selx$height %in% max(selx$height,na.rm=T)) &selx$n.overlaps==1,"IND"])
    
    #### modify output
    area_part <- sum(oo[,"area"]/length(nammax))
    
    for (j in 1:length(nammax)){
      out[out$IND==nammax[j],"area_total"] <- out[out$IND==nammax[j],"area_total"]+ area_part
    }
  }
}

out_fin <- out[!(grepl("\\|", out$IND)),]

### Save final result
write.csv(out_fin, paste0("./" ,ID[k],"_crowns.csv"),row.names = F)

# }, error=function(e){writeLines(as.character(e),errfile)})
# close(errfile)
#   }




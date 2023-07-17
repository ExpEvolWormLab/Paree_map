#####################################
### Ploting Against Rockman's map ###
#####################################
library('ggplot2')
library(data.table)
# load rockman maps and put relative genetic and physical distance

source("utils.R")
recmaps = as.data.frame(fread("recmaps.csv"))
cemap = as.data.frame(fread("rockman_recombinationmaps.csv"))

recmaps = subset(recmaps, intercross == "consensus" & rec=="wt")



### Maps centered on rockman maps

cemap = do.call(rbind, lapply(split(cemap, cemap$chrom), function(x){ x$relgenetic = x$genetic/max(x$genetic); x}))

mapsrel = do.call(rbind, lapply(split(recmaps, paste0(recmaps$chrom, recmaps$cross)), function(tar){
  
  #tar = split(recmaps, paste0(recmaps$chrom, recmaps$cross))[[5]]
  chr = tar$chrom[1]
  tar = tar[order(tar$pos),]
  
  
  #relative
  rockman_genetic = approx(subset(cemap, chrom==chr)$cpos, subset(cemap, chrom==chr)$relgenetic, xout = c(min(tar$pos), max(tar$pos)), rule=2)$y
  explength = diff(rockman_genetic)
  a = explength/diff(c(min(tar$genetic), max(tar$genetic)))
  relgen = tar$genetic*a
  b = rockman_genetic[2] - max(relgen)
  tar$relgenetic= relgen + b
  
  #plot(tar$pos, tar$relgenetic)
  
  #centered
  centerpos = suppressWarnings(approx(subset(cemap, chrom==chr)$relgenetic, subset(cemap, chrom==chr)$cpos, xout = 0.5, rule=2)$y)
  centergen = approx(tar$pos, tar$relgenetic, xout = centerpos, rule=2)$y
  tar$centergenetic=tar$relgenetic + (0.5 - centergen)
  
  # If chromosome V, not marker in the center. do not center it and just keep rel genetic
  if(chr == 5){
    tar$centergenetic=tar$relgenetic
  }
  
  #plot(tar$pos, tar$centergenetic)
    
  
  tar
  
}))


# Figure S6
ggplot(mapsrel, aes(pos/1e6, centergenetic))+
  geom_point(data=cemap, aes(cpos/1e6, relgenetic), color='grey', size=1.5)+
  geom_point(size=1.5)+facet_grid(numtorom(chrom)~.)+
  facet_wrap(~numtorom(chrom), ncol = 2)+theme_perso()+
  ylab("Relative genetic position") + xlab("Physical position (Mb)")+
  theme(strip.text.y.right = element_text(angle = 0))+
  theme(axis.text = element_text(size=13, color="black"),
        axis.title = element_text(size=18),
        strip.text = element_text(size=15),
        legend.position = "none")




# Correlation
unlist(lapply(split(mapsrel, paste0(mapsrel$cross, mapsrel$chrom)), function(x){
  
  chr = x$chrom[1]
  
  plot(x$pos, x$genetic)
  
  #Inetrpolate Rockman's genetic position for the physical positon of our SNV
  rockmangenetic = approx(subset(cemap, chrom==chr)$cpos, subset(cemap, chrom==chr)$relgenetic, xout = x$pos, rule=2)$y
  
  x = cor(rockmangenetic, x$genetic)
  
  names(x)=numtorom(chr)
  
  x
  
}))


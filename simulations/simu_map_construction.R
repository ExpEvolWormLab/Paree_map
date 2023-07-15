path = './Data/Linkage_map/clean/'
simpath = paste0(path, "simulation/")

library(qtl)
library(data.table)


#load(paste0(path, "/rqtl_mapthis_clean.Rdata"))
#allmarkers = unlist(lapply(mapthis[["geno"]], function(chr){names(chr[["map"]])}))
#rm(mapthis)
#fileid = "_riails_simu_drift_1to250.Rdata"
#crosses = LETTERS[1:8]

#sims.allcross = list()
#for(cross in crosses){
#  load(paste0(simpath, cross, fileid))
#  sims.allcross[[which(crosses==cross)]] = riails.simulation[sub]
#}

#names(sims.allcross) = crosses

######################################
######### FUNCTIONS ##################
######################################

expand.missing.markers = function(lines, allmarkers){
  lines = lines[match(allmarkers, rownames(lines)),]
  rownames(lines) = allmarkers
  lines
}



matrix.to.cross.qtl = function(geno.matrix){
  # Get phenotype which is just the id of the lines here
  pheno = colnames(geno.matrix)
  pheno = data.frame(pheno) #df format is important
  #pheno[,1] = as.factor(pheno[,1]) 
  
  # Pass for 0:1 allele format to 1:2
  geno.matrix[geno.matrix==1]=2
  geno.matrix[geno.matrix==0]=1
  
  # Get chromosome vector which come from rownames here
  chrom = data.table::tstrsplit(rownames(geno.matrix), '_')[[1]]
  
  #split by chromosomes
  geno.list = split(as.data.frame(geno.matrix), chrom)
  
  #Get the geno object at the R/qtl format
  geno = lapply(1:length(geno.list), function(gchr){
    gchr = geno.list[[gchr]]
    data = t(gchr)
    rownames(data)=NULL
    map = as.numeric(data.table::tstrsplit(rownames(gchr), '_')[[2]])
    names(map)=rownames(gchr)
    gchr = list(data=data, map=map)
    class(gchr)="A"
    gchr
  })
  
  names(geno) = 1:length(geno.list)
  class(geno[[1]])='A'
  
  # Get the geno and pheno object is one riself cross object at the good format
  mp=list()
  mp[["geno"]]=geno
  mp[["pheno"]]= pheno
  names(mp[["pheno"]])="NA."
  class(mp) = c("riself", "cross")
  
  mp
}


Do.map = function(mp){
  nt.bymar <- ntyped(mp, "mar")
  todrop <- unique(c(names(nt.bymar[nt.bymar < 1])))
  mp <- drop.markers(mp, todrop)
  
  map <- est.map(mp, map.function="haldane", maxit=50000)
  map = unlist(map)
  markers = tstrsplit( names(map), ".", fixed = TRUE)[[2]]
  map = data.frame(marker=markers, genetic = map)
  map$chrom =as.numeric(tstrsplit(map$marker, '_')[[1]])
  map$pos =as.numeric(tstrsplit(map$marker, '_')[[2]])
  for(i in unique(map$chrom)){map$genetic[map$chrom==i] = map$genetic[map$chrom==i] - min(map$genetic[map$chrom==i])}
  map
}


#########################################
##### All maps for one simulation #######
#########################################


Do.maps.simu = function(nsim, sims.allcross, allmarkers,cross.map = c('AB', 'CD', 'EF', 'GH', 'ABCD', 'EFGH')){
  crosses =  names(sims.allcross)
  riails.sim = dplyr::bind_cols(lapply(crosses, function(cross){
    xx = sims.allcross[[cross]]
    xx = xx[[nsim]]
    xx = expand.missing.markers(xx, allmarkers = allmarkers)
    colnames(xx) = paste0(cross, '_', colnames(xx))
    xx
  }))
  
  riails.sim = matrix.to.cross.qtl(riails.sim)
  
  simmaps=NULL
  for(cross in cross.map){
    # Subset the qtl object for this cross
    xcross =  data.table::tstrsplit(riails.sim[["pheno"]][["NA."]], '_')[[1]] %in% unlist(strsplit(cross, ""))
    mpcross = subset(riails.sim, ind = which(xcross))
    #Construct the map
    segdis = suppressWarnings(geno.table(mpcross))
    segdis = segdis[(segdis$AA + segdis$BB) > 0, ]
    simmap = Do.map(mpcross)
    simmap$ref.count = segdis$AA
    simmap$all.count = segdis$AA + segdis$BB
    simmap$cross=cross
    simmap$nsimulation = nsim
    simmaps = rbind(simmaps, simmap)
  }
  
  simmaps
}


#test = Do.maps.simu(1, sims.allcross, allmarkers)

#ggplot(subset(simmaps, cross %in% c('ABCD', "EFGH")), aes(pos/1e6, genetic))+geom_point()+facet_grid(chrom~cross)


######################################
#### Do Map for simulations ##########
######################################

X.BD = T # only keeping B and D cross

type="sel" #drift or sel or flat

load(paste0(path, "/rqtl_mapthis_clean.Rdata"))
allmarkers = unlist(lapply(mapthis[["geno"]], function(chr){names(chr[["map"]])}))


ff = list.files(paste0(simpath, type, "_simulation/"))

crosses = LETTERS[1:8]
cross.map = c('ABCD', 'EFGH')

if(type == 'flat'){
  crosses = LETTERS[5:8] # mutant only
  cross.map = c('EFGH')} 

if(X.BD){
  crosses = c("B","D") # mutant only
  cross.map = c('BD')
  allmarkers = names(mapthis[["geno"]][[6]][["map"]])
}

rm(mapthis)


start = seq(1,1000, 125)
end = (seq(1,1000, 125)+124)

#length(start)
for(nn in 1:length(start)){
  print(nn)
  print("Begin at:")
  print(Sys.time())
  #Select the file id for this batch (125 per 125)
  simix = start[nn]:end[nn]
  #simix=1:10
  fftar = ff[grepl(paste0("_",min(simix), "to"),ff) | grepl(paste0("to",max(simix), "."),ff)]
  #fftar = ff[grepl("1to10",ff)]
  fileid = paste0("_riails", tstrsplit(fftar, '_riails')[[2]][1])
  
  #Load the files
  sims.allcross = list()
  for(cross in crosses){
    load(paste0(simpath,type, "_simulation/", cross, fileid))
    sims.allcross[[which(crosses==cross)]] = riails.simulation[as.character(simix)]
  }
  
  names(sims.allcross) = crosses
  
  if(X.BD){ 
    
    for(i in 1:length(sims.allcross)){
      sims.allcross[[i]] = lapply(sims.allcross[[i]], function(x){x[which(rownames(x) %in% allmarkers),]})
    } 
    
    } 
  
  # Simulate maps
  
  simulated.maps = dplyr::bind_rows(parallel::mclapply(simix, mc.cores = 4, function(nsim){
    Do.maps.simu(nsim=as.character(nsim), sims.allcross=sims.allcross, allmarkers=allmarkers, cross.map=cross.map)
  }))
  
  simulated.maps$type = type
  
  #save the data
  
  if(!X.BD){
    xfile = paste0(simpath, "simulated_maps/","simu_maps_",type ,sprintf("_%1.0fto%1.0f", min(simix), max(simix)), ".Rdata")
  }else{
    xfile = paste0(simpath, "simulated_maps/","X_BD_simu_maps_",type ,sprintf("_%1.0fto%1.0f", min(simix), max(simix)), ".Rdata")
  }
  
  save(simulated.maps, file = xfile)
}






#ggplot(subset(simulated.maps, nsimulation==3), aes(pos,genetic))+geom_point()+facet_wrap(~chrom)


################################################
#### Combine X BD with simulated maps ##########
################################################

spath = paste0(path, "simulation/simulated_maps/")

ff_X_BD = list.files(spath)
ff_X_BD = ff_X_BD[grepl("X_BD_simu_maps",ff_X_BD)]

ffall = list.files(paste0(spath, "without_X_correction/"))

ids = data.table::tstrsplit(ff_X_BD, "maps_")[[2]]
ids = data.table::tstrsplit(ids, ".Rdata")[[1]]

for(id in ids){
  
  fX = ff_X_BD[grepl(id,ff_X_BD)]
  fall = ffall[grepl(id,ffall)]
  
  load(paste0(spath, fX))
  simulated.maps.X = simulated.maps
  load(paste0(spath,"without_X_correction/", fall))
  
  simulated.maps.X$cross = 'ABCD'
  
  simulated.maps = rbind(subset(simulated.maps, !(chrom==6 & cross=="ABCD")),
                         simulated.maps.X)
  
  save( simulated.maps, file = paste0(spath,fall))
  
}








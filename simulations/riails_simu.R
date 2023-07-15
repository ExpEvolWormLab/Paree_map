##########################################
######### WITH SELECTION #################
##########################################


path = './Paree_map/'
simpath = paste0(path, "simulations/")

######################################
######### PARAMETERS #################
######################################

# The functions below need some parameters in the environment
# Parameters are uploaded below
# Some parameter vary between cross

## Different parameters are:
# ngen.intercross ; ngen.selfing : the n generation for the intercross and selfing phase
# nplates: the number of starting lines
# nextinction: number of lines lost during riails construction => lines extinct randomly during simulation
# failed.seq: number of lines failed during sequencing => supressed after simulation

# chromosomes: vector of chromosomes c(1,1,1, ..., 2,2,2, ..., 3, 3,..., 6)
# num.loci: the total number of markers
# recombination.rate: vector of rr (size = num.loci) from the consensus genetic maps
# => rr serve as weight for the placement of the one CO per chrom

# sex.chromosomes: vector for sex.chrom (non recombining in male XO)
# sex.locus: the locus of the sex.chrom used for sex determination
# male.allele: the allele of the sex.locus which determine male

# markers: the name of the markers



######################################
######### FUNCTIONS ##################
######################################

# Simulation function upload
source(paste0("path","simulation_functions.R"))

# Need to modify some function :

###########################################################
######### Force sex ratio in the progeny ##################

# First: we need precise number of males and females in progeny
# So we need to force the sex of the progeny
# 1) Force one haplotype to be selected (The one containing the male locus)
make.gamete.disto <- function(sim,id, distortion.locus, distortion.alleles, select=T) {
  indiv = select.ind(sim,id)
  recomb.vector = sample.breakpoints(sim,id)
  recomb.mat = breakpoints.tomatrix(recomb.vector)
  gamete=indiv$genotype[recomb.mat]
  #Inverse segregation if needed matrix to choose one haplotype
  present=!(gamete[distortion.locus] %in% distortion.alleles)
  invert=ifelse(select, !(gamete[distortion.locus] %in% distortion.alleles), 
                (gamete[distortion.locus] %in% distortion.alleles))
  if(invert) recomb.mat[,2] = ifelse(recomb.mat[,2] == 1, 2, 1)
  gamete=indiv$genotype[recomb.mat]
  gamete
}

# Use the distorted make.gamete function in make.offspring
# select = T => force in male; select = F => force in female

#make.offspring.forcesex <- function(sim,id.mother, id.father, sex.locus, male.allele, select) {
#  # Makes an individual out of two parents. 
#  genotype <- cbind(make.gamete(sim,id.mother), make.gamete.disto(sim,id.father, sex.locus, male.allele, select))
#  sex = sample(c(select.ind(sim,id.mother)$sex, select.ind(sim,id.father)$sex),1)
#  list(
#    genotype  = genotype, 
#    sex = sex,
#    fitness   = 1
#  )
#}


## SMALL MODIFICATION TO INCORPORATE SELECTION from make.offspring.forcesex function 

make.offspring.viability <- function(sim,id.mother, id.father, sex.locus,forcesex=F,male.allele=NULL, select=NULL, id.subpop) {
  
  #If force sex = T : force to select (select=T) or reject(select=F) the Y chromosome (male.allele, sex.locus)
  alive = 0
  # We create an offspring, which have some chance to die in function of fitness
  # If die, recreate an offspring
  # If possibility of offspring fitness is really low => long so not to hard selection
  
  while(alive==0){
    # Makes an individual out of two parents. 
    
    if(forcesex==T){ 
      gam.father = make.gamete.disto(sim,id.father, sex.locus, male.allele, select)
    } else {
      gam.father = make.gamete(sim,id.father) }
    
    genotype <- cbind(g1=make.gamete(sim,id.mother), g2=gam.father)
    
    sex = sample(c(select.ind(sim,id.mother)$sex, select.ind(sim,id.father)$sex),1)
    
    s.vector = select.selcoef(sim,id=c(which(names(sim)==id.subpop),0))
    viability.fitness=geno.to.fitness(genotype, s.vector)
    
    #print(viability.fitness)
    
    # if viability.fitness is >= 1 => always viable
    if(runif(1) < viability.fitness){alive = 1}else{alive = 0} #Chance to die 
  }
  
  list(
    genotype  = genotype, 
    sex = sex,
    fitness   = 1 #we keep this fitness at 1 = no competition between parents
  )
}





#########################################################################
### C.elegans recombination: One obligate CO + supression on XO #########

sample.breakpoints = function(sim,id, chrom=chromosomes, sex.chrom=sex.chromosome){
  rr.vec = select.rr(sim, id)
  sex=select.ind(sim,id)$sex
  breakpoints = unlist(lapply(unique(chrom), function(chr){
    chr.pos = which(chrom==chr)
    chr.pos = chr.pos[-length(chr.pos)]
    bk = sample(chr.pos, 1, prob=rr.vec[chr.pos]) #sample obligate CO inside the chrom
    bk = c(bk, max(chr.pos)+1) # recombination 
    bk #Vector of 100% recombination intra and inter chromosoomes
  }))
  
  #50% chance to be non recombinant (100cM => 50cM maps)
  breakpoints = breakpoints[sample(c(T,F), length(breakpoints), replace=T)] 
    
  breakpoints = breakpoints[breakpoints<length(rr.vec)+1] #for the last loci, need to supress the breakpoint if it arose
  if(sex=="male") breakpoints = breakpoints[!(breakpoints %in% sex.chrom[-length(sex.chrom)])] #supress recombination in male
  breakpoints=c(0,breakpoints, length(rr.vec)+1 ) # c(0, breakpoints, num.loci) is the format for the following function
  breakpoints
}



####################################################################
###### Modified reproduction call for previous function ############

### Intercross phase
intercross = function(sim, sex.locus=NA, male.allele=NA, sex.progeny){
  
  repro.subpop = function(sim,id.subpop, sex.locus, male.allele, sex){
    parents = chooseparents.dioecious(sim,id.subpop)
    offspring=make.offspring.viability(sim, parents[[1]], parents[[2]], forcesex = T,
                                       sex.locus = sex.locus, male.allele = male.allele, 
                                       select=ifelse(sex=="male", T, F),id.subpop=id.subpop)
    offspring = update.sex(offspring, sex.locus = sex.locus, male.allele = male.allele)
    offspring
  }
  
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  for(subpop in subpops.ids){
    #print(subpop)
    newsubpop = lapply(sex.progeny, function(sx){repro.subpop(sim, subpop, sex.locus=sex.locus, male.allele=male.allele,sex=sx)})
    names(newsubpop) = paste0("i", 1:length(sex.progeny))
    sim[[subpop]] = newsubpop
  }
  
  return(sim)
}



### Obigate selfing phase
selfing = function(sim){
  
  repro.subpop = function(sim,id.subpop){
    nsubpop=which(names(sim)==id.subpop)
    parent = c(nsubpop,1) # Pure selfing of one individual each time
    make.offspring.viability(sim, parent, parent, forcesex=F, id.subpop = id.subpop)
    
  }
  
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  for(subpop in subpops.ids){
    newsubpop = list(repro.subpop(sim, subpop)) #only one individual
    names(newsubpop) = paste0("i", 1)
    sim[[subpop]] = newsubpop
  }
  
  return(sim)
}



######################################################
###### Random intercross migration scheme ############

#Random "pair" mating with equal contributions (RPMEC), [Rockman & Kruglyak, 2008]
# Except here = not pair as there is two males
# Keep the same number of subpopulation by default but possibility for extinction
# The female (hermaphrodite) randomly migrate
# The two males randomly migrate in the same

RPMEC.migration.scheme = function(sim, extinction = 0){
  migration.scheme = empty.migration.scheme(sim)
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  sexes = unlist(lapply(1:length(subpops.ids), function(x){
    subpop = subpops.ids[x]
    sapply(sim[[subpop]], "[[", "sex")
  }))
  
  persist = !(migration.scheme$original.pop %in% sample(unique(migration.scheme$original.pop), extinction))
  migration.scheme = migration.scheme[persist,]
  sexes = sexes[persist]
  subpop = unique(migration.scheme$original.pop)
  nsubpop = length(unique(migration.scheme$original.pop))
  mig.males = data.frame(original.pop=subpop, newpop=sample(1:nsubpop,nsubpop))
  mig.females = data.frame(original.pop=subpop, newpop=sample(1:nsubpop, nsubpop))
  migration.scheme = rbind(merge(migration.scheme[sexes=="male",1:2],mig.males), merge(migration.scheme[sexes=="female",1:2],mig.females))
  migration.scheme = migration.scheme[order(migration.scheme$original.pop, migration.scheme$ind),]
  migration.scheme
}

#################################################################################
### Keep only female (hermaphrodite) migraryion scheme for selfing transition ###

keep.female.scheme = function(sim){
  migration.scheme = empty.migration.scheme(sim)
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  sexes = unlist(lapply(1:length(subpops.ids), function(x){
    subpop = subpops.ids[x]
    sapply(sim[[subpop]], "[[", "sex")
  }))
  
  migration.scheme = migration.scheme[sexes=="female",]
  migration.scheme$newpop = migration.scheme$original.pop #keep metapopulation structire
  migration.scheme
}

################################
### Some other function used ###

### Randomly extinct n subpop (used during selfing phase) to mimick the loss of lineages
subpop.extinction = function(sim, nextinction){
  ext = sample(1:(length(sim) - 3), nextinction)
  sim = sim[-ext]
  #important for some function to have subpop going from p1=>pn without missing
  names(sim)[1:(length(sim) - 3)] = paste0("p", 1:(length(sim) - 3))
  sim
}

### return the haplotypes of the individual of the simulation
### homozygous 1:1 = 1 ; homozygous 0:0 = 0 ; heterozygous = NA 
return.haplotypes = function(sim){
  onepop = empty.migration.scheme(sim)
  onepop$newpop = 1
  sim = migration(sim, onepop, verbose=F)
  haplotypes = dplyr::bind_cols(lapply(sim[["p1"]], function(x){apply(x$genotype, 1, sum)}))
  haplotypes[haplotypes==1] = NA  #heterozygous = NA => treated like this in R/qtl
  haplotypes[haplotypes==2] = 1
  return(haplotypes)
}


### Sequencing noise 
sequencing.noise=function(haplotype.matrix, missing.matrix){haplotype.matrix*missing.matrix}

# haplotype matrix have some rare NA value to rare heterozygous loci
# If they do not overlap with NA value in the missing.matrix,
# resulting matrix from sequencing.noise will have slighly more NA than missing.matrix
# heterozygous loci are very rare, so I assume it should be ok
# below, some draft for a function to switch columns so NA values in both matrix overlap


### Sample NA values in the haplotypes matrix = unsequenced markers
### Need the depth
### Rem: we sample per loci 
### Maybe a beter way would be to use a "Missing marker matrix" to mimick the same distribution of missing values in the riails

#sequencing.noise = function(haplotypes, depth){
#  alreadymissing = apply(haplotypes, 1, function(x) sum(is.na(x)))
#  nhaplotypes = ncol(haplotypes)
#  #number of loci to turn in NA (unsequenced markers) to get same depth
#  nunsequenced = nhaplotypes - depth - alreadymissing
#  haplotypes=dplyr::bind_rows(lapply(1:nrow(haplotypes), function(l){
#    geno.l=haplotypes[l,]
#    miss = sample((1:nhaplotypes)[!is.na(geno.l)], nunsequenced)
#    geno.l[miss]=NA
#    geno.l
#  }))
#  
#  haplotypes
#}


# Some draft for a function that would change the order of haplotype.matrix columns
# so NA values in missing.matrix and haplotype.matrix would overlap
# => thus, the final number of NA values is the same as in missing.matrix and not a bit over 
#haplotype.matrix = riails.sim()

#x = which(is.na(haplotype.matrix))
#cols = ceiling(x/nrow(haplotype.matrix))
#rows = x - ((cols -1)*nrow(haplotype.matrix))

###each element n of the list = which(is.na()) for col n
#missrow = apply(missing.matrix, 2, function(x){which(is.na(x))})

###return the column of the missing matrix which satisfy the "NA structure" of the haplotype
###If haplotype.matrix[,1] have NA value at rows 3,5,6 :
### => we want to have corresponding missing.matrix colum which have NA at at least rows 3,5,6
#coresp = lapply(unique(cols), function(cc){
#  rr = rows[which( cols ==cc)]
#  which(unlist(lapply(missrow, function(mm){ sum(rr %in% mm) })) == length(rr))
#})
#names(coresp)=unique(cols)


##############################################
######### RANDOM SELECTED LOCO CHOICE ########


# Function designed for negative viability selection
# minfitness = the min fitness possible if all selected loci are homozygous 1
sample.sel.coef = function(n.sel.loci.per.chr, minfitness, chromosomes){
  
  s = -1*((1 - minfitness)/ (2*n.sel.loci.per.chr*length(unique(chromosomes))))
  
  selected.loci = unlist(lapply(split(chromosomes, chromosomes), function(chr){
    
    unwhich <- function (which, dim = max(which)) {
      y <- array(logical(length(which)), dim = dim)
      y[which] <- TRUE
      y}
    
    sel.loci = sample(1:length(chr), n.sel.loci.per.chr)
    unwhich(sel.loci, dim = length(chr))
  }))
  
  sel.coef = rep(0, length(selected.loci))
  sel.coef[selected.loci] = s
  sel.coef
}


######################################
######### SIMULATION FUNCTION ########


riails.sim.sel = function(s){
  # First create a F2 population of 3 individuals containng all the genotypes possible
  geno.matrix = cbind(rep(1, num.loci), rep(0, num.loci))[,rep(c(1,2),3)] # Three fully heterozygous individual
  # Add the male chromsome = 2 => = "null" in C.elegans 
  geno.matrix[sex.chromosome, c(3,6)]= male.allele # one female (herm) and Two males differing only for the X chrom
  
  sim = init.simul(list(geno.matrix), list(c("female", "male", "male")))
  sim = set.recombinationrate(sim, recombination.rate)
  sim = set.mutationrate(sim, 0)
  sim = set.selection(sim, s)
  
  
  # Sample extinction events
  extinction.events = sample(1:(ngen.intercross + ngen.selfing), nextinction, replace=T) 
  #Put in form of a vector of length = ngeneraton and each element = nextinction event for this generation
  extinction.events = table(extinction.events)
  extinction.events = cbind(generation = as.numeric(names(extinction.events)), nextinction = extinction.events)
  nonextinct = (1:(ngen.intercross + ngen.selfing))[!( 1:(ngen.intercross + ngen.selfing) %in% extinction.events[,1] )]
  nonextinct = cbind(generation = nonextinct, nextinction = rep(0, length(nonextinct)))
  extinction.events = rbind(extinction.events,nonextinct)
  extinction.events = extinction.events[order(extinction.events[,1]),2]
  names(extinction.events) = 1:(ngen.intercross + ngen.selfing)
  
  # Reproduction F2 => F3
  # nplates = 55
  # Needed progeny = 55 females and 55*2 males => force sex
  sim = intercross(sim, sex.locus=sex.locus, male.allele=male.allele, sex.progeny = c(rep("male", nplates*2), rep("female", nplates)))
  
  # We then dispacth one herm and two males in 55 "metapopulations" = 55 plates
  dispatch.F3 = empty.migration.scheme(sim)
  sexes = sapply(sim[["p1"]], "[[", "sex")
  dispatch.F3$newpop[sexes == "female"] = sample(1:nplates, nplates)
  dispatch.F3$newpop[sexes == "male"] = sample(rep(1:nplates, 2), 2*nplates)
  
  sim = migration(sim, dispatch.F3, verbose=F)
  
  # Random intercross until F5
  for(gg in 1:ngen.intercross){
    sim = intercross(sim, sex.locus=sex.locus, male.allele=male.allele, sex.progeny = c("male", "male", "female"))
    random.mig = RPMEC.migration.scheme(sim, extinction = extinction.events[gg])
    sim = migration(sim, random.mig, verbose=F) 
  }
  
  
  # Selfing until F13
  
  herm  = keep.female.scheme(sim)
  sim = migration(sim, herm, verbose=F) #Keep only female (=hermpaphrodite)
  for(gg in ((1:ngen.selfing) + ngen.intercross)){
    sim=selfing(sim)
    #Extinction 
    if(extinction.events[gg]>0){sim=subpop.extinction(sim, extinction.events[gg])}
  }
  
  return.haplotypes(sim)
}




##############################################################
################## SIMULATIONS ###############################
##############################################################

allcrosses = LETTERS[1:8]
total.fitness.loss = 0.9 #0 if drift, >0 if negative selection
flat = F

if(flat){
  allcrosses = LETTERS[5:8] # mutant only
  load("./C_elegans_ressources/domain_boundaries.RData")
  endchrom=aggregate(pos~chrom, boundaries, max)
  } 

nsim = 1000
batch = 4

for(cross in allcrosses){
  # LOAD PARAMETERS
  print(cross)
  load(paste0(simpath, 'parameters/', cross, "_parameters_simu.Rdata"))
  
  if(flat){
    ppos = data.table::tstrsplit(markers, '_')[[2]]
    ppos[c(1,which(diff(chromosomes)==1) + 1)]=0
    ppos[c(which(diff(chromosomes)==1), length(chromosomes))]=endchrom$pos
    ppos = as.numeric(ppos)
    # totally flat landscape 
    # recombination probability = interval physical size
    recombination.rate = diff(ppos)
    recombination.rate[recombination.rate<0] = 0 #between chrom, number doesn't matter, not taken into account
    }
  
  for(n in 1:batch){
    print("Begin at:")
    print(Sys.time())
    
    ix = ((n-1)*(nsim/batch) +1):(n * (nsim/batch))
    
    riails.simulation=parallel::mclapply(ix, mc.cores = 4,function(ii){
    #riails.simulation=lapply(ix,function(ii){
      ## Simulation 
      sel.coef=sample.sel.coef(1, 1-total.fitness.loss, chromosomes)
      riails = riails.sim.sel(s=sel.coef)
      
      ## delete the markers at extemities added for the simulation:
      truemarkers = !(grepl("start",markers) |  grepl("end",markers))
      sel.coef = sel.coef[truemarkers]
      riails = riails[truemarkers,]
      truemarkers = markers[truemarkers]
      
      ## Mimick lines lost at sequencing step:
      
      if(failed.seq>0){riails = riails[, -sample(ncol(riails), failed.seq)]}
      
      ## Add missing markers in sequence data
      riails=sequencing.noise(riails, missing.matrix)
      rownames(riails)=truemarkers
      riails
      
      #af = apply(riails, 1, sum, na.rm=T)/apply(riails, 1, function(x) sum(!is.na(x)))
      #af = data.frame(markers = rownames(riails), freq = af)
      #af$pos = as.numeric( data.table::tstrsplit(af$markers, '_')[[2]] )
      #af$chrom = as.numeric( data.table::tstrsplit(af$markers, '_')[[1]] )
      #af$s = sel.coef
      #mean(af$freq)
      #ggplot(af, aes(pos/1e6, freq, color=as.factor(sel.coef)))+geom_point()+facet_wrap(.~chrom)+geom_hline(yintercept = 0.5, color="red")+ylim(0,1)
    })
    
    names(riails.simulation)=ix
    
    save(riails.simulation, file = paste0(simpath,"sel_simulation/", cross, sprintf("_riails_simu_sel_%1.0fto%1.0f", min(ix), max(ix)), ".Rdata"))
    print("End at:")
    print(Sys.time())
  }
}


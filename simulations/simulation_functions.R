###########################################################
####### Initialize the simulation #########################
###########################################################

# Create one individual
init.individual <- function(genotype, sex=NA) {
  # Generates a random individual for the starting population
  list(
    genotype  = genotype, 
    sex = sex,
    fitness   = 1
  )
}

# Create randomly a genotype matrix for a (sub)population
# Not necessary if we have a known genotype matrix
# genotype matrix = loci * haplotype
# diploid individuals= ind1 = col 1 & 2; ind2 = col 3 & 4;  etc.

init.genotype.matrix = function(pop.size, num.loci, alt.freq){
  #pop.size = the population size
  #num.loci = the number of loci
  #alt.freq = the alternative allele frequency (biallelic)
  #if alt.freq is a number => allele.freq for all loci
  #if alt.freq is a vector of size == num.loci, each loci have a different initial allele frequency
  
  if(length(alt.freq)==1) alt.freq = rep(alt.freq, num.loci)
  
  geno.matrix = NULL
  for(l in 1:num.loci){
    freq = alt.freq[l]
    alleles=rbinom(pop.size*2,1,freq)
    geno.matrix = rbind(geno.matrix, alleles)
  }
  
  geno.matrix = as.matrix(geno.matrix)
  colnames(geno.matrix) = paste0(rep("ind", pop.size), rep(1:pop.size, each=2), 
                                 rep('_hap', pop.size), rep(c(1,2),pop.size))
  
  rownames(geno.matrix) = 1:num.loci
  return(geno.matrix)
}

#geno.matrix=init.genotype.matrix(4,10,0.2)


# Transform a geno.matrix at the subpopulation format
init.subpopulation <- function(geno.matrix, sex.vector=NA) {
  # Generates the initial population
  # geno.matrix = list if the inital haplotype
  #sex.vector = vector of sex of size = pop.size, or =NA
  pop.size = ncol(geno.matrix)/2
  num.loci = nrow(geno.matrix)
  if(length(sex.vector)==1) sex.vector = rep(sex.vector, pop.size)
  
  pop = lapply(1:pop.size, function(i) {
    ind = init.individual(genotype=geno.matrix[,c(i*2 -1, i*2)], 
                          sex = sex.vector[i])
    ind
  })
  
  names(pop) = paste0(rep("i",pop.size), 1:pop.size)
  
  return(pop)
}

# Initialize the simulation
# => several population possible (metapopulations)
init.simul = function(list.geno.matrix, list.sex.vector){
  lgm = length(list.geno.matrix)
  lsv = length(list.sex.vector)
  if(lgm != lsv) print("Error: list.geno.matrix and list.sex.vector must have same length")
  sim = list()
  for(m in 1:lgm){
    sim[[m]] = init.subpopulation(list.geno.matrix[[m]], list.sex.vector[[m]])
  }
  
  names(sim) = paste0(rep("p", lgm), 1:lgm)
  
  num.loci = nrow(sim[["p1"]][["i1"]]$genotype)
  
  param = list(rr = matrix(NA, ncol=lgm, nrow=num.loci-1),
               sel.coef = matrix(NA, ncol=lgm, nrow=num.loci),
               mut.rate = matrix(NA, ncol=lgm, nrow=num.loci))
  
  sim = c(sim, param)
  
  return(sim)
}


#Take simu + new rr and return simu with new rr
# if rr is a number, same rr everywhere
#if rr is a vector of length num.loci -1, each element is the rr between the loci and the next one
# If a matrix of dim: nrow = num.loci & ncol = n subpop, different recombination rate for each subpop
set.recombinationrate = function(sim, recombinationrates){
  n = nrow(sim[["rr"]])
  nn = length(sim[["rr"]])
  if(length(recombinationrates) == 1){
    sim[["rr"]] = matrix(recombinationrates, ncol=ncol(sim[["rr"]]), nrow=nrow(sim[["rr"]]))}else{
      
      if(length(recombinationrates) == n){
        sim[["rr"]] = as.matrix(matrix(recombinationrates, ncol=1)[,rep(1,length(sim)-3)]) }
      
      if(length(recombinationrates) == nn & nn!=nn){sim[["rr"]] = recombinationrates}
    }
  return(sim)
}



set.mutationrate = function(sim, mutationrates){
  n = nrow(sim[["mut.rate"]])
  nn = length(sim[["mut.rate"]])
  
  if(length(mutationrates) == 1){
    sim[["mut.rate"]] = matrix(mutationrates, ncol=ncol(sim[["mut.rate"]]), nrow=nrow(sim[["mut.rate"]]))}else{
      
      if(length(mutationrates) == n){
        sim[["mut.rate"]] = as.matrix(matrix(mutationrates, ncol=1)[,rep(1,length(sim)-3)]) }
      
      if(length(mutationrates) == nn & nn!=nn){sim[["mut.rate"]] = mutationrates}
    }
  return(sim)
}


set.selection = function(sim, selectioncoef){
  n = nrow(sim[["sel.coef"]])
  nn = length(sim[["sel.rate"]])
  
  if(length(selectioncoef) == 1){
    sim[["sel.coef"]] = matrix(selectioncoef, ncol=ncol(sim[["sel.coef"]]), nrow=nrow(sim[["sel.coef"]]))}else{
      
      if(length(selectioncoef) == n){
        sim[["sel.coef"]] = as.matrix(matrix(selectioncoef, ncol=1)[,rep(1,length(sim)-3)]) }
      
      if(length(selectioncoef) == nn & nn!=nn){sim[["sel.coef"]] = selectioncoef}
    }
  return(sim)
}


#Simple additive fitness without dominance
# if negative => 0 = lethal
# Warning: Here does no count allele different from 0 and 1
geno.to.fitness = function(genotype, s.vector){
  genotype[!(genotype %in% c(0,1))]=0
  genotype=apply(genotype, 1, sum)
  fitness = genotype*s.vector
  fitness=1+sum(fitness)
  if(fitness<0) fitness=0
  fitness
}

# This way of updating fitness is quite inefficent and should be used on a newly created sim
# If possible, geno.to.fitness function should be added during the reproduction call (make.offspring function)
# However, this can be useful to udpate it like this
# Notably if there is migration event and we want to update fitness after migration
update.sim.fitness = function(sim){
  allids = sim.ids(sim)
  
  for(rid in 1:nrow(allids)){
    id = unlist(allids[rid,])
    indiv = select.ind(sim, id)
    genotype = indiv$genotype
    s.vector = select.selcoef(sim,id)
    fitness=geno.to.fitness(genotype, s.vector)
    sim[[id[1]]][[id[2]]]$fitness = fitness
  }
  
  sim
}

# If we are explicit about sex, this will be needed
# sex.locus is the row of the genotype matrix that contain the sex determining locus
update.sex=function(indiv, sex.locus=NA, male.allele=NA){
  indiv$sex=ifelse(male.allele %in% indiv$genotype[sex.locus,], "male", "female") 
  indiv
}


#sim = init.simul(list(geno.matrix), list(c("male", "male", "female", "female")))
#sim = set.recombinationrate(sim, 0.1)
#sim = set.mutationrate(sim, 0) #Mutation not implemented in function yet
#sim = set.selection(sim, 0)  #selection not implemented in function yet


###########################################################
####### Generation ########################################
###########################################################

#id is the identification of an individual = c(p,i); p being the n of the subpop and i being the n of the ind
#id = c(1,1) is the first individual of the first subpopulation
select.ind = function(sim,id){sim[[paste0("p", id[1])]][[paste0("i", id[2])]]} #return one specific indicvidual
select.rr = function(sim,id){sim[["rr"]][,id[1]]} #return rr vector from one indiviual
select.mutrate = function(sim,id){sim[["mut.rate"]][,id[1]]} # ...
select.selcoef = function(sim,id){sim[["sel.coef"]][,id[1]]}

sim.ids = function(sim){
  nsubpop=sum(!names(sim) %in% c("rr", "sel.coef", "mut.rate"))
  dplyr::bind_rows(lapply(1:nsubpop, function(pop){
    nind=length(sim[[pop]])
    data.frame(p=pop, i = 1:nind)}))
}


##################################
## MAKE GAMETE: Recombination ####

#function to sample brekpoints of recombination knowing id
# brekpoints = c(0, x, x, final.locus)
sample.breakpoints = function(sim,id){
  rr.vec = select.rr(sim, id)
  breakpoints = rbinom(length(rr.vec),1, rr.vec)
  breakpoints=c(0,which(breakpoints==1), length(breakpoints)+1 )
  breakpoints
}

### => This is the function you want to modifiy if you want to play with recombination

## Example of modification: supression of recimbination on the male sex chromosomes

#sex.chromosome = 1:9
#sample.breakpoints = function(sim,id){
#  rr.vec = select.rr(sim, id)
#  sex=select.ind(sim,id)$sex
#  breakpoints = rbinom(length(rr.vec),1, rr.vec)
#  if(sex=="male") breakpoints[sex.chromosome]=0
#  breakpoints=c(0,which(breakpoints==1), length(breakpoints)+1 )
#  breakpoints
#}

## Example of modification: one obligate CO per chromosome + supression on male sex chromosome 

#chromosomes = c(rep(1,4), rep(2,6))
#sex.chromosome = 1:4

#sample.breakpoints = function(sim,id){
#  rr.vec = select.rr(sim, id)
#  sex=select.ind(sim,id)$sex
#  breakpoints = unlist(sapply(unique(chromosomes), function(chr){
#    chr.pos = which(chromosomes==chr)
#    chr.pos = chr.pos[-length(chr.pos)]
#    bk = sample(chr.pos, 1, prob=rr.vec[chr.pos])
#    if(sample(c(T,F), 1)) bk = c(bk, max(chr.pos)+1)
#    bk
#  }))
#  breakpoints = breakpoints[breakpoints<length(rr.vec)+1]
#  if(sex=="male") breakpoints = breakpoints[!(breakpoints %in% sex.chromosome[-length(sex.chromosome)])]
#  breakpoints=c(0,breakpoints, length(rr.vec)+1 )
#  breakpoints
#}


#Function to turn breakpoints into recombination matrix
# ex: cbind(1:10, c(rep(1,5), rep(2,5)))
# => Take the five first locus of haplotype 1 and the 5 last of haplotype 2 to create gamete = recomb event in the middle
breakpoints.tomatrix = function(breakpoints){
  breakpoints = diff(breakpoints)
  nbk = length(breakpoints)
  break.haplo=rep(unlist(sample(list(c(1,2), c(2,1)), 1)), ceiling(nbk/2 ))[1:nbk]
  recomb.matrix = unlist(lapply(1:nbk, function(i){rep(break.haplo[i],breakpoints[i])}))
  recomb.matrix = cbind(1:length(recomb.matrix),recomb.matrix)
  return(recomb.matrix)
}


# Function to create gamete from one individual id
make.gamete <- function(sim,id) {
  indiv = select.ind(sim,id)
  recomb.vector = sample.breakpoints(sim,id)
  recomb.mat = breakpoints.tomatrix(recomb.vector)
  indiv$genotype[recomb.mat]
}



##We can modify the function to force segregation distortion
##If so, we also need to modify make.offspring function below
#make.gamete.disto <- function(sim,id, distortion.locus, distortion.alleles, select=T) {
#  indiv = select.ind(sim,id)
#  recomb.vector = sample.breakpoints(sim,id)
#  recomb.mat = breakpoints.tomatrix(recomb.vector)
#  gamete=indiv$genotype[recomb.mat]
#  #Inverse segregation if needed matrix to choose one haplotype
#  present=!(gamete[distortion.locus] %in% distortion.alleles)
#  invert=ifelse(select, !(gamete[distortion.locus] %in% distortion.alleles), 
#                (gamete[distortion.locus] %in% distortion.alleles))
#  if(invert) recomb.mat[,2] = ifelse(recomb.mat[,2] == 1, 2, 1)
#  gamete=indiv$genotype[recomb.mat]
#  gamete
#}


######################################
## MAKE offspring: choose parents ####

# Here sex is not determined by an explicit locus and in sampled from parent sexes
# But the sex can be updated later
make.offspring <- function(sim,id.mother, id.father) {
  # Makes an individual out of two parents. 
  genotype <- cbind(make.gamete(sim,id.mother), make.gamete(sim,id.father))
  sex = sample(c(select.ind(sim,id.mother)$sex, select.ind(sim,id.father)$sex),1)
  list(
    genotype  = genotype, 
    sex = sex,
    fitness   = 1
  )
}

#make.offspring(sim, c(1,1), c(1,2))

## If we have implemented segregation distortion in make.gamete, we also need to modify this function
## In this example, the segregation distortion was implemented to force the sex (distortion of the Y chromosome)
#make.offspring.forcesex <- function(sim,id.mother, id.father, sex.locus, sex.allele, select) {
#  # Makes an individual out of two parents. 
#  genotype <- cbind(make.gamete(sim,id.mother), make.gamete.disto(sim,id.father, sex.locus, sex.allele, select))
#  sex = sample(c(select.ind(sim,id.mother)$sex, select.ind(sim,id.father)$sex),1)
#  list(
#    genotype  = genotype, 
#    sex = sex,
#    fitness   = 1
#  )
#}


# make. offspring function need the id of parents
# Choose parents:

#Hermpaphroditic population with possibility for selfing
chooseparents.herm <- function(sim,id.subpop, selfing) {
  # Returns the next generation
  population = sim[[id.subpop]]
  nind = length(population)
  nsubpop=which(names(sim)==id.subpop)
  fitnesses <- sapply(population, "[[", "fitness")
  mother=unlist(sample(1:nind, 1, prob=fitnesses), recursive=FALSE)
  if(selfing<runif(1,0,1)){
    father=unlist(sample(1:nind, 1, prob=fitnesses), recursive=FALSE)
  }else{father=mother}
  list(c(nsubpop, mother),c(nsubpop, father))
}

#Dioecious population: male mate with female
# This function could be modify to androdioecious quite easily
chooseparents.dioecious <- function(sim,id.subpop, selfing=0) {
  # Returns the next generation
  population = sim[[id.subpop]]
  sexes  <- sapply(population, "[[", "sex")
  nind = length(population)
  nsubpop=which(names(sim)==id.subpop)
  fitnesses <- sapply(population, "[[", "fitness")
  males = which(sexes=="male")
  females = which(sexes=="female")
  if(length(females)==1){mother = (1:nind)[females]}else{
    mother=sample((1:nind)[females], 1, prob=fitnesses[females])}
  
  if(length(males)==1){father = (1:nind)[males]}else{
    father=sample((1:nind)[males], 1, prob=fitnesses[males])}
  
  list(c(nsubpop, mother),c(nsubpop, father))
}



# Basic reproduction call which allow to choose to mating system
# Call on the simulation and return the nex generation
# It call all the above function expect the initializing ones

reproduction = function(sim, subpop.size, dioecious=F, selfing=NA, explicit.sex.determination=F, sex.locus=NA, male.allele=NA){
  
  if(dioecious){chooseparents = chooseparents.dioecious }else{chooseparents = chooseparents.herm}
  repro.subpop = function(sim,id.subpop,selfing,explicit.sex.determination=explicit.sex.determination,sex.locus=sex.locus, male.allele=male.allele){
    parents = chooseparents(sim,id.subpop,selfing)
    offspring=make.offspring(sim, parents[[1]], parents[[2]])
    if(explicit.sex.determination) update.sex(offspring, sex.locus = sex.locus, male.allele = male.allele)
    offspring
  }
  
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  for(subpop in subpops.ids){
    newsubpop = replicate(subpop.size,repro.subpop(sim, subpop, selfing), simplify = F)
    names(newsubpop) = paste0("i", 1:subpop.size)
    sim[[subpop]] = newsubpop
  }
  
  return(sim)
}


###########################################
####### Migration #########################
###########################################

# After reporduction we would want to do migration but ne necessay
# Here a very flexible migration function
# In migration scheme (see below), the new population for each individual of the population is specified
# if it the number of subpopulation change, by default the parameter are set as in previous "p1" for all the subpop
# You can then change it with set.xxx function if needed
# if all individual migration in p2 & p3 but none in p1 => The new subopulation will be p1 and p2 and not p2 and p3

# WARNING: Doing explicit migration like this is NON Wright-Fisher simulation
# If we want to do WF simulation, we need to implement the migration in dunring the choosing parents
# => 10% migration = 90% sampling parent in "p1" and 10% to sample it in "p2"

migration = function(sim, migration.scheme, verbose=T){
  newpopstructure=unique(migration.scheme[,3])
  newpopstructure = order(newpopstructure)
  nsubpop = length(newpopstructure)
  
  sim.aftermigration = list()
  for(m in 1:nsubpop){
    newsubpop = newpopstructure[m]
    newsubpop=migration.scheme[migration.scheme[,3]==newsubpop,1:2]
    newsubpop = lapply(1:nrow(newsubpop), function(x){
      id = newsubpop[x,]
      select.ind(sim,id)
    })
    
    names(newsubpop) = paste0("i", 1:length(newsubpop))
    sim.aftermigration[[m]] = newsubpop
  }
  
  names(sim.aftermigration) = paste0("p", 1:length(sim.aftermigration))
  sim.aftermigration = c(sim.aftermigration, sim[c("rr", "sel.coef", "mut.rate")])
  
  if(nsubpop!=length(unique(migration.scheme[,1]))){
    sim.aftermigration[["rr"]] =as.matrix(matrix(sim.aftermigration[["rr"]][,1], ncol=1)[,rep(1,nsubpop)])
    sim.aftermigration[["sel.coef"]] =as.matrix(matrix(sim.aftermigration[["sel.coef"]][,1], ncol=1)[,rep(1,nsubpop)])
    sim.aftermigration[["mut.rate"]] =as.matrix(matrix(sim.aftermigration[["mut.rate"]][,1], ncol=1)[,rep(1,nsubpop)])
    if(verbose){print("WARNING: rr, sel.coef and mut.rate set as p1")}
  }
  
  return(sim.aftermigration)
}

# Example of migration scheme
# Empty for new population:
empty.migration.scheme = function(sim){
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  migration.scheme = dplyr::bind_rows(lapply(1:length(subpops.ids), function(x){
    subpop = subpops.ids[x]
    nind = length(sim[[subpop]])
    data.frame(original.pop = rep(x, nind), ind = as.numeric(1:nind))
  }))
  migration.scheme$newpop=NA
  migration.scheme
}

#Random migration scheme
# This is not an exact random mating population
# Some individual will be more related than in random mating
# Especially for small subpopulation size: ex:
# For N = 2 dioecious, 2 individual will share the exact same parent which can be unlikely in pure radom mating
# For big population size however, it should be quite identical to random mating, I guess
# This function is more for illustration purpose of the migration scheme
random.migration.scheme = function(sim){
  subpops.ids = names(sim)[!names(sim) %in% c("rr", "sel.coef", "mut.rate")]
  migration.scheme = dplyr::bind_rows(lapply(1:length(subpops.ids), function(x){
    subpop = subpops.ids[x]
    nind = length(sim[[subpop]])
    data.frame(original.pop = rep(x, nind), ind = as.numeric(1:nind))
  }))
  migration.scheme$newpop = sample(migration.scheme$original.pop, nrow(migration.scheme),replace=T)
  migration.scheme
}




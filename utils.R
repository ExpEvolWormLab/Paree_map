################################
### Maps Analysis Functions ####
################################

################################
### Utility ####################


#Assign recombination domains

assigndomain = function(newpos,newchrom,pos.vec,chrom.vec,domain.vec){
  
  #newpos = maps$pos
  #newchrom = maps$chrom
  #chrom.vec = cemap$chrom
  #domain.vec = cemap$domain
  #pos.vec = cemap$cpos
  
  domains = NULL
  for(chr in unique(newchrom)){
    dom = c(domain.vec[chrom.vec==chr], rep(NA,length(newpos[newchrom==chr])))
    pos = c(pos.vec[chrom.vec==chr], newpos[newchrom==chr])
    dom = dom[order(pos)]
    pos = pos[order(pos)]
    unkdom = which(is.na(dom))
    for( u in unkdom){if(u == 1){dom[u]='tip'}else{dom[u] = dom[(u-1)]}}
    dp = data.frame(pos=pos[unkdom], dom = dom[unkdom])
    dp$chrom=chr
    domains=rbind(domains,dp)
  }
  
  domains = domains[!duplicated(paste0(domains$chrom,"_", domains$pos)),]
  return(domains)
}



assign.cross.info = function(data){
  data$background[data$cross %in% c('AB', 'EF', "A", "B", "E", 'F')] = 1
  data$background[data$cross %in% c('CD', 'GH', "C", "D", "G", "H")] = 2
  data$background[data$cross %in% c('ABCD', 'EFGH')] = "all"
  data$rec[data$cross %in% c('AB', 'CD', 'ABCD',"A", "B", "C", "D")] = "wt"
  data$rec[data$cross %in% c('EF','GH', 'EFGH',"E", 'F',"G", "H")] = "mut"
  data$rec = factor(data$rec, levels = c("wt","mut"))
  return(data)
}

# Change chromosomes id Roman <-> num

numtorom = function(vector_chromid){
  vector_chromid = as.character(as.roman(vector_chromid))
  vector_chromid[vector_chromid=="VI"]="X"
  vector_chromid = factor(vector_chromid, levels = c("I","II","III","IV","V","X"))
  return(vector_chromid)
}

romtonum = function(vector_chromid){
  vector_chromid[vector_chromid=="I"]=1
  vector_chromid[vector_chromid=="II"]=2
  vector_chromid[vector_chromid=="III"]=3
  vector_chromid[vector_chromid=="IV"]=4
  vector_chromid[vector_chromid=="V"]=5
  vector_chromid[vector_chromid=="X"]=6
  return(as.numeric(vector_chromid))
}



# ggplot
theme_perso=function(){
  require(ggplot2)
  theme_minimal()+
    theme(
      strip.text.x = element_text(color="black",face="bold"),
      strip.text.y = element_text(color="black",face="bold"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.4))
}

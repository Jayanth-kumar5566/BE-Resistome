library(phyloseq)
genus_org<-readRDS("./../Data/R_input_files/Ivan_AMR_analysis_2020/CAMEB2/analysis/KAIJU/ps_Genus.RData")
species_org<-readRDS("./../Data/R_input_files/Ivan_AMR_analysis_2020/CAMEB2/analysis/KAIJU/ps_Species.RData")

genus<-as.data.frame(otu_table(genus_org))
species<-as.data.frame(otu_table(species_org))

#Filtering
filter<-function(x){
  x_tmp<-(x/rowSums(x))*100
  ind<-colSums(x_tmp>1)>0.05*251
  res<-x[,ind]
  return((res/rowSums(res))*100)
}

genus<-filter(genus)
species<-filter(species)


genus[,"SC_AMR_alt"]<-factor(sample_data(genus_org)$SC_AMR_alt)
species[,"SC_AMR_alt"]<-factor(sample_data(species_org)$SC_AMR_alt)


source("~/tmp/TRIKAFTA/correlation_network.R")

data<-species
for (i in levels(data$SC_AMR_alt)){
  sub_dat<-subset(data,data$SC_AMR_alt==i,select=-SC_AMR_alt)
  res<-correlation_network(sub_dat)
  saveRDS(res,paste0("./../Data/R_output_files/",i,"_adj.rds"))
}


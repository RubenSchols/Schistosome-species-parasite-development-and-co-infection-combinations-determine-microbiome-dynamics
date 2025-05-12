####### Install packages ####### 
#comes from: https://github.com/viennadelnat/microbiome-pesticide-selection/blob/main/Rcode.r

#Installing Phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("decontam")

#Install package that imports QIIME artifacts into phyloseq - requires devtools first 
install.packages("devtools")
install.packages("usethis")
library(usethis)
library(devtools)
devtools::install_github("jbisanz/qiime2R")   

## Installer Code for Packages Below If You Don't Have Them 
install.packages(c("vegan","microbiome","RColorBrewer","ggplot2","gplots","cowplot"))
install.packages(c("labdsv","tidyverse","ggfittext","dplyr","stats","ggpubr"))
install.packages(c("lme4","car","effects","emmeans","MuMIn","lmerTest","lattice","afex","DCA","picante"))

library(BiocManager)
BiocManager::install("microbiome")

####### Load packages ####### 

library(phyloseq);library(decontam);library(vegan);library(qiime2R); 

library(ggplot2);library(RColorBrewer);library(ggfittext);library(ggpubr);library(ggrepel)

library(dplyr); library("writexl");library(micrUBIfuns);library(moments);library(nlme);library(multcomp);library("QsRutils");library(microbiome)
library(reshape);library(matrixStats);library(lsmeans);library(pscl);library(readxl);library("DESeq2")

library("pheatmap");library("rlang")

library(lme4);library(car);library(effects);library(emmeans);library(lmerTest);library(afex);library(stats);library(picante);library(afex);library(MASS)
#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
#https://microbiome.github.io/tutorials/

sessionInfo()

####### Read data files ####### 

#Session --> Set Working Directory --> To Source File Location
setwd("~directory")

temporal <- readRDS("~directory/phyloseq_temporal.rds")

####### Cleaning and Filtering #######

## Decontam prevalence, 1 contaminating ASVs. PCR neg is seen as negatives here.
temporal@sam_data[["Negative"]] <- as.logical(temporal@sam_data[["Negative"]])
contamdf.prev01 <- isContaminant(temporal, method="prevalence",neg=temporal@sam_data[["Negative"]], threshold=0.05)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)
#print the contaminant information to file
badtaxa <- row.names(as.data.frame(subset(contamdf.prev01, contaminant == TRUE)))
alltaxa <- taxa_names(temporal)
# Exclude contaminant taxa from alltaxa to create goodtaxa
goodtaxa <- alltaxa[!alltaxa %in% badtaxa]
#remove the contaminants from the dataset
temporal <- prune_taxa(goodtaxa, temporal)

#Note,"Chloroplast" and "Mitochondria" ASVs in Genus --> Remove these ASVs (179)
temporal <- subset_taxa(temporal, 
                        !Genus %in% c("Mitochondria", "Chloroplast", "Lineage_IV", "LWQ8", "vadinHA49", 
                                      "WD2101_soil_group", "WPS-2", "Blfdi19", "CL500-29_marine_group", 
                                      "NS11-12_marine_group", "NS9_marine_group", 
                                      "SAR324_clade(Marine_group_B)", "PeM15","mle1-27", "Pla3_lineage", 
                                      "Yoonia-Loktanella", "SM2D12", "OM60(NOR5)_clade", "SM1A07", 
                                      "[Eubacterium]_brachy_group"))
temporal <-prune_taxa(taxa_sums(temporal)>0,temporal)

#Replace "uncultured*", "metagenome", "microbial*" , "unidentified*",.. as NA values in Species
NAlistSpecies= c("uncultured*","*like*","*_str.", ".*[0-9].*","mouse_gut","*_sp.", "filamentous*", "metagenome", "microbial*", "unidentified*", "iron-reducing*", "bacterium_enrichment*", "bacterium_RS5A")
for (pattern in NAlistSpecies) {tax_table(temporal)[, "Species"] <- gsub(tax_table(temporal)[, "Species"], pattern = pattern, replacement =  NA)}

#Note, "uncultured" in Genus, Family and Order; also "Unknown_Family" in Family; 
#Combined with all the names that have numbers etc in them
#NA values --> Replace these with nothing but keep ASVs!

NAlist= c("Candidatus_*", ".*_.*", ".*-.*", "uncultured","Unknown_.*", ".*[0-9].*", "filamentous_*")
for (pattern in NAlist) {tax_table(temporal)[, "Genus"] <- gsub(tax_table(temporal)[, "Genus"], pattern = pattern, replacement =  NA)}
for (pattern in NAlist) {tax_table(temporal)[, "Family"] <- gsub(tax_table(temporal)[, "Family"], pattern = pattern, replacement =  NA)}
for (pattern in NAlist) {tax_table(temporal)[, "Order"] <- gsub(tax_table(temporal)[, "Order"], pattern = pattern, replacement =  NA)}
for (pattern in NAlist) {tax_table(temporal)[, "Class"] <- gsub(tax_table(temporal)[, "Class"], pattern = pattern, replacement =  NA)}
for (pattern in NAlist) {tax_table(temporal)[, "Phylum"] <- gsub(tax_table(temporal)[, "Phylum"], pattern = pattern, replacement =  NA)}

## Remove non-bacterial ASVs (0)
temporal <- subset_taxa(temporal,!Kingdom %in% c("Unassigned"))
temporal <-prune_taxa(taxa_sums(temporal)>0,temporal)
DCDF=as.data.frame(tax_table(temporal))

#Create phylogenetic tree to assess contaminant positions or strange outliers.
tree1 = phy_tree(temporal)
ape::write.tree(tree1, "~directory/tree1.tree")

# if NA for order, family or genus level, make sure that all lower taxonomic levels are also NA. Otherwise you have species information but nothing higher so species is most likely wrong
#order level
for (i in 1:nrow(tax_table(temporal)[, "Order"])){
  if (is.na(tax_table(temporal)[, "Order"][i,1]) == FALSE) {
    NULL
  }
  else {
    tax_table(temporal)[, "Family"][i,1] <- NA
    tax_table(temporal)[, "Genus"][i,1] <- NA
    tax_table(temporal)[, "Species"][i,1] <- NA
  }
}
#family level
for (i in 1:nrow(tax_table(temporal)[, "Family"])){
  if (is.na(tax_table(temporal)[, "Family"][i,1]) == FALSE) {
    NULL
  }
  else {
    tax_table(temporal)[, "Genus"][i,1] <- NA
    tax_table(temporal)[, "Species"][i,1] <- NA
  }
}
#genus level
for (i in 1:nrow(tax_table(temporal)[, "Genus"])){
  if (is.na(tax_table(temporal)[, "Genus"][i,1]) == FALSE) {
    NULL
  }
  else {
    tax_table(temporal)[, "Species"][i,1] <- NA
  }
}

# for loop with if statement to check if genus and species levels are in agreement. If not it replaces the species level with NA because the silva database is focused up to genus level.
for (i in 1:nrow(tax_table(temporal)[, "Genus"])){
  if (is.na(tax_table(temporal)[, "Species"][i,1])){
    NULL
  }  else if (tax_table(temporal)[, "Genus"][i,1] == sapply(strsplit(tax_table(temporal)[, "Species"][i,1], split='_', fixed=TRUE), function(x) (x[1]))){
    NULL
  }   else {tax_table(temporal)[, "Species"][i,1] <- NA
  }
}

#Save temporal phyloseq object to save all the manual edits
saveRDS(temporal, "~directory/phyloseq_temporal_cleaned.rds") 
temporal <- readRDS("~directory/phyloseq_temporal_cleaned.rds")

####### rarefaction #######
##rarefaction curve after removing ASVs with less than 4 representatives per sample
temporal_tripletons <- temporal
# Loop over all samples in the phyloseq object
for (sample_id in sample_names(temporal_tripletons)) {
  # Get the abundance data for the current sample
  sample_abundance <- as(otu_table(temporal_tripletons), "matrix")[, sample_id]
  # Find ASVs detected less than 4 times in the current sample
  less_than_4 <- which(sample_abundance < 4)
  # Set the abundance of ASVs detected less than 4 times to 0
  sample_abundance[less_than_4] <- 0
  # Update the abundance data in the phyloseq object
  otu_table(temporal_tripletons)[, sample_id] <- sample_abundance
}
samples <- prune_taxa(taxa_sums(temporal_tripletons)>0,temporal_tripletons)

#Determine what the minimum coverage is 
sort(sample_sums(samples))
tab_samples <- otu_table(samples)
class(tab_samples) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab_samples <- t(tab_samples) # transpose observations to rows
#Plot rarefaction curve (species discovery curve) in vegan to see if plateau is reached for each sample
rarecurve(tab_samples, step=500, xlab = "# reads", ylab = "ASVs", rngseed = 711, label = FALSE)
abline(v=10417, col="blue", lty=2)
abline(v=8306, col="red", lty=2)
abline(v=3555, col="orange", lty=2)

# Function to rarefy multiple times and compute an average OTU table. 
rarefy_multiple_times <- function(physeq, iterations = 100, seed = 711,min_depth) {
  set.seed(seed)  # Set seed for reproducibility
  # Extract OTU table
  otu_matrix <- as(otu_table(physeq), "matrix")
  # Get minimum sample depth
  if (is.null(min_depth)) {
    min_depth <- min(sample_sums(physeq))}
  # Initialize an empty list to store rarefied OTU tables
  rarefied_otus <- vector("list", iterations)
  # Perform rarefaction multiple times
  for (i in seq_len(iterations)) {
    print(i)
    set.seed(seed + i)  # Change seed in every iteration for randomness
    rarefied_physeq <- rarefy_even_depth(physeq, sample.size = min_depth, verbose = FALSE, replace = FALSE)
    # Convert to matrix and ensure all taxa are aligned
    rarefied_otus[[i]] <- as(otu_table(rarefied_physeq), "matrix")
  }
  return(rarefied_otus)
}
combine_avg_multiple <- function(df_list) {
  # Get all unique row names across all data frames
  all_rows <- unique(unlist(lapply(df_list, rownames)))
  
  # Align all data frames to include all row names, filling missing values with zeros
  aligned_dfs <- lapply(df_list, function(df) {
    # Create a new data frame with all row names and zero-filled values
    df_aligned <- as.data.frame(matrix(0, nrow = length(all_rows), ncol = ncol(df), 
                                       dimnames = list(all_rows, colnames(df))))
    # Fill in existing values
    df_aligned[rownames(df), ] <- df
    return(df_aligned)
  })
  
  # Sum all aligned data frames and divide by the number of data frames to compute the average
  averaged_df <- Reduce(`+`, aligned_dfs) / length(df_list)
  
  return(averaged_df)
}

ps_rare <- rarefy_multiple_times(samples, iterations = 100, seed = 711,min_depth = 10417) #seed is here a single number, but in the function 1 is added in every loop to the set.seed to not have identical seeds; change min_depth to the necessary read depth
averaged_otus <- combine_avg_multiple(ps_rare)
averaged_otus <- as.matrix(averaged_otus)
averaged_otu_table <- otu_table(averaged_otus, taxa_are_rows = TRUE)
physeq_averaged <- phyloseq(averaged_otu_table, sample_data(samples), tax_table(samples),phy_tree(samples))
ps_rare <- physeq_averaged
otu_table(ps_rare) <- round(otu_table(ps_rare))
#check whether all ASVs are observed in the final dataset
ps_rare <- prune_taxa(taxa_sums(ps_rare) > 0, ps_rare) 
ps_rare
samples

#removing the last mock community
ps_final <- subset_samples(ps_rare, Code != "GC139407_G12")
ps_final <- prune_taxa(taxa_sums(ps_final) > 0, ps_final) 

saveRDS(ps_final, "~/directory/phyloseq_temporal_100averageRarefied.rds")
temporal_trip_RF <- readRDS("~/directory/phyloseq_temporal_100averageRarefied.rds")

saveRDS(samples, "~/directory/phyloseq_temporal_trip.rds")
samples <- readRDS("~/directory/phyloseq_temporal_trip.rds")



######### 3.1.1. Comparing miracidia exposed snail samples to control snails ##############
## prepare dataset
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
sample_data(testmem)
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
physeqFilter2RF <- testmem

#Determine Shannon index 
Shannon<-estimate_richness(physeqFilter2RF,measures=c("Shannon"))
#Determine Faith's phylogenetic index
OTUphyseqthesis=as.data.frame(physeqFilter2RF@otu_table)
TREEphyseqthesis=physeqFilter2RF@phy_tree 
TREEphyseqthesis #Rooted tree
FaithPD=pd(t(OTUphyseqthesis), TREEphyseqthesis,include.root=T)
#Add Richness Onto Our Metadata
Shannon$sampleID<-rownames(Shannon)
FaithPD$sampleID<-rownames(FaithPD)
alphaRF_1=sample_data(physeqFilter2RF)
alphaRF_1$sampleID<-rownames(alphaRF_1)
alphaRF_1<-left_join(Shannon,alphaRF_1,"sampleID")
alphaRF_1<-left_join(FaithPD,alphaRF_1,"sampleID")
alphaRF_1$daysPE <- as.factor(alphaRF_1$daysPE)
alphaRF_1$Sub_exp <- as.factor(alphaRF_1$Sub_exp)

#Statistics Shannon index
set_sum_contrasts()
fit=glm(formula = Shannon ~ Sub_exp , family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = Shannon ~ Sub_exp , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs # empty, so no outliers here
influenceIndexPlot(fit2,vars="Cook")

#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))

#The direction of skewness is given by the sign of the skewness coefficient:
#A zero means no skewness at all (normal distribution).
#A negative value means the distribution is negatively skewed.
#A positive value means the distribution is positively skewed.
#The skewness coefficient can be computed using the moments R packages
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/
skewness(alphaRF_1$Shannon, na.rm = TRUE) # -1.329376
#solutions for negatively skewed data:
#log transform
alphaRF_1$Shannonlog <- log10(max(alphaRF_1$Shannon+1) - alphaRF_1$Shannon)
skewness(alphaRF_1$Shannonlog, na.rm = TRUE) # -1.329376
#square root transform
alphaRF_1$Shannonsqrt <- sqrt(max(alphaRF_1$Shannon+1) - alphaRF_1$Shannon)
skewness(alphaRF_1$Shannonsqrt, na.rm = TRUE) # 1.035139
#inverse transform
alphaRF_1$Shannoninv <- 1/(max(alphaRF_1$Shannon+1) - alphaRF_1$Shannon)
skewness(alphaRF_1$Shannoninv, na.rm = TRUE) # -0.1916467

fit3=glm(formula = Shannoninv ~ Sub_exp , family = Gamma(link = log), data = alphaRF_1)
fit4=glm(formula = Shannoninv ~ Sub_exp , family = gaussian, data = alphaRF_1)
fit5=glm(formula = Shannonsqrt ~ Sub_exp , family = Gamma(link = log), data = alphaRF_1)
fit6=glm(formula = Shannonsqrt ~ Sub_exp , family = gaussian, data = alphaRF_1)
fit8=glm(formula = Shannonlog ~ Sub_exp , family = gaussian, data = alphaRF_1)
AIC(fit, fit2, fit3, fit4, fit5, fit6, fit8)

#testing normality
shapiro.test(resid(fit8))
histogram(resid(fit8))
# testing for heteroscedasticity
leveneTest(resid(fit8), alphaRF_1$Sub_exp)
leveneTest(alphaRF_1$Shannoninv, alphaRF_1$Sub_exp)
#Check of rule of thumb that ratio of maximum variance to minimum variance should not exceed 5
max(by(alphaRF_1$Shannonlog,alphaRF_1$Sub_exp,sd))^2/min(by(alphaRF_1$Shannonlog,alphaRF_1$Sub_exp,sd))^2

##normal but heteroscedastic data. Focus on the log transformed data but with heteroscedasticity accounted for.

#1) ANOVA with unequal variances:
oneway.test(Shannonlog ~ Sub_exp, data = alphaRF_1, var.equal = FALSE) ## Welsh ANOVA if heteroscedastic
#2) model heteroscedasticity:
mod.gls = gls(Shannonlog ~ Sub_exp, data=alphaRF_1,weights=varIdent(form= ~ 1 | Sub_exp))
anova(mod.gls)
summary(mod.gls)
summary(fit2)
# create ome helper functions to help multcomp for the pairwise comparissons
model.matrix.gls <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)
model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)
terms.gls <- function(object, ...)
  terms(model.frame(object),...)
# run the pairwise comparissons
mod.gls.mc = glht(mod.gls, linfct = mcp(Sub_exp = "Tukey"))
summary(mod.gls.mc)

SEplotShannon <- summary(emmeans(fit2, ~Sub_exp, type = "response"))

### Suppl. Figure 1A
#Plot shannon diversity
ggplot(SEplotShannon, aes(x = Sub_exp, y = emmean)) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5)) +
  #Set the shape of the points manually
  geom_point(aes(fill = Sub_exp), size = 5, 
             position = position_dodge(0), color=c("#08c0c5","#08c0c5","#08c0c5","#08c0c5", "#f87d74", "#f87d74")) +
  #Set x-label to 'Selection treatment'
  xlab("sample vs control") +
  #Set y-label to 'Shannon index'
  ylab(expression("untransformed Shannon")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(25, 2, 2, 30,"pt"), hjust = 1.1, vjust = 1.4, angle = 45,colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20,  vjust =1.9),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20,  hjust =0.8),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )

### Suppl. Figure 1B
SEplotShannon <- summary(emmeans(mod.gls, ~Sub_exp, type = "response"))
#Plot shannon diversity
ggplot(SEplotShannon, aes(x = Sub_exp, y = emmean)) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5)) +
  #Set the shape of the points manually
  geom_point(aes(fill = Sub_exp), size = 5, 
             position = position_dodge(0), color=c("#08c0c5","#08c0c5","#08c0c5","#08c0c5", "#f87d74", "#f87d74")) +
  #Set x-label to 'Selection treatment'
  xlab("sample vs control") +
  #Set y-label to 'Shannon index'
  ylab(expression("log Shannon")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(25, 2, 2, 30,"pt"), hjust = 1.1, vjust = 1.4, angle = 45,colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20,  vjust =1.9),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20,  hjust =0.8),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )



#Statistics PD index
set_sum_contrasts()
fit=glm(formula = PD ~ Sub_exp , family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = PD ~ Sub_exp , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs # empty, so no outliers here
influenceIndexPlot(fit2,vars="Cook")

# testing for heteroscedasticity
leveneTest(resid(fit2), alphaRF_1$Sub_exp)
leveneTest(alphaRF_1$PD, alphaRF_1$Sub_exp)
#Check of rule of thumb that ratio of maximum variance to minimum variance should not exceed 5
max(by(alphaRF_1$PD,alphaRF_1$Sub_exp,sd))^2/min(by(alphaRF_1$PD,alphaRF_1$Sub_exp,sd))^2

#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))

##normal but heteroscedastic data. Focus on the untransformed data but with heteroscedasticity accounted for.

#1) ANOVA with unequal variances:
oneway.test(PD ~ Sub_exp, data = alphaRF_1, var.equal = FALSE) ## Welsh ANOVA if heteroscedastic
car::Anova(fit2,type="II",white.adjust=TRUE) # I don't have an interaction factor so type II should be more appropriate here.
#2) model heteroscedasticity:
mod.gls = gls(PD ~ Sub_exp, data=alphaRF_1,weights=varIdent(form= ~ 1 | Sub_exp))
anova(mod.gls)
summary(mod.gls)
summary(fit2)

# run the pairwise comparissons
mod.gls.mc = glht(mod.gls, linfct = mcp(Sub_exp = "Tukey"))
summary(mod.gls.mc)


### Suppl. Figure 1C
SEplotPD <- summary(emmeans(mod.gls, ~Sub_exp, type = "response"))
#Plot PD diversity
ggplot(SEplotPD, aes(x = Sub_exp, y = emmean)) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5)) +
  #Set the shape of the points manually
  geom_point(aes(fill = Sub_exp), size = 5, 
             position = position_dodge(0), color=c("#08c0c5","#08c0c5","#08c0c5","#08c0c5", "#f87d74", "#f87d74")) +
  #Set x-label to 'Selection treatment'
  xlab("sample vs control") +
  #Set y-label to 'PD index'
  ylab(expression("Faith's PD")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(25, 2, 2, 30,"pt"), hjust = 1.1, vjust = 1.4, angle = 45,colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20,  vjust =1.9),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20,  hjust =0.8),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )




temporal_trip_RF_beta <- testmem
temporal_trip_RF_beta@phy_tree
#This error seems to occur when a phylogenetic tree ends up with an edge matrix with an odd number of rows
#One way to solve this problen is transforming all multichotomies into a series of dichotomies with one (or several) branch(es) of length zero.
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
#reroot the tree; especially important for the unifrac calculations
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@sam_data[["daysPE"]] <- as.factor(temporal_trip_RF_beta@sam_data[["daysPE"]])

temporal_trip_RF_beta <- testmem
temporal_trip_RF_beta@phy_tree
#This error seems to occur when a phylogenetic tree ends up with an edge matrix with an odd number of rows
#One way to solve this problen is transforming all multichotomies into a series of dichotomies with one (or several) branch(es) of length zero.
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
#reroot the tree; especially impoortant for the unifrac calculations
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@sam_data[["daysPE"]] <- as.factor(temporal_trip_RF_beta@sam_data[["daysPE"]])
set.seed=711
ordu_sa_ra = ordinate(temporal_trip_RF_beta ,method="NMDS",k=3, "bray", trymax= 400) #increased trymax to 1000 so that the model can converge. 20 was not enough
OrdiPlot_sa_ra = plot_ordination(temporal_trip_RF_beta, ordu_sa_ra, color="type", shape = "Sub_exp")
### Suppl. Figure 1D
OrdiPlot_sa_ra + geom_point(size=4, stroke=1.5, alpha=0.95) + coord_fixed(1/1)+theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  axis.line = element_line(color = "black")
)

df = as(sample_data(temporal_trip_RF_beta), "data.frame")
dist = phyloseq::distance(temporal_trip_RF_beta, "bray")
#betadisper tests homogeneity of dispersion among groups (dyasPE in this case), which is a condition (assumption) for adonis
physeq.disper <- betadisper(dist, df$type)
permutest(physeq.disper, pairwise = TRUE) 
#P <0.05 --> not ok; differences in composition within sampletypes (heterogeneous dispersion)

#Weighted unifrac distance
### Suppl. Figure 1E
wUF.ordu = ordinate(temporal_trip_RF_beta, method="NMDS", distance="unifrac", weighted=TRUE) #does not converge
plot_ordination(temporal_trip_RF_beta, wUF.ordu, color="type", shape = "Sub_exp") + theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  axis.line = element_line(color = "black")
)+  ggtitle("Weighted UniFrac")+ geom_point(size = 4)

#unweighted unifrac distance
### Suppl. Figure 1F
uwUF.ordu = ordinate(temporal_trip_RF_beta, method="NMDS", distance="unifrac", weighted=FALSE) #does not converge
plot_ordination(temporal_trip_RF_beta, uwUF.ordu, color="type", shape = "Sub_exp") +  theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  axis.line = element_line(color = "black")
)+   ggtitle("unWeighted UniFrac")+ geom_point(size = 4)

### Core microbiome
physeqFilter2RF <- testmem
#ASV level
x=0.001 # detection value
z=2 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
###prepare dataset
sample_data(physeqFilter2RF)
y=as.numeric(nrow(sample_data(physeqFilter2RF))-z)/nrow(sample_data(physeqFilter2RF)) # *100 and /100 are required to have the trunc function working.
# For this section of the data the rounding corresponds to one extra sample allowed to be below the 0.001 threshold. (99% of 212 samples)
pseq.rel <- microbiome::transform(physeqFilter2RF, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance))
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(physeqFilter2RF, rownames(tax_table(physeqFilter2RF)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)

######### 3.1.2. The temporal pattern throughout schistosome development (Figure 1B) #########  
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "1mir")
testmem <- subset_samples(testmem, daysPE != "2") # remove the first time point from this dataset to avoid aliased coeficients
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
physeqFilter2RF <- testmem
#Determine Shannon index 
Shannon<-estimate_richness(physeqFilter2RF,measures=c("Shannon", "Fisher", "Simpson", "ACE"))
#Determine Faith's phylogenetic index
OTUphyseqthesis=as.data.frame(physeqFilter2RF@otu_table)
TREEphyseqthesis=physeqFilter2RF@phy_tree 
TREEphyseqthesis #Rooted tree
FaithPD=pd(t(OTUphyseqthesis), TREEphyseqthesis,include.root=T)
#Add Richness Onto Our Metadata
Shannon$sampleID<-rownames(Shannon)
Shannon$sampleID<-gsub("X","",as.character(Shannon$sampleID))
FaithPD$sampleID<-rownames(FaithPD)
alphaRF_1=sample_data(physeqFilter2RF)
alphaRF_1$sampleID<-rownames(alphaRF_1)
alphaRF_1<-left_join(Shannon,alphaRF_1,"sampleID")
alphaRF_1<-left_join(FaithPD,alphaRF_1,"sampleID")
alphaRF_1$daysPE <- as.factor(alphaRF_1$daysPE)

### Suppl. Figure 2
#Statistics Shannon's index
set_sum_contrasts()
fit=glm(formula = Shannon ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = Shannon ~ daysPE*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
alphaRF_2=alphaRF_1[-outl,]
fit=glm(formula = Shannon ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_2)
fit2=glm(formula = Shannon ~ daysPE*infected , family = gaussian, data = alphaRF_2)
AIC(fit, fit2)
outlierTest(fit2) 
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs # empty, so no outliers here
influenceIndexPlot(fit2,vars="Cook")
# testing for heteroscedasticity
leveneTest(resid(fit2), interaction(alphaRF_2$daysPE,alphaRF_2$infected))
#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))
##normal homoscedastic data. 
Anova(fit2, type=3)
SEplotShannon <- summary(emmeans(fit2, ~daysPE*infected, type = "response"))
#Plot shannon diversity
Plot_Shannon <- ggplot(SEplotShannon, aes(x = daysPE, y = emmean)) +
  #Set white/black instead of red/green as color
  scale_fill_manual( values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = daysPE), shape = c(22,22,22,22,22,22, 21, 21, 21, 21, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("daysPE") +
  #Set y-label to 'Shannon index'
  ylab(expression("Shannon's index")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", hjust=1),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon
#pairwise comparissons
lsmeans(fit2, pairwise ~ infected| daysPE)

### Figure 2A
#### Statistics Faith's phylogenetic index
set_sum_contrasts()
fit=glm(formula = PD ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = PD ~ daysPE*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
alphaRF_2=alphaRF_1[-outl,]
fit=glm(formula = PD ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_2)
fit2=glm(formula = PD ~ daysPE*infected , family = gaussian, data = alphaRF_2)
AIC(fit, fit2)
outlierTest(fit2) 
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs 
influenceIndexPlot(fit2,vars="Cook")
# testing for heteroscedasticity
leveneTest(resid(fit2), interaction(alphaRF_2$daysPE,alphaRF_2$infected))
#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))
##normal homoscedastic data. 
Anova(fit2, type=3)
SEplotlmFaithsPD <- summary(emmeans(fit2, ~daysPE*infected, type = "response"))
#Plot PD diversity
Plot_lmFaithsPD <- ggplot(SEplotlmFaithsPD, aes(x = daysPE, y = emmean)) +
  #Set white/black instead of red/green as color
  scale_fill_manual(values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = daysPE), shape = c(22,22,22,22,22,22, 21, 21, 21, 21, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("daysPE") +
  #Set y-label to 'Shannon index'
  ylab(expression("lmFaithsPD")) +
  #Set x-ticks labels
  #scale_x_discrete(labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_lmFaithsPD
#pairwise comparissons
lsmeans(fit2, pairwise ~ infected| daysPE)


### Suppl. Figure 3A
#### Statistics Shannon index
set_sum_contrasts() 
Shannon<-glm(Shannon ~ prop_Schisto*daysPE, family = Gamma(link = log), data = alphaRF_1)
summary(Shannon)
shapiro.test(residuals(Shannon))
outlierTest(Shannon) 
outl=as.numeric(names(which(outlierTest(Shannon)$bonf.p<0.05)))
alphaRF_2 = alphaRF_1[-outl,] 
influenceIndexPlot(Shannon,vars=c("Studentized","Bonf"))
Shannon<-glm(Shannon ~ prop_Schisto*daysPE, family = Gamma(link = log), data = alphaRF_2)
Anova(Shannon, type=3)
shapiro.test(residuals(Shannon)) 
Anova(Shannon, type=3)
my_exp <- as.character(expression('McFadden pseudo-R²'))
R2 <- round(pR2(Shannon)[4],3) 
R2_sq <-paste(my_exp, R2, sep= ": ")
ggplot(Shannon, aes(x = prop_Schisto, y = Shannon, color = daysPE)) +
  geom_point() +  
  scale_color_manual( values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) + 
  labs(x = "Proportion of trematode reads", y = "Shannon", title = "Effect of Days PE and Proportion of Schisto on Shannon") +
  theme_classic() +
    annotate(
    "label",
    x = 0.2, y = 1.9,
    label = R2_sq,
    label.size = 0.5,
    label.r = unit(0.15, "lines")
  )


### Suppl. Figure 3B
#### Statistics Faith's phylogenetic index
set_sum_contrasts() 
lmFaithsPD<-glm(PD ~ prop_Schisto*daysPE, family = Gamma(link = log), data = alphaRF_1)
shapiro.test(residuals(lmFaithsPD)) 
outlierTest(lmFaithsPD) 
outl=as.numeric(names(which(outlierTest(lmFaithsPD)$bonf.p<0.05)))
alphaRF_2 = alphaRF_1[-outl,] 
influenceIndexPlot(lmFaithsPD,vars=c("Studentized","Bonf"))
lmFaithsPD<-glm(PD ~ prop_Schisto*daysPE, family = Gamma(link = log), data = alphaRF_2)
Anova(lmFaithsPD, type=3)
shapiro.test(residuals(lmFaithsPD)) #
my_exp <- as.character(expression('McFadden pseudo-R²'))
R2 <- round(pR2(lmFaithsPD)[4],3) 
R2_sq <-paste(my_exp, R2, sep= ": ")
ggplot(lmFaithsPD, aes(x = prop_Schisto, y = PD, color = daysPE)) +
  geom_point() +  
  scale_color_manual( values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) + 
  labs(x = "Proportion of trematode reads", y = "Faiths PD", title = "Effect of Days PE and Proportion of Schisto on Faiths PD") +
  theme_classic() +
  annotate(
    "label",
    x = 0.2, y = 18,
    label = R2_sq,
    label.size = 0.5,
    label.r = unit(0.15, "lines")
  )

### Suppl. Figure 3C
alphaRF_1_30 <- subset(alphaRF_1, daysPE == 30)
rownames(alphaRF_1_30) <- NULL
Shannon<-glm(Shannon ~ prop_Schisto, family = Gamma(link = log), data = alphaRF_1_30)
Shannon2<-glm(Shannon ~ prop_Schisto, family = gaussian, data = alphaRF_1_30)
AIC(Shannon, Shannon2)
outlierTest(Shannon2) 
outl=as.numeric(names(which(outlierTest(Shannon2)$bonf.p<0.05)))
influenceIndexPlot(Shannon2,vars=c("Studentized","Bonf"))
alphaRF_2=alphaRF_1_30[-outl,]
Shannon2<-glm(Shannon ~ prop_Schisto, family = gaussian, data = alphaRF_2)
shapiro.test(resid(Shannon2))
histogram(resid(Shannon2))
Anova(Shannon2, type=2)
summary(Shannon2)
my_exp <- as.character(expression('McFadden pseudo-R²'))
R2 <- round(pR2(Shannon2)['McFadden'],3)
R2_sq <-paste(my_exp, R2, sep= ":")
ggplot(Shannon2, aes(x = prop_Schisto, y = Shannon, color = "#85221c")) +
  geom_point() +  
  scale_color_manual(values = c("#85221c")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) + 
  labs(x = "Proportion of trematode reads", y = "Shannon", title = "Effect of Days PE and Proportion of Schisto on Shannon") +
  theme_classic() +
  theme (legend.position = "none")+
  annotate( "label",
    x = 0.2, y = 2,
    label = R2_sq,
    label.size = 0.5,
    label.r = unit(0.15, "lines")
  )

### Suppl. Figure 3D
#### Statistics lmFaithsPD index
lmFaithsPD<-glm(PD ~ prop_Schisto, family = Gamma(link = log), data = alphaRF_1_30)
lmFaithsPD2<-glm(PD ~ prop_Schisto, family = gaussian, data = alphaRF_1_30)
AIC(lmFaithsPD, lmFaithsPD2)
outlierTest(lmFaithsPD2) 
outl=as.numeric(names(which(outlierTest(lmFaithsPD2)$bonf.p<0.05)))
influenceIndexPlot(lmFaithsPD2,vars=c("Studentized","Bonf"))
alphaRF_2=alphaRF_1_30[-outl,]
lmFaithsPD2<-glm(PD ~ prop_Schisto, family = gaussian, data = alphaRF_2)
summary(lmFaithsPD2)
Anova(lmFaithsPD2, type=2)
my_exp <- as.character(expression('McFadden pseudo-R²'))
R2 <- round(pR2(lmFaithsPD2)['McFadden'],3)
R2_sq <-paste(my_exp, R2, sep= ":")
ggplot(lmFaithsPD2, aes(x = prop_Schisto, y = PD, color = "#85221c")) +
  geom_point() +  
  scale_color_manual(values = c("#85221c")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) + 
  labs(x = "Proportion of trematode reads", y = "Faiths PD", title = "Effect of Days PE and Proportion of Schisto on Faiths PD") +
  theme_classic() +
  theme (legend.position = "none")+
  annotate( "label",
            x = 0.2, y = 20,
            label = R2_sq,
            label.size = 0.5,
            label.r = unit(0.15, "lines")
  )


### Suppl. Figure 4A
#Statistics Shannon's index
set_sum_contrasts()
fit=glm(formula = Shannon ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = Shannon ~ daysPE*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs # empty, so no outliers here
influenceIndexPlot(fit2,vars="Cook")
# testing for heteroscedasticity
leveneTest(resid(fit2), interaction(alphaRF_1$daysPE,alphaRF_1$infected))
#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))
##normal homoscedastic data. 
Anova(fit2, type=3)
SEplotShannon <- summary(emmeans(fit2, ~daysPE*infected, type = "response"))
#Plot shannon diversity
Plot_Shannon <- ggplot(SEplotShannon, aes(x = daysPE, y = emmean)) +
  #Set white/black instead of red/green as color
  scale_fill_manual( values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = daysPE), shape = c(22,22,22,22,22,22, 21, 21, 21, 21, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("daysPE") +
  #Set y-label to 'Shannon index'
  ylab(expression("Shannon")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", hjust=1),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon
#pairwise comparissons
lsmeans(fit2, pairwise ~ infected| daysPE)

### Suppl. Figure 4B
#### Statistics Faith's phylogenetic index
set_sum_contrasts()
fit=glm(formula = PD ~ daysPE*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = PD ~ daysPE*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs 
influenceIndexPlot(fit2,vars="Cook")
# testing for heteroscedasticity
leveneTest(resid(fit2), interaction(alphaRF_1$daysPE,alphaRF_1$infected))
#testing normality
shapiro.test(resid(fit2))
histogram(resid(fit2))
##normal homoscedastic data. 
Anova(fit2, type=3)
SEplotlmFaithsPD <- summary(emmeans(fit2, ~daysPE*infected, type = "response"))
#Plot PD diversity
Plot_lmFaithsPD <- ggplot(SEplotlmFaithsPD, aes(x = daysPE, y = emmean)) +
  #Set white/black instead of red/green as color
  scale_fill_manual(values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = daysPE), shape = c(22,22,22,22,22,22, 21, 21, 21, 21, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("daysPE") +
  #Set y-label to 'Shannon index'
  ylab(expression("FaithsPD")) +
  #Set x-ticks labels
  #scale_x_discrete(labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_lmFaithsPD
#pairwise comparissons
lsmeans(fit2, pairwise ~ infected| daysPE)


### Suppl. Figure 4C
alphaRF_1_30 <- subset(alphaRF_1, daysPE == 30)
rownames(alphaRF_1_30) <- NULL
Shannon<-glm(Shannon ~ prop_Schisto, family = Gamma(link = log), data = alphaRF_1_30)
Shannon2<-glm(Shannon ~ prop_Schisto, family = gaussian, data = alphaRF_1_30)
AIC(Shannon, Shannon2)
# check for outliers
outlierTest(Shannon2) 
outl=as.numeric(names(which(outlierTest(Shannon2)$bonf.p<0.05)))
influenceIndexPlot(Shannon2,vars=c("Studentized","Bonf"))
alphaRF_1_30$outlier <- "no"
alphaRF_1_30$outlier[outl] <- "yes"
#testing normality
shapiro.test(resid(Shannon2))
histogram(resid(Shannon2))
Anova(Shannon2, type=2)
summary(Shannon2)
#Adj-R2 label
my_exp <- as.character(expression('McFadden pseudo-R'^2))
R2 <- round(pR2(Shannon2)['McFadden'],3)
#combine R2 value with expression to label the plot
R2_sq <-paste(my_exp, R2, sep= ":")
ggplot(alphaRF_1_30, aes(x = prop_Schisto, y = Shannon)) +
  geom_point(aes(color = outlier)) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) +  # Set the same color for se with lower alpha
  labs(x = "Proportion of trematode reads", y = "Shannon", title = "Effect of Days PE and Proportion of Schisto on Shannon") +
  theme_classic() +
  geom_label(aes(x=0.25,y=1.85, label = R2_sq,size = 15), parse = TRUE)+ 
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt"))
  )

### Suppl. Figure 4D
#### Statistics lmFaithsPD index
lmFaithsPD<-glm(PD ~ prop_Schisto, family = Gamma(link = log), data = alphaRF_1_30)
lmFaithsPD2<-glm(PD ~ prop_Schisto, family = gaussian, data = alphaRF_1_30)
AIC(lmFaithsPD, lmFaithsPD2)
# check for outliers
outlierTest(lmFaithsPD2) 
outl=as.numeric(names(which(outlierTest(lmFaithsPD2)$bonf.p<0.05)))
influenceIndexPlot(lmFaithsPD2,vars=c("Studentized","Bonf"))
alphaRF_1_30$outlier <- "no"
alphaRF_1_30$outlier[outl] <- "yes"
summary(lmFaithsPD2)
Anova(lmFaithsPD2, type=2)
#Adj-R2 label
my_exp <- as.character(expression('McFadden pseudo-R'^2))
R2 <- round(pR2(lmFaithsPD2)['McFadden'],3)
#combine R2 value with expression to label the plot
R2_sq <-paste(my_exp, R2, sep= ":")
ggplot(alphaRF_1_30, aes(x = prop_Schisto, y = PD)) +
  geom_point(aes(color = outlier)) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05) +  # Set the same color for se with lower alpha
  labs(x = "Proportion of trematode reads", y = "Faiths PD", title = "Effect of Days PE and Proportion of Schisto on FaithsPD") +
  theme_classic() +
  geom_label(aes(x=0.25,y=20, label = R2_sq,size = 15), parse = TRUE)+ 
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt"))
  )

### Suppl. Figure 5
qPCR <- read_excel("~/directory/qPCR.xlsx")
# Assuming your dataset is named 'qPCR_data'
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "1mir")
testmem <- subset_samples(testmem, daysPE != "2")
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
phydat <- sample_data(testmem)
dataframe1 <- as.data.frame(qPCR)
dataframe2 <- as.data.frame(phydat)
dataframe2$Sample <- row.names.data.frame(dataframe2) 
# Filter dataframe1 to only include samples that are in dataframe2
filtered_dataframe1 <- dataframe1 %>%  filter(Sample %in% dataframe2$Sample)
merged_data <- inner_join(dataframe1, dataframe2, by = c("Sample" = "Sample"))
merged_data
min(merged_data$Mean)
max(merged_data$Mean)
mean(merged_data$Mean)
sd(merged_data$Mean)
merged_data$daysPE <- as.factor(merged_data$daysPE)
fit=glm(formula = Mean ~ daysPE*infected, family = Gamma(link = log), data = merged_data) 
fit2=glm(formula = Mean ~ daysPE*infected, family = gaussian, data = merged_data) 
AIC(fit, fit2)
stepAIC(fit)
# testing for heteroscedasticity
leveneTest(resid(fit), interaction(merged_data$daysPE,merged_data$infected))
#testing normality
shapiro.test(resid(fit))
histogram(resid(fit))
Anova(fit, type=3)
outlierTest(fit) 
outl=as.numeric(names(which(outlierTest(fit)$bonf.p<0.05)))
outl
influenceIndexPlot(fit,vars=c("Studentized","Bonf"))
# check for influential observations
cd=cooks.distance(fit)
inflobs=which(cd>1) 
inflobs # empty, so no outliers here
influenceIndexPlot(fit,vars="Cook")
SEplotbactload <- summary(emmeans(fit, ~daysPE*infected, type = "response"))
#Plot shannon diversity
Plot_bactload <- ggplot(SEplotbactload, aes(x = daysPE, y = response)) +
  #Set white/black instead of red/green as color
  scale_fill_manual( values = c("#f4deac","#eab676", "#D0ADAD","#b16a40","#85221c","#643519")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = daysPE), shape = c(22,22,22,22,22,22, 21, 21, 21, 21, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("daysPE") +
  #Set y-label to 'Shannon index'
  ylab(expression("Bacterial load (ng)")) +
  #Set x-ticks labels
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", hjust=1),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_bactload
lsmeans(fit, pairwise ~ infected| daysPE)

#Figure 2B
TRISS_trip_RF_gen <- subset_samples(temporal_trip_RF, Sub_exp == "1mir" )
TRISS_trip_RF_gen <-prune_taxa(taxa_sums(TRISS_trip_RF_gen)>0,TRISS_trip_RF_gen)
sample_data(TRISS_trip_RF_gen)
OTU1 = as(otu_table(TRISS_trip_RF_gen), "matrix")
OTU1=t(OTU1)
df = as(sample_data(TRISS_trip_RF_gen), "matrix")
df=data.frame(df)
df$daysPE <- as.factor(df$daysPE)
null_tRDA = capscale(OTU1 ~ 1, distance = "bray")
tbRDA <- capscale(OTU1 ~ df$infected + df$daysPE + df$infected:df$daysPE, data = df, distance = "bray")
finalmodel<- ordistep(null_tRDA, scope=formula(tbRDA), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
#RDA-1 & 2 label creation
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA3 <- as.numeric(round(constrained_eig[3]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
db_RDA3 <-paste("db-RDA 3: ", db_RDA3, "%", sep= "")
#Adj-R2 label
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
#combine R2 value with expression to label the plot
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE==" 2"),aes(x=CAP1,y=CAP2,color="#119DA4",shape = as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE 2"),aes(x=CAP1,y=CAP2, color="#119DA4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE 2"),aes(x=CAP1,y=CAP2-0.02, color="#119DA4"),size=4,label="2")+ 
  geom_point(data=points_CAP%>%filter(daysPE==" 6"),aes(x=CAP1,y=CAP2,color="#13505B",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE 6"),aes(x=CAP1,y=CAP2,color="#13505B"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE 6"),aes(x=CAP1,y=CAP2-0.02,color="#13505B"),size=4,label="6")+
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color="#A18276",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2,color="#A18276"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2-0.02,color="#A18276"),size=4,label="10")+
  geom_point(data=points_CAP%>%filter(daysPE=="16"),aes(x=CAP1,y=CAP2,color="black",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2,color="black"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2-0.02,color="black"),size=4,label="16")+
  geom_point(data=points_CAP%>%filter(daysPE=="20"),aes(x=CAP1,y=CAP2,color="darkgrey",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2,color="darkgrey"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2-0.02,color="darkgrey"),size=4,label="20")+
  geom_point(data=points_CAP%>%filter(daysPE=="30"),aes(x=CAP1,y=CAP2,color="brown4",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2,color="brown4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2-0.02,color="brown4"),size=4,label="30")+
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color="#FFD999",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2,color="#FFD999"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2-0.02,color="#FFD999"),size=4,label="40")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.11,y=-0.3, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values=c("#ffe8b8","#f4deac","#eab676", "#643519","#D0ADAD","#85221c","#b16a40"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1

#Suppl. Figure 6A
### weighted UniFrac
temporal_trip_RF_beta <- TRISS_trip_RF_gen
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
unifrac_weighted <- UniFrac(temporal_trip_RF_beta, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
cap_w <- capscale(unifrac_weighted ~ df$infected + df$daysPE + df$infected:df$daysPE, data = df)
null_w <- capscale(unifrac_weighted ~ 1)
finalmodel<- ordistep(null_w, scope=formula(cap_w), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE=="2"),aes(x=CAP1,y=CAP2,color="#119DA4",shape = as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE2"),aes(x=CAP1,y=CAP2, color="#119DA4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE2"),aes(x=CAP1,y=CAP2-0.02, color="#119DA4"),size=4,label="2")+ 
  geom_point(data=points_CAP%>%filter(daysPE=="6"),aes(x=CAP1,y=CAP2,color="#13505B",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE6"),aes(x=CAP1,y=CAP2,color="#13505B"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE6"),aes(x=CAP1,y=CAP2-0.02,color="#13505B"),size=4,label="6")+
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color="#A18276",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2,color="#A18276"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2-0.02,color="#A18276"),size=4,label="10")+
  geom_point(data=points_CAP%>%filter(daysPE=="16"),aes(x=CAP1,y=CAP2,color="black",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2,color="black"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2-0.02,color="black"),size=4,label="16")+
  geom_point(data=points_CAP%>%filter(daysPE=="20"),aes(x=CAP1,y=CAP2,color="darkgrey",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2,color="darkgrey"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2-0.02,color="darkgrey"),size=4,label="20")+
  geom_point(data=points_CAP%>%filter(daysPE=="30"),aes(x=CAP1,y=CAP2,color="brown4",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2,color="brown4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2-0.02,color="brown4"),size=4,label="30")+
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color="#FFD999",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2,color="#FFD999"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2-0.02,color="#FFD999"),size=4,label="40")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.2,y=-0.35, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values=c("#ffe8b8","#f4deac","#eab676", "#643519","#D0ADAD","#85221c","#b16a40"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1


### Suppl. Figure 6B
#unweighted unifrac
temporal_trip_RF_beta <- TRISS_trip_RF_gen
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)
unifrac_unweighted <- UniFrac(temporal_trip_RF_beta, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
cap_unw <- capscale(unifrac_unweighted ~ df$infected + df$daysPE + df$infected:df$daysPE, data = df)
null_unw <- capscale(unifrac_unweighted ~ 1)
finalmodel<- ordistep(null_unw, scope=formula(cap_unw), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE=="2"),aes(x=CAP1,y=CAP2,color="#119DA4",shape = as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE2"),aes(x=CAP1,y=CAP2, color="#119DA4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE2"),aes(x=CAP1,y=CAP2-0.02, color="#119DA4"),size=4,label="2")+ 
  geom_point(data=points_CAP%>%filter(daysPE=="6"),aes(x=CAP1,y=CAP2,color="#13505B",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE6"),aes(x=CAP1,y=CAP2,color="#13505B"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE6"),aes(x=CAP1,y=CAP2-0.02,color="#13505B"),size=4,label="6")+
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color="#A18276",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2,color="#A18276"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE10"),aes(x=CAP1,y=CAP2-0.02,color="#A18276"),size=4,label="10")+
  geom_point(data=points_CAP%>%filter(daysPE=="16"),aes(x=CAP1,y=CAP2,color="black",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2,color="black"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE16"),aes(x=CAP1,y=CAP2-0.02,color="black"),size=4,label="16")+
  geom_point(data=points_CAP%>%filter(daysPE=="20"),aes(x=CAP1,y=CAP2,color="darkgrey",shape = as.factor(infected)),size=3, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2,color="darkgrey"),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE20"),aes(x=CAP1,y=CAP2-0.02,color="darkgrey"),size=4,label="20")+
  geom_point(data=points_CAP%>%filter(daysPE=="30"),aes(x=CAP1,y=CAP2,color="brown4",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2,color="brown4"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE30"),aes(x=CAP1,y=CAP2-0.02,color="brown4"),size=4,label="30")+
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color="#FFD999",shape = as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2,color="#FFD999"),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="daysPE40"),aes(x=CAP1,y=CAP2-0.02,color="#FFD999"),size=4,label="40")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.2,y=-0.3, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values=c("#ffe8b8","#f4deac","#eab676", "#643519","#D0ADAD","#85221c","#b16a40"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1

### Figure 2C
beta_violin_plot <- function(physeq, method) {
  
  # Load required packages
  require("phyloseq") # for phylogenetic sequence data management
  require("ggplot2")  # for plotting
  
  # Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq, method = method)
  beta_div_dist <- as.matrix(beta_div_dist)
  
  # Coerce distance matrix into a tidy data frame
  dist_df <- as.data.frame(as.table(beta_div_dist))
  names(dist_df) <- c("Sample1", "Sample2", "BetaDiv")
  dist_df$Infected1 <- sample_data(physeq)$infected[dist_df$Sample1]
  dist_df$Infected2 <- sample_data(physeq)$infected[dist_df$Sample2]
  dist_df$DaysPE <- sample_data(physeq)$daysPE[dist_df$Sample1] # Assuming daysPE is the same for both samples
  dist_df$DaysPE2 <- sample_data(physeq)$daysPE[dist_df$Sample2]
  
  # Exclude self-comparisons (where Sample1 == Sample2)
  dist_df <- dist_df[dist_df$Sample1 != dist_df$Sample2, ]
  
  # Filter to only include comparisons within the same infection status
  dist_df <- dist_df[dist_df$Infected1 == dist_df$Infected2, ]
  # Filter to only include comparisons within the same day status
  dist_df <- dist_df[dist_df$DaysPE == dist_df$DaysPE2, ]
  dist_df$Group <- paste(dist_df$DaysPE, dist_df$Infected1)  # combine daysPE and infected status for grouping
  
  # Create a new column with sorted Sample1 and Sample2
  dist_df$sorted_samples <- apply(dist_df[, c("Sample1", "Sample2")], 1, function(x) paste(sort(x), collapse = "-"))
  
  # Remove duplicates based on the sorted samples
  dist_df <- dist_df[!duplicated(dist_df$sorted_samples), ]
  
  # Remove the sorted_samples column
  dist_df$sorted_samples <- NULL
  
  # Define the order of the days and infection status
  days_order <- c("2", "6", "10", "16", "20", "30", "40")
  infection_order <- c("0", "1")
  ordered_groups <- outer(days_order, infection_order, paste)
  
  # Define the color palette
  color_palette <- setNames(c("#ffe8b8", "#ffe8b8", "#f4deac", "#f4deac", "#eab676", "#eab676", "#D0ADAD", 
                                       "#D0ADAD", "#b16a40", "#b16a40", "#85221c", "#85221c", "#643519", "#643519"),
                                       c("2 0" , "2 1" , "6 0" , "6 1" , "10 0" ,"10 1", "16 0" ,"16 1" ,"20 0" ,"20 1",
                                         "30 0" ,"30 1" ,"40 0" ,"40 1"))
  
  # Define comparisons
  comparisons <- list(c("6 0" , "6 1"), c("10 0", "10 1"), c("16 0", "16 1"), c("20 0", "20 1"), c("30 0", "30 1"), c("40 0", "40 1"))
  
  comparisons <- Filter(function(x) {
    all(x %in% unique(dist_df$Group) & table(dist_df$Group %in% x) > 1) & all(sapply(x, function(g) sum(dist_df$Group == g) > 1))
  }, comparisons)
  
  # Create a ggplot2 violin plot
  plot_violin <- ggplot(dist_df, aes(x = factor(Group, levels = as.vector(t(ordered_groups))), y = BetaDiv, fill = Group)) +
    geom_violin(trim = TRUE) +
    geom_jitter( width = 0.1, alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = color_palette) +
    labs(title = "Beta Diversity Violin Plot", y = "Bray-Curtis Dissimilarity", x = "Day and Infection Status") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
    stat_compare_means(comparisons = comparisons, label.y = c(0.09, 0.05, 0.01, 0.85, 0.89, 0.93), tip.length = 0.0, method = "wilcox.test", p.adjust.method = "bonferroni",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1, exact = FALSE), 
                                          symbols = c("****", "***", "**", "*", "ns")), hide.ns = TRUE) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          text = element_text(size = 12), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Save df and boxplot into a list
  list_Out <- list("data" = dist_df, "plot" = plot_violin)
  
  return(list_Out)
}
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "1mir")
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
temporal_trip_RF_beta <- testmem
temporal_trip_RF_beta@sam_data[["daysPE"]] <- as.factor(temporal_trip_RF_beta@sam_data[["daysPE"]])
temporal_trip_RF_beta@sam_data[["infected"]] <- as.factor(temporal_trip_RF_beta@sam_data[["infected"]])
set.seed(711)
plot2 <- beta_violin_plot(temporal_trip_RF_beta, method = "bray")
plot2$plot  


### Core microbiome
##### uninfected
#ASV level
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
testmem <- subset_samples(temporal_trip_RF, Sub_exp == "1mir" )
testmem <- subset_samples(testmem, infected == "0" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance))
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
rm(core.taxa.standard, core.abundance, testmem.subset)

##### infected
#ASV level
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
testmem <- subset_samples(temporal_trip_RF, Sub_exp == "1mir" )
testmem <- subset_samples(testmem, infected == "1" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance))
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
rm(core.taxa.standard, core.abundance, testmem.subset)


##### uninfected
#ASV level
x=0.001 # detection value
z=0 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
testmem <- subset_samples(temporal_trip_RF, Sub_exp == "1mir" )
testmem <- subset_samples(testmem, daysPE == "30" )
testmem <- subset_samples(testmem, infected == "0" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance))
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
rm(core.taxa.standard, core.abundance, testmem.subset)

##### infected
#ASV level
x=0.001 # detection value
z=0 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
testmem <- subset_samples(temporal_trip_RF, Sub_exp == "1mir" )
testmem <- subset_samples(testmem, daysPE == "30" )
testmem <- subset_samples(testmem, infected == "1" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance))
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
rm(core.taxa.standard, core.abundance, testmem.subset)


### Suppl Figure 7
set.seed(711)
temporal_trip <- readRDS("~/directory/phyloseq_temporal_merged2_trip.rds")
temporal_trip <- subset_samples(temporal_trip, experiment == "temporal")
temporal_trip <-prune_taxa(taxa_sums(temporal_trip)>0,temporal_trip)
sa <- temporal_trip
for (i in 1:ntaxa(sa)){
  if (nchar(tax_table(sa)[i,7], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,7], sep = "")
  }
  else if (nchar(tax_table(sa)[i,6], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,6], "_sp.", sep = "")
  }
  else if (nchar(tax_table(sa)[i,5], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,5], sep = "")
  }
  else if (nchar(tax_table(sa)[i,4], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,4], sep = "")
  }
  else if (nchar(tax_table(sa)[i,3], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,3], sep = "")
  }
  else if (nchar(tax_table(sa)[i,2], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,2], sep = "")
  }
  else {
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,1], sep = "")
  }
}
sa_sub <- sa
perform_deseq_analysis_ASV <- function(physeq_obj) {
  
  # DESeq2 analysis
  tryCatch({
    Class_glom_sa_Genus   <- tax_glom(physeq_obj, taxrank="Genus")
    deseq_obj = phyloseq_to_deseq2(Class_glom_sa_Genus, ~ infected)
    deseq_obj$Sub_exp <- relevel(deseq_obj$infected, ref = "0")
    deseq_obj <- DESeq(deseq_obj, test = "Wald", fitType = "local")
    res_obj = results(deseq_obj, cooksCutoff = FALSE)
    alpha = 0.05
    sig_res_obj = res_obj[which(res_obj$padj < alpha), ]
    
    if (nrow(sig_res_obj) == 0) {
      cat("No significant results to plot.\n")
      return(invisible(NULL))
    }
    sig_res_obj = cbind(as(sig_res_obj, "data.frame"), as(tax_table(Class_glom_sa_Genus)[rownames(sig_res_obj), ], "matrix"))
    
    # Create ggplot visualizations
    p <- ggplot(sig_res_obj, aes_string(x = "Family", y = "log2FoldChange", color = "Phylum")) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
      ggtitle(paste("DESeq2 Results at", "ASV", "level", "ref = uninfected"))
    print(p)
    # Save the ggplot as SVG
    plot_filename <- paste("DESeq2 Results at", "ref = uninfected","ASV", "level.svg")
    ggsave(plot_filename, plot = p, device = "svg", width = 10, height = 8)
    
  }, error = function(e) {
    cat("Error during DESeq2 execution:", e$message, "\nContinuing to next input value...\n")
    return(invisible(NULL))
  })
  # Set custom colors for DaysPE and infection
  custom_colors <- list(
    daysPE = c("2" = "#e9d9b5ff", "6"="#f3deabff", "10" = "#e8b575ff", "16"= "#cfababff", "20"="#b0683fff", "30"="#852119ff", "40" = "#633317ff"),
    infected = c("0" = "lightgrey", "1" = "black")
  )
  
  # Heatmap generation with error handling
  if (exists("sig_res_obj") && nrow(sig_res_obj) > 1) {
    tryCatch({
      select <- rownames(sig_res_obj)
      nt <- normTransform(deseq_obj)
      log2_norm_counts <- assay(nt)[select, ]
      df <- as.data.frame(colData(deseq_obj)[, c("daysPE", "infected")])
      title <- paste("log2(counts + 1) at", "ASV", "level", "ref = uninfected")
      heatmap_filename <- paste0("heatmap_at_","refco_uninfected" ,"ASV", "_level.svg")
      p2 <- pheatmap::pheatmap(log2_norm_counts, annotation_col = df, main = title, annotation_colors = custom_colors)
      ggsave(heatmap_filename, plot = p2, device = "svg", width = 10, height = 8)
    }, error = function(e) {
      cat("Error in generating heatmap :", e$message, "\n")
      dev.off()  # Ensure the device is closed even in case of error
    })
  } else {
    if (exists("sig_res_obj")) {
      cat("No significant differences to display in heatmap at ASV level.\n")
    } else {
      cat("No significant results object found for heatmap generation.\n")
    }
  }
  
  return(invisible(NULL))  # Ensure function returns invisibly under all conditions
}
sa_sub@sam_data[["daysPE"]] <- as.factor(sa_sub@sam_data[["daysPE"]])
sa_sub@sam_data[["infected"]] <- as.factor(sa_sub@sam_data[["infected"]])
sa_sub_2 <- subset_samples(sa_sub, Sub_exp == "1mir")
perform_deseq_analysis_ASV(sa_sub_2)


######### 3.1.3. The effect of parasite population and species in single infections on the microbiome (Figure 1C) #########  
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "10mir")
sample_data(testmem)
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
physeqFilter2RF <- testmem
Shannon<-estimate_richness(physeqFilter2RF,measures=c("Shannon"))
OTUphyseqthesis=as.data.frame(physeqFilter2RF@otu_table)
TREEphyseqthesis=physeqFilter2RF@phy_tree 
TREEphyseqthesis
FaithPD=pd(t(OTUphyseqthesis), TREEphyseqthesis,include.root=T)
Shannon$sampleID<-rownames(Shannon)
FaithPD$sampleID<-rownames(FaithPD)
alphaRF_1=sample_data(physeqFilter2RF)
alphaRF_1$sampleID<-rownames(alphaRF_1)
alphaRF_1<-left_join(Shannon,alphaRF_1,"sampleID")
alphaRF_1<-left_join(FaithPD,alphaRF_1,"sampleID")
alphaRF_1$daysPE <- as.factor(alphaRF_1$daysPE)

### Figure 3 & Suppl. Figure 8A 
#Statistics Shannon index
set_sum_contrasts()
fit=glm(formula = Shannon ~ X10mira*RD.PCR*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = Shannon ~ X10mira*RD.PCR*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
fit2 <-stepAIC(fit2)
shapiro.test(resid(fit2))
histogram(resid(fit2))
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
alphaRF_1[outl,]
alphaRF_2=alphaRF_1[-outl,]
fit2=glm(formula = Shannon ~ RD.PCR, family = gaussian, data = alphaRF_2)
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs
influenceIndexPlot(fit2,vars="Cook")
leveneTest(resid(fit2), alphaRF_2$RD.PCR)
max(by(alphaRF_2$Shannon,alphaRF_2$RD.PCR,sd))^2/min(by(alphaRF_2$Shannon,alphaRF_2$RD.PCR,sd))^2
shapiro.test(resid(fit2))
histogram(resid(fit2))
skewness(alphaRF_2$Shannon, na.rm = TRUE) 
alphaRF_2$Shannonlog <- log10(max(alphaRF_2$Shannon+1) - alphaRF_2$Shannon)
skewness(alphaRF_2$Shannonlog, na.rm = TRUE) 
alphaRF_2$Shannonsqrt <- sqrt(max(alphaRF_2$Shannon+1) - alphaRF_2$Shannon)
skewness(alphaRF_2$Shannonsqrt, na.rm = TRUE) 
alphaRF_2$Shannoninv <- 1/(max(alphaRF_2$Shannon+1) - alphaRF_2$Shannon)
skewness(alphaRF_2$Shannoninv, na.rm = TRUE) 
fit3=glm(formula = Shannoninv ~ RD.PCR , family = Gamma(link = log), data = alphaRF_2)
fit4=glm(formula = Shannoninv ~ RD.PCR , family = gaussian, data = alphaRF_2)
fit5=glm(formula = Shannonsqrt ~ RD.PCR , family = Gamma(link = log), data = alphaRF_2)
fit6=glm(formula = Shannonsqrt ~ RD.PCR , family = gaussian, data = alphaRF_2)
fit8=glm(formula = Shannonlog ~ RD.PCR , family = gaussian, data = alphaRF_2)
AIC(fit2, fit3, fit4, fit5, fit6, fit8)
shapiro.test(resid(fit4)) #Choose fit 4 since it meets normality, fit 8 does not.
histogram(resid(fit4))
leveneTest(resid(fit4), alphaRF_2$RD.PCR)
##normal but heteroscedastic data. Focus on the inverse transformed data but with heteroscedasticity accounted for.
oneway.test(Shannoninv ~ RD.PCR, data = alphaRF_2, var.equal = FALSE) ## Welsh ANOVA if heteroscedastic
mod.gls = gls(Shannoninv ~ RD.PCR, data=alphaRF_2,weights=varIdent(form= ~ 1 | RD.PCR))
anova(mod.gls)
summary(mod.gls)
summary(fit4)
model.matrix.gls <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)
model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)
terms.gls <- function(object, ...)
  terms(model.frame(object),...)
mod.gls.mc = glht(mod.gls, linfct = mcp(RD.PCR = "Tukey"))
summary(mod.gls.mc)

### Suppl. Figure 8A
SEplotShannon <- summary(emmeans(fit2, ~RD.PCR, type = "response"))
Plot_Shannon <- ggplot(SEplotShannon, aes(x = RD.PCR, y = emmean)) +
  scale_fill_manual(values = c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5.6, 
             position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5, 
             position = position_dodge(.5), color=c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  xlab("Parasite infection") +
  ylab(expression("untransformed Shannon")) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 17, margin = margin(2, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 17, margin = margin(20, 10, 25, 30, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon

### Figure 3A
#including the heteroscedasticity
SEplotShannon <- summary(emmeans(mod.gls, ~RD.PCR, type = "response"))
Plot_Shannon <- ggplot(SEplotShannon, aes(x = RD.PCR, y = emmean)) +
  scale_fill_manual(values = c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5.6, 
             position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5, 
             position = position_dodge(.5), color=c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  xlab("Parasite infection") +
  ylab(expression("1/(max(Shannon+1)
    - Shannon)")) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 17, margin = margin(2, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 17, margin = margin(20, 10, 25, 30, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon

### Suppl. Figure 8B & 8C
#Statistics Faith's phylogenetic index
set_sum_contrasts()
fit=glm(formula = PD ~ X10mira*RD.PCR*infected, family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = PD ~ X10mira*RD.PCR*infected , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
fit2 <-stepAIC(fit2)
shapiro.test(resid(fit2))
histogram(resid(fit2))
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
alphaRF_2=alphaRF_1[-outl,]
fit2=glm(formula = PD ~ RD.PCR, family = gaussian, data = alphaRF_2)
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs
influenceIndexPlot(fit2,vars="Cook")
leveneTest(resid(fit2), alphaRF_2$RD.PCR)
max(by(alphaRF_2$PD,alphaRF_2$RD.PCR,sd))^2/min(by(alphaRF_2$PD,alphaRF_2$RD.PCR,sd))^2
shapiro.test(resid(fit2))
histogram(resid(fit2))
skewness(alphaRF_2$PD, na.rm = TRUE) 
alphaRF_2$PDlog <- log10(max(alphaRF_2$PD+1) - alphaRF_2$PD)
skewness(alphaRF_2$PDlog, na.rm = TRUE) # -1.329376
alphaRF_2$PDsqrt <- sqrt(max(alphaRF_2$PD+1) - alphaRF_2$PD)
skewness(alphaRF_2$PDsqrt, na.rm = TRUE) # 1.035139
alphaRF_2$PDinv <- 1/(max(alphaRF_2$PD+1) - alphaRF_2$PD)
skewness(alphaRF_2$PDinv, na.rm = TRUE) # -0.1916467
fit3=glm(formula = PDinv ~ RD.PCR , family = Gamma(link = log), data = alphaRF_2)
fit4=glm(formula = PDinv ~ RD.PCR , family = gaussian, data = alphaRF_2)
fit5=glm(formula = PDsqrt ~ RD.PCR , family = Gamma(link = log), data = alphaRF_2)
fit6=glm(formula = PDsqrt ~ RD.PCR , family = gaussian, data = alphaRF_2)
fit8=glm(formula = PDlog ~ RD.PCR , family = gaussian, data = alphaRF_2)
AIC(fit3, fit4, fit5, fit6, fit8)
shapiro.test(resid(fit3))
histogram(resid(fit3))
leveneTest(resid(fit3), alphaRF_2$RD.PCR)
Anova(fit3, type=3)
lsmeans(fit3, pairwise ~  RD.PCR)
fit=glm(formula = PD ~ RD.PCR, family = Gamma(link = log), data = alphaRF_2)
SEplotPD <- summary(emmeans(fit, ~RD.PCR, type = "response"))
Plot_PD <- ggplot(SEplotPD, aes(x = RD.PCR, y = response)) +
  scale_fill_manual(values = c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5.6, 
             position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5, 
             position = position_dodge(.5), color=c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  xlab("Parasite infection") +
  ylab(expression("untransformed Faiths PD")) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 17, margin = margin(2, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 17, margin = margin(20, 10, 25, 30, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_PD

SEplotPD <- summary(emmeans(fit3, ~RD.PCR, type = "response"))
Plot_PD <- ggplot(SEplotPD, aes(x = RD.PCR, y = response)) +
  scale_fill_manual(values = c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5.6, 
             position = position_dodge(.5), color="black") +
  geom_point(aes(fill = RD.PCR), shape = c(22, 16, 16, 16), size = 5, 
             position = position_dodge(.5), color=c("#119DA4", "#A18276", "#5A7E0B", "#13505B")) +
  xlab("Parasite infection") +
  ylab(expression("1/(max(Faiths PD+1) - Faiths PD)")) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black", angle = 45, hjust=1),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 17, margin = margin(2, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 17, margin = margin(20, 10, 25, 30, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_PD



### Figure 3B
set.seed(711)
temporal_trip_RF_beta <- subset_samples(temporal_trip_RF, Sub_exp == "10mir")
temporal_trip_RF_beta <-prune_taxa(taxa_sums(temporal_trip_RF_beta)>0,temporal_trip_RF_beta)
sample_data(temporal_trip_RF_beta)
OTU1 = as(otu_table(temporal_trip_RF_beta), "matrix")
OTU1=t(OTU1)
df = as(sample_data(temporal_trip_RF_beta), "matrix")
df=data.frame(df)
null_tRDA = capscale(OTU1 ~ 1, distance = "bray")
tbRDA <- capscale(OTU1 ~  df$infected*df$X10mira, data = df, distance = "bray")
anova(tbRDA, permutations = 10000)
finalmodel<- ordistep(null_tRDA, scope=formula(tbRDA), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=ifelse(segments_CAP$var=="df$infected", "infected", "X10mira")
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(X10mira=="B"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1,y=CAP2, color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1-0.02,y=CAP2+0.02, color=as.factor(cluster)),size=4,label="SmBre")+ 
  geom_point(data=points_CAP%>%filter(X10mira=="L"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="SmLE")+
  geom_point(data=points_CAP%>%filter(X10mira=="R"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=15, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="uninfected")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="infected")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=0.11,y=-0.3, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("X10miraB" = "#A18276", "B" = "#A18276", "X10miraL"= "#5A7E0B", "L"="#5A7E0B", "X10miraR"="#13505B", "R"="#13505B", "infected1"= "#60564a", "infected0"= "#60564a"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1



### Suppl. Figure 9A
### weighted UniFrac
set.seed(711)
temporal_trip_RF_beta <- subset_samples(temporal_trip_RF, Sub_exp == "10mir")
temporal_trip_RF_beta <-prune_taxa(taxa_sums(temporal_trip_RF_beta)>0,temporal_trip_RF_beta)
sample_data(temporal_trip_RF_beta)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)
unifrac_weighted <- UniFrac(temporal_trip_RF_beta, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
cap_w <- capscale(unifrac_weighted ~ df$infected*df$X10mira, data = df)
anova(cap_w, permutations = 10000)
null_w <- capscale(unifrac_weighted ~ 1)
finalmodel<- ordistep(null_w, scope=formula(cap_w), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=ifelse(segments_CAP$var=="df$infected", "infected", "X10mira")
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(X10mira=="B"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1,y=CAP2, color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1-0.02,y=CAP2+0.02, color=as.factor(cluster)),size=4,label="SmBre")+ 
  geom_point(data=points_CAP%>%filter(X10mira=="L"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="SmLE")+
  geom_point(data=points_CAP%>%filter(X10mira=="R"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=15, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="uninfected")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="infected")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.2,y=-0.5, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("X10miraB" = "#A18276", "B" = "#A18276", "X10miraL"= "#5A7E0B", "L"="#5A7E0B", "X10miraR"="#13505B", "R"="#13505B", "infected1"= "#60564a", "infected0"= "#60564a"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1


### Suppl. Figure 9B
### unweighted UniFrac
unifrac_unweighted <- UniFrac(temporal_trip_RF_beta, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
cap_unw <- capscale(unifrac_unweighted ~ df$infected*df$X10mira, data = df)
anova(cap_unw, permutations = 10000)
null_unw <- capscale(unifrac_unweighted ~ 1)
finalmodel<- ordistep(null_unw, scope=formula(cap_unw), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=ifelse(segments_CAP$var=="df$infected", "infected", "X10mira")
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(X10mira=="B"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) + 
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1,y=CAP2, color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraB"),aes(x=CAP1-0.02,y=CAP2+0.02, color=as.factor(cluster)),size=4,label="SmBre")+ 
  geom_point(data=points_CAP%>%filter(X10mira=="L"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraL"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="SmLE")+
  geom_point(data=points_CAP%>%filter(X10mira=="R"),aes(x=CAP1,y=CAP2,color=as.factor(X10mira),shape=as.factor(infected)),size=3) +
  geom_point(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="X10miraR"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=15, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected0"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="uninfected")+
  geom_point(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=3.5,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="infected1"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="infected")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.09,y=-0.3, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("X10miraB" = "#A18276", "B" = "#A18276", "X10miraL"= "#5A7E0B", "L"="#5A7E0B", "X10miraR"="#13505B", "R"="#13505B", "infected1"= "#60564a", "infected0"= "#60564a"))+
  scale_shape_manual(values = c("0" = 15, "1" = 16))
plot1



### Figure 3C
# Within group Bray-Curtis dissimilarity
beta_violin_plot <- function(physeq, method) {
  
  # Packages
  require("phyloseq") # for phylogenetic sequence data management
  require("ggplot2")  # for plotting
  require("ggpubr")   # for stat_compare_means()
  
  # Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq, method = method)
  beta_div_dist <- as.matrix(beta_div_dist)
  
  # Coerce distance matrix into a tidy data frame
  dist_df <- as.data.frame(as.table(beta_div_dist))
  names(dist_df) <- c("Sample1", "Sample2", "BetaDiv")
  dist_df$Group1 <- sample_data(physeq)$RD.PCR[dist_df$Sample1]
  dist_df$Group2 <- sample_data(physeq)$RD.PCR[dist_df$Sample2]
  
  # Exclude self-comparisons (where Sample1 == Sample2)
  dist_df <- dist_df[dist_df$Sample1 != dist_df$Sample2, ]
  
  # Filter to only include comparisons within the same group
  dist_df <- dist_df[dist_df$Group1 == dist_df$Group2, ]
  dist_df$Group <- dist_df$Group1  # simplify to one group column
  
  # Create a new column with sorted Sample1 and Sample2
  dist_df$sorted_samples <- apply(dist_df[, c("Sample1", "Sample2")], 1, function(x) paste(sort(x), collapse = "-"))
  
  # Remove duplicates based on the sorted samples
  dist_df <- dist_df[!duplicated(dist_df$sorted_samples), ]
  
  # Remove the sorted_samples column
  dist_df$sorted_samples <- NULL
  
  # Define comparisons for statistical tests
  comparisons <- list(c("None", "SmBRE"), c("None", "SmLE"), c("None", "Sr"), 
                      c("SmBRE", "SmLE"), c("SmBRE", "Sr"), c("SmLE", "Sr"))
  
  comparisons <- Filter(function(x) {
    all(x %in% unique(dist_df$Group) & table(dist_df$Group %in% x) > 1) & all(sapply(x, function(g) sum(dist_df$Group == g) > 1))
  }, comparisons)
  
  # Create a ggplot2 violin plot
  plot_violin <- ggplot(dist_df, aes(x = Group, y = BetaDiv, fill = Group)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values=c("None"="#119DA4", "SmBRE"="#A18276", "SmLE"="#5A7E0B", "Sr"="#13505B")) +
    labs(title = "Beta Diversity Violin Plot", y = "Bray-Curtis dissimilarity", x = "Infection Diagnostic") +
    stat_compare_means(comparisons = comparisons, label.y = c(0.09, 0.05, 0.01, 0.93, 0.95, 0.9), 
                       tip.length=0.0, method = "wilcox.test", p.adjust.method = "bonferroni",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "ns")), hide.ns = TRUE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          text = element_text(size = 20), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Save df and violin plot into a list
  list_Out <- list("data" = dist_df, "plot" = plot_violin)
  
  return(list_Out)
}
temporal_trip_RF_beta <- subset_samples(temporal_trip_RF, experiment == "temporal")
temporal_trip_RF_beta <- subset_samples(temporal_trip_RF_beta, Sub_exp == "10mir")
sample_data(temporal_trip_RF_beta)
temporal_trip_RF_beta <-prune_taxa(taxa_sums(temporal_trip_RF_beta)>0,temporal_trip_RF_beta)
temporal_trip_RF_beta@sam_data[["Sub_exp"]] <- as.factor(temporal_trip_RF_beta@sam_data[["Sub_exp"]])
temporal_trip_RF_beta <-prune_taxa(taxa_sums(temporal_trip_RF_beta)>0,temporal_trip_RF_beta)
set.seed(711)
plot2 <- beta_violin_plot(temporal_trip_RF_beta, method = "bray")
plot2$plot  



### Core microbiome
TRISS_trip_RF_gen <- subset_samples(temporal_trip_RF, Sub_exp == "10mir")
TRISS_trip_RF_gen <-prune_taxa(taxa_sums(TRISS_trip_RF_gen)>0,TRISS_trip_RF_gen)
sample_data(TRISS_trip_RF_gen)
#ASV level
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
#all samples
testmem<-TRISS_trip_RF_gen
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) 
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
#uninfected
testmem <- subset_samples(TRISS_trip_RF_gen, RD.PCR == "None" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) 
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
#Sr
testmem <- subset_samples(TRISS_trip_RF_gen, RD.PCR == "Sr" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
#SmBRE
testmem <- subset_samples(TRISS_trip_RF_gen, RD.PCR == "SmBRE" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)
#SmLE
testmem <- subset_samples(TRISS_trip_RF_gen, RD.PCR == "SmLE" )
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) # *100 and /100 are required to have the trunc function working.
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)



######### 3.1.4. The effect of parasite population and species under co-infection on the microbiome (Figure 1D-E)####### 
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
sample_data(testmem)
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
physeqFilter2RF <- testmem
Shannon<-estimate_richness(physeqFilter2RF,measures=c("Shannon"))
OTUphyseqthesis=as.data.frame(physeqFilter2RF@otu_table)
TREEphyseqthesis=physeqFilter2RF@phy_tree 
FaithPD=pd(t(OTUphyseqthesis), TREEphyseqthesis,include.root=T)
Shannon$sampleID<-rownames(Shannon)
FaithPD$sampleID<-rownames(FaithPD)
alphaRF_1=sample_data(physeqFilter2RF)
alphaRF_1$sampleID<-rownames(alphaRF_1)
alphaRF_1<-left_join(Shannon,alphaRF_1,"sampleID")
alphaRF_1<-left_join(FaithPD,alphaRF_1,"sampleID")
alphaRF_1$daysPE <- as.factor(alphaRF_1$daysPE)

###Suppl. Figure 11
set_sum_contrasts()
fit=glm(formula = Shannon ~ RD.PCR*infected*daysPE*Sub_exp , family = Gamma(link = log), data = alphaRF_1)
fit2=glm(formula = Shannon ~ RD.PCR*infected*daysPE*Sub_exp , family = gaussian, data = alphaRF_1)
AIC(fit, fit2)
fit2 <-stepAIC(fit2)
outlierTest(fit2) 
outl=as.numeric(names(which(outlierTest(fit2)$bonf.p<0.05)))
outl
influenceIndexPlot(fit2,vars=c("Studentized","Bonf"))
cd=cooks.distance(fit2)
inflobs=which(cd>1) 
inflobs
influenceIndexPlot(fit2,vars="Cook")
leveneTest(resid(fit2), interaction(alphaRF_1$daysPE, alphaRF_1$Sub_exp))
max(by(alphaRF_1$Shannon,interaction(alphaRF_1$daysPE, alphaRF_1$Sub_exp),sd))^2/min(by(alphaRF_1$Shannon,interaction(alphaRF_1$daysPE, alphaRF_1$Sub_exp),sd))^2
shapiro.test(resid(fit2))
histogram(resid(fit2))
##normal homoscedastic data
Anova(fit2, type=3)
SEplotShannon <- summary(emmeans(fit2, ~daysPE + Sub_exp + daysPE:Sub_exp, type = "response"))
Plot_Shannon <- ggplot(SEplotShannon, aes(x = daysPE, y = emmean)) +
  scale_fill_manual(values = c("#A18276","#13505B")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = Sub_exp), shape = c(22), size = 5, 
             position = position_dodge(0), color="black") +
  xlab("DaysPE") +
  ylab(expression("Shannon")) +
  theme_classic() +
  scale_colour_manual(name="Legend",
                      labels=c("CoinfLeBr", "CoinfSrSm"),
                      values=c("white", "grey"))+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 16, margin = margin(25, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 16, margin = margin(2, 10, 4, 2, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon
lsmeans(fit, pairwise ~  Sub_exp|daysPE)
lsmeans(fit, pairwise ~  daysPE|Sub_exp)

### Figure 4A
set_sum_contrasts() 
lmFaithsPD<-glm(PD ~ daysPE*Sub_exp, family = Gamma(link = log), data = alphaRF_1)
fit2<-glm(PD ~ daysPE*Sub_exp, family = gaussian, data = alphaRF_1)
AIC(lmFaithsPD,fit2)
lmFaithsPD <- stepAIC(lmFaithsPD)
Anova(lmFaithsPD, type=3)
shapiro.test(residuals(lmFaithsPD))
hist(resid(lmFaithsPD)) 
qqnorm(resid(lmFaithsPD))
qqline(resid(lmFaithsPD))
outlierTest(lmFaithsPD) 
outl=as.numeric(names(which(outlierTest(lmFaithsPD)$bonf.p<0.05)))
outl 
cd=cooks.distance(lmFaithsPD)
inflobs=which(cd>1) 
inflobs 
influenceIndexPlot(lmFaithsPD,vars="Cook")
leveneTest(resid(lmFaithsPD), interaction(alphaRF_1$daysPE, alphaRF_1$Sub_exp))
#Homoscedastic and normal data
Anova(lmFaithsPD, type=3)
SEplotlmFaithsPD <- summary(emmeans(lmFaithsPD, ~daysPE + Sub_exp + daysPE:Sub_exp, type = "response"))
Plot_lmFaithsPD <- ggplot(SEplotlmFaithsPD, aes(x = daysPE, y = response)) +
  scale_fill_manual(values = c("#A18276","#13505B")) +
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  geom_point(aes(fill = Sub_exp), shape = c(22), size = 5, 
             position = position_dodge(0), color="black") +
  xlab("DaysPE") +
  ylab(expression("FaithsPD")) +
  theme_classic() +
  scale_colour_manual(name="Legend",
                      labels=c("CoinfLeBr", "CoinfSrSm"),
                      values=c("white", "grey"))+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    axis.title.x = element_text(size = 15, margin = margin(25, 2, 2, 2, "pt")),
    axis.title.y = element_text(size = 15, margin = margin(2, 30, 2, 2, "pt")),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_lmFaithsPD
lsmeans(lmFaithsPD, pairwise ~  Sub_exp|daysPE)
lsmeans(lmFaithsPD, pairwise ~  daysPE|Sub_exp)



###Figure 4B
TRISS_trip_RF_gen <- subset_samples(temporal_trip_RF, experiment == "temporal")
TRISS_trip_RF_gen <- subset_samples(TRISS_trip_RF_gen, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
TRISS_trip_RF_gen <-prune_taxa(taxa_sums(TRISS_trip_RF_gen)>0,TRISS_trip_RF_gen)
sample_data(TRISS_trip_RF_gen)
OTU1 = as(otu_table(TRISS_trip_RF_gen), "matrix")
OTU1=t(OTU1)
df = as(sample_data(TRISS_trip_RF_gen), "matrix")
df=data.frame(df)
df$daysPE <- as.factor(df$daysPE)
null_tRDA = capscale(OTU1 ~ 1, distance = "bray")
tbRDA <- capscale(OTU1 ~  df$daysPE*df$Sub_exp*df$infected, data = df, distance = "bray")
anova(tbRDA, permutations = 10000)
finalmodel<- ordistep(null_tRDA, scope=formula(tbRDA), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
test <- points_CAP%>%filter(daysPE==" 2"&Sub_exp=="co_infSrSm")
centroid_x_2_co_infSrSm <- mean(test$CAP1)
centroid_y_2_co_infSrSm <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_2_co_infSrSm, centroid_y_2_co_infSrSm,"na","na","na","2_co_infSrSm"))
test <- points_CAP%>%filter(daysPE==" 2"&Sub_exp=="co_infLeBr")
centroid_x_10_co_infLeBr <- mean(test$CAP1)
centroid_y_10_co_infLeBr <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_10_co_infLeBr, centroid_y_10_co_infLeBr,"na","na","na","2_co_infLeBr"))
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infSrSm")
centroid_x_10_co_infSrSm <- mean(test$CAP1)
centroid_y_10_co_infSrSm <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_10_co_infSrSm, centroid_y_10_co_infSrSm,"na","na","na","10_co_infSrSm"))
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infLeBr")
centroid_x_10_co_infLeBr <- mean(test$CAP1)
centroid_y_10_co_infLeBr <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_10_co_infLeBr, centroid_y_10_co_infLeBr,"na","na","na","10_co_infLeBr"))
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infSrSm")
centroid_x_40_co_infSrSm <- mean(test$CAP1)
centroid_y_40_co_infSrSm <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_40_co_infSrSm, centroid_y_40_co_infSrSm,"na","na","na","40_co_infSrSm"))
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infLeBr")
centroid_x_40_co_infLeBr <- mean(test$CAP1)
centroid_y_40_co_infLeBr <- mean(test$CAP2)
centroids_CAP <- rbind(centroids_CAP, list(centroid_x_40_co_infLeBr, centroid_y_40_co_infLeBr,"na","na","na","40_co_infLeBr"))
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE==" 2"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=4, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label= "40daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="40daysPE SmLE-Sr")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.12,y=-0.35, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("co_infLeBr" = "#A18276", "co_infSrSm" = "#13505B", "2_co_infLeBr"= "#A18276", "2_co_infSrSm"="#13505B", "10_co_infLeBr"="#A18276", "10_co_infSrSm"="#13505B", "40_co_infLeBr"="#A18276", "40_co_infSrSm"="#13505B"))+
  scale_shape_manual(values = c(" 2" = 17, "10" = 16, "40"= 18))
plot1



### Suppl. figure 12
TRISS_trip_RF_gen <- subset_samples(temporal_trip_RF, experiment == "temporal")
TRISS_trip_RF_gen <- subset_samples(TRISS_trip_RF_gen, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
TRISS_trip_RF_gen <-prune_taxa(taxa_sums(TRISS_trip_RF_gen)>0,TRISS_trip_RF_gen)
sample_data(TRISS_trip_RF_gen)
temporal_trip_RF_beta <- TRISS_trip_RF_gen
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta@phy_tree <- ape::multi2di(temporal_trip_RF_beta@phy_tree)
temporal_trip_RF_beta@phy_tree
temporal_trip_RF_beta <- root_phyloseq_tree(temporal_trip_RF_beta)

# Suppl. Figure 12A
# weighted unifrac
set.seed(711)
unifrac_weighted <- UniFrac(temporal_trip_RF_beta, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
cap_w <- capscale(unifrac_weighted ~ df$daysPE*df$Sub_exp*df$infected, data = df)
anova(cap_w, permutations = 10000)
null_w <- capscale(unifrac_weighted ~ 1)
finalmodel<- ordistep(null_w, scope=formula(cap_w), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
test <- points_CAP%>%filter(daysPE=="2"&Sub_exp=="co_infSrSm")
centroid_x_2_co_infSrSm <- mean(test$CAP1)
centroid_y_2_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_2_co_infSrSm,
  CAP2 = centroid_y_2_co_infSrSm,
  CAP3 = NA,
  cluster = "2_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="2"&Sub_exp=="co_infLeBr")
centroid_x_2_co_infLeBr <- mean(test$CAP1)
centroid_y_2_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_2_co_infLeBr,
  CAP2 = centroid_y_2_co_infLeBr,
  CAP3 = NA,
  cluster = "2_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infSrSm")
centroid_x_10_co_infSrSm <- mean(test$CAP1)
centroid_y_10_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_10_co_infSrSm,
  CAP2 = centroid_y_10_co_infSrSm,
  CAP3 = NA,
  cluster = "10_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infLeBr")
centroid_x_10_co_infLeBr <- mean(test$CAP1)
centroid_y_10_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_10_co_infLeBr,
  CAP2 = centroid_y_10_co_infLeBr,
  CAP3 = NA,
  cluster = "10_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infSrSm")
centroid_x_40_co_infSrSm <- mean(test$CAP1)
centroid_y_40_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_40_co_infSrSm,
  CAP2 = centroid_y_40_co_infSrSm,
  CAP3 = NA,
  cluster = "40_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infLeBr")
centroid_x_40_co_infLeBr <- mean(test$CAP1)
centroid_y_40_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_40_co_infLeBr,
  CAP2 = centroid_y_40_co_infLeBr,
  CAP3 = NA,
  cluster = "40_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE=="2"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=4, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label= "40daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="40daysPE SmLE-Sr")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=0.2,y=-0.45, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("co_infLeBr" = "#A18276", "co_infSrSm" = "#13505B", "2_co_infLeBr"= "#A18276", "2_co_infSrSm"="#13505B", "10_co_infLeBr"="#A18276", "10_co_infSrSm"="#13505B", "40_co_infLeBr"="#A18276", "40_co_infSrSm"="#13505B"))+
  scale_shape_manual(values = c("2" = 17, "10" = 16, "40"= 18))
plot1



#Suppl. Figure 12B
#unweighted unifrac
unifrac_unweighted <- UniFrac(temporal_trip_RF_beta, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)
df <- as(sample_data(temporal_trip_RF_beta), "data.frame")
df$daysPE <- as.factor(df$daysPE)
cap_unw <- capscale(unifrac_unweighted ~ df$daysPE*df$Sub_exp*df$infected, data = df)
anova(cap_unw, permutations = 10000)
null_unw <- capscale(unifrac_unweighted ~ 1)
finalmodel<- ordistep(null_unw, scope=formula(cap_unw), direction = "both", permutations = 10000)
finalmodel$anova
vif.cca(finalmodel) 
finalmodel
anova(finalmodel, permutations = 10000)
anova.cca(finalmodel, by="axis", step=1000) 
RsquareAdj (finalmodel)
ordiplot(finalmodel)
constrained_eig <- finalmodel$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- finalmodel$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
segments_CAP=data.frame(finalmodel$CCA$biplot)
segments_CAP$var=rownames(segments_CAP)
segments_CAP$var=gsub("df\\$", "",segments_CAP$var)
centroids_CAP=data.frame(finalmodel$CCA$centroids)
centroids_CAP$cluster=rownames(centroids_CAP)
test <-centroids_CAP$cluster
centroids_CAP$cluster=gsub("df\\$", "",centroids_CAP$cluster)
points_CAP=data.frame(finalmodel$CCA$wa)
points_CAP$id=rownames(points_CAP)
df$id=rownames(df)
points_CAP=left_join(points_CAP,df)
db_RDA1 <- as.numeric(round(constrained_eig[1]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA2 <- as.numeric(round(constrained_eig[2]/RsquareAdj(finalmodel)$r.squared,2))
db_RDA1 <-paste("db-RDA 1: ", db_RDA1, "%", sep= "")
db_RDA2 <-paste("db-RDA 2: ", db_RDA2, "%", sep= "")
my_exp <- as.character(expression('adjusted-R'^2))
R2 <- round(RsquareAdj(finalmodel)$adj.r.squared,2) 
R2_sq <-paste(my_exp, R2, sep= ":")
test <- points_CAP%>%filter(daysPE=="2"&Sub_exp=="co_infSrSm")
centroid_x_2_co_infSrSm <- mean(test$CAP1)
centroid_y_2_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_2_co_infSrSm,
  CAP2 = centroid_y_2_co_infSrSm,
  CAP3 = NA,
  cluster = "2_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="2"&Sub_exp=="co_infLeBr")
centroid_x_2_co_infLeBr <- mean(test$CAP1)
centroid_y_2_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_2_co_infLeBr,
  CAP2 = centroid_y_2_co_infLeBr,
  CAP3 = NA,
  cluster = "2_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infSrSm")
centroid_x_10_co_infSrSm <- mean(test$CAP1)
centroid_y_10_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_10_co_infSrSm,
  CAP2 = centroid_y_10_co_infSrSm,
  CAP3 = NA,
  cluster = "10_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="10"&Sub_exp=="co_infLeBr")
centroid_x_10_co_infLeBr <- mean(test$CAP1)
centroid_y_10_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_10_co_infLeBr,
  CAP2 = centroid_y_10_co_infLeBr,
  CAP3 = NA,
  cluster = "10_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infSrSm")
centroid_x_40_co_infSrSm <- mean(test$CAP1)
centroid_y_40_co_infSrSm <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_40_co_infSrSm,
  CAP2 = centroid_y_40_co_infSrSm,
  CAP3 = NA,
  cluster = "40_co_infSrSm")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
test <- points_CAP%>%filter(daysPE=="40"&Sub_exp=="co_infLeBr")
centroid_x_40_co_infLeBr <- mean(test$CAP1)
centroid_y_40_co_infLeBr <- mean(test$CAP2)
new_row <- tibble(
  CAP1 = centroid_x_40_co_infLeBr,
  CAP2 = centroid_y_40_co_infLeBr,
  CAP3 = NA,
  cluster = "40_co_infLeBr")
centroids_CAP <- bind_rows(centroids_CAP, new_row)
plot1 <-ggplot()+ ylab(db_RDA2) + xlab(db_RDA1)+ 
  geom_point(data=points_CAP%>%filter(daysPE=="2"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="10"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=3, stroke = 1.2) +
  geom_point(data=points_CAP%>%filter(daysPE=="40"),aes(x=CAP1,y=CAP2,color=as.factor(Sub_exp),shape=as.factor(daysPE)),size=4, stroke = 1.2) +
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=17, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="2_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="2daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=4,shape=16, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="10_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="10daysPE SmLE-Sr")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infLeBr"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label= "40daysPE SmLE-SmBRE")+
  geom_point(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2,color=as.factor(cluster)),size=5,shape=18, stroke = 1.2)+ 
  geom_text(data=centroids_CAP%>%filter(cluster=="40_co_infSrSm"),aes(x=CAP1,y=CAP2-0.02,color=as.factor(cluster)),size=4,label="40daysPE SmLE-Sr")+
  theme_classic()+ geom_hline(yintercept = 0, lty = 3, size=0.3) +
  geom_label(aes(x=-0.12,y=-0.35, label = R2_sq), parse = TRUE)+ 
  geom_vline(xintercept = 0, lty = 3, size=0.3)+
  theme(panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  scale_color_manual(values = c("co_infLeBr" = "#A18276", "co_infSrSm" = "#13505B", "2_co_infLeBr"= "#A18276", "2_co_infSrSm"="#13505B", "10_co_infLeBr"="#A18276", "10_co_infSrSm"="#13505B", "40_co_infLeBr"="#A18276", "40_co_infSrSm"="#13505B"))+
  scale_shape_manual(values = c("2" = 17, "10" = 16, "40"= 18))
plot1



###Figure 4C
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
sample_data(testmem)
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
temporal_trip_RF_beta <- testmem
temporal_trip_RF_beta@sam_data[["daysPE"]] <- as.factor(temporal_trip_RF_beta@sam_data[["daysPE"]])
temporal_trip_RF_beta@sam_data[["infected"]] <- as.factor(temporal_trip_RF_beta@sam_data[["infected"]])
temporal_trip_RF_beta@sam_data[["RD.PCR"]] <- as.factor(temporal_trip_RF_beta@sam_data[["RD.PCR"]])
temporal_trip_RF_beta <-prune_taxa(taxa_sums(temporal_trip_RF_beta)>0,temporal_trip_RF_beta)
beta_violin_plot <- function(physeq, method) {
  
  # Load required packages
  require("phyloseq") # for phylogenetic sequence data management
  require("ggplot2")  # for plotting
  
  # Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq, method = method)
  beta_div_dist <- as.matrix(beta_div_dist)
  
  # Coerce distance matrix into a tidy data frame
  dist_df <- as.data.frame(as.table(beta_div_dist))
  names(dist_df) <- c("Sample1", "Sample2", "BetaDiv")
  dist_df$Sub_exp1 <- sample_data(physeq)$Sub_exp[dist_df$Sample1]
  dist_df$Sub_exp2 <- sample_data(physeq)$Sub_exp[dist_df$Sample2]
  dist_df$DaysPE <- sample_data(physeq)$daysPE[dist_df$Sample1] # Assuming daysPE is the same for both samples
  dist_df$DaysPE2 <- sample_data(physeq)$daysPE[dist_df$Sample2]
  
  # Exclude self-comparisons (where Sample1 == Sample2)
  dist_df <- dist_df[dist_df$Sample1 != dist_df$Sample2, ]
  
  # Filter to only include comparisons within the same Sub_exp
  dist_df <- dist_df[dist_df$Sub_exp1 == dist_df$Sub_exp2, ]
  # Filter to only include comparisons within the same day status
  dist_df <- dist_df[dist_df$DaysPE == dist_df$DaysPE2, ]
  dist_df$Group <- paste(dist_df$DaysPE, dist_df$Sub_exp1)  # combine daysPE and Sub_exp for grouping
  
  # Create a new column with sorted Sample1 and Sample2
  dist_df$sorted_samples <- apply(dist_df[, c("Sample1", "Sample2")], 1, function(x) paste(sort(x), collapse = "-"))
  
  # Remove duplicates based on the sorted samples
  dist_df <- dist_df[!duplicated(dist_df$sorted_samples), ]
  
  # Remove the sorted_samples column
  dist_df$sorted_samples <- NULL
  
  # Define the order of the days and Sub_exp
  days_order <- c("2", "6", "10", "16", "20", "30", "40")
  sub_exp_order <- c("co_infLeBr", "co_infSrSm")
  ordered_groups <- outer(days_order, sub_exp_order, paste)
  
  # Define the color palette
  color_palette <- setNames(c("#A18276", "#13505B", "#A18276", "#13505B", "#A18276", "#13505B"),
                            c("2 co_infLeBr", "2 co_infSrSm", "10 co_infLeBr", "10 co_infSrSm", "40 co_infLeBr", "40 co_infSrSm"))
  
  # Define comparisons
  comparisons <- list( c("10 co_infSrSm", "40 co_infSrSm"),c("40 co_infLeBr", "40 co_infSrSm"), c("2 co_infSrSm", "10 co_infSrSm"), c("2 co_infLeBr", "10 co_infLeBr"), c("10 co_infLeBr", "40 co_infLeBr"))
  
  comparisons <- Filter(function(x) {
    all(x %in% unique(dist_df$Group) & table(dist_df$Group %in% x) > 1) & all(sapply(x, function(g) sum(dist_df$Group == g) > 1))
  }, comparisons)
  
  # Create a ggplot2 violin plot
  plot_violin <- ggplot(dist_df, aes(x = factor(Group, levels = as.vector(t(ordered_groups))), y = BetaDiv, fill = Group)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = color_palette) +
    labs(title = "Beta Diversity Violin Plot", y = "Bray-Curtis Dissimilarity", x = "Day and Sub Exp") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
    stat_compare_means(comparisons = comparisons, label.y = c(0.09, 0.05, 0.01, 0.85, 0.89, 0.93), tip.length = 0.0, method = "wilcox.test", p.adjust.method = "bonferroni",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "ns")), hide.ns = TRUE) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          text = element_text(size = 12), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Save df and violin plot into a list
  list_Out <- list("data" = dist_df, "plot" = plot_violin)
  
  return(list_Out)
}
set.seed(711)
plot2 <- beta_violin_plot(temporal_trip_RF_beta, method = "bray")
plot2$plot  



### core microbiome
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
#ASV level
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
sample_data(testmem)
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) 
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance), useNames=FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)

#Co-exposure SmLE-SmBRE
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "co_infLeBr")
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) 
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance), useNames = FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)

#Co-exposure SmLE-Sr
testmem <- subset_samples(temporal_trip_RF, experiment == "temporal")
testmem <- subset_samples(testmem, Sub_exp == "co_infSrSm")
testmem <-prune_taxa(taxa_sums(testmem)>0,testmem)
x=0.001 # detection value
z=1 #number of samples were the core microbiome is allowed missing
# prevalence value determined below according to the number of samples -z
y=as.numeric(nrow(sample_data(testmem))-z)/nrow(sample_data(testmem)) 
pseq.rel <- microbiome::transform(testmem, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE)
core.abundance <- sample_sums(core(pseq.rel, detection = x, prevalence = y, include.lowest=TRUE))
colMeans(as.matrix(core.abundance))
colSds(as.matrix(core.abundance),useNames=FALSE)
min(as.matrix(core.abundance))
max(as.matrix(core.abundance))
sort(as.matrix(core.abundance))
testmem.subset <- subset_taxa(testmem, rownames(tax_table(testmem)) %in% unique(core.taxa.standard))
tax_table(testmem.subset)



###Suppl. Figure 13
set.seed(711)
temporal_trip <- readRDS("~/directory/phyloseq_temporal_merged2_trip.rds")
temporal_trip <- subset_samples(temporal_trip, experiment == "temporal")
temporal_trip <-prune_taxa(taxa_sums(temporal_trip)>0,temporal_trip)
sa <- temporal_trip
for (i in 1:ntaxa(sa)){
  if (nchar(tax_table(sa)[i,7], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,7], sep = "")
  }
  else if (nchar(tax_table(sa)[i,6], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,6], "_sp.", sep = "")
  }
  else if (nchar(tax_table(sa)[i,5], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,5], sep = "")
  }
  else if (nchar(tax_table(sa)[i,4], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,4], sep = "")
  }
  else if (nchar(tax_table(sa)[i,3], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,3], sep = "")
  }
  else if (nchar(tax_table(sa)[i,2], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,2], sep = "")
  }
  else {
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,1], sep = "")
  }
}
sa_sub <- sa
setwd("~/directory/folder of choice")
perform_deseq_analysis_ASV <- function(physeq_obj) {
  
  # DESeq2 analysis
  tryCatch({
    deseq_obj = phyloseq_to_deseq2(physeq_obj, ~ Sub_exp)
    deseq_obj$Sub_exp <- relevel(deseq_obj$Sub_exp, ref = "co_infLeBr")
    deseq_obj <- DESeq(deseq_obj, test = "Wald", fitType = "local")
    res_obj = results(deseq_obj, cooksCutoff = FALSE)
    alpha = 0.05
    sig_res_obj = res_obj[which(res_obj$padj < alpha), ]
    
    if (nrow(sig_res_obj) == 0) {
      cat("No significant results to plot.\n")
      return(invisible(NULL))
    }
    sig_res_obj = cbind(as(sig_res_obj, "data.frame"), as(tax_table(physeq_obj)[rownames(sig_res_obj), ], "matrix"))
    
    # Create ggplot visualizations
    p <- ggplot(sig_res_obj, aes_string(x = "Family", y = "log2FoldChange", color = "Phylum")) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
      ggtitle(paste("DESeq2 Results at", "ASV", "level", "ref = co_infLeBr"))
    print(p)
    # Save the ggplot as SVG
    plot_filename <- paste("DESeq2 Results at", "ref = co_infLeBr","ASV", "level.svg")
    ggsave(plot_filename, plot = p, device = "svg", width = 10, height = 8)
    
  }, error = function(e) {
    cat("Error during DESeq2 execution:", e$message, "\nContinuing to next input value...\n")
    return(invisible(NULL))
  })
  # Set custom colors for DaysPE and Sub_exp
  custom_colors <- list(
    daysPE = c("2" = "#ffe8b8", "10" = "#eab676", "40" = "#643519"),
    Sub_exp = c("co_infLeBr" = "#A18276", "co_infSrSm" = "#13505B")
  )
  
  # Heatmap generation with error handling
  if (exists("sig_res_obj") && nrow(sig_res_obj) > 1) {
    tryCatch({
      select <- rownames(sig_res_obj)
      nt <- normTransform(deseq_obj)
      log2_norm_counts <- assay(nt)[select, ]
      df <- as.data.frame(colData(deseq_obj)[, c("daysPE", "Sub_exp")])
      title <- paste("log2(counts + 1) at", "ASV", "level", "ref = co_infLeBr")
      heatmap_filename <- paste0("heatmap_at_","refco_infLeBr" ,"ASV", "_level.svg")
      p2 <- pheatmap::pheatmap(log2_norm_counts, annotation_col = df, main = title, annotation_colors = custom_colors) # Colors for Sub_exp
      ggsave(heatmap_filename, plot = p2, device = "svg", width = 10, height = 8)
    }, error = function(e) {
      cat("Error in generating heatmap :", e$message, "\n")
      dev.off()  # Ensure the device is closed even in case of error
    })
  } else {
    if (exists("sig_res_obj")) {
      cat("No significant differences to display in heatmap at ASV level.\n")
    } else {
      cat("No significant results object found for heatmap generation.\n")
    }
  }
  
  return(invisible(NULL))  # Ensure function returns invisibly under all conditions
}
sa_sub@sam_data[["daysPE"]] <- as.factor(sa_sub@sam_data[["daysPE"]])
sa_sub_2 <- subset_samples(sa_sub, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
perform_deseq_analysis_ASV(sa_sub_2)



###Suppl. Figure 14
set.seed(711)
temporal_trip <- readRDS("~/directory/phyloseq_temporal_merged2_trip.rds")
temporal_trip <- subset_samples(temporal_trip, experiment == "temporal")
temporal_trip <-prune_taxa(taxa_sums(temporal_trip)>0,temporal_trip)
sa <- temporal_trip
for (i in 1:ntaxa(sa)){
  if (nchar(tax_table(sa)[i,7], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,7], sep = "")
  }
  else if (nchar(tax_table(sa)[i,6], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,6], "_sp.", sep = "")
  }
  else if (nchar(tax_table(sa)[i,5], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,5], sep = "")
  }
  else if (nchar(tax_table(sa)[i,4], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,4], sep = "")
  }
  else if (nchar(tax_table(sa)[i,3], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,3], sep = "")
  }
  else if (nchar(tax_table(sa)[i,2], keepNA = FALSE)>3){
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,2], sep = "")
  }
  else {
    taxa_names(sa)[i] = paste("ASV_", i, "_", tax_table(sa)[i,1], sep = "")
  }
}
sa_sub <- sa
setwd("~/directory/12may-coinfection")
perform_deseq_analysis_ASV <- function(physeq_obj) {
  
  # DESeq2 analysis
  tryCatch({
    Class_glom_sa_Genus   <- tax_glom(physeq_obj, taxrank="Genus")
    deseq_obj = phyloseq_to_deseq2(Class_glom_sa_Genus, ~ Sub_exp)
    deseq_obj$Sub_exp <- relevel(deseq_obj$Sub_exp, ref = "co_infLeBr")
    deseq_obj <- DESeq(deseq_obj, test = "Wald", fitType = "local")
    res_obj = results(deseq_obj, cooksCutoff = FALSE)
    alpha = 0.05
    sig_res_obj = res_obj[which(res_obj$padj < alpha), ]
    
    if (nrow(sig_res_obj) == 0) {
      cat("No significant results to plot.\n")
      return(invisible(NULL))
    }
    sig_res_obj = cbind(as(sig_res_obj, "data.frame"), as(tax_table(Class_glom_sa_Genus)[rownames(sig_res_obj), ], "matrix"))
    
    # Create ggplot visualizations
    p <- ggplot(sig_res_obj, aes_string(x = "Family", y = "log2FoldChange", color = "Phylum")) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
      ggtitle(paste("DESeq2 Results at", "genus", "level", "ref = co_infLeBr"))
    print(p)
    # Save the ggplot as SVG
    plot_filename <- paste("DESeq2 Results at", "ref = co_infLeBr","genus", "level.svg")
    ggsave(plot_filename, plot = p, device = "svg", width = 10, height = 8)
    
  }, error = function(e) {
    cat("Error during DESeq2 execution:", e$message, "\nContinuing to next input value...\n")
    return(invisible(NULL))
  })
  # Set custom colors for DaysPE and Sub_exp
  custom_colors <- list(
    daysPE = c("2" = "#ffe8b8", "10" = "#eab676", "40" = "#643519"),
    Sub_exp = c("co_infLeBr" = "#A18276", "co_infSrSm" = "#13505B"),
    RD.PCR =c("None" = "#FFFFFF",         # white
              "SmBRE" = "#CCCCCC",        # light grey
              "SmLE" = "#999999",         # medium-light grey
              "Sr" = "#666666",           # medium grey
              "SmBre+SmLE" = "#333333",   # dark grey
              "Sr+SmLE" = "#000000")       # black
  )
  
  # Heatmap generation with error handling
  if (exists("sig_res_obj") && nrow(sig_res_obj) > 1) {
    tryCatch({
      select <- rownames(sig_res_obj)
      nt <- normTransform(deseq_obj)
      log2_norm_counts <- assay(nt)[select, ]
      df <- as.data.frame(colData(deseq_obj)[, c("daysPE", "Sub_exp", "RD.PCR")])
      title <- paste("log2(counts + 1) at", "genus", "level", "ref = co_infLeBr")
      heatmap_filename <- paste0("heatmap_at_","refco_infLeBr" ,"genus", "_level.svg")
      p2 <- pheatmap::pheatmap(log2_norm_counts, annotation_col = df, main = title, annotation_colors = custom_colors) # Colors for Sub_exp
      ggsave(heatmap_filename, plot = p2, device = "svg", width = 10, height = 8)
    }, error = function(e) {
      cat("Error in generating heatmap :", e$message, "\n")
      dev.off()  # Ensure the device is closed even in case of error
    })
  } else {
    if (exists("sig_res_obj")) {
      cat("No significant differences to display in heatmap at genus level.\n")
    } else {
      cat("No significant results object found for heatmap generation.\n")
    }
  }
  
  return(invisible(NULL))  # Ensure function returns invisibly under all conditions
}
sa_sub@sam_data[["daysPE"]] <- as.factor(sa_sub@sam_data[["daysPE"]])
sa_sub@sam_data[["RD.PCR"]] <- as.factor(sa_sub@sam_data[["RD.PCR"]])
sa_sub_2 <- subset_samples(sa_sub, Sub_exp == "co_infLeBr" | Sub_exp == "co_infSrSm")
perform_deseq_analysis_ASV(sa_sub_2)


This repository contains the metadata associated with the manuscript "Schistosome species, parasite development, and co-infection combinations determine microbiome dynamics in the snail _Biomphalaria glabrata_" authored by Ruben Schols, Cyril Hammoud, Karen Bisschop, Isabel Vanoverberghe, Tine Huyse and Ellen Decaestecker

**Abstract **

Background
Schistosomiasis is a snail-borne disease affecting over 200 million people worldwide. Despite dedicated control efforts and effective diagnostic tools, schistosomiasis remains prevalent. Novel and sustainable control measures are urgently needed. Bacteria might offer such a solution as links between bacteria, disease resistance and transmission potential of intermediate hosts have been established in other systems. To better understand the tripartite interaction potentially driving snail-schistosome compatibility patterns, microbial communities must be investigated throughout and across various parasite exposure conditions. Therefore, we studied _Biomphalaria glabrata_ snails exposed to a high- and low-shedder population of _Schistosoma mansoni_ and _Schistosoma rodhaini_ in single and co-exposure experiments. Snails were sacrificed at different time points post-exposure and their bacterial communities and trematode (co-)infection status were determined through metabarcoding tools.
Results
Snails infected by low- and high-shedder _S. mansoni_ populations were more likely to have bacterial community dysbiosis than those infected by _S. rodhaini_ but this was also affected by miracidial load. Moreover, the single-infection hierarchical effect on the bacterial component of the microbiome is not maintained under co-infection with _S. rodhaini_, which appears to stabilize the snail’s bacterial profile even after being outcompeted by high-shedder _S. mansoni_. Finally, alpha diversity differed significantly between infected and uninfected snails around the onset period of shedding at 30 days post-miracidial exposure. 
Conclusion
The timing of this bacterial shift suggests an intricate parasite-snail interaction around key parasite development moments. Future studies investigating the tripartite interaction are advised to consider the effect of outcompeted or prepatent infections on the snail’s microbiome.


**File structure **

Raw sequencing data can be downloaded from ENA under bioproject PRJEB77058: https://www.ebi.ac.uk/ena/browser/view/PRJEB77058

The associated metadata can be found in the file: Datafiles ENA upload Temporal signal paper.xlsx

The assocaited bacterial load (ng) as determined by qPCR can by found in the file: qPCR.xlsx

The raw phyloseq object which has not been cleaned can be found here: phyloseq_temporal.rds

The  phyloseq object which has  been cleaned by removing contaminants: phyloseq_temporal_cleaned.rds

The phyloseq object which has been cleaned by removing up to and including tripletons per sample can be found here: phyloseq_temporal_trip.rds

The phyloseq object which has been cleaned by removing up to and including tripletons per sample and has been rarefied and averaged across 100 rarefaction events can be found here: phyloseq_temporal_100averageRarefied.rds

The Rscript used to produce all figures and analysis of the manuscript can be found here: Rscript_microbiome_dynamics_in_the_snail_Biomphalaria_glabrata.R
  Please adapt the directory to your specifications.

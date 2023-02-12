## install CRAN packages
lapply(c("BiocManager","RcppParallel","R.utils","clue","gtools",
         "targets","tarchetypes","visNetwork",
         "tidyverse","igraph","RColorBrewer",
         "Radviz","cytofan","hilbertSimilarity",
         "tsne","lme4","emmeans"),
       function(pck) {
           if(!requireNamespace(pck, quietly = TRUE)) {
               install.packages(pck)
           }
       })

## install BioConductor packages
BiocManager::install(c("CytoML","flowWorkspace",
                       "FlowSOM","flowViz","flowStats"))

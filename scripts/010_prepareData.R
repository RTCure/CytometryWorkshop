#
# Prepare data for cytometry workflow
#
################################################################

################# libraries & parameters #######################

library(CytoML)
library(flowWorkspace)
library(ncdfFlow)
library(flowCore)
library(openCyto)
library(ggcyto)
library(tidyverse)

## parameters
input.dir <- "input"
process.dir <- "process"
output.dir <- "output"

exp.id <- "workshop" # specify which folder contains the experiment

if(!file.exists(file.path(output.dir,exp.id))) {
    dir.create(file.path(output.dir,exp.id),
               recursive = TRUE,
               showWarnings = FALSE)
}

params_ref_markers <- "(cells|dead|CD45|CD3)"
params_pheno_markers <- "CD127|CD4|CD8|TBET|CD25|gdTCR|CD45RA|CD27|FOXP3|CD56"
params_func_markers <- "CCR6|IL4|CD38|IL22|CD28|CXCR4|CCR4|GMCSF|TGFb|IL17a|CXCR3|CD161|IL-10|CCR7|IL21|CTLA4|CXCR5|IL17F|HLA-DR|IFNg|TNF"

params_gate_filename <- "gates_workshop.xml"
params_cur_pop <- "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3"

######################## load data #############################

# unlink(file.path(process.dir,exp.id,"live"),
#        recursive = TRUE)

raw.files <- dir(file.path(input.dir,"raw"),
                 pattern="fcs$")
raw.files <- raw.files[str_detect(raw.files,"THC")]

rawSet <- cytobank_to_gatingset(file.path(input.dir,"raw",
                                          params_gate_filename),
                                file.path(input.dir,"raw",
                                          raw.files))

##################### prepare markers ##########################

as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
          rownames = "id")

.markers <- as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
                      rownames = "id") %>% 
    mutate(desc_orig = desc) %>% 
    separate(desc,
             c("isotope","desc"),
             sep = "_") %>% 
    mutate(type = case_when(
        str_detect(desc,params_ref_markers) ~ "reference",
        str_detect(desc,params_pheno_markers) ~ "phenotypic",
        str_detect(desc,params_func_markers) ~ "functional",
        TRUE ~ "CyTOF"),
    desc = if_else(is.na(desc),
                   name,
                   desc),
    desc = make.names(desc),
    type = factor(type,
                  levels = c("reference","phenotypic","functional","CyTOF")))

################### prepare annotations ########################
tibble(name = sampleNames(rawSet))

.annots <- tibble(name = sampleNames(rawSet)) %>% 
    mutate(Experiment = str_extract(name,"[0-9]{6}"),
           Experiment = factor(Experiment),
           Stimulation = str_extract(name,"(rest)|(stim)"),
           Stimulation = factor(Stimulation,
                                levels = c("rest","stim")),
           Individual = str_extract(name,"SK[0-9]{3}"),
           Individual = factor(Individual))

################## apply transformations #######################

asinhTrans <- asinhtGml2_trans(T=5.8760059682190064,
                               M=0.43429448190325176,
                               A=0.0)

trans <- gh_get_transformations(rawSet[[1]])

trans_channels <- .markers %>% 
    dplyr::filter(!name %in% names(trans),
                  type %in% c("phenotypic","functional")) %>% 
    distinct(name)

(trans_channels <- trans_channels[,"name",drop=TRUE])

rawSet <- transform(rawSet,
                    transformerList(from=trans_channels,
                                    trans=asinhTrans))

##################### add extra gates ##########################

gs_get_pop_paths(rawSet)

gs_add_gating_method(rawSet, 
                     gating_method = "mindensity", 
                     dims = "145Nd_CD4", 
                     parent = "CD45:CD3",
                     pop = "+",
                     alias = "CD4+")

gs_add_gating_method(rawSet, 
                     gating_method = "mindensity", 
                     dims = "146Nd_CD8", 
                     parent = "CD45:CD3", 
                     pop = "+",
                     alias = "CD8+")

ggcyto(rawSet,
       aes(x = CD4,
           y = CD8),
       subset = "CD45:CD3")+
    geom_hex()+
    geom_gate("CD4+")+
    geom_gate("CD8+")

## add a TNF gate to CD4 and CD8
ggcyto(rawSet,
       aes(x = TNF,
           y = CD8),
       subset = "CD8+")+
    geom_hex()

g <- lapply(rawSet,function(gh) {
    mindensity(gh_pop_get_data(gh,
                               "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3/CD8+"),
               "Lu175Di", filterId = "TNFstim")
})

lapply(str_which(names(g),"stim"),function(i) {
    gs_pop_add(rawSet[[c(i-1,i)]],
               g[[i]],
               parent = "CD8+")
})

# gs_pop_add(rawSet,
#            any.gate,
#            parent = "CD8+")

ggcyto(rawSet,
       aes(x = TNF,
           y = CD8),
       subset = "CD8+")+
    geom_hex()+
    geom_gate("TNFstim")

g <- lapply(rawSet,function(gh) {
    mindensity(gh_pop_get_data(gh,
                               "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3/CD4+"),
               "Lu175Di", filterId = "TNFstim")
})

lapply(str_which(names(g),"stim"),function(i) {
    gs_pop_add(rawSet[[c(i-1,i)]],
               g[[i]],
               parent = "CD4+")
})

## add a IFNg gate to CD4 and CD8
g <- lapply(rawSet,function(gh) {
    mindensity(gh_pop_get_data(gh,
                               "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3/CD8+"),
               "Yb174Di", filterId = "IFNstim")
})

lapply(str_which(names(g),"stim"),function(i) {
    gs_pop_add(rawSet[[c(i-1,i)]],
               g[[i]],
               parent = "CD8+")
})

g <- lapply(rawSet,function(gh) {
    mindensity(gh_pop_get_data(gh,
                               "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3/CD4+"),
               "Yb174Di", filterId = "IFNstim")
})

lapply(str_which(names(g),"stim"),function(i) {
    gs_pop_add(rawSet[[c(i-1,i)]],
               g[[i]],
               parent = "CD4+")
})

## make sure all gates are up-to-date
recompute(rawSet)

save_gs(rawSet,
        path = file.path(process.dir,exp.id,"live"),
        backend_opt = "copy")

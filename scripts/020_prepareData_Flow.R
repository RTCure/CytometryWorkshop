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

params_flow_filename <- "your_workspace.wsp"
params_cur_pop <- "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3"

######################## load data #############################

unlink(file.path(process.dir,exp.id,"live"),
       recursive = TRUE)

fj <- CytoML::open_flowjo_xml(params_flow_filename)
rawSet <- flowjo_to_gatingset(fj)

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

################## apply transformations #######################

trans <- gh_get_transformations(rawSet[[1]])

trans_channels <- .markers %>% 
    dplyr::filter(!name %in% names(trans),
                  type %in% c("phenotypic","functional")) %>% 
    distinct(name)

(trans_channels <- trans_channels[,"name",drop=TRUE])

if(length(trans_channels)>0) {
    trans.list <- estimateLogicle(rawSet[[1]],
                                  trans_channels)
    rawSet <- transform(rawSet,
                        trans.list)
}

## make sure everything is up-to-date
recompute(rawSet)

save_gs(rawSet,
        path = file.path(process.dir,exp.id,"live"),
        backend_opt = "copy")

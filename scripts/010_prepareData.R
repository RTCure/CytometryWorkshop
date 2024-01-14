#
# Prepare data for cytometry workflow
#
################################################################

## ---- libraries & parameters ----

library(CytoML)
library(flowWorkspace)
library(ncdfFlow)
library(flowCore)
library(openCyto)
library(ggcyto)
library(tidyverse)

## parameters
input.dir <- file.path("input","KI")
process.dir <- "process"
output.dir <- "output"

exp.id <- "workshop" # specify which folder contains the experiment

if(!file.exists(file.path(output.dir,exp.id))) {
    dir.create(file.path(output.dir,exp.id),
               recursive = TRUE,
               showWarnings = FALSE)
}

if(!file.exists(file.path(process.dir,exp.id,"live"))) {
    dir.create(file.path(process.dir,exp.id),
               recursive = TRUE,
               showWarnings = FALSE)
}

## ---- load data ----

params_flow_filename <- file.path("5hPMA_Iono","30-Oct-2023MGS.wsp")

fj <- CytoML::open_flowjo_xml(file.path(input.dir,params_flow_filename))
rawSet <- flowjo_to_gatingset(fj,
                              name = 1)

## make sure everything is up-to-date
recompute(rawSet)

plot(rawSet)

## ---- review markers ----

as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
          rownames = "id") %>% 
    View()

## ---- review transformations ----

trans <- gh_get_transformations(rawSet[[1]])

names(trans)
lapply(trans,attr,"type")

as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
          rownames = "id") %>% 
    filter(!name %in% names(trans))

(trans_channels <- c()) # add extra markers that need to be transformed here

if(length(trans_channels)>0) {
    trans.list <- estimateLogicle(rawSet[[1]],
                                  trans_channels)
    rawSet <- transform(rawSet,
                        trans.list)
}

## ---- review annotations ----

pData(rawSet) %>% 
    View()

## ---- review gates & channels ----

gs_get_pop_paths(rawSet)

autoplot(rawSet[[1]],
         gs_get_pop_paths(rawSet)[-1])

autoplot(rawSet[[1]],
         str_subset(gs_get_pop_paths(rawSet),"CD8"))+
    geom_stats(type = c("gate_name","percent"))

ggcyto(rawSet[1],
       aes(x = CD3),
       subset = "live")+
    geom_density()

ggcyto(rawSet[c(1,2)],
       aes(x = CD4,
           y = gd),
       subset = "CD3 subset, FSC-A")+
    geom_hex()+
    geom_gate("CD4")+
    geom_stats("CD4",
               type = c("gate_name","percent"))+
    geom_gate("gd")+
    geom_stats("gd",
               type = c("gate_name","percent"))+
    geom_gate("CD8 ")+
    geom_stats("CD8 ",
               type = c("gate_name","percent"))

## ---- save GatingSet ----

## make sure everything is up-to-date
recompute(rawSet)

unlink(file.path(process.dir,exp.id,"live"),
       recursive = TRUE)

save_gs(rawSet,
        path = file.path(process.dir,exp.id,"live"),
        backend_opt = "copy")

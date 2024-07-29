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
input.dir <- file.path("input","HELIOS share")
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

(fcsFiles <- dir(input.dir,pattern = "fcs",
                 recursive = TRUE,
                 full.names = TRUE))

fs <- read.ncdfFlowSet(fcsFiles)

rawSet <- GatingSet(fs)

## ---- review markers ----

as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
          rownames = "id") %>% 
    View()

## ---- review transformations ----

asinhTrans <- asinhtGml2_trans(T=5*sinh(1),
                               M=1/log(10),
                               A=0.0)

(trans_channels <- as_tibble(pData(parameters(gs_cyto_data(rawSet)[[1]])),
                             rownames = "id") %>% 
        filter(str_detect(desc,"_")) %>% 
        pull(name))

if(length(trans_channels)>0) {
    trans.list <- transformerList(trans_channels,
                                  asinhTrans)
    rawSet <- transform(rawSet,
                        trans.list)
}

## make sure everything is up-to-date
recompute(rawSet)

## ---- review gates & channels ----

gs_get_pop_paths(rawSet)

ggcyto(rawSet,
       aes(x = CD3,
           y = CD14),
       subset = "root")+
    geom_hex()

ggcyto(rawSet,
       aes(x = L_D),
       subset = "root")+
    geom_histogram()

## add some basic gates

ggcyto(rawSet,
       aes(x = Event_length,
           y = DNA1),
       subset = "root")+
    geom_hex()

g <- rectangleGate(filterId = "Singlets",
                   "Event_length" = c(20,80))

gs_pop_add(rawSet,g,parent = "root")
recompute(rawSet)

ggcyto(rawSet,
       aes(x = DNA1,
           y = DNA2),
       subset = "Singlets")+
    geom_hex()+
    geom_gate("DNA")

g <- rectangleGate(filterId = "DNA",
                   "Ir191Di" = c(5,8),
                   "Ir193Di" = c(6,9))

gs_pop_add(rawSet,g,parent = "Singlets")
recompute(rawSet)

ggcyto(rawSet[[1]],
       aes(x = L_D,
           y = CD45),
       subset = "DNA")+
    geom_hex()

g <- rectangleGate(filterId = "Live singlets",
                   "Rh103Di" = c(0,1))

gs_pop_add(rawSet,g,parent = "DNA")
recompute(rawSet)

g <- rectangleGate(filterId = "LymphoMono",
                   "Y89Di" = c(1,3),
                   "Eu151Di" = c(0,5))

gs_pop_add(rawSet,g,parent = "Live singlets")
recompute(rawSet)

ggcyto(rawSet,
       aes(x = CD3,
           y = CD56),
       subset = "LymphoMono")+
    geom_hex()

g1 <- rectangleGate(filterId = "CD3+",
                    "Pr141Di" = c(2,6))

g1neg <- rectangleGate(filterId = "CD3-",
                       "Pr141Di" = c(0,2))

g2 <- rectangleGate(filterId = "CD14+",
                    "Eu151Di" = c(1.5,4))

g2neg <- rectangleGate(filterId = "CD14-",
                       "Eu151Di" = c(0,1.5))

g3 <- rectangleGate(filterId = "CD19+",
                    "Nd142Di" = c(2,6))

g3neg <- rectangleGate(filterId = "CD19-",
                       "Nd142Di" = c(0,2))

g4 <- rectangleGate(filterId = "CD56+",
                    "Dy163Di" = c(1.5,7))

g4neg <- rectangleGate(filterId = "CD56-",
                       "Dy163Di" = c(0,1.5))

tcells <- g2neg * g3neg
tcells@filterId <- "T cells"

monos <- g2
monos@filterId <- "Monocytes"

bcells <- g3 * g4neg
bcells@filterId <- "B cells"

nkcells <- g4 * g3neg
nkcells@filterId <- "NK cells"

nktcells <- g4
nktcells@filterId <- "NKT cells"

gs_pop_add(rawSet,
           g1,
           parent = "LymphoMono")

gs_pop_add(rawSet,
           g1neg,
           parent = "LymphoMono")

gs_pop_add(rawSet,
           tcells,
           parent = "CD3+")

gs_pop_add(rawSet,
           nktcells,
           parent = "CD3+")

gs_pop_add(rawSet,
           monos,
           parent = "CD3-")

gs_pop_add(rawSet,
           bcells,
           parent = "CD3-")

gs_pop_add(rawSet,
           nkcells,
           parent = "CD3-")

plot(rawSet)
recompute(rawSet)

ggcyto(rawSet,
       aes(x = CD8a,
           y = TCRgd),
       subset = "T cells")+
    geom_hex()

g5 <- rectangleGate(filterId = "CD4+",
                    "Yb176Di" = c(2,6))

g5neg <- rectangleGate(filterId = "CD4-",
                       "Yb176Di" = c(0,2))

g6 <- rectangleGate(filterId = "CD8+",
                    "Nd146Di" = c(2,7))

g6neg <- rectangleGate(filterId = "CD8-",
                       "Nd146Di" = c(0,2))

g7 <- rectangleGate(filterId = "TCRgd+",
                    "Sm152Di" = c(1.5,5))

g7neg <- rectangleGate(filterId = "TCRgd-",
                       "Sm152Di" = c(0,1.5))

cd4s <- g5 * g6neg
cd4s@filterId <- "CD4 T cells"

cd8s <- g5neg * g6
cd8s@filterId <- "CD8 T cells"

gs_pop_add(rawSet,
           g7,
           parent = "T cells")

gs_pop_add(rawSet,
           g7neg,
           parent = "T cells")

gs_pop_add(rawSet,
           cd4s,
           parent = "TCRgd-")

gs_pop_add(rawSet,
           cd8s,
           parent = "TCRgd-")

plot(rawSet)
recompute(rawSet)

ggcyto(rawSet,
       aes(x = CD56,
           y = CD8a),
       subset = "NK cells")+
    geom_hex()

dims <- rectangleGate(filterId = "CD56dim NK cells",
                      "Dy163Di" = c(1.5,3))

brights <- rectangleGate(filterId = "CD56Br NK cells",
                         "Dy163Di" = c(3,7))

gs_pop_add(rawSet,
           dims,
           parent = "NK cells")

gs_pop_add(rawSet,
           brights,
           parent = "NK cells")

plot(rawSet)
recompute(rawSet)

ggcyto(rawSet,
       aes(x = CD45,
           y = CD34),
       subset = "Live singlets")+
    geom_hex()

autoplot(rawSet[[1]],y="Y89Di")

## ---- save GatingSet ----

## make sure everything is up-to-date
recompute(rawSet)

unlink(file.path(process.dir,exp.id,"live"),
       recursive = TRUE)

save_gs(rawSet,
        path = file.path(process.dir,exp.id,"live"),
        backend_opt = "copy")

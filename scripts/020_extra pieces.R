## ---- if extra transformations are required ----

(trans_channels <- c()) # add extra markers that need to be transformed here

if(length(trans_channels)>0) {
    trans.list <- estimateLogicle(rawSet[[1]],
                                  trans_channels)
    rawSet <- transform(rawSet,
                        trans.list)
}

## ---- if extra gates are required ----

# gs_add_gating_method(rawSet, 
#                      gating_method = "mindensity", 
#                      dims = "145Nd_CD4", 
#                      parent = "CD45:CD3",
#                      pop = "+",
#                      alias = "CD4+")
# 
# gs_remove_gating_method(rawSet)
# 
# g <- lapply(rawSet,function(gh) {
#     mindensity(gh_pop_get_data(gh,
#                                "/Bead Removal/Intact cells/Singlets/Viable/CD45+/CD45:CD3/CD8+"),
#                "Lu175Di", filterId = "TNFstim")
# })
# 
# gs_pop_add(rawSet,
#            g,
#            parent = "live")
# 
# gs_pop_remove(rawSet,
#               "TNFstim")


## ---- Where are the gd T cells in Radviz? ----
gate_gd <- lapply(rawSet,gh_pop_get_indices,
                  "/Lymphocytes+Monocytes/Single Cells/live/CD3 subset, FSC-A/gd")
gate_gd <- unlist(gate_gd)

gate_live <- lapply(rawSet,gh_pop_get_indices,
                    params_population)
gate_live <- unlist(gate_live)

sum(gate_live)

live_gate_gd <- gate_gd[gate_live]

plot(pheno.rv)+
    geom_density2d(data = . %>% 
                       group_by(name) %>% 
                       sample_frac(0.1),
                   aes(color = Condition))+
    geom_density2d(data = . %>% 
                       mutate(Population = live_gate_gd) %>% 
                       dplyr::filter(Population),
                   color = "grey30")+
    facet_grid(cols = vars(Condition))

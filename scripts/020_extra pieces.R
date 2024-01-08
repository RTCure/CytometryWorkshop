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

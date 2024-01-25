## ---- adding extra automated gates ----
#
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

## ---- adding extra manual gates ----
# 
# ggcyto(rawSet[c(1,2)],
#        aes(x = TNF),
#        subset = "CD4")+
#     geom_density()+
#     geom_vline(xintercept = 85,
#                col = 2, lty = 2)
# 
# g <- rectangleGate(filterId = "TNF+",
#                    "Comp-BV510-A" = c(85,Inf))
# 
# gs_pop_add(rawSet,g,parent = "CD4")
# 
# plot(rawSet)
# 
# ggcyto(rawSet[c(1,2)],
#        aes(x = TNF),
#        subset = "CD4")+
#     geom_density()+
#     geom_gate("TNF+")
# 
# recompute(rawSet)

library(targets)
library(dplyr)

tar_load(live.df)
tar_load(live.fsom.clusters)
tar_load(params_fsom_channels)

live.mat <- live.df %>% 
    select(all_of(params_fsom_channels)) %>% 
    as.matrix()

do_BIC_SPSS <- function(mat,i) {
    # X is data matrix, N objects x P quantitative variables.
    # Y is column of length N designating cluster membership; clusters 1, 2,..., K.
    # 1. Compute 1 x K row Nc showing number of objects in each cluster.
    # 2. Compute P x K matrix Vc containing variances by clusters.
    # Use denominator "n", not "n-1", to compute those, because there may be clusters with just one object.
    # 3. Compute P x 1 column containing variances for the whole sample. Use "n-1" denominator.
    # Then propagate the column to get P x K matrix V.
    # 4. Compute log-likelihood LL, 1 x K row. LL = -Nc &* csum( ln(Vc + V)/2 ),
    # where "&*" means usual, elementwise multiplication;
    # "csum" means sum of elements within columns.
    # 5. Compute BIC value. BIC = -2 * rsum(LL) + 2*K*P * ln(N),
    # where "rsum" means sum of elements within row.
    # https://stats.stackexchange.com/questions/55147/compute-bic-clustering-criterion-to-validate-clusters-after-k-means/55160#55160
    
    N <- nrow(mat)
    P <- ncol(mat)
    
    Nc <- table(i)
    K <- length(Nc)
    
    Vc <- by(mat,i, function(x) apply(x,2,var))
    Vc <- do.call("rbind",Vc)
    
    V <- apply(mat,2,var)
    
    LL <- log(sweep(Vc,2,V,"+"))
    LL <- rowSums(LL/2)
    
    LL <- -Nc * LL
    return(-2 * sum(LL) + 2*K*P * log(N))
}

do_BIC_SPSS(live.mat,live.fsom.clusters) # this one should be smaller
do_BIC_SPSS(live.mat,sample(live.fsom.clusters)) # than this one which is basically random!


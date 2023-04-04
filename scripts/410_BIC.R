library(targets)
library(dplyr)

do_BIC <- function(mat,i) {
    # after https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
    .dst <- function(x) {
        .center <- colMeans(x)
        x <- sweep(x,2,.center)
        return(sum(x^2))
    }
    m <- ncol(mat)
    n <- nrow(mat)
    k <- length(unique(i))
    D <- by(mat,i,.dst)
    D <- sum(D)
    return(D+log(n)*m*k)
}

tar_load(live.df)
tar_load(live.fsom.clusters)
tar_load(params_fsom_channels)

live.mat <- live.df %>% 
    select(all_of(params_fsom_channels)) %>% 
    as.matrix()

do_BIC(live.mat,live.fsom.clusters) # this one should be smaller
do_BIC(live.mat,sample(live.fsom.clusters)) # than this one which is basically random!

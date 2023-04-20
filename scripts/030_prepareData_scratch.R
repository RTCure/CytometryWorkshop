###################### Loading arbitrary data ######################
# author : ayann@its.jnj.com

###################### flowFrame ######################

mat <- matrix(rnorm(1500),ncol = 15)
colnames(mat) <- LETTERS[1:15]
flowFrame(mat)

###################### flowSet ######################

ffs <- lapply(1:10,function(i) {
    mat <- matrix(rnorm(1500),ncol = 15)
    colnames(mat) <- LETTERS[1:15]
    flowFrame(mat)
})
names(ffs) <- colors()[1:10]

fs <- flowSet(ffs)
fs
pfs <- pData(fs)
pfs$Condition <- factor(sample(c("treated","untreated"),
                               10,
                               replace = TRUE),
                        levels = c("untreated","treated"))
pData(fs) <- pfs
pData(fs)

###################### GatingSet ######################

gs <- GatingSet(fs)
pData(gs)

ggcyto(gs,
       aes(x = A,
           y = B),
       subset = "root")+
    geom_hex()

# row names = genes
# colnames = tumor specimens
# Matrix is binary alteration matrix
alterationRank <- function(dataset){
  gene.names = rownames(dataset)
  alteration.freq = data.frame()
  logical <- (dataset > 0)
  for (i in 1:nrow(logical)) {
    alteration.freq[i,1] = gene.names[i]
    alteration.freq[i,2] = sum(logical[i,], na.rm = TRUE)
    alteration.freq[i,3] = sum(logical[i,] == FALSE, na.rm = TRUE)
    alteration.freq[i,4] = alteration.freq[i,2] / (alteration.freq[i,2] + alteration.freq[i,3])
  }
  alteration.freq <- alteration.freq[order(alteration.freq[,4], decreasing = TRUE),]
  colnames(alteration.freq) = c("Hugo_Symbol", "total altered", 
                       "total unaltered", "alteration freq")
  alteration.freq$`Gene rank` <- 1:nrow(alteration.freq)
  return(alteration.freq)
}
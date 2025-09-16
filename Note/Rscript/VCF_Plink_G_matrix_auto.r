# install.packages("vcfR")
# install.packages("bigsnpr")

library(vcfR)
library(bigsnpr)

# Core G function
calc_G <- function(W) {
  n_ind <- nrow(W)
  p <- colSums(W) / (2 * n_ind)
  P <- matrix(rep(2 * p, each = n_ind), nrow = n_ind)
  M <- W - P
  k <- 2 * sum(p * (1 - p))
  G <- (M %*% t(M)) / k
  rownames(G) <- paste0("Ind", 1:n_ind)
  colnames(G) <- paste0("Ind", 1:n_ind)
  return(G)
}

# Unified function
calc_G_auto <- function(input) {
  if (grepl("\\.vcf$", input, ignore.case = TRUE)) {
    # --- VCF ---
    vcf <- read.vcfR(input)
    gt <- extract.gt(vcf)
    W <- apply(gt, 2, function(col) {
      sapply(col, function(g) {
        if (g %in% c("0/0","0|0")) return(0)
        if (g %in% c("0/1","1/0","0|1","1|0")) return(1)
        if (g %in% c("1/1","1|1")) return(2)
        return(NA)
      })
    })
    for (j in 1:ncol(W)) {
      W[is.na(W[,j]), j] <- mean(W[,j], na.rm=TRUE)
    }
    return(calc_G(W))
  } else if (grepl("\\.bed$", input, ignore.case = TRUE)) {
    # --- PLINK ---
    obj.bed <- snp_attach(snp_readBed(input))
    Gmat <- obj.bed$genotypes[]
    W <- apply(Gmat, 2, function(x) { x[is.na(x)] <- mean(x, na.rm=TRUE); x })
    return(calc_G(W))
  } else {
    stop("Input must be a VCF or PLINK .bed file.")
  }
}

# Example usage:
# G1 <- calc_G_auto("genotypes.vcf")
# G2 <- calc_G_auto("genotypes.bed")

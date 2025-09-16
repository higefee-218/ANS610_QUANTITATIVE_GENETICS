## --- Example: Build the Genomic Relationship Matrix (VanRaden, 2008) ---

## Raw genotype matrix W: rows = individuals, cols = markers
## Genotypes coded as 0, 1, 2 (copies of reference allele)
W <- matrix(c(
  1, 2, 0,
  1, 1, 1,
  2, 2, 0,
  0, 1, 2
), nrow = 4, byrow = TRUE)

rownames(W) <- paste0("Ind", 1:4)
colnames(W) <- paste0("M", 1:3)

cat("Raw genotype matrix W:\n")
print(W)

## Step 1. Calculate allele frequencies
n_ind <- nrow(W)
p <- colSums(W) / (2 * n_ind)   # vector of allele frequencies
cat("\nAllele frequencies p:\n")
print(p)

## Step 2. Centered genotype matrix M = W - 2p
P <- matrix(rep(2 * p, each = n_ind), nrow = n_ind)  # expected values
M <- W - P
cat("\nCentered genotype matrix M:\n")
print(M)

## Step 3. Scaling factor k = 2 * sum(p * (1 - p))
k <- 2 * sum(p * (1 - p))
cat("\nScaling factor k:\n")
print(k)

## Step 4. Genomic Relationship Matrix G = (M M') / k
G <- (M %*% t(M)) / k
cat("\nGenomic Relationship Matrix G:\n")
print(round(G, 3))

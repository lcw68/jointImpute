#' Microbiome data input
#'
#' Reads a file in list format and eliminate invalid taxa and samples (preprocessing)
#'
#' @param file the name of the file which the data are to be read from.
#' @param eliminate_taxa_threshold threshold for deciding unexpressed level to be elminated
#'
#'
#' @return A list of microbiome data after preprocessing
#'
#' @examples
#'
#'
#' humann_2 <- read_microbiome_data("data.bactRPK.marginal.DRNA.ZOE2.rds",eliminate_taxa_threshold=0.1)
#'
#'
#'
#' @export
#'
read_microbiome_data <- function(file, eliminate_taxa_threshold=0.0){
  data_bact <- readRDS(file)
  otu <- data_bact[['otu']]
  meta_data <- data_bact[['meta']]
  #eliminate invalid samples
  invalid_sample_idx <- which(meta_data$exclusion!='included')
  meta_data <- meta_data[-c(invalid_sample_idx),]
  otu <- otu[,-c(invalid_sample_idx),]
  otu_dna <- otu[,,1]
  otu_rna <- otu[,,2]
  #eliminate invalid genes(unexpressed in each samples)
  threshold = dim(meta_data)[1]* eliminate_taxa_threshold
  invalid_gene_idx <- which(rowSums(otu_dna!=0)<threshold | rowSums(otu_rna!=0)<threshold)
  if(length(invalid_gene_idx)>0){
    otu <- otu[-c(invalid_gene_idx),,]
    otu_dna <- otu[,,1]
    otu_rna <- otu[,,2]
  }
  return (list(otu_dna=otu_dna, otu_rna=otu_rna,meta_data=meta_data))
}


#' Ratio Density plot
#'
#' Plot the log ratio of four peak group
#'
#' @param otu_dna DNA OTU table
#' @param otu_rna RNA OTU table
#' @param with_zero indicator whether we should included the double zero peak
#'
#'
#' @return A list of Density plot and double zero proportion
#'
#' @examples
#'
#'
#' humann_2 <- read_microbiome_data("data.bactRPK.marginal.DRNA.ZOE2.rds",eliminate_taxa_threshold=0.1)
#' plot_ratio_density(humann_2[[1]],human_2[[2]])
#'
#' @importFrom ggplot2 ggplot geom_density scale_fill_discrete scale_color_discrete labs
#' @export
#'
plot_ratio_density <- function(otu_dna, otu_rna, with_zero=TRUE){

  epsilon_dna <- min(otu_dna[which(otu_dna>0)])
  epsilon_rna <- min(otu_rna[which(otu_rna>0)])

  ratio <- as.vector(log((otu_rna+epsilon_rna)/(otu_dna+epsilon_dna)))
  id_00 <- which(otu_dna==0 & otu_rna==0)
  id_10 <- which(otu_dna!=0 & otu_rna==0)
  id_01 <- which(otu_dna==0 & otu_rna!=0)
  id_11 <- which(otu_dna!=0 & otu_rna!=0)

  class <- rep(0, length(otu_rna))
  class[id_00] <- 1
  class[id_10] <- 2
  class[id_01] <- 3
  class[id_11] <- 4
  class_fac <- factor(x=class, levels=c(1,2,3,4),labels=c("(0,0)","(+,0)","(0,+)","(+,+)"))

  df <- data.frame(ratio, class=class_fac)
  if(with_zero)
    density_plot <- ggplot(df, aes(ratio, after_stat(count), fill=class, color=class)) + geom_density(
      adjust=1/10, alpha=0.5) + scale_fill_discrete(name = "class (DNA,RNA)") +  scale_color_discrete(name = "class (DNA,RNA)")+labs(x="log-ratio", y='count', title='log-ratio plot with (0,0) situation')
  else
    density_plot <- ggplot(subset(df, class %in% c("(+,0)","(0,+)","(+,+)")), aes(ratio, after_stat(count), fill=class, color=class)) + geom_density(
      adjust=1/10, alpha=0.5) + scale_fill_discrete(name = "class (DNA,RNA)") +  scale_color_discrete(name = "class (DNA,RNA)")+labs(x="log-ratio", y='count', title='log-ratio plot without (0,0) situation')
  zero_ratio <- length(id_00)/length(otu_rna)
  return(list(density_plot, zero_ratio))
}

#' Random Drop-out Matrix Generation
#'
#' For simulating the false zero, use this function to generate random drop-out matrix
#'
#' @param reads reads that decide the dropping probabilities
#' @param dims Dimension of the drop matrix
#' @param ratio proportional constant
#'
#'
#' @return A 0-1 drop-out matrix
#'
#' @examples
#'
#' generate_drop_matrix(meta$reads,dims = c(20,50))
#'
#'
#' @export
#'
generate_drop_matrix <- function(reads, dims, ratio=0.01){
  drop_prob <- exp(-ratio*reads^2)
  drop_matrix <- matrix(0, dims[1], dims[2])
  for(i in 1:dims[2]){
    while(sum(drop_matrix[,i]!=0)==0)
      drop_matrix[,i] <- rbinom(dims[1], 1, 1-drop_prob[i])
  }
  return(drop_matrix)
}

#' Mean Correlation in log transformed
#'
#'
#'
#' @param a matrix
#' @param b matrix with same dimension of a
#'
#'
#' @return mean correlation
#'
#' @examples
#'
#' log_corr(matrix(rnorm(10),2,5),matrix(runif(10),2,5))
#'
#'
#' @noRd
#'
log_corr <- function(a, b){
  return(mean(sapply(1:nrow(a), function(i){cor(log(a[i,]+1), log(b[i,]+1))})))
}

#' Mean Square Error in log transformed
#'
#'
#'
#' @param a matrix
#' @param b matrix with same dimension of a
#'
#'
#' @return mean MSE
#'
#' @examples
#'
#' log_mse(matrix(rnorm(10),2,5),matrix(runif(10),2,5))
#'
#'
#' @noRd
#'
log_mse <- function(a, b){
  return(mean(sapply(1:nrow(a), function(i){(log(a[i,]+1)-log(b[i,]+1))^2})))
}


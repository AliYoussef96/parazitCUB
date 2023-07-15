#' Relative Synonymous Codon Usage (RSCU) PCA
#'
#' PCA for the Relative Synonymous Codon Usage (RSCU) of DNA sequence.
#'
#' @usage rscu.pca(virus.list.rscu, host.rscu, codons.exclude = c("ATG", "TAA", "TAG", "TGA", "TGG"))
#'
#' @param virus.list.rscu  A list of viruses' RSCU values
#' @param host.rscu  A list of host RSCU values
#' @param codons.exclude exclude codons like start and termination codons
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' fasta.host <- read.host(host.fasta)
#'
#' # Calculate RSCU
#' rscu.virus <- RSCU.values(fasta.virus)
#' rscu.host <- RSCU.values(fasta.host)
#'
#' rscu.pca(rscu.virus, rscu.host)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'




rscu.pca <- function(virus.list.rscu, host.rscu, codons.exclude = c("ATG", "TAA", "TAG", "TGA", "TGG")){

  rsci.df <- data.frame(none = rep("none",64))
  for (i in 1:length(virus.list.rscu)){
    i.df <- virus.list.rscu[[i]]
    i.df <- data.frame(apply(i.df, 2, mean))
    colnames(i.df) <- names(virus.list.rscu)[i]
    rsci.df <- cbind(rsci.df, i.df )
  }

  rsci.df$none <- NULL

  i.df <- data.frame(apply(host.rscu[[1]], 2, mean))
  colnames(i.df) <- "Host"
  rsci.df <- cbind(rsci.df,i.df)
  rsci.df <- rsci.df[!row.names(rsci.df) %in%  tolower(codons.exclude),]

  ##PCA
  res.pca <- prcomp(t(rsci.df), scale.= T, center = T)

  pca_data <- data.frame(sample=rownames(res.pca$x),#1 column with sample ids
                         x=res.pca$x[,1],#2 columns for the X and Y
                         y=res.pca$x[,2])#coordinates for each sample

  pca.var <- res.pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

  pca_data$group <- ifelse(pca_data$sample == "Host" , "Host", "Virus")

  pca.plot <- ggplot(data=pca_data,aes(x = x , y = y, color= group)) +

    geom_point() +

    xlab(paste("PC1: ",pca.var.per[1],"%",sep="")) +

    ylab(paste("PC2: ",pca.var.per[2],"%",sep="")) +
    theme_bw()

  return(pca.plot)
}

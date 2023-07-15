#' Relative Synonymous Codon Usage (RSCU) cluster using factoextra package
#'
#' PCA for the Relative Synonymous Codon Usage (RSCU) of DNA sequence.
#'
#' @usage rscu.cluster(virus.list.rscu, host.rscu, FUNcluster = "kmeans",
#'          k = 3, hc_metric = "euclidean" , hc_method = "ward.D2", rank = 4,
#'          codons.exclude = c("ATG", "TAA", "TAG", "TGA", "TGG"))
#'
#' @param virus.list.rscu  A list of viruses' RSCU values
#' @param host.rscu  A list of host RSCU values
#' @param FUNcluster a clustering function including "kmeans", "pam", "clara", "fanny", "hclust", "agnes" and "diana". Abbreviation is allowed.
#' @param K the number of clusters to be generated. If NULL, the gap statistic is used to estimate the appropriate number of clusters. In the case of kmeans, k can be either the number of clusters, or a set of initial (distinct) cluster centers.
#' @param hc_metric character string specifying the metric to be used for calculating dissimilarities between observations. Allowed values are those accepted by the function dist() [including "euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski"] and correlation based distance measures ["pearson", "spearman" or "kendall"]. Used only when FUNcluster is a hierarchical clustering function such as one of "hclust", "agnes" or "diana".
#' @param hc_method the agglomeration method to be used (?hclust): "ward.D", "ward.D2", "single", "complete", "average", ...
#' @param rank a number specifying the number of PCAs used in the clustering the data, i.e., maximal number of principal components to be used.
#' @param codons.exclude exclude codons like start and termination codons
#'
#'
#'
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @import factoextra
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
#' rscu.cluster(rscu.virus, rscu.host)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

rscu.cluster <- function(virus.list.rscu, host.rscu, FUNcluster = "kmeans",
                         k = 3, hc_metric = "euclidean" , hc_method = "ward.D2", rank = 4,
                         codons.exclude = c("ATG", "TAA", "TAG", "TGA", "TGG")){

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
  pca <- prcomp(t(rsci.df),center= F, scale.= TRUE, rank. = rank)
  results <- pca$x

  km1<-eclust(results, FUNcluster = FUNcluster, hc_metric= hc_metric, k= k,
              hc_method = hc_method)
  cluster.plot <- fviz_cluster(km1, labelsize = 10, main = "") + theme_bw()

  return(list(cluster.plot, data.frame(km1[["cluster"]])))
}

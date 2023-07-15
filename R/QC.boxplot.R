#' QC boxplot
#'
#' boxplot for the amino acid sequence length as a QC
#'
#' @usage QC.boxplot(list.virus)
#'
#' @param list.virus  A list of viruses' DNA returned from read.virus function
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' QC.boxplot(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'


QC.boxplot <- function(list.virus){

  calc.len <- function(seq){

    seq <- as.character(seq)
    return(floor(str_count(seq)/3))
  }

  df.len.all <- data.frame()
  for(i.fasta in 1:length(list.virus)){
    dna <- list.virus[[i.fasta]]
    df.len <- data.frame(sapply(dna,calc.len))
    df.len.all <- rbind(df.len.all, df.len)
  }
  colnames(df.len.all) <- "sequence length"

  pp <- ggplot(df.len.all, aes(y = `sequence length`))+
    geom_boxplot(outlier.colour = "red") + theme_classic()


  return(pp)


}


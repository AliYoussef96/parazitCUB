#' QC.cutoff
#'
#' remove the genes with maximum and minimum threshold amino acid sequence length
#'
#' @usage QC.cutoff(list.virus, cut.off.up, cut.off.down)
#'
#' @param list.virus  A list of viruses' DNA returned from read.virus function
#' @param cut.off.up  The maximum threshold of amino acid sequence length
#' @param cut.off.down The minimum threshold of amino acid sequence length
#' 
#' @return modified list.virus
#'
#' @import seqinr
#' 
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' fasta.virus.mod <- QC.cutoff(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'


QC.cutoff <- function(list.virus, cut.off.up, cut.off.down){
  
  calc.len <- function(seq){
    
    seq <- as.character(seq)
    return(floor(str_count(seq)/3))
  }
  
  
  for(i.fasta in 1:length(list.virus)){
    dna <- list.virus[[i.fasta]]
    df.len <- data.frame(sapply(dna,calc.len))
    df.len$ids <- row.names(df.len)
    df.len <- df.len[df.len$sapply.dna..calc.len. <= cut.off.up,]
    df.len <- df.len[df.len$sapply.dna..calc.len. >= cut.off.down,]

    list.virus[[i.fasta]] <- dna[df.len$ids]
    
  }
  
  return(list.virus)
  
}
  
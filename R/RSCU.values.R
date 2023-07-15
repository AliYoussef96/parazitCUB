#' Relative Synonymous Codon Usage (RSCU)
#'
#' Measure the Relative Synonymous Codon Usage (RSCU) of DNA sequence.
#'
#' For more information about ENc \href{https://academic.oup.com/nar/article-abstract/14/13/5125/1143812?redirectedFrom=fulltext}{Sharp et al., 1986}.
#'
#' @usage RSCU.values(df.fasta)
#'
#' @param df.fasta  a data frame with seq_name and its DNA sequence.
#'
#' @return A data.frame containing the computed RSCU values for each codon for each DNA sequences within df.fasta.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate RSCU 
#' rscu.virus <- RSCU.values(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

RSCU.values <- function(df.fasta) {
  
  rscu.calc <- function(sequence){
    rscu <- uco(s2c( as.character(sequence)), index = "rscu", NA.rscu = 0)
    return(rscu)
  }
  
  rscu.list <- list()
  
  for(i.fasta in 1:length(df.fasta)){
    dna <- df.fasta[[i.fasta]]
    
    rscu <- data.frame(t(sapply(dna,rscu.calc)))
    
    tryCatch({
      rscu.list[dna@virus.name] <- list(rscu)},
    error = function(e){
      rscu.list[dna@host.name] <<- list(rscu)}
    )
    
    
  }
  
  return(rscu.list)
}

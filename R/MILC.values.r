#' Measure Independent of Length and Composition  (MILC).
#'
#' For more information about MILC \href{https://pubmed.ncbi.nlm.nih.gov/16029499/}{Supek & Vlahovicek, 2005}.
#'
#' @usage MILC.values <- function(df.fasta, genetic.code = "1", threshold = 0, alt.init = TRUE , stop.rm = FALSE, filtering = "none", host = "none")
#'
#' @param df.fasta   A list returned from read.virus function or read.host function
#' @param genetic.code  A single string that uniquely identifies the genetic code to extract. Should be one of the values in the id or name2 columns of GENETIC_CODE_TABLE.
#' @param threshold Optional numeric, specifying sequence length, in codons, used for filtering
#' @param stop.rm Logical, whether to remove stop codons. Default is FALSE.
#' @param alt.init Logical, whether to use alternative initiation codons. Default is TRUE.
#' @param filtering Character vector, one of c("none", "soft", "hard"). Specifies whether sequences shorther than some threshold value of length (in codons), len.threshold, should be excluded from calculations. If "none" (default), length of sequences is not checked, if "soft", a warrning is printed if there are shorter sequences, and if "hard", these sequences are excluded from calculation.
#' @param host MILC could be calculated without refrence set host = "none". If host = to a list of dataframe of host genes returned from read.host function it will be calculated using the host genes as a refrence set. 
#'
#' @return A list of data.frame(s) containing the computed MILC values for each DNA sequences.
#'
#' @import coRdon
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate enc version Novembre (2002)
#' milc <- MILC.values(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'


MILC.values <- function(df.fasta, genetic.code = "1", 
                        threshold = 0, alt.init = TRUE , stop.rm = FALSE ,
                        filtering = "none", host = "none" ) {
  
  milc.list <- list()
  
  for(i.fasta in 1:length(df.fasta)){
    dna <- df.fasta[[i.fasta]]
    
    cT <- codonTable(dna)
    
    if( host == "none"){
    
    milc <- MILC(cT,id_or_name2 = genetic.code,
                     alt.init = alt.init, stop.rm = stop.rm,
                     filtering = filtering, len.threshold = threshold,
                 ribosomal = FALSE, self = TRUE)
    
    milc <- data.frame(milc, row.names = names(dna))
    
    colnames(milc) <- "MILC"
    
    }else{
      
      ct.host <- codonTable(host[[1]])
      
      milc <- MILC(cT,id_or_name2 = genetic.code,
                   alt.init = alt.init, stop.rm = stop.rm,
                   filtering = filtering, len.threshold = threshold,
                   ribosomal = FALSE, self = FALSE, subsets  = list(ct.host))
      
      milc <- data.frame(milc, row.names = names(dna))
      
      colnames(milc) <- "MILCvsHost"
      
      
    }
    
    
    tryCatch({
      milc.list[dna@virus.name] <- list(milc)},
      error = function(e){
        milc.list[dna@host.name] <<- list(milc)}
    )
    
  }
  
  return(milc.list)
  
}

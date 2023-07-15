#' Effective Number of Codons (ENc).
#'
#' Measure the Effective Number of Codons (ENc) of DNA sequence. Using version Wright (1990).
#'
#' For more information about ENc \href{https://pubmed.ncbi.nlm.nih.gov/2110097/}{Wright, 1990}.
#'
#' @usage ENc.values.old <- function(df.fasta, genetic.code = "1", threshold = 0, alt.init = TRUE , stop.rm = FALSE , filtering = "none")
#'
#' @param df.fasta   A list returned from read.virus function or read.host function
#' @param genetic.code  A single string that uniquely identifies the genetic code to extract. Should be one of the values in the id or name2 columns of GENETIC_CODE_TABLE.
#' @param threshold Optional numeric, specifying sequence length, in codons, used for filtering
#' @param stop.rm Logical, whether to remove stop codons. Default is FALSE.
#' @param alt.init Logical, whether to use alternative initiation codons. Default is TRUE.
#' @param filtering Character vector, one of c("none", "soft", "hard"). Specifies whether sequences shorther than some threshold value of length (in codons), len.threshold, should be excluded from calculations. If "none" (default), length of sequences is not checked, if "soft", a warrning is printed if there are shorter sequences, and if "hard", these sequences are excluded from calculation.
#' 
#' @return A list of data.frame(s) containing the computed ENc values for each DNA sequences.
#'
#' @import coRdon
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate enc version Wright (1990)
#' enc.wright <- ENc.values.old(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

ENc.values.old <- function(df.fasta, genetic.code = "1", threshold = 0, 
                           alt.init = TRUE , stop.rm = FALSE, filtering = "none") {
  
  
  enc.list <- list()
  
  for(i.fasta in 1:length(df.fasta)){
    dna <- df.fasta[[i.fasta]]
    
    cT <- codonTable(dna)
    
    enc <- ENC(cT,
                    id_or_name2 = genetic.code,
                    alt.init = alt.init, stop.rm = stop.rm,
                    filtering = filtering, len.threshold = threshold)
    
    enc <- data.frame(enc, row.names = names(dna))
    
    
    
    tryCatch({
      enc.list[dna@virus.name] <- list(enc)},
      error = function(e){
        enc.list[dna@host.name] <<- list(enc)}
    )
    
  }
  
  return(enc.list)
  
}

#####################################
#############################
#####################################


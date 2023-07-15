#' Statistical dinucleotide over- and underrepresentation (codon model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of codons.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.codon(df.virus,permutations=100,exact_numbers = FALSE)
#'
#' @param df.virus  a list returned from read.virus function
#' @param permutations  the number of permutations for the z-score computation.  Default is 100
#' @param exact_numbers if TRUE exact analytical calculation will be used.  Default is FALSE
#'
#' @return A list of data.frame(s) containing the computed statistic for each dinucleotide in all DNA sequences within df.virus.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate zscore using (base model)
#' base <- dinuc.codon(fasta.virus, permutations = 100)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}

dinuc.codon <- function(df.virus, permutations = 100, exact_numbers = FALSE) {
  
  dinuc.codon.calc <- function(sequence){
    
    sequence <- tolower(sequence)
    
    codon.result <- zscore(s2c( as.character(sequence)), simulations = permutations, modele = "codon", 
                          exact = exact_numbers)
    
    return(codon.result)
  }
  
  dinuc.codon.list <- list()
  
  for(i.fasta in 1:length(df.virus)){
    dna <- df.virus[[i.fasta]]
    
    dinuc.codon <- data.frame(sapply(dna, dinuc.codon.calc))
    
    dinuc.codon.list[dna@virus.name] <- list(dinuc.codon)
  }
  
  return(dinuc.codon.list)
  
}




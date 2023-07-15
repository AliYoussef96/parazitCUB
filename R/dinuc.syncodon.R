#' Statistical dinucleotide over- and underrepresentation (syncodon model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of synonymous codons.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.syncodon(df.virus,permutations=100,exact_numbers = FALSE)
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
#' base <- dinuc.syncodon(fasta.virus, permutations = 100)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}


dinuc.syncodon <- function(df.virus, permutations = 100, exact_numbers = FALSE) {
 
  dinuc.syncodon.calc <- function(sequence){
    
    sequence <- tolower(sequence)
    
    syncodon.result <- zscore(s2c( as.character(sequence)), simulations = permutations, modele = "syncodon", 
                           exact = exact_numbers)
    
    return(syncodon.result)
  }
  
  dinuc.syncodon.list <- list()
  
  for(i.fasta in 1:length(df.virus)){
    dna <- df.virus[[i.fasta]]
    
    dinuc.syncodon <- data.frame(sapply(dna, dinuc.syncodon.calc))
    
    dinuc.syncodon.list[dna@virus.name] <- list(dinuc.syncodon)
  }
  
  return(dinuc.syncodon.list)
  
}
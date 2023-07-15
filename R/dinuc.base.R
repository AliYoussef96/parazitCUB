#' Statistical dinucleotide over- and underrepresentation (base model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of all bases in the sequence.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.base(df.virus,permutations=100,exact_numbers = FALSE)
#'
#' @param df.virus  a list returned from read.virus function
#' @param permutations  the number of permutations for the z-score computation.  Default is 100
#' @param exact_numbers if TRUE exact analytical calculation will be used.  Default is FALSE
#'
#' @return A list of data.frame(s) containing the computed statistic for each dinucleotide in all DNA sequences.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate zscore using (base model)
#' base <- dinuc.base(fasta.virus, permutations = 100)
#' }
#' @export
#' 
#' @note If there are NAs values in the data it will by imputed by  1/5 the lowest value in the data.
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}


dinuc.base <- function(df.virus, permutations = 100, exact_numbers = FALSE) {
  
  dinuc.base.calc <- function(sequence){
    
    sequence <- tolower(sequence)
    
    base.result <- zscore(s2c( as.character(sequence)), simulations = permutations, modele = "base", 
           exact = exact_numbers)
    
    return(base.result)
  }
  
  dinuc.base.list <- list()
  
  for(i.fasta in 1:length(df.virus)){
    dna <- df.virus[[i.fasta]]
    
    dinuc.base <- data.frame(sapply(dna, dinuc.base.calc))
    
    dinuc.base.list[dna@virus.name] <- list(dinuc.base)
  }
  
  return(dinuc.base.list)
  
}
  
  
 

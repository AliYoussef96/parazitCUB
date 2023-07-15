#' GC content
#'
#' Calculates overall GC content as well as GC at first, second, and third codon positions.
#'
#' @usage GC.content(df.virus)
#'
#' @param df.virus  a list returned from read.virus function
#'
#' @return A list of data.frame(s) with overall GC content as well as GC at first, second, and third codon positions of all DNA sequence.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate GC content
#' gc.df <- GC.content(fasta.virus)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}

GC.content <- function(df.virus) {
  
  GC.overall <- function(sequence){
    GC <- GC(s2c( as.character(sequence)) )
    return(GC)
  }
  
  GC1.calc <- function(sequence){
    GC1 <- GCpos(s2c( as.character(sequence)) , "1")
    return(GC1)
  }

  GC2.calc <- function(sequence){
    GC2 <- GCpos(s2c( as.character(sequence)) , "2")
    return(GC2)
  }
  
  GC3.calc <- function(sequence){
    GC3 <- GCpos(s2c( as.character(sequence)) , "3")
    return(GC3)
  }
  
  GC.list <- list()
  
  for(i.fasta in 1:length(df.virus)){
    dna <- df.virus[[i.fasta]]
    
    GC <- data.frame(sapply(dna,GC.overall))
    colnames(GC) <- "GC.overall"
    
    GC1 <- data.frame(sapply(dna,GC1.calc))
    colnames(GC1) <- "GC1"
    
    GC2 <- data.frame(sapply(dna,GC2.calc))
    colnames(GC2) <- "GC2"
    
    GC3 <- data.frame(sapply(dna,GC3.calc))
    colnames(GC3) <- "GC3"
    
    GC.result <- cbind.data.frame(GC,GC1,GC2, GC3)
    
    GC.list[dna@virus.name] <- list(GC.result)
    
  }
  

  
  return(GC.list)
  
  
}



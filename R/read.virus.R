#' Read fasta formate and convert it to data frame
#'
#' @usage read.virus(virus.fastas, sep = "|")
#'
#' @param virus.fastas  a Vector of directory path to the virus fasta file(s).
#' @param sep for long sequence name separate the fasta header with?
#'
#' @return A list with DNAStringSet.
#'
#' @import stringr
#' @import coRdon
#' @import stringr
#' @importFrom  Biostrings readDNAStringSet
#'
#' @examples
#' \dontrun{
#' fastas.virus <- read.virus(c("File1.fasta", "File2.fasta", "File3.fasta") , sep = "|")
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

read.virus <- function(virus.fastas, sep  = "|"){
  list.virus <- ""
  count <- 0
  for(i.fasta in virus.fastas){
    count <- count  + 1
    
    virus.file <- readSet(file = i.fasta)
    
    attributes(virus.file)$virus.name <- str_remove_all(basename(i.fasta), fixed(".fasta"))
    
    sep.names <- data.frame(str_split_fixed(names(virus.file) , fixed(sep), 2 ))
    
    names(virus.file) <- sep.names$X1
    
    list.virus[count] <- list(virus.file)
  }
  return(list.virus)
}

#' Read fasta formate and convert it to data frame
#'
#' @usage read.virus(host.fasta , sep = "|")
#'
#' @param host.fasta  a directory path to the host fasta file.
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
#' fastas.host <- read.host(c("File1.fasta", "File2.fasta", "File3.fasta"), sep = "|")
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

read.host <- function(host.fasta, sep = "|"){
  
  if (length(host.fasta) != 1){
    stop("Only one host could be used for the analysis!!")
  }

  list.host <- ""
  count <- 0
  for(i.fasta in host.fasta){
    count <- count  + 1
    
    host.file <- readSet(file = i.fasta)
    
    attributes(host.file)$host.name <- str_remove_all(basename(i.fasta), fixed(".fasta"))
    
    sep.names <- data.frame(str_split_fixed(names(host.file) , fixed(sep), 2 ))
    
    names(host.file) <- sep.names$X1
    
    list.host[count] <- list(host.file)
  }
  return(list.host)
}

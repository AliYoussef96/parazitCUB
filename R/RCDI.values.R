#' Relative Codon Deoptimization Index (RCDI)
#'
#' Measure the Relative Codon Deoptimization Index (RCDI) of DNA sequence.
#'
#' For more information about RCDI \href{https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-87}{Puigbò et al., 2010}
#'
#' @usage RCDI.calc(list.virus, list.host , rscu.host, enc.host, set.length = 5)
#'
#' @param list.virus  a list of data frame(s) with virus seq_name and its DNA sequence.
#' @param list.host  a list of data frame with host seq_name and its DNA sequence.
#' @param rscu.host   a list of data frame of a hosts' RSCU values.
#' @param enc.host   a list data frame of a hosts' ENc values.
#' @param set.len  a number represents a percent that will be used as reference genes from the total host genes.
#'
#' @return A list contains data.frame containing the computed RCDI values for each DNA sequences of each virus.
#'
#' @importFrom  Biostrings DNAStringSet
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' df.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' df.host <- read.host(c("host.fasta") , sep = " [locus")
#' # Calculate RSCU 
#' rscu.virus <- RSCU.values(fasta.virus)
#' rscu.host <- RSCU.values(df.host)
#' enc.novembre.host <- ENc.values.old(df.host)
#' RCDI.calc(df.virus , df.host, rscu.host, enc.novembre.host, set.length = 100)
#' 
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'
#'


RCDI.calc <- function(list.virus, list.host , rscu.host, enc.host, set.length = 5){
  
  RCDI.values <- function(fasta.virus, fasta.host , rscu.host , 
                          enc.host, set.len = 5) {
    

    fasta.virus <- data.frame(fasta.virus)
    fasta.host <- data.frame(fasta.host)
    
    colnames(fasta.virus)[1] <- "sequence"
    colnames(fasta.host)[1] <- "sequence"
    
    fasta.virus$seq_name <- row.names(fasta.virus)
    fasta.host$seq_name <- row.names(fasta.host)
    
    
    colnames(enc.host) <- "ENc"
    enc.host$gene.name <- row.names(enc.host)
    
    newENc <- enc.host[order(enc.host$ENc), ]
    set.len <- length(newENc$gene.name) * (set.len / 100)
    gene.set <- newENc$gene.name[1:set.len]
    gene.set <- fasta.host[fasta.host$seq_name %in% gene.set, ]
    
    #rscu.virus <- RSCU.values(fasta.virus)
    
    #rscu.ref <- RSCU.values(gene.set)
    
    #rscu.host$codons <- row.names(rscu.host)
    rscu.ref <- rscu.host[row.names(rscu.host) %in% gene.set$seq_name,]
    
    #rscu.ref <- data.frame(t(rscu.ref))
    #rscu.ref <- rscu.ref[!row.names(rscu.ref) %in% "codons",]
    
    df.rscu.ref <- data.frame(apply(rscu.ref, 2, mean, na.rm = TRUE))
    colnames(df.rscu.ref) <-  "rscu.ref"
    
    # df.rscu.ref <- data.frame()
    # length <- 1:length(rscu.ref)
    # for (i_mean in length) {
    #   df.rscu <- NULL
    #   means <- mean( as.numeric(rscu.ref[[i_mean]]), na.rm = TRUE)
    #   codon <- colnames(rscu.ref)
    #   codon <- codon[i_mean]
    #   df.rscu <- data.frame(codon = codon, rscu.ref = means)
    #   df.rscu.ref <- rbind(df.rscu.ref, df.rscu)
    # }
    
    codon.name <- row.names(df.rscu.ref)
    #df.rscu.ref <- as.data.frame(t(df.rscu.ref))
    #colnames(df.rscu.ref) <- codon.name
    #df.rscu.ref <- df.rscu.ref[-c(1), ]
    df.rscu.ref <- data.frame(t(df.rscu.ref))
    
    
    
    RCDI.df <- data.frame()
    length <- 1:length(fasta.virus$seq_name)
    
    for (i_seq in length) {
      sequence <- as.character(fasta.virus$sequence[[i_seq]])
      firstframe <- function(sequence) {
        sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
        return(sequence)
      }
      sequence <- firstframe(sequence)
      seq_name <- as.character(fasta.virus$seq_name[[i_seq]])
      
      rscu <- uco(s2c(sequence),
                  index = "rscu",
                  as.data.frame = FALSE, NA.rscu = 0
      )
      
      rscu <- as.data.frame(t(rscu))
      rownames(rscu) <- "rscu"
      
      
      dna <- DNAStringSet(c(sequence, "NNN"))
      count <- codonTable(dna)
      count <- as.data.frame(count@counts[1, ])
      count <- as.data.frame(t(count))
      colnames(count) <- tolower(colnames(count))
      rownames(count) <- "count"
      
      CiFa <- rbind(rscu, count)
      CiFa.CiFh <- rbind(CiFa, df.rscu.ref)
      CiFa.CiFh <- as.data.frame(t(CiFa.CiFh))
      N <- as.numeric(floor(nchar(sequence) / 3))
      
      
      #CiFa.CiFh$rscu <- as.numeric(CiFa.CiFh$rscu)
      #CiFa.CiFh$rscu.ref <- as.numeric(CiFa.CiFh$rscu.ref)
      #CiFa.CiFh$count <- as.numeric(CiFa.CiFh$count)
      
      CiFa.CiFh$RCDI <- ((CiFa.CiFh$rscu / CiFa.CiFh$rscu.ref) * CiFa.CiFh$count) / N
      RCDI <- sum(CiFa.CiFh$RCDI, na.rm = T)
      
      df <- NULL
      df <- data.frame(gene.name = seq_name, RCDI = RCDI)
      RCDI.df <- rbind(RCDI.df, df)
    }
    return(RCDI.df)
  }
  
  
  ##################################
  ##########################
  ##################################
  
  RCDI.list.result <- list()
  
  list.host <- list.host[[1]]
  rscu.host <- rscu.host[[1]]
  enc.host <- enc.host[[1]]
  
  for(i.fasta in 1:length(list.virus)){
    dna <- list.virus[[i.fasta]]

    rcdi.result <- RCDI.values(fasta.virus = dna, 
                               fasta.host = list.host,
                               rscu.host = rscu.host, 
                               enc.host = enc.host, 
                               set.len = set.length)
  
    RCDI.list.result[dna@virus.name] <- list(rcdi.result)
    
  }
  
  return(RCDI.list.result)
}
  
  
  

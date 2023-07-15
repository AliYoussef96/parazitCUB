#' Similarity Index (SiD)
#'
#' Measure the Similarity Index (SiD) between a virus and its host codon usage.
#'
#' For more information about SiD \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077239}{Zhou et al., 2013}.
#'
#'
#' @usage SiD.list(RSCU.virus,RSCU.host)
#'
#' @param RSCU.virus  a list of data.frame(s) with RSCU a virus codon values.
#' @param RSCU.host  a list of a data.frame with RSCU a host codon values.
#'
#' @return A data.frame of numeric represent a SiD value.
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' fasta.host <- read.host(c("host.fasta"))
#' # Calculate RSCU 
#' rscu.virus <- RSCU.values(fasta.virus)
#' rscu.host <- RSCU.values(fasta.host)
#' # Calculate Sid
#' SiD.list(rscu.virus, rscu.host)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'


SiD.list <- function(RSCU.virus, RSCU.host){
  
  
  
  SiD.value <- function(rscu.host, rscu.virus) {
    
    df.rscu.host <- data.frame()
    length <- 1:length(rscu.host)
    for (i_mean in length) {
      df.rscu <- NULL
      means <- mean(rscu.host[,i_mean], na.rm = TRUE)
      codon <- colnames(rscu.host)
      codon <- codon[i_mean]
      df.rscu <- data.frame(codon = codon, rscu.host = means)
      df.rscu.host <- rbind(df.rscu.host, df.rscu)
    }
    
    df.rscu.virus <- data.frame()
    length <- 1:length(rscu.virus)
    for (i_mean in length) {
      df.rscu <- NULL
      means <- mean(rscu.virus[[i_mean]], na.rm = TRUE)
      codon <- colnames(rscu.virus)
      codon <- codon[i_mean]
      df.rscu <- data.frame(codon = codon, rscu.virus = means)
      df.rscu.virus <- rbind(df.rscu.virus, df.rscu)
    }
    
    rscu.df.all <- merge(df.rscu.host, df.rscu.virus, by = "codon")
    rscu.df.all$rscu.all <- rscu.df.all$rscu.host * rscu.df.all$rscu.virus
    
    up <- sum(rscu.df.all$rscu.all)
    down <- sqrt((sum(rscu.df.all$rscu.host)^2) * (sum(rscu.df.all$rscu.virus)^2))
    R.a.b <- up / down
    
    D.a.b <- (1 - R.a.b) / 2
    
    return(D.a.b)
  }
  
  #######################
  ################
  ######################
  
  #RSCU.host <- data.frame(t(RSCU.host[[1]]))
  RSCU.host <- data.frame(RSCU.host[[1]])
  
  all.sid.result <- data.frame()
  
  for(i.fasta in 1:length(RSCU.virus)){
    #df.rscu.virus <- data.frame(t(RSCU.virus[[i.fasta]]))
    
    df.rscu.virus <- data.frame(RSCU.virus[[i.fasta]])
    
    
    sid.result <- SiD.value(RSCU.host , df.rscu.virus)
    
    sid.df <- data.frame("sid" = sid.result )
    
    row.names(sid.df) <- names(RSCU.virus[i.fasta])
    
    all.sid.result <- rbind(all.sid.result , sid.df)
    
    
  }
  
  return(all.sid.result)
  
}

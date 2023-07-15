#' Parity rule 2 (PR2) plot
#'
#' Make a Parity rule 2 (PR2) plot, where the AT-bias [A3/(A3 +T3)] at the third codon position of the four-codon amino acids of entire genes is the ordinate and the GC-bias [G3/(G3 +C3)] is the abscissa. The center of the plot, where both coordinates are 0.5, is where A = U and G = C (PR2), with no bias between the influence of the mutation and selection rates.
#'
#' For more information about PR2 plot \href{https://www.tandfonline.com/doi/full/10.1038/emi.2016.106}{Butt et al., 2016}.
#'
#' @usage PR2.plot.group(df.list, fasta.df, fold4 = TRUE)
#'
#' @param fasta.df  a DNAStringSet object contain only one data frame.
#' @param fold4 if TRUE only 4-fold amino acids will be included in the analysis.
#' @param df.list a list returned from read.virus function
#' 
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' df <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta") , sep = "_")
#'
#' PR2.plot.group(df, fold4 = TRUE)
#' }
#' @export
#'
#' @author AAli Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'

PR2.plot.group <- function(df.list , fold4 = TRUE) {

  
  pr2.analysis <- function(fasta.df, v.name){
  
  fasta.df <- data.frame(fasta.df)
  fasta.df$seq_name <- names(fasta.df)
  colnames(fasta.df)[1] <- "sequence"
  
  codons.Fold4 <- tolower( c("GCA",
                             "GCC",
                             "GCG",
                             "GCT",
                             "CGA",
                             "CGC",
                             "CGG",
                             "CGT",
                             "GGA",
                             "GGC",
                             "GGG",
                             "GGT",
                             "CTA",
                             "CTC",
                             "CTG",
                             "CTT",
                             "CCA",
                             "CCC",
                             "CCG",
                             "CCT",
                             "TCA",
                             "TCC",
                             "TCG",
                             "TCT",
                             "ACA",
                             "ACC",
                             "ACG",
                             "ACT",
                             "GTA",
                             "GTC",
                             "GTG",
                             "GTT") ) 
  
  freq.nt.all <- data.frame()
  length <- 1:length(fasta.df$seq_name)
  
  for (i_seq in length) {
    sequence <- tolower(as.character(fasta.df$sequence[[i_seq]]))
    seq_name <- as.character(fasta.df$seq_name[[i_seq]])
    
    sequence <- s2c(sequence)
    
    if (fold4 == TRUE){
      l <- 3 * ((length(sequence) - 0)%/%3)
      c1 <- seq(from = 0 + 1, to = 0 + l, by = 3)
      sequence.new <- c()
      for(i in c1){
        ii <- i + 2
        codons <- sequence[i:ii]
        codons <- paste0(codons , collapse = "")
        if (codons %in% codons.Fold4){
          sequence.new <- c(sequence.new , codons)
        }
      }
      sequence.new <- paste0(sequence.new, collapse = "")
      sequence <- s2c(sequence.new)
    }
    
    freq.nt <- count(sequence, wordsize = 1, by = 3, start = 2)
    freq.nt <- as.data.frame(freq.nt)
    col.name <- freq.nt$Var1
    freq.nt <- as.data.frame(t(as.data.frame(freq.nt)))
    colnames(freq.nt) <- col.name
    freq.nt <- freq.nt[-c(1), ]
    rownames(freq.nt) <- seq_name
    freq.nt.all <- rbind(freq.nt.all, freq.nt)
  }
  
  freq.nt.all$a <- as.numeric(freq.nt.all$a)
  freq.nt.all$t <- as.numeric(freq.nt.all$t)
  freq.nt.all$g <- as.numeric(freq.nt.all$g)
  freq.nt.all$c <- as.numeric(freq.nt.all$c)
  
  A3T3 <- NULL
  G3C3 <- NULL
  freq.nt.all$A3T3 <- freq.nt.all$a / (freq.nt.all$a + freq.nt.all$t)
  freq.nt.all$G3C3 <- freq.nt.all$g / (freq.nt.all$g + freq.nt.all$c)
  
  freq.nt.all.mean <- freq.nt.all[,5:6]
  freq.nt.all.mean.df <- data.frame(a3t3.mean = mean(freq.nt.all.mean$A3T3, na.rm = TRUE) , 
                                    g3c3.mean = mean(freq.nt.all.mean$G3C3, na.rm = TRUE))
  
  row.names(freq.nt.all.mean.df) <- v.name
  return(freq.nt.all.mean.df)
  }
  
  pr2.df.all <- data.frame()
  for(i.fasta in 1:length(df.list)){
    dna <- df.list[[i.fasta]]
    v.name <- dna@virus.name
    pr2.df <- pr2.analysis(dna, v.name = v.name)
    pr2.df.all <- rbind(pr2.df.all, pr2.df)
  }
  
  pr2.df.all$virus.name <- row.names(pr2.df.all)
  
  
  plot <- ggplot(pr2.df.all, aes(x = a3t3.mean, y = g3c3.mean, color= virus.name)) + 
    geom_point(size = 4, alpha = 0.8) +
    ylab("A3/(A3 + T3)") + xlab("G3/(G3 + C3)") + 
    ylim(0, 1) + xlim(0, 1) + theme_classic(base_size = 20) +
    geom_hline(yintercept = 0.5, color = "red", size = 1.2 , alpha = 0.5) + 
    geom_vline(xintercept = 0.5, color = "red", size = 1.2, alpha = 0.5)
  
  return(plot)
}


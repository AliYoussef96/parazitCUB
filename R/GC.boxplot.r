#' GC box plot
#'
#' Make a GC box plot
#'
#' @usage GC.boxplot(fasta.list)
#'
#' @param fasta.list  A list of GC content from GC.content function.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#' \dontrun{
#' 
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate GC content
#' gc.df <- GC.content(fasta.virus)
#' GC.boxplot(gc.df)
#' 
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#' 
#'  
#'  



GC.boxplot <- function(GC.list){
  
  all.GC <- data.frame()
  for(i.fasta in 1:length(GC.list)){
    gc.df <- data.frame(GC.list[[i.fasta]])
  
    gc.df$virus <- rep(names(GC.list[i.fasta]), nrow(gc.df))
    
    gc.df <- melt(gc.df)
    
    all.GC <- rbind(all.GC, gc.df)
  }
  
  all.GC$value <- all.GC$value * 100
  
  plot <- ggplot(all.GC , aes(x  = virus, y = value , fill = variable)) + 
    
    geom_boxplot(color="black", alpha = 0.8) + 
    
    scale_fill_brewer(palette="Set1") +

    theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                           colour ="black"),legend.position="top") + 

    theme(panel.background = element_rect(fill = "white",
                                            colour = "white")) + 
    
    theme(axis.text.x = element_text(colour="black"), 
           axis.text.y = element_text(colour="black"), 
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
          text=element_text(size=16, face = "bold", colour = "black")) +
     
    ylab("%") + xlab("Virus(es)") 
  
  
  return(plot)
}


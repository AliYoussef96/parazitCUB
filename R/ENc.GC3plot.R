#' ENc-GC3 scatterplot.
#'
#' Make an ENc-GC3 scatterplot. Where the y-axis represents the ENc values and the x-axis represents the GC3 content.
#' The red fitting line shows the expected ENc values when codon usage bias affected solely by GC3.
#'
#' For more information about ENc-GC3 plot \href{https://www.tandfonline.com/doi/full/10.1038/emi.2016.106}{Butt et al., 2016}.
#'
#' @usage ENc.GC3plot(enc.df, gc.df)
#'
#' @param enc.df  a data frame with ENc values.
#' @param gc.df  a data frame with GC3 values.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' df <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta") , sep = "_")
#'
#' GC.list <- GC.content(df)
#' enc.novembre <- ENc.values.new(df)
#' 
#' ENc.GC3plot(enc.novembre[["test1"]] , GC.list[["test1"]])
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}

ENc.GC3plot <- function(enc.df, gc.df) {
  
  
  colnames(enc.df) <- "ENc"
  
  x <- NULL
  eq <- function(x) {
    2 + x + (29 / (x^2 + (1 - x)^2))
  }
  plot <- ggplot() + geom_point(data = enc.df, aes(x = gc.df$GC3, y = ENc)) +
    stat_function(fun = eq, geom = "line", color = "red", size = 1, data = data.frame(x = c(seq(0, 1, 0.001))), aes(x)) +
    theme_classic() + xlab("GC3") + ylab("ENc")

  return(plot)
}

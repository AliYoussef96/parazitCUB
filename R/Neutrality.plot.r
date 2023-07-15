#' Neutrality plot plot
#'
#' Neutrality plots of GC12 in function of GC3 was generated to assess the extent of mutational pressure on the usage of codons A correlation among GC12 and GC3 indicates mutational forces.
#'
#' For more information about Neutrality plot \href{https://doi.org/10.1016/j.meegid.2020.104471}{Nambou.K and Anakpa.M, 2020}.
#'
#' @usage Neutrality.plot(gc.df, size = 5)
#'
#' @param gc.df  A dataframe contains the GC12 and GC3 calculated using the GC.content() function.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import ggpmisc
#' @import ggpubr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' df <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' GC.list <- GC.content(df)
#' Neutrality.plot(GC.list[["test1"]], size = 5)
#' }
#'
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#'


Neutrality.plot <- function(gc.df, size = 5){

  gc.df$GC12 <- (gc.df$GC1 + gc.df$GC2) / 2

  gc.df <- gc.df * 100

  gc.df <- round(gc.df,1)

  my.formula <- y ~ x


  plot <- ggplot(gc.df, aes(x=GC3, y=GC12)) +
    geom_point(size = 3) +
    geom_smooth(method=lm, se=F, colour="red" , size = 1.3) +
    theme(text=element_text(size=21)) +
    ylab("GC12%") + xlab("GC3%") +
    theme_minimal() + theme_classic() +

    stat_poly_eq(formula = my.formula,
                 aes(label = paste( after_stat(eq.label) , after_stat(rr.label) , sep = "*\", \"*") ) ,
                 parse = TRUE  , label.y = (max(gc.df$GC12 + 65) / 100), label.x = (min(gc.df$GC3 -16) / 100), size = size) +

    theme(text=element_text(size=20 ))  +
    stat_cor(label.y = max(gc.df$GC12), label.x = min(gc.df$GC3 + 1)  , size =size) #this means at 35th unit in the y axis, the r squared and p value will be shown

  return(plot)
}

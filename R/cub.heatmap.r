#' CUB Heatmap
#'
#' Make a GC box plot
#'
#' @usage cub.heatmap <- function(cub.list , path.save = "./" , cellwidth = 45 , 
#'                   cellheight = 30 , fontsize = 25 , height = 22 , 
#'                   width = 30 , clustering_distance_rows = "correlation" , 
#'                   clustering_distance_cols = "correlation" , 
#'                   clustering_method = "ward.D")
#'
#' @param cub.list  A list of dinucleotide or codon usage.
#' @param path.save A file path where to save the picture.
#' @param cellwidth An individual cell width in points. If left as NA, then the values depend on the size of plotting window. Default is 45
#' @param cellheight An individual cell height in points. If left as NA, then the values depend on the size of plotting window. Default is 30
#' @param fontsize A base fontsize for the plot. Default is 25
#' @param height A manual option for determining the output file height in inches. Default is 22
#' @param width manual option for determining the output file width in inches. Default is 30
#' @param clustering_distance_rows distance measure used in clustering rows. Possible values are "correlation" for Pearson correlation and all the distances supported by dist, such as "euclidean", etc. If the value is none of the above it is assumed that a distance matrix is provided. Default is euclidean
#' @param clustering_distance_cols distance measure used in clustering columns. Possible values the same as for clustering_distance_rows. Default is euclidean
#' @param clustering_method clustering method used. Accepts the same values as hclust. Default is ward.D
#' @param draw Draw Heatmap for dinucleotide then draw should be = "di" Default. if RSCU then draw should be = "rscu"
#'
#' @return Invisibly a pheatmap object
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @examples
#' \dontrun{
#' 
#' # read DNA from fasta file
#' fasta.virus <- read.virus(c("test1.fasta", "test2.fasta", "test3.fasta"))
#' # Calculate zscore using (base model)
#' base <- dinuc.base(fasta.virus, permutations = 100)
#' cub.heatmap(base)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{aliali.mostafa99@gmail.com} , Salma Bayoumi {salma.ismail.hamed@gmail.com}
#' 
#'  
#'  
#'  


cub.heatmap <- function(cub.list , path.save = "./heat_map" , cellwidth = 45 , 
                        cellheight = 30 , fontsize = 25 , height = 22 , 
                        width = 30 , clustering_distance_rows = "euclidean" , 
                        clustering_distance_cols = "euclidean" , cluster_rows = FALSE,
                        clustering_method = "ward.D", draw = "di", aver = TRUE){
  
  cub.df.all <- data.frame(none = rep("none" , 16))
  cub.group.all <- data.frame()
  
  for(i.fasta in 1:length(cub.list)){
    cub.df <- data.frame(cub.list[[i.fasta]])
    
    if(draw == "rscu"){
      cub.df <- data.frame(t(cub.df))
    }
    
    cub.df[cub.df == "NaN"] <- NA
    
    cub.df[is.na(cub.df)] <- min(cub.df , na.rm = T) * (1/5)
    
    if(draw == "rscu"){
      if (aver == TRUE){
        
        cub.df <- data.frame(apply(cub.df, 1 , mean, na.rm= TRUE), row.names= row.names(cub.df))
        
        colnames(cub.df) <- names(cub.list)[i.fasta]
        
        cub.group <- data.frame(group = names(cub.list)[i.fasta], row.names = names(cub.list)[i.fasta])
      }else{
      colnames(cub.df) <- paste0(row.names(cub.list[[i.fasta]]) , "_" , names(cub.list[i.fasta]))
      
      cub.group <- data.frame(group = rep(names(cub.list[i.fasta]),nrow(cub.list[[i.fasta]])),
                              
                              row.names = paste0(row.names(cub.list[[i.fasta]]) , "_" , names(cub.list[i.fasta])))
      }
      
    }else{
      
      
      if (aver == TRUE){
        
        cub.df <- data.frame(apply(cub.df, 1 , mean, na.rm= TRUE), row.names= row.names(cub.df))
        
        colnames(cub.df) <- names(cub.list)[i.fasta]
        
        cub.group <- data.frame(group = names(cub.list)[i.fasta], row.names = names(cub.list)[i.fasta])
      }else{
    colnames(cub.df) <- paste0(names(cub.list[[i.fasta]]) , "_" , names(cub.list[i.fasta]))
    cub.group <- data.frame(group = rep(names(cub.list[i.fasta]),length(names(cub.list[[i.fasta]]))),
                            row.names = paste0(names(cub.list[[i.fasta]]) , "_" , names(cub.list[i.fasta])))
      }
    }
    
    

    cub.df.all <- cbind(cub.df.all, cub.df)
    cub.group.all <- rbind(cub.group.all, cub.group)
    
  }
  
  cub.df.all$none <- NULL
  
  row.names(cub.df.all) <- toupper(row.names(cub.df.all))
  
  if (draw == "di"){
    
  cub.df.all <- data.frame(t(cub.df.all))
  
  
  plot.heat <- pheatmap(cub.df.all , border_color = "black" , 
                        cellwidth = cellwidth, cellheight = cellheight, 
                        clustering_distance_rows = clustering_distance_rows, 
                        clustering_distance_cols = clustering_distance_cols , 
                        clustering_method = clustering_method,
                        cluster_rows = cluster_rows,
                        fontsize = fontsize,  color =  rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
                        angle_col = 90 , na_col = "#DDDDDD",
                        cluster_cols = T ,
                        height =  height, width = width , 
                        filename = paste0( path.save ,".di.jpeg"),
                        annotation_row = cub.group.all,
                        annotation_names_row = F,
                        annotation_colors = NA )
  
  }else if(draw == "rscu"){

    plot.heat <- pheatmap(cub.df.all , border_color = "black" , 
                          cellwidth = cellwidth, cellheight = cellheight, 
                          clustering_distance_rows = clustering_distance_rows, 
                          clustering_distance_cols = clustering_distance_cols , 
                          clustering_method = clustering_method,
                          fontsize = fontsize,  color =  rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
                          angle_col = 90 , na_col = "#DDDDDD",
                          cluster_cols = T , cluster_rows = cluster_rows,
                          height =  height, width = width , 
                          filename = paste0( path.save ,".rscu.jpeg"),
                          annotation_col = cub.group.all,
                          annotation_names_row = F,
                          annotation_names_col = F,
                          annotation_colors = NA )
  }
  
  

  
  return(plot.heat)
}

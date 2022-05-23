library(ComplexHeatmap)
library(circlize)

## Draws a standard proportion heatmap
propHeatmap = function(proportions, mixture.labels=NA, upper=1, file="~/proportion_heatmap.png"){
  ## Create an annotation for our heatmap
  if(length(mixture.labels)>1){
    row_anno = HeatmapAnnotation(
        mixture.type = mixture.labels,
        which="row")
  }else{
    row_anno = NULL
  }

  ## Another way to define annotations
  col_anno = HeatmapAnnotation(
      cell.type = colnames(proportions),
      which="column")

  ## Plot results
  #png(file, width=16, height=16, units='in', res=300)
  heatmap = Heatmap(proportions,
                    top_annotation = col_anno,
                    col = colorRamp2(c(0, upper), c('white', 'red')),
                    cluster_rows = T,
                    cluster_columns = T,
                    show_row_names = F,
                    show_column_names = T,
                    show_row_dend = T,
                    show_column_dend = F,
                    column_title = "Cell type proportions",
                    row_title = "Mixture profiles",
                    na_col = 'white',
                    column_names_gp = gpar(fontsize = 24),
                    column_title_gp = gpar(fontsize = 24),
                    row_title_gp = gpar(fontsize = 24),
                    column_names_max_height = unit(20, "cm"),
                    row_dend_width = unit(3, "cm"),
                    row_km=4,
                    show_heatmap_legend=T,
                    border=T) + row_anno
  draw(heatmap)
}

# Tutorial: Alignment of C57BL/6 hemoteopietic stem cells (HSCs) from young and old mice from Kowalcyzk et al.

This tutorial provides a guided alignment for two groups of cells from [Kowalczyk et al, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4665007/). In this experiment, single cell RNA (scRNA) sequencing profiles were generated from HSCs in young (2-3 mo) and old (>22 mo) C57BL/6 mice. Age related expression programs make a joint analysis of three isolated cell types long-term (LT), short-term (ST) and  multipotent progenitors (MPPs). In this tutorial we demonstrate the unsupervised alignment strategy of `scAlign` described in [Johansen et al, 2018](https://www.biorxiv.org/content/10.1101/504944v2) along with typical analysis utilizing the aligned dataset, and show how `scAlign` can identify and match cell types across age without using the labels as input.

## Alignment goals
The following is a walkthrough of a typical alignment problem for `scAlign` and has been designed to provide an overview of data preprocessing, alignment and finally analysis in our joint embedding space. Here, our primary goals include:

1. Learning a low-dimensional cell state space in which cells group by function and type, regardless of condition (age).
2. Computing a single cell paired differential expression map from paired cell projections.

## Data setup
The gene count matrices used for this tutorial are hosted [here](https://github.com/quon-titative-biology/examples/blob/master/scAlign_paired_alignment/kowalcyzk_gene_counts.rda).

First, we perform a typical scRNA preprocessing step using the `Seurat` package. Then, reduce to the top 3,000 highly variable genes from both datasets to improve convergence and reduce computation time.


```R
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)

## User paths
working.dir = "." #where our data file, kowalcyzk_gene_counts.rda is located
results.dir = "." #where the output should be stored

## Load in data
load(url('https://github.com/quon-titative-biology/examples/raw/master/scAlign_paired_alignment/kowalcyzk_gene_counts.rda'))

## Extract age and cell type labels
cell_age = unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 1))
cell_type = gsub('HSC', '', unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 2)))

## Separate young and old data
young_data = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "young")]
old_data   = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "old")]

## Set up young mouse Seurat object
youngMouseSeuratObj <- CreateSeuratObject(counts = young_data, project = "MOUSE_AGE", min.cells = 0)
youngMouseSeuratObj <- NormalizeData(youngMouseSeuratObj)
youngMouseSeuratObj <- ScaleData(youngMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)

## Set up old mouse Seurat object
oldMouseSeuratObj <- CreateSeuratObject(counts = old_data, project = "MOUSE_AGE", min.cells = 0)
oldMouseSeuratObj <- NormalizeData(oldMouseSeuratObj)
oldMouseSeuratObj <- ScaleData(oldMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)

## Gene selection
youngMouseSeuratObj <- FindVariableFeatures(youngMouseSeuratObj, do.plot = F, nFeature=3000)
oldMouseSeuratObj <- FindVariableFeatures(oldMouseSeuratObj, do.plot = F, nFeature=3000,)
genes.use = Reduce(intersect, list(VariableFeatures(youngMouseSeuratObj),
                                   VariableFeatures(oldMouseSeuratObj),
                                   rownames(youngMouseSeuratObj),
                                   rownames(oldMouseSeuratObj)))
```

## scAlign setup
The general design of `scAlign` makes it agnostic to the input RNA-seq data representation. Thus, the input data can either be
gene-level counts, transformations of those gene level counts or a preliminary step of dimensionality reduction such
as canonical correlates or principal component scores. Here we create the scAlign object from the previously defined
`Seurat` objects and perform both PCA and CCA on the unaligned data.

```R
## Create paired dataset SCE objects to pass into scAlignCreateObject
youngMouseSCE <- SingleCellExperiment(
    assays = list(counts = youngMouseSeuratObj@assays$RNA@counts[genes.use,],
                  logcounts  = youngMouseSeuratObj@assays$RNA@data[genes.use,],
                  scale.data = youngMouseSeuratObj@assays$RNA@scale.data[genes.use,])
)

oldMouseSCE <- SingleCellExperiment(
  assays = list(counts = oldMouseSeuratObj@assays$RNA@counts[genes.use,],
                logcounts  = oldMouseSeuratObj@assays$RNA@data[genes.use,],
                scale.data = oldMouseSeuratObj@assays$RNA@scale.data[genes.use,])
)
```
We now build the scAlign SCE object and compute PCs and/or CCs using Seurat for the assay defined by `data.use`. It is assumed that `data.use`, which is being used for the initial step of dimensionality reduction, is properly normalized and scaled. 
Resulting combined matrices will always be ordered based on the sce.objects list order.

```R
scAlignHSC = scAlignCreateObject(sce.objects = list("YOUNG"=youngMouseSCE, "OLD"=oldMouseSCE),
                                 labels = list(cell_type[which(cell_age == "young")], cell_type[which(cell_age == "old")]),
                                 data.use="scale.data",
                                 pca.reduce = TRUE,
                                 pcs.compute = 50,
                                 cca.reduce = TRUE,
                                 ccs.compute = 15,
                                 project.name = "scAlign_Kowalcyzk_HSC")
```

## Alignment of young and old HSCs
Now we align the young and old cpopulations for multiple input types which are specified by `encoder.data`. `scAlign` returns a 
low-dimensional joint embedding space where the effect of age is removed allowing us to use the complete dataset for downstream analyses such as clustering or differential expression. For the gene level input we also run the decoder procedure which projects each cell into logcount space for both conditions to perform paired single cell differential expressional.

```R
## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                    encoder.data="scale.data",
                    decoder.data="logcounts",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=TRUE,
                    log.dir=file.path(results.dir, 'models','gene_input'),
                    device="CPU")

## Additional run of scAlign with PCA, the early.stopping heuristic terminates the training procedure too early with PCs as input so it is disabled.
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=15000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE),
                    encoder.data="PCA",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=FALSE,
                    log.dir=file.path(results.dir, 'models','pca_input'),
                    device="CPU")

## Additional run of scAlign with CCA
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=5000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=TRUE),
                    encoder.data="CCA",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=FALSE,
                    log.dir=file.path(results.dir, 'models','cca_input'),
                    device="CPU")

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by input labels
set.seed(5678)
gene_plot = PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", perplexity=30)
pca_plot = PlotTSNE(scAlignHSC, "ALIGNED-PCA", title="scAlign-PCA", perplexity=30)
cca_plot = PlotTSNE(scAlignHSC, "ALIGNED-CCA", title="scAlign-CCA", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot, pca_plot, cca_plot, legend, nrow = 1, layout_matrix=cbind(1,1,1,2,2,2,3,3,3,4))

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by dataset.
set.seed(5678)
gene_plot = PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", perplexity=30)
pca_plot = PlotTSNE(scAlignHSC, "ALIGNED-PCA", title="scAlign-PCA", cols=c("red","blue"), labels.use="group.by", perplexity=30)
cca_plot = PlotTSNE(scAlignHSC, "ALIGNED-CCA", title="scAlign-CCA", cols=c("red","blue"), labels.use="group.by", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot, pca_plot, cca_plot, legend, nrow = 1, layout_matrix=cbind(1,1,1,2,2,2,3,3,3,4))
```
![Cell_label](https://github.com/quon-titative-biology/examples/blob/master/scAlign_paired_alignment/figures/combined_plot_alignment_label.png)
![Cell_batch](https://github.com/quon-titative-biology/examples/blob/master/scAlign_paired_alignment/figures/combined_plot_alignment_stim.png)

## Paired differential expression of young and old cells.
Since we have run the decoder for "ALIGNED-GENE" alignment we can also investigate the paired differential of cells with respect to the young and old conditions. The projections are saved as "YOUNG2OLD" and "OLD2YOUNG" in `scAlignHSC` reducededDims. "YOUNG2OLD" indicates the projection of young cells into the expression space of old cells. As a reminder the combined matrices are always in the order of the `sce.objects` list passed to `scAlignCreateObject`.

```R 
library(ComplexHeatmap)
library(circlize)

all_data_proj_2_old = reducedDim(scAlignHSC, "YOUNG2OLD")
all_data_proj_2_young = reducedDim(scAlignHSC, "OLD2YOUNG")
scDExpr = all_data_proj_2_old - all_data_proj_2_young
colnames(scDExpr) <- genes.use

## Identify most variant DE genes
gene_var <- apply(scDExpr,2,var)
high_var_genes <- sort(gene_var, decreasing=T, index.return=T)$ix[1:40]

## Subset to top variable genes
scDExpr = scDExpr[,high_var_genes]

## Cluster each cell type individually
hclust_LT <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "LT"),], method = "euclidean"), method="ward.D2"))
hclust_ST <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "ST"),], method = "euclidean"), method="ward.D2"))
hclust_MPP <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "MPP"),], method = "euclidean"), method="ward.D2"))

clustering <- merge(hclust_ST, hclust_LT)
clustering <- merge(clustering, hclust_MPP)

## Reorder based on cluster groups
reorder = c(which(cell_type == "LT"), which(cell_type == "ST"), which(cell_type == "MPP"))

## Colors for ComplexHeatmap annotation
cell_type_colors = c("red", "green", "blue")
names(cell_type_colors) = unique(cell_type)

ha = rowAnnotation(df = data.frame(type=cell_type[reorder]),
    						   col = list(type = cell_type_colors),
                   na_col = "white")

png(file.path(results.dir,"scDExpr_heatmap.png"), width=28, height=12, units='in', res=150, bg="transparent")
Heatmap(as.matrix(scDExpr[reorder,]),
        col = colorRamp2(c(-0.75, 0, 0.75), c('blue', 'white', 'red')),
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(1, "cm"),
  	    cluster_columns = TRUE,
        cluster_rows = clustering,
        row_dend_reorder=FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        na_col = 'white',
        column_names_gp = gpar(fontsize = 24),
        split=3,
        gap = unit(0.3, "cm")) + ha
dev.off()
```

![PairedDE](https://github.com/quon-titative-biology/examples/blob/master/scAlign_paired_alignment/figures/scDExpr_heatmap.png)

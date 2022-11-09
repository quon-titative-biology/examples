# scProjection benchmarking on spatial transcriptomics dataset (MERFISH)

Benchmarking using data from Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region https://www.science.org/doi/full/10.1126/science.aau5324


```
conda create -y --name scProjection -c conda-forge -c bioconda python=3.8.13
conda activate scProjection
```

## Setup dependencies
```
conda install numpy scipy scikit-learn pandas tensorflow tensorflow-probability
pip3 install scProjection
```


## benchmarking
### spatialDWLS
```
# https://rubd.github.io/Giotto_site/index.html
#library(devtools)  # if not installed: install.packages('devtools')
#library(remotes)  # if not installed: install.packages('remotes')
#remotes::install_github("RubD/Giotto", lib=.libPaths()[2])  # only for Hongru's env
library(Giotto) # spatialDWLS
# merfish
# https://github.com/RubD/spatial-datasets/tree/master/data/2018_merFISH_science_hypo_preoptic
# https://rubd.github.io/Giotto_site/articles/mouse_merFISH_preoptic_region_200909.html
```



## Tangram [tangram_environment.yml]
```
name: tangram
channels:
-   conda-forge
dependencies:
-   python=3.8
-   pip
-   pip:
    -   squidpy>=1.1.0
    -   tangram-sc==0.4.0
```
```
conda env create -f tangram_environment.yml
conda activate tangram
conda install pytorch==1.12.0 torchvision torchaudio cudatoolkit=11.6 -c pytorch -c conda-forge
```

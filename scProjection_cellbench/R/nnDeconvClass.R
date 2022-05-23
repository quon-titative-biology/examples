nnDeconvClass = setClass("nnDeconvResults", slots=c(deconv_data="list", deconv_emb="list", deconv_logp="list", deconv_var="list", deconv_inv_data="list",
                                                    deconv_data_post="list", deconv_emb_post="list", deconv_logp_post="list", deconv_var_post="list", deconv_inv_data_post="list",
                                                    deconv_samples_emb="list", deconv_samples_rec="list",
                                                    component_metrics="list", mixture_metrics="list",
                                                    proportions="list", weights="list", celltypes="character", genes="character", singleCellNames="character", mixtureCellNames="character",
                                                    flags="list", train_index="array", test_index="array"))
## This function converts the python class from a deconvolution model to the R equivalent
convertDeconv = function(pyClassBase, gene.ids, sc.ids, mixture.ids, low.mem=F){
  nnDeconvResults = nnDeconvClass()
  pyClass = pyClassBase[["deconvResults"]]
  for(name in names(pyClass)){
    if(.hasSlot(nnDeconvResults, name)){
      print(paste0("Recording: ", name))
      if(name %in% c("deconv_data", "deconv_emb", "deconv_logp", "deconv_var", "deconv_inv_data") & low.mem == F){
        slot(nnDeconvResults, name)[["component"]] = list(train = pyClass[[name]][['component']][['train']], test = pyClass[[name]][['component']][['test']]);
        slot(nnDeconvResults, name)[["purified"]]  = list(train = pyClass[[name]][['purified']][['train']], test = pyClass[[name]][['purified']][['test']]);
      }else if(name %in% c("deconv_data_post", "deconv_emb_post", "deconv_logp_post", "deconv_var_post", "deconv_inv_data_post") & low.mem == F){
        slot(nnDeconvResults, name)[["component"]] = list(train = pyClass[[name]][['component']][['train']], test = pyClass[[name]][['component']][['test']]);
        slot(nnDeconvResults, name)[["purified"]]  = list(train = pyClass[[name]][['purified']][['train']], test = pyClass[[name]][['purified']][['test']]);
      }else if(name == "component_metrics"){
        for(celltype in pyClassBase[["celltypes"]]){
          slot(nnDeconvResults, name)[[celltype]][["train"]] = data.frame(step = unlist(pyClass[[name]][[celltype]]$step[['train']]),
                                                                          loss = unlist(pyClass[[name]][[celltype]]$loss[['train']]),
                                                                          mse  = unlist(pyClass[[name]][[celltype]]$mse[['train']]),
                                                                          logp = unlist(pyClass[[name]][[celltype]]$log_probability[['train']]),
                                                                          kl   = unlist(pyClass[[name]][[celltype]]$kl_divergence[['train']]), stringsAsFactors=F);
          slot(nnDeconvResults, name)[[celltype]][["test"]] = data.frame(step = unlist(pyClass[[name]][[celltype]]$step[['test']]),
                                                                         loss = unlist(pyClass[[name]][[celltype]]$loss[['test']]),
                                                                         mse  = unlist(pyClass[[name]][[celltype]]$mse[['test']]),
                                                                         logp = unlist(pyClass[[name]][[celltype]]$log_probability[['test']]),
                                                                         kl   = unlist(pyClass[[name]][[celltype]]$kl_divergence[['test']]), stringsAsFactors=F);
        }
      }else if(name == 'deconv_samples_emb' & low.mem == F){
        slot(nnDeconvResults, name)[["prebatch"]] = pyClass[[name]][['prebatch']];
        slot(nnDeconvResults, name)[["postbatch"]] = pyClass[[name]][['postbatch']];
      }else if(name == 'deconv_samples_rec' & low.mem == F){
        slot(nnDeconvResults, name)[["prebatch"]] = pyClass[[name]][['prebatch']];
        slot(nnDeconvResults, name)[["postbatch"]] = pyClass[[name]][['postbatch']];
      }else if(name == "mixture_metrics"){
        slot(nnDeconvResults, name)[["train"]] = data.frame(loss = unlist(pyClass[[name]]$loss[['train']]),
                                                            mse  = unlist(pyClass[[name]]$mse[['train']]), stringsAsFactors=F);
        slot(nnDeconvResults, name)[["test"]]  = data.frame(loss = unlist(pyClass[[name]]$loss[['test']]),
                                                            mse  = unlist(pyClass[[name]]$mse[['test']]), stringsAsFactors=F);
      }else if(name %in% c("proportions", "weights")){
        res.list = list()
        if(length(pyClass[[name]]) > 0){
          for(step in names(pyClass[[name]])){
            proportion = pyClass[[name]][[step]]
            colnames(proportion) = pyClassBase[["celltypes"]]
            res.list[[step]] = proportion
          }
          res.list[["final"]] = proportion
          slot(nnDeconvResults, name) = res.list
        }
      }else if(name == "celltypes"){
        slot(nnDeconvResults, name) = as.character(pyClassBase[["celltypes"]])
      }else if(name %in% c("train_index", "test_index")){
        slot(nnDeconvResults, name) = pyClassBase[[name]]+1
      }
    }else{
      print(paste0("Not recorded: ", name))
    }
  }
  if(!is.null(pyClassBase[["flags"]])){slot(nnDeconvResults, "flags") = pyClassBase[["flags"]]}
  slot(nnDeconvResults, "genes") = as.character(gene.ids)
  slot(nnDeconvResults, "singleCellNames") = as.character(sc.ids)
  slot(nnDeconvResults, "mixtureCellNames") = as.character(mixture.ids)
  return(nnDeconvResults)
}

## ----setup---------------------------------------------------------------
library(SCORNET)

## ------------------------------------------------------------------------
surelda_run_phenorm <- with(
  simdata, sureLDA(X, ICD, NLP, HU, filter, nEmpty = 10))

## ---- message=FALSE------------------------------------------------------
surelda_scores_phenorm_aucs <- sapply(1:ncol(simdata$filter),function(k){
  pROC::auc(simdata$Y[,k],surelda_run_phenorm$scores[,k])
})

## ---- message=FALSE------------------------------------------------------
surelda_ensemble_phenorm_aucs <- sapply(1:ncol(simdata$filter),function(k){
  auc(simdata$Y[,k],surelda_run_phenorm$ensemble[,k])
})

## ------------------------------------------------------------------------
surelda_result_combined <- rbind(surelda_scores_phenorm_aucs,surelda_ensemble_phenorm_aucs)
rownames(surelda_result_combined) <- c('sureLDA Scores','sureLDA Probs')
print(surelda_result_combined)

## ------------------------------------------------------------------------
surelda_prediction <- with(simdata,
                           sureLDA(X, ICD, NLP, HU, filter, prior = surelda_run_phenorm$prior, nEmpty = 10,
                                   weight = surelda_run_phenorm$weight, phi = surelda_run_phenorm$phi))

## ---- message=FALSE------------------------------------------------------
surelda_scores_prediction_aucs <- sapply(1:ncol(simdata$filter),function(k){
  auc(simdata$Y[,k],surelda_prediction$scores[,k])
})

## ---- message=FALSE------------------------------------------------------
surelda_ensemble_prediction_aucs <- sapply(1:ncol(simdata$filter),function(k){
  auc(simdata$Y[,k],surelda_prediction$ensemble[,k])
})

## ------------------------------------------------------------------------
surelda_prediction_result_combined <- rbind(surelda_scores_prediction_aucs,surelda_ensemble_prediction_aucs)
rownames(surelda_prediction_result_combined) <- c('sureLDA Scores','sureLDA Probs')
print(surelda_prediction_result_combined)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  surelda_run_map <- with(simdata, sureLDA(X, ICD, NLP, HU, filter, prior = 'MAP'))
#  surelda_scores_map <- surelda_run_phenorm$scores
#  surelda_ensemble_map <- surelda_run_phenorm$ensemble

## ------------------------------------------------------------------------
proc.time()


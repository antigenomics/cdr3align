#!/usr/bin/env Rscript

library(pROC)
library(dplyr)
library(reshape2)
library(ggplot2)


compute_roc <- function(tcr_distances, random_tcr_distances, nlim, chain, ant, nbrdist_percentiles, datainfo){
  epitopes = as.character(unlist(unique(tcr_distances[ant])))
  
  names(tcr_distances)[names(tcr_distances)==ant] <- 'epitope.work'
  
  tcr_sums = tcr_distances %>% group_by(epitope.work) %>% summarise(tcr_summ = length(epitope.work)) %>% as.data.frame()
  total_sum = sum(tcr_sums$tcr_summ)
  
  
  for (antig in epitopes) {
    numobserv = as.integer(tcr_sums[tcr_sums['epitope.work']==antig][2])
    if (numobserv < nlim) next
    tcr_distances[antig] = as.numeric(tcr_distances['epitope.work'] == antig)
    random_tcr_distances[antig] = 0
    for (nbrdist in nbrdist_percentiles){
      roccolname <- paste(antig, chain, paste0('nbrdist', nbrdist), sep='_')
      rocwtdcolname <- paste(antig, chain, 'wtd', paste0('nbrdist', nbrdist), sep='_')
      
      rocdataframe <- data.frame(same_ant=tcr_distances[[antig]], dist=as.numeric(unlist(tcr_distances[[roccolname]])))
      
      if (length(epitopes) > 1) {
        inforoc <- roc(tcr_distances[[antig]], as.numeric(unlist(tcr_distances[[roccolname]])))
        infowtdroc <- roc(tcr_distances[[antig]], as.numeric(unlist(tcr_distances[[rocwtdcolname]])))
        
        #P1
        datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('others'), 
                                               auc=c(inforoc$auc[1]), npositives=c(numobserv), nnegatives=c(total_sum-numobserv)))#, roc=c(inforoc)))
        
        #P2
        
        datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('others_wtd'), 
                                               auc=c(infowtdroc$auc[1]), npositives=c(numobserv), nnegatives=c(total_sum-numobserv), roc=c(infowtdroc)))
      }
      
      #P3
      random_tcr = mutate(random_tcr_distances, !!quo_name(antig):=0)
      tcrdist = filter(tcr_distances, epitope.work == !!antig) %>% select(!!antig, !!roccolname)
      random_dist = rbind(tcrdist, select(random_tcr, !!roccolname, !!antig))
      randomroc <- roc(random_dist[[antig]], as.numeric(unlist(random_dist[[roccolname]])))
      
      datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('random'), 
                                             auc=c(randomroc$auc[1]), npositives=c(numobserv), nnegatives=c(nrow(random_tcr)), roc=c(randomroc)))
      
      #P4
      tcr_wtd_dist = filter(tcr_distances, epitope.work == !!antig) %>% select(!!antig, !!rocwtdcolname)
      random_wtd_dist = rbind(tcr_wtd_dist, select(random_tcr, !!rocwtdcolname, !!antig))
      randomwtdroc <- roc(random_wtd_dist[[antig]], as.numeric(unlist(random_wtd_dist[[rocwtdcolname]])))
      
      datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('random_wtd'), 
                                             auc=c(randomwtdroc$auc[1]), npositives=c(numobserv), nnegatives=c(nrow(random_tcr)), roc=c(randomwtdroc)))
      
    }
  }
}

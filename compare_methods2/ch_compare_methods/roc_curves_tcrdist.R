#!/usr/bin/env Rscript

library(pROC)
library(dplyr)
library(reshape2)
library(ggplot2)
library(optparse)
#../../imgt_work/dists/vdjdb_imgted_filtered_P.T_nbrdists.tsv
#../../workofothers/tcr-dist-master-upd/dists/human_pairseqs_v1_parsed_seqs_probs_mq20_clones_nbrdists.tsv

option_list = list(
  make_option(c("--tcr_distances"), type="character", default="../../imgt_work/dists/vdjdb_imgted_filtered_HomoSapiens_AB_nbrdists.tsv", 
              help="input distances path"),
  make_option(c("--random_tcr_distances"), type="character", 
              default='../../imgt_work/random_tcr_distances/vdjdb_imgted_filtered_HomoSapiens_AB_random_nbrdists.tsv', 
              help="input random distances path"),
  make_option(c("--organism"), type="character", default="HomoSapiens"),
  make_option(c("--nbrdist_percentiles"), type="character", default="5;10;25"),
  make_option(c("--epitope_col"), type="character", default="epitope.seq"),
  make_option(c("--output_dir"), type="character", default="../ch_compare_methods/compare/HomoSapiens_filtered_default"),
  make_option(c("--output_file"), type="character", default="filtered"),
  make_option(c("--script.dir"), type="character"),
  make_option(c("--chains"), type="character", default="AB"),
  make_option(c("--nlim"), type="integer", default=30)
)

parser = 
  parse_args(OptionParser(option_list=option_list))
nbrdist_percentiles = as.numeric(unlist(strsplit(parser$nbrdist_percentiles, ";")))

if (is.null(parser$script.dir)) {
  script.dir <- dirname(sys.frame(1)$ofile)
} else {
  script.dir <- parser$script.dir
}
script.dir
tcr_distances = read.table(file.path(script.dir, parser$tcr_distances), sep='\t', header=TRUE)
random_tcr_distances = read.table(file.path(script.dir, parser$random_tcr_distances), sep='\t', header=TRUE)

nlim = parser$nlim

chain = parser$chains

ant = parser$epitope_col

epitopes = as.character(unlist(unique(tcr_distances[ant])))

names(tcr_distances)[names(tcr_distances)==ant] <- 'epitope.work'

tcr_sums = tcr_distances %>% group_by(epitope.work) %>% summarise(tcr_summ = length(epitope.work)) %>% as.data.frame()
total_sum = sum(tcr_sums$tcr_summ)

datainfo = data.frame()

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
                                             auc=c(inforoc$auc[1]), npositives=c(numobserv), nnegatives=c(total_sum-numobserv)))
      png(filename=file.path(script.dir, parser$output_dir, paste0(roccolname, '.png')))
      plot(inforoc)
      text(0.5, 0, paste0('AUC = ', round(inforoc$auc[1], 3)))
      dev.off()
      
      #P2
      
      datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('others_wtd'), 
                                             auc=c(infowtdroc$auc[1]), npositives=c(numobserv), nnegatives=c(total_sum-numobserv)))
      png(filename=file.path(script.dir, parser$output_dir, paste0(rocwtdcolname, '.png')))
      plot(infowtdroc)
      text(0.5, 0, paste0('AUC = ', round(infowtdroc$auc[1], 3)))
      dev.off()
    }
    
    #P3
    random_tcr = mutate(random_tcr_distances, !!quo_name(antig):=0)
    tcrdist = filter(tcr_distances, epitope.work == !!antig) %>% select(!!antig, !!roccolname)
    random_dist = rbind(tcrdist, select(random_tcr, !!roccolname, !!antig))
    randomroc <- roc(random_dist[[antig]], as.numeric(unlist(random_dist[[roccolname]])))
    
    datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('random'), 
                                           auc=c(randomroc$auc[1]), npositives=c(numobserv), nnegatives=c(nrow(random_tcr))))
    png(filename=file.path(script.dir, parser$output_dir, paste0(roccolname, '_random.png')))
    plot(randomroc)
    text(0.5, 0, paste0('AUC = ', round(randomroc$auc[1], 3)))
    dev.off()
    
    #P4
    tcr_wtd_dist = filter(tcr_distances, epitope.work == !!antig) %>% select(!!antig, !!rocwtdcolname)
    random_wtd_dist = rbind(tcr_wtd_dist, select(random_tcr, !!rocwtdcolname, !!antig))
    randomwtdroc <- roc(random_wtd_dist[[antig]], as.numeric(unlist(random_wtd_dist[[rocwtdcolname]])))
    
    datainfo =  rbind(datainfo, data.frame(epitope=c(antig), chains=c(chain), nbrdist=c(nbrdist), roctype=c('random_wtd'), 
                                           auc=c(randomwtdroc$auc[1]), npositives=c(numobserv), nnegatives=c(nrow(random_tcr))))
    png(filename=file.path(script.dir, parser$output_dir, paste0(rocwtdcolname, '_random.png')))
    plot(randomwtdroc)
    text(0.5, 0, paste0('AUC = ', round(randomwtdroc$auc[1], 3)))
    dev.off()
  }
}


write.table(datainfo, file=file.path(script.dir, parser$output_dir, paste0(parser$output_file, '.datainfo.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
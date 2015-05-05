##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###############
##
## Demostration figures for variant calling algorithm validation review
##
## written by: Nate Olson
## NIST
## March 2015
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###############

## Figures for SNP method evaluation review
library(stringr)
library(reshape2)
library(ggplot2)
library(ROCR)
library(plyr)
library(dplyr)
library(tidyr)

## Load data
snp_calls <- read.csv("snp_demo_data.csv")

## predictions -----------------------------------------------------------------
predict.UG <- as.list(snp_calls[,seq(1,ncol(snp_calls)-1,2)])
predict.HC <- as.list(snp_calls[,seq(2,ncol(snp_calls)-1,2)])
labels <- matrix(rep(snp_calls$ref_call,ncol(snp_calls)-1),
                 nrow=nrow(snp_calls),ncol=(ncol(snp_calls)-1)/2,byrow=F)
pred.UG <- prediction(predict.UG, labels)
pred.HC <- prediction(predict.HC, labels)

## ROCR S4 to DF function ------------------------------------------------------
perfS4Todf <- function(perf){
  perf_df <- data.frame()
  for(i in 1:length(perf@x.values)){
    df <- data.frame(x.value= perf@x.values[[i]], y.value = perf@y.values[[i]]) 
    df$name <- i
    perf_df <- rbind(perf_df, df)
  }
  perf_df$x.value[is.infinite(perf_df$x.value)] <- 5000
  return(perf_df)
}

perf_merge <- function(ds1, ds2, name1, name2, metric, metric.2 = NULL){
  if(is.null(metric.2)){
    perf.1 <- performance(ds1, metric)
    perf.2 <- performance(ds2, metric) 
  }else{
    perf.1 <- performance(ds1, metric, metric.2)
    perf.2 <- performance(ds2, metric, metric.2)
  }
  #dataset 1
  perf.1.df <- perfS4Todf(perf.1)
  perf.1.df$class <- name1
  perf.1.df$metric <- metric
  
  #dataset 2
  perf.2.df <- perfS4Todf(perf.2)
  perf.2.df$class <- name2
  perf.2.df$metric <- metric
  
  return(rbind(perf.1.df,perf.2.df))  
}


## Calculating Metrics ---------------------------------------------------------
metrics <- as.list(c("tpr","fpr"))
single_metrics_df <- ldply(.data = metrics,
                           ds1 = pred.HC, ds2 = pred.UG,
                           name1 = "A", name2 = "B",
                           .fun = perf_merge)
# Haplotype caller - anonymized to A and Unified Genotyper anonymized to B
single_metrics_df$metric_name[single_metrics_df$metric == "tpr"] <- "True Positive Rate"
single_metrics_df$metric_name[single_metrics_df$metric == "fpr"] <- "False Positive Rate"


## static metrics --------------------------------------------------------------

static_metrics_df <- single_metrics_df[single_metrics_df$x.value == 3195.79,]
static_metrics_df$x.value <- ifelse(static_metrics_df$class == "A", 100000, 200000)
static_metrics_df$table_type <- "Static"
single_metrics_df$table_type <- "Dynamic"

single_metrics_df <- rbind(single_metrics_df, static_metrics_df)

## Fig 4 -----------------------------------------------------------------------


ggplot() + 
    geom_vline(mapping = aes(xintercept = 3195.79), 
               data = single_metrics_df[single_metrics_df$table_type == "Dynamic",],
               linetype = 2) + 
    geom_smooth(data = single_metrics_df[single_metrics_df$table_type == "Dynamic",],
                aes(x = x.value, y = y.value, color = class), family = "binomial") + 
    geom_boxplot(data = single_metrics_df[single_metrics_df$table_type == "Static",],
                 aes(x = x.value, y = y.value, color = class)) + 
    scale_x_continuous(breaks = c(0:5,100,200)*1000, 
                       labels = c(0,1000,2000,3000,4000,5000,
                                  "A","B")) +
    facet_grid(metric_name~table_type, scale = "free_x")  +
    labs(y =  "Metric Value", color = "Method") + 
    theme_bw() + theme(axis.title.x = element_blank(), legend.position=c(0.95,0.85))
ggsave("single_metric.png", height = 6, width = 8)




## Fig 5 -----------------------------------------------------------------------
twoMetricSummary <- ddply(static_metrics_df, .(class, metric), summarize, metric_mean = mean(y.value), 
                          lci = quantile(y.value, probs = c(0.025)),
                          uci = quantile(y.value, probs = c(0.975)))
metric_means <- dcast(twoMetricSummary, class~metric, value.var = "metric_mean")
metric_uci <- dcast(twoMetricSummary, class~metric, value.var = "uci")
metric_lci <- dcast(twoMetricSummary, class~metric, value.var = "lci")

met_summary <- join_all(list(metric_means, metric_uci, metric_lci), by = "class")
colnames(met_summary) <- c("Method", "fpr.mean","tpr.mean","fpr.uci", "tpr.uci","fpr.lci", "tpr.lci")
ggplot(met_summary) + geom_point(aes(x = fpr.mean, y = tpr.mean, color = Method, shape = Method), size = 4) +
  geom_errorbar(aes(x = fpr.mean, ymax = tpr.uci,ymin = tpr.lci, color = Method), width = 0) +
  geom_errorbarh(aes(x = fpr.mean, y = tpr.mean, xmax = fpr.uci,xmin = fpr.lci, color = Method), height = 0) +
  ylim(0,1)+xlim(0,1) +
  theme_bw() +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme(legend.position=c(0.9,0.2))
ggsave("scatter_static.png", height = 5.5, width = 4.5)



## Fig 6------------------------------------------------------------------------
### generating dataframe
dynamic_roc <- perf_merge(pred.HC, pred.UG,"A","B", "tpr","fpr")  %>% arrange(name, class, metric, x.value, y.value)
dynamic_roc$cl_name <- str_c(dynamic_roc$name, dynamic_roc$class, sep = "")

### 2-metric dynamic plot
ggplot(dynamic_roc) + 
    geom_path(aes(x = x.value,y = y.value, color = class, group= cl_name)) +
    labs(x= "False Positive Rate", y = "True Positive Rate", color = "Method") + 
#     facet_wrap(~class) + 
    theme_bw() + theme(legend.position = c(0.9,0.2))
ggsave("roc_dynamic.png", height = 6, width = 6)
## SNP review figures with uncertianty
library(stringr)
library(reshape2)
library(ggplot2)
library(ROCR)
library(plyr)

## Load data
snp_calls <- read.csv("snp_demo_data.csv")

## predictions ----------------------------------------------------------------------------------------------------
predict.UG <- as.list(snp_calls[,seq(1,ncol(snp_calls)-1,2)])
predict.HC <- as.list(snp_calls[,seq(2,ncol(snp_calls)-1,2)])
labels <- matrix(rep(snp_calls$ref_call,ncol(snp_calls)-1),nrow=nrow(snp_calls),ncol=(ncol(snp_calls)-1)/2,byrow=F)
pred.UG <- prediction(predict.UG, labels)
pred.HC <- prediction(predict.HC, labels)

## ROCR S4 to DF function ----------------------------------------------------------------------------------------
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


## Fig 2 ------------------------------------------------------------------------------------------------
metrics <- as.list(c("sens"))
single_metrics_df <- ldply(.data = metrics,ds1 = pred.HC, ds2 = pred.UG,name1 = "A", name2 = "B",.fun = perf_merge)
single_metrics_df$metric_name[single_metrics_df$metric == "sens"] <- "Sensitivity"


### Need to work out better method for finding on value for each dataset
static_metrics_df <- single_metrics_df[single_metrics_df$x.value == 3195.79,]
static_metrics_df$x.value <- ifelse(static_metrics_df$class == "A", 100000, 200000)
static_metrics_df$table_type <- "Static"
single_metrics_df$table_type <- "Dynamic"

single_metrics_df <- rbind(single_metrics_df, static_metrics_df)

ggplot() + 
  geom_vline(mapping = aes(xintercept = 3195.79), 
             data = single_metrics_df[single_metrics_df$table_type == "Dynamic",],
             linetype = 2) + 
  geom_smooth(data = single_metrics_df[single_metrics_df$table_type == "Dynamic",],
              aes(x = x.value, y = y.value, color = class), family = "binomial", show_guide = F) + 
  geom_boxplot(data = single_metrics_df[single_metrics_df$table_type == "Static",],
               aes(x = x.value, y = y.value, color = class), show_guide = F) + 
  scale_x_continuous(breaks = c(0:5,100,200)*1000, labels = c(0,20,40,60,80,100,"A","B")) +
  facet_wrap(~table_type, scale = "free_x")  +
  labs(y =  "Sensitivity", color = "Algorithm") + 
  theme_bw() + theme(axis.title.x = element_blank())
ggsave("~/Desktop/single_metric.png", height = 5, width = 8)
ggsave("~/Desktop/single_metric_small.png", height = 2.5, width = 4)

## Fig 3 ------------------------------------------------------------------------------------------------------------
# measure="tpr", x.measure="fpr"
dynamic_roc <- perf_merge(pred.HC, pred.UG,"HC","UG", "tpr","fpr")
ggplot(dynamic_roc) + 
  geom_smooth(aes(x = x.value,y = y.value, color = class), family = "binomial") +
  labs(x= "False Positive Rate", y = "True Positive Rate", color = "Algorithm") + 
  theme_bw()
ggsave("roc_metric.png", height = 6, width = 8)

## Need to work out code for static roc scatter plot
#two metric static scatter plot
# twoMetric <- single_metrics_df[single_metrics_df$x.value == 3195.79 & single_metrics_df$metric %in% c("Specificity", "Sensitivity"),]
# twoMetricSummary <- ddply(twoMetric, .(class, metric), summarize, metric_mean = mean(y.value), 
#                           lci = quantile(y.value, probs = c(0.025)),
#                           uci = quantile(y.value, probs = c(0.975)))
# metric_means <- dcast(twoMetricSummary, class~metric, value.var = "metric_mean")
# metric_uci <- dcast(twoMetricSummary, class~metric, value.var = "uci")
# metric_lci <- dcast(twoMetricSummary, class~metric, value.var = "lci")
# 
# met_summary <- join_all(list(metric_means, metric_uci, metric_lci), by = "class")
# colnames(met_summary) <- c("Algorithm", "sens.mean","spec.mean","sens.uci", "spec.uci","sens.lci", "spec.lci")
# ggplot(met_summary) + geom_point(aes(x = sens.mean, y = spec.mean, color = Algorithm, shape = Algorithm), size = 4) +
#   geom_errorbar(aes(x = sens.mean, ymax = spec.uci,ymin = spec.lci, color = Algorithm), width = 0) +
#   geom_errorbarh(aes(x = sens.mean, y = spec.mean, xmax = sens.uci,xmin = sens.lci, color = Algorithm), height = 0) +
#   ylim(0,1)+xlim(0,1) +
#   theme_bw() +
#   labs(x = "Sensitivity", y = "Specificity") +
#   theme(legend.position = "bottom", legend.position = "horizontal")
# ggsave("scatter_static.png", height = 5.5, width = 4.5)

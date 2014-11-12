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
labels <- matrix(rep(snp_calls$ref_call,ncol(snp_calls)-1),
                 nrow=nrow(snp_calls),ncol=(ncol(snp_calls)-1)/2,byrow=F)
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
metrics <- as.list(c("tpr","fpr"))
single_metrics_df <- ldply(.data = metrics,ds1 = pred.HC, ds2 = pred.UG,name1 = "Haplotype Caller", name2 = "Unified Genotyper",.fun = perf_merge)
single_metrics_df$metric_name[single_metrics_df$metric == "tpr"] <- "True Positive Rate"
single_metrics_df$metric_name[single_metrics_df$metric == "fpr"] <- "False Positive Rate"


# static
#static_metrics_df <- single_metrics_df[single_metrics_df$x.value == 3195.79,]
single_metrics_df$diff <- abs(single_metrics_df$x.value - 3195.79)
static_filter <- ddply(single_metrics_df,.(name,class), summarize, diff = min(diff))
static_filter <- join(static_filter, static_metrics_df, type = "left")

static_filter$x.value <- ifelse(static_filter$class == "Haplotype Caller", 100000, 200000)
static_filter$table_type <- "Static"

single_metrics_df$table_type <- "Dynamic"

single_metrics_df <- rbind(single_metrics_df, static_filter)


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
                                "Haplotype\nCaller","Unified\nGenotyper")) +
  facet_grid(metric_name~table_type, scale = "free_x")  +
  labs(y =  "Metric Value", color = "Algorithm") + 
  theme_bw() + theme(axis.title.x = element_blank())
ggsave("single_metric.png", height = 6, width = 8)



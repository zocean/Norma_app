#!~/R/bin/Rscript

library(ggplot2)
library(scales)
library(reshape2)

pdf_file = paste("figures", "cumulative_housekeeping.pdf", sep='/')

# load the annotation result as a data.frame
table <- read.table(file="results/20kb_anno.txt", header = TRUE, sep='\t', comment.char = "")

# check data
#head(table)
#summary(table)

# select rows with no NA values in samples
filter_row <- complete.cases(table[, c("HCT116", "HFF", "K562")])
filter_table <- table[filter_row, ]
# check how many rows are kept
nrow(filter_table)
nrow(table)

# convert the TSA score into _ercentile
# pay attention to ties, by default, rank function uses 'average', which means that return the average ranking of ties
filter_table$HCT116_per <- rank(filter_table$HCT116)/nrow(filter_table)
filter_table$HFF_per <- rank(filter_table$HFF)/nrow(filter_table)
filter_table$K562_per <- rank(filter_table$K562)/nrow(filter_table)
# check percentile
#summary(filter_table)

# calculate the average of the percentile
filter_table$metric_mean_per <- apply(filter_table[, c("HCT116_per", "HFF_per", "K562_per")], 1, mean)
# calculate the standard divation
filter_table$metric_sd_per <- apply(filter_table[, c("HCT116_per", "HFF_per", "K562_per")], 1, sd)
# calculate the coefficient of variation (CV)
filter_table$metric_cv_per <- filter_table$metric_mean_per/filter_table$metric_sd_per

# function
cal_cumu_count <- function(table, metric_col, hk_col, non_hk_col) {
    "
    idea is to sort metric_col from small to large
    then count how many rows in the gene_col is not NA
    for example,
    metric_col  gene_col    cumulative_count
    0           NA          0
    0.1         gene1       1
    0.2         NA          0
    0.3         gene2       2
    - table
    - metric_col
    - hk_col
    - non_hk_col
    "
    x <- table[order(table[, metric_col]), c(metric_col, hk_col, non_hk_col)]
    x$cum_hk <- cumsum(!is.na(x[, hk_col]))
    x$cum_non_hk <- cumsum(!is.na(x[, non_hk_col]))
    return(x)
}

# calculate the cumulative number of housekeeping genes relative to the mean percentage
result_table <- cal_cumu_count(filter_table, 'metric_mean_per', 'hk', 'non_hk')
result_table$hk_per <- result_table$cum_hk/max(result_table$cum_hk)
result_table$non_hk_per <- result_table$cum_non_hk/max(result_table$cum_non_hk)
# melt the data.frame
per_table_mean <- melt(result_table, id.vars = "metric_mean_per", measure.vars = c("hk_per", "non_hk_per"), variable.name = "gene", value.name = "per")

# calculate the cumulative number of housekeeping genes relative to the cv percentage
result_table <- cal_cumu_count(filter_table, 'metric_cv_per', 'hk', 'non_hk')
result_table$hk_per <- result_table$cum_hk/max(result_table$cum_hk)
result_table$non_hk_per <- result_table$cum_non_hk/max(result_table$cum_non_hk)
# melt the data.frame
per_table_cv <- melt(result_table, id.vars = "metric_cv_per", measure.vars = c("hk_per", "non_hk_per"), variable.name = "gene", value.name = "per")

# plot
plot_per_mean <- ggplot(per_table_mean, aes(x = metric_mean_per, y = per, color = gene)) +
  geom_line(size = 1.5) + 
  xlab("TSA scores percentile (mean)") +
  ylab("Cumulative percentage") +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(name="Gene",values=c("hk_per"="#e41a1c", "non_hk_per"="#377eb8")) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(2))) +
  theme(panel.border = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size=rel(1.5),margin=margin(0,10,0,0))) +
  theme(axis.title.x = element_text(size=rel(1.5),margin=margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=rel(1.5), color="black")) +
  theme(axis.line = element_line(color="black"))

plot_per_cv <- ggplot(per_table_cv, aes(x = metric_cv_per, y = per, color = gene)) +
  geom_line(size = 1.5) + 
  xlab("TSA scores percentile (CV)") +
  ylab("Cumulative percentage") +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(name="Gene",values=c("hk_per"="#e41a1c", "non_hk_per"="#377eb8")) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(2))) +
  theme(panel.border = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size=rel(1.5),margin=margin(0,10,0,0))) +
  theme(axis.title.x = element_text(size=rel(1.5),margin=margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=rel(1.5), color="black")) +
  theme(axis.line = element_line(color="black"))

pdf(file=pdf_file, width=5, height=5)
print(plot_per_mean)
dev.off()

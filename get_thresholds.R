
dis <- read.table("ref_msh_dis_clu.tsv", sep='\t', header=TRUE, comment.char='@')

between_dis <- dis[dis$ref_cluster != dis$met_cluster, ]
within_dis <- dis[dis$ref_cluster == dis$met_cluster, ]

within_median <- unlist(by(within_dis$distance, within_dis$ref_cluster, median))
within_min <- unlist(by(within_dis$distance, within_dis$ref_cluster, min))
within_max <- unlist(by(within_dis$distance, within_dis$ref_cluster, max))

between_median <- unlist(by(between_dis$distance, between_dis$ref_cluster, median))
between_min <- unlist(by(between_dis$distance, between_dis$ref_cluster, min))
between_max <- unlist(by(between_dis$distance, between_dis$ref_cluster, max))

clu <- read.table("ref_clu_comp.tsv", sep='\t', header=TRUE, comment.char='@')

thr_prop_exp <- 0.2
thr_prop_min <- 0.2

threshold <- rep(0, nrow(clu))
names(threshold) <- clu[, 1]

within.max.fake <- rep(0, nrow(clu))
names(within.max.fake) <- clu[, 1]
within.max.fake[names(within_max)] <- within_max

threshold <- within.max.fake*(1 + thr_prop_exp)
threshold[threshold > between_min] <- within.max.fake[threshold > between_min]

t_min <- median(between_median*thr_prop_min)
t_med <- median(threshold[threshold > 0])
t_comp <- ifelse(t_med < t_min, t_min, t_med)

t_summary <- cbind(clu[, 1:2], threshold, within.max.fake, rep(median(within_median), nrow(clu)), rep(median(between_median), nrow(clu)))
colnames(t_summary) <- c("cluster", "n", "threshold", "dis_same_max", "dis_same_med_all", "dis_diff_med_all")

t_summary$threshold <- ifelse(t_summary$threshold < t_comp, t_comp, t_summary$threshold)

write.table(t_summary, file = "ref_clu_thr.tsv", sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

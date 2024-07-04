#' @import pbapply
#' @import kohonen
#' @import reshape2
#' @import circlize
#' @import ggplot2
#' @import patchwork
#' @import ggsci
#' @import vegan
#' @import ggpubr
#' @import permute
#' @import lattice
#' @import dplyr
#' @import rstatix
#' @import forcats
#' @importFrom cowplot plot_grid
#' @importFrom caret rfe rfeControl rfFuncs trainControl twoClassSummary train
#' @importFrom ggh4x force_panelsizes
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom linkET mantel_test qcorrplot geom_square geom_mark geom_couple correlate nice_curvature color_pal
#' @importFrom pROC roc ci.auc ggroc ci.se
#' @import tibble
#' @import viridis
#' @import scales

#' @export
X_som <- function(gated_fcs = NULL, fcs_data = NULL, m = 1e3, n_hclust = 300,
                  n_cells_sub = 3e5, cluster_to_show = NULL, out_path = './') {

  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }

  cat("Drawing random sample\n")

  sub_sample <- downsample(
    gated_fcs, # list of matrices from fcs, gated
    n = n_cells_sub, # total no of cells to sub-sample from all samples; CHANGE ME!
    samples = names(gated_fcs) # run all samples; CHANGE ME if needed
  )

  # as list (match with format of the list with sample data)
  sub_sample <- list(name = "random sample",
                     data = sub_sample)

  # compute SOM on combined data subset
  cat("\tDONE.\nTraining SOM\n")
  sub_sample<- compute_som(
    sub_sample, n_cells = nrow(sub_sample$data) / m
  ) # train SOM, average c. 0.05% of cells per node (i.e. 1/2000)
  cat("\tDONE.\nSaving\n")

  ## extracting the SOM cluster information
  cohonen_information<-as.data.frame(sub_sample[["som"]][["codes"]][[1]])
  assign("cohonen_information", cohonen_information, envir = .GlobalEnv)
  write.csv(cohonen_information, paste0(out_path, "/cluster_information.csv"))

  ########
  min_cells <- 2e5
  if (any(!sapply(gated_fcs, function(x) nrow(x$data)) >= min_cells)) {
    warning(
      paste0("\nDropping samples:",
             paste(names(gated_fcs)[!sapply(gated_fcs,
                                            function(x) nrow(x$data)) >= min_cells],
                   collapse = ", "),
             "low cell count!\n"))
  }
  gated_fcs <- gated_fcs[sapply(gated_fcs, function(x) nrow(x$data)) >= min_cells]
  n_subset_large <- 1e10

  # map to trained SOM
  SOM_fcs <- pbapply::pblapply(gated_fcs, map_som,
                      trained = sub_sample$som,
                      n_subset = n_subset_large)

  cat("Annotating Clusters\n")
  cluster_number<-as.matrix(as.list(1:1024))
  colnames(cluster_number)<-as.character(1024)
  SOM_fcs<-
    pbapply::pblapply(SOM_fcs, assign_clusters, clusters = cluster_number)
  cat("\tDONE.\n")

  cat("Counting Observations\n")
  SOM_fcs <- pbapply::pblapply(SOM_fcs, count_observations,
                      clusters = colnames(cluster_number))

  count_tables <- get_counts(SOM_fcs)
  tmp_path <- file.path(out_path, "/count_tables/")
  if (!dir.exists(tmp_path)) dir.create(tmp_path, recursive = TRUE)
  save(count_tables,
       file = file.path(tmp_path, "count_table_SOM.save"),
       compress = TRUE)

  original <- as.data.frame(t(count_tables[[1]]))
  count_table_2025 <- original / rowSums(original)


  if (sum(sapply(count_table_2025, min) < 1e-4) > 0) {
    warning(paste0(sum(sapply(count_table_2025, min)  < 1e-4),
                   ' clusters contains cells less than 0.01%!'))
  }

  # Hcluster for cohonen
  hc <- hclust(dist(cohonen_information), method = "ward.D2")
  # cut the hcluster by n
  hc_label <- cutree(hc, n_hclust)

  selected_rows <- c('V20','V102')

  new_colnames <- sapply(colnames(cohonen_information), function(x) {
    if (x == "FSC PAR") {
      return("FSC.PAR")
    } else {
      return(sub("\\.[^.]+$", "", x) %>% sub(' ', '\\.',.) %>% sub('-', '\\.', .))
    }
  })

  bins <- cohonen_information %>%
    mutate_all(~ 10^((4 * .x) / 65000)) %>%
    mutate(hclust = hc_label) %>%
    `colnames<-`(sub(' ', '\\.', sub("\\.[^.]+$", "", colnames(.)))) %>%
    .[row.names(cohonen_information) %in% selected_rows, ]



  plot_data <- exprs(fcs_data) %>%
    as.data.frame() %>%
    mutate_at(vars(-'classes'), ~ 10^((4 * .x) / 65000)) %>%
    `colnames<-`(c(new_colnames, "classes"))

  gplots <- function(dat, x, y, bins, selected_rows_plot){
    dat %>%
      round(2) %>%
      ggplot(aes_string(x = x, y = y)) +
      geom_hex(bins = 100) +
      geom_point(data = dat[dat$classes %in% cluster_to_show, ],
                 aes_string(x=x, y=y),
                 color = "grey20", alpha = 0.3, size =0.5) +
      scale_fill_viridis(discrete = F, trans = 'log') +
      scale_x_log10(breaks = c(0, 10, 100, 1000, 10000),
                    labels = trans_format("log10", scales::math_format(10^.x)),
                    limits = c(0.9, 11000)) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    labels = trans_format("log10", scales::math_format(10^.x)),
                    limits = c(0.9, 11000)) +
      annotate('text', x = bins[, x], y = bins[, y], label = cluster_to_show, color = 'white',size = 3) +
      theme_bw() +
      theme(text = element_text(family = 'Times'))
  }

  g1 <- gplots(plot_data, "FSC.PAR", 'SSC', bins, selected_rows_plot)
  g2 <- gplots(plot_data, 'FSC.PAR', 'Hoechst.Red', bins, selected_rows_plot)
  g3 <- gplots(plot_data, 'FITC', 'APC', bins, selected_rows_plot)
  g4 <- gplots(plot_data, 'Pe.TR', 'BV650', bins, selected_rows_plot)

  ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2, common.legend = T)

  write.csv(count_table_2025, paste0(out_path, "/Ig_SOM_1024.csv"))
  assign("original", original, envir = .GlobalEnv)
  assign("count_table", count_table_2025, envir = .GlobalEnv)
}

#' @export
X_stat <- function(data = NULL, meta_data = NULL, group_col = NULL,
                   test_type = 'wilcox', cutoff = 0.05,
                   correction = 'none', out_path = './') {

  cat(paste0(test_type, " test\n"))
  if (length(unlist(meta_data[group_col])) != nrow(data)) {
    stop("Error: The lengths of group_labels and number of samples are not equal.")
  }

  if (test_type == 'wilcoxon') {
    pvalue <- sapply(data, function(t) {
      tryCatch({wilcox.test(t ~ as.vector(unlist(meta_data[group_col])))$p.value}, warning = function(e){1})})
  } else if (test_type == 'kruskal') {
    pvalue <- sapply(data, function(t) {
      tryCatch({kruskal.test(t ~ as.vector(unlist(meta_data[group_col])))$p.value}, warning = function(e){1})})
  } else if (test_type == 'anova') {
    pvalue <- sapply(data, function(t) {
      tryCatch({anova(aov(t ~ as.vector(unlist(meta_data[group_col]))))$`Pr(>F)`[1]}, warning = function(e){1})})
  } else {
    stop("Error: test_type should be one of 'wilcoxon', 'kruskal', 'anova'.")
  }

  if (correction == 'fdr') {
    pvalue <- p.adjust(pvalue, method = 'fdr')
  } else if (correction == 'bonferroni') {
    pvalue <- p.adjust(pvalue, method = 'bonferroni')
  } else if (correction == 'BH') {
    pvalue <- p.adjust(pvalue, method = 'BH')
  } else if (correction != 'none') {
    stop('Error: correction should be one of "none", "fdr", "bonferroni", "BH".')
  }
  pvalue_data = as.data.frame(pvalue) %>% replace(is.na(.), 1)
  significant_data = data[, pvalue < cutoff]

  assign("pvalue_data", pvalue_data, envir = .GlobalEnv)
  assign("significant_data", significant_data, envir = .GlobalEnv)

  write.csv(pvalue_data, paste0(out_path, "/Ig_SOM_1024_pvalue.csv"))
  write.csv(significant_data, paste0(out_path, "/Ig_SOM_1024_significant.csv"))
  cat("\tDONE.\n")
}

#' @export
X_circle <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     out_path = './', width = 8, height = 8,
                     point_colors = c("red", "blue"),
                     cell_colors = c("blue", "white", "red")) {
  cat("Circled heatmap\n")
  group_labels <- as.vector(unlist(meta_data[group_col]))
  mean_table <- t(sapply(data, function(t) tapply(t, group_labels, mean)))
  rownames(mean_table) <- gsub('V', 'Bin', rownames(mean_table))
  min_val <- min(mean_table)
  max_val <- max(mean_table)
  mean_val <- (min_val + max_val) / 2
  col_fun = colorRamp2(c(min_val, mean_val, max_val), cell_colors)

  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })

  pdf(paste0(out_path, '/circle.pdf'), height = height, width = width)
  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })
  invisible(dev.off())
  cat("\tDONE.\n")
}

#' @export
X_violin <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     pvalue_data = NULL, cluster = 1, out_path = './',
                     colors = c('#E41A1C', '#377EB8'),
                     width = 4, height = 4) {
  cat("Violin plot\n")
  group_labels <- as.vector(unlist(meta_data[group_col]))
  clusters <- as.numeric(gsub('V', '', colnames(data)))
  cluster_data <- data[, clusters == cluster]
  p_value <- round(pvalue_data[clusters == cluster, ], 3)
  p <- ggplot(mapping = aes(x = group_labels, y = cluster_data)) +
    geom_violin(aes(fill = group_labels), alpha = 0.5) +
    geom_boxplot(aes(fill = group_labels), width = 0.1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text(aes(label = ifelse(p_value < 0.05, paste0('P = ', p_value, '*'), paste0('P = ', p_value))),
              x = 1.5, y = max(cluster_data) * 0.95, size = 4) +
    scale_fill_manual(values = colors) +
    labs(title = paste0('Cluster ', cluster),
         x = 'Group',
         y = 'Relative Anundance',
         fill = 'Group')
  print(p)
  pdf(paste0(out_path, '/violin.pdf'), height = height, width = width)
  print(p)
  invisible(dev.off())
  cat("\tDONE.\n")
}

#' @export
X_beta <- function(data, out_path = './', test = 'wilcox.test', meta_data = NULL,
                   group_name = NULL, colors = c('#E41A1C', '#377EB8'),
                   width = 5, height = 5) {

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL
  groups<-as.data.frame(cbind(sample=paste('sample',rownames(df),sep = ''),group=df$Health.State))
  write.table(groups, paste0(out_path, '/.group.txt'), row.names = F,sep = '\t',quote = F)
  df$Health.State<-groups$sample
  rownames(df)<-df$Health.State
  dataT<-df[,-1]

  dist <- vegdist(dataT, method="bray")
  dist <- as.matrix(dist)
  adist<- as.dist(dist)

  options(stringsAsFactors=F)

  sd <- groups
  rownames(sd) <- as.character(sd[,1])
  sd[,2] <- as.character(sd[,2])

  dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]

  pc_num <-c(1,2)
  pc_x <- pc_num[1]
  pc_y <- pc_num[2]

  pcoa <- cmdscale(dist, k=3, eig=TRUE)
  pc12 <- pcoa$points[,pc_num]
  pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits = 2)
  pc12 <- as.data.frame(pc12)
  colnames(pc12) <- c("pc_x","pc_y")
  pc12['sample'] <- rownames(pc12)
  colnames(sd)[1:2] <- c("sample","group")
  sd$group<-factor(sd$group,levels=sd$group[!duplicated(sd$group)])
  pc12 <- merge(pc12,sd,by="sample")
  pc12$group<-factor(pc12$group,levels=levels(sd$group))

  p<- ggscatter(pc12, x = "pc_x", y = "pc_y",
               color = "group", shape = "group", linewidth=3,
               ellipse = FALSE, conf.int.level = 0.95,
               #palette = colors,
               mean.point = TRUE,
               star.plot = TRUE,star.plot.lty = 3) +
    #ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    #xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "solid")+
    geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "solid")+
    scale_color_manual(values = colors) +

    theme(axis.title = element_blank(),
          legend.position = "top",legend.title = element_blank(),
          panel.border = element_rect(color = "black",linewidth  = 1.0, fill = NA),
          text = element_text(size=12)) +theme(plot.margin = unit(c(0,0,0,0),'cm'))

  mycols=pal_npg("nrc")(10)[1:length(levels(sd$group))]

  ADONIS<-suppressMessages(adonis2(dist~sd$group))

  TEST<-ADONIS$`Pr(>F)`[1]
  R2adonis<-round(ADONIS$R2[1],digits = 3)
  sink(paste0(out_path, '/adonis.txt'))
  print(ADONIS)
  sink()

  p<-p+ggtitle(label =expr(paste(bold(Bray)," ",bold(Curtis),bold(","),
                                 bold(Adnois:R^2),bold('='),bold(!!R2adonis),
                                 bold(","),bold(P),bold('='),!!TEST))) +
    theme(title = element_text(size = 10))

  cp <- combn(levels(pc12$group),2)
  comp <- list()
  for(i in 1:ncol(cp)){
    comp[[i]] <- cp[,i]
  }

  pl<-ggboxplot(pc12, x="group", y="pc_y", fill = "group", palette = colors) +
    stat_compare_means(comparisons = comp, label = "p.signif",method="wilcox.test")+
    ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    theme(panel.border = element_rect(color = "black",linewidth = 1.0,fill = NA),
          #axis.text.y = element_blank(),
          #axis.ticks.y= element_blank(),
          #axis.title.y = element_blank(),
          legend.position = "none",
          axis.text = element_text(face='bold'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12,angle = 60,hjust = 1,face='bold'),
          axis.title.y = element_text(size = 15,face='bold'))+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'))

  pt<-ggboxplot(pc12, x="group", y="pc_x", fill = "group", palette = colors) + coord_flip() +
    stat_compare_means(comparisons = comp, label = "p.signif",method = test) +
    scale_x_discrete(limits = rev(levels(pc12$group)))+
    ylab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          #axis.title.x = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 12, angle = 0,face='bold'),
          axis.title.y=element_blank(),
          axis.title.x = element_text(size = 15,face='bold'))+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'))

  p0 <- ggplot() + theme(panel.background = element_blank(),
                         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"))

  p_final<- pl+p+p0+pt+plot_layout(ncol = 2,nrow = 2,heights = c(4,1),widths = c(1,4))+
    plot_annotation(theme = theme(plot.margin = margin()))

  ggsave(paste0(out_path, '/pcoa.pdf') ,width = width, height = height)
  return(p_final)
  cat("\tDONE.\n")
}

#' @export
X_fs <- function(data = NULL, out_path = './',
                 test = 'wilcox.test',
                 meta_data = NULL,
                 group_name = NULL,
                 nfolds_cv = 5,
                 top_n_features = 20,
                 rfe_size = 10,
                 ref_group = NULL,
                 colors = c('#E41A1C', '#377EB8'),
                 seed = 123,
                 width = 5, height = 5) {

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL
  groups<-as.data.frame(cbind(sample=paste('sample',rownames(df),sep = ''),group=df$Health.State))
  write.table(groups, paste0(out_path, '/.group.txt'), row.names = F,sep = '\t',quote = F)
  df$Health.State<-groups$sample
  rownames(df)<-df$Health.State
  df<-df[,-1]
  df<-as.data.frame(t(df))
  df2<-as.data.frame(t(df))
  df2$group<-factor(groups$group)

  set.seed(seed)
  control <- rfeControl(functions=rfFuncs, method="cv", number = nfolds_cv)
  # run the RFE algorithm
  results <- rfe(df2[,1:(ncol(df2)-1)], df2[,ncol(df2)], sizes=c(1:rfe_size), rfeControl=control)
  # summarize the results
  print(results)
  # list the chosen features
  predictors(results)
  p1_out <- plot(results, type=c("g", "o"), main = 'Feature Selection',col = '#377EB8', lwd = 2)

  best_selection<-results$optVariables
  print(paste0('the number of feature is ',length(best_selection)))

  df2<-df2[,c(best_selection,'group')]
  write.table(df2, paste0(out_path, '/figure1.txt'),row.names = F,sep = '\t',quote = F)

  varimp_data <- data.frame(feature = row.names(varImp(results))[1:top_n_features],
                            importance = varImp(results)[1:top_n_features, 1])

  p2_out <- ggplot(data = drop_na(varimp_data),
         aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") +
    geom_text(aes(label = round(importance, 2)), vjust=1.6, color="white", size=4) +
    theme_bw() + theme(legend.position = "none") +
    labs(title = "Variable Importance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  # effect size
  esize <- df2 %>%
    #select(varimp_data$feature, group) %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    group_by(variable) %>%
    suppressMessages() %>%
    cohens_d(value ~ group, conf.level = 0.95, ci = T, ref.group = ref_group) %>%
    arrange(desc(effsize))

  # figure 2
  fig2 <- esize %>%
    ggplot(aes(y = reorder(variable, effsize), x = effsize, fill = variable)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = .2) +
    geom_point(aes(x = effsize), size = 4, color = 'black') +
    geom_point(aes(x = effsize, color = ifelse(effsize > 0, 'positive', 'negative')), size = 3) +
    theme_minimal() + theme(legend.position = "none") +
    scale_color_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    labs(title = "Effect Size", x = "Cohen's d", y = '') +
    #remove y text
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    force_panelsizes(rows = 0.5, cols = 0.5)

  # figure1
  fig1 <- df2 %>%
    #elect(varimp_data$feature, group) %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    ggplot(aes(x = value, y = fct_rev(factor(variable, esize$variable)), fill = group)) +
    geom_boxplot(position = 'dodge') +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "log10(abundance)", y = "Variable") +
    scale_fill_manual(values = colors) +
    force_panelsizes(rows = 0.5, cols = 0.5)

  # figure3 barplor of varimp
  fig3 <- merge(varimp_data, esize, by.x = 'feature', by.y = 'variable') %>%
    ggplot(aes(x = importance, y = fct_rev(factor(feature, esize$variable)))) +
    geom_bar(stat="identity", color = 'black', aes(fill = ifelse(effsize > 0, 'positive', 'negative'))) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = "Variable Importance", y = '') +
    scale_fill_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    #remove y text
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    force_panelsizes(rows = 0.5, cols = 0.5)

  p3_out <- plot_grid(fig1, fig2, fig3, ncol = 3, align = 'h', axis = 't')

  pdf(paste0(out_path, '/feature_exploration.pdf'), width = width, height = height)
  p1_out
  p2_out
  p3_out
  dev.off()

  show_plots <- function(plots) {
    for (plot in plots) {
      print(plot)
      readline(prompt="Press [enter] to see the next plot")
    }
  }

  show_plots(list(p1_out, p2_out, p3_out))
  cat("\tDONE.\n")
}

#' @export
X_heatmap <- function(data = NULL, cohonen_information = NULL,
                      out_path = './',
                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                      width = 10, height = 10,
                      scale = 'row', cluster_rows = F, cluster_cols = F,
                      display_numbers = T, boarder_color = 'grey60',
                      legend = T) {

  df <- cohonen_information[rownames(cohonen_information) %in% names(data), ]
  p <- pheatmap(df, scale = scale, cluster_rows = cluster_rows,
           cluster_cols = cluster_cols, display_numbers = display_numbers,
           border_color = boarder_color, color = color, legend = legend)
  pdf(paste0(out_path, '/heatmap.pdf'), width = width, height = height)
  print(p)
  dev.off()
  print(p)
  cat("\tDONE.\n")
}

#' @export
X_mantel <- function(data = NULL, meta_data = NULL, clinical_cols = NULL,
                     demographic_cols = NULL,
                     colors=RColorBrewer::brewer.pal(11, "RdBu"),
                     out_path = './', width = 8, height = 8) {
  mantel <- suppressMessages(mantel_test(meta_data, data,
                        spec_select = list(cli = clinical_cols,
                                           dem = demographic_cols)) %>%
    mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))))
  #> `mantel_test()` using 'bray' dist method for 'spec'.
  #> `mantel_test()` using 'euclidean' dist method for 'env'.

  p <- qcorrplot(correlate(data), type = "lower", diag = FALSE) +
    geom_square() +
    geom_mark(sep='\n', size = 1.8, sig_level = c(0.05, 0.01, 0.001),
              sig_thres = 0.05, color="white")+
    geom_couple(aes(colour = pd, size = rd),
                data = mantel,
                curvature = nice_curvature()) +
    scale_fill_gradientn(colours = colors) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))

  #pdf
  pdf(paste0(out_path, '/mantel.pdf'), width = width, height = height)
  print(p)
  dev.off()

  print(p)
  cat("\tDONE.\n")
}

#' @export
X_ml <- function(data = NULL, meta_data = NULL,
                 group_name = 'Group', out_path = './',
                 reference_level = 'A',
                 width = 5, height = 5) {
  rownames(data) <- NULL
  groups <- as.data.frame(cbind(sample=paste('sample',rownames(data),sep = ''), group=meta_data[group_name]))
  write.table(groups,paste0(out_path, 'ml_group.txt'),row.names = F,sep = '\t',quote = F)
  data <- cbind(data, Health.State = groups$sample)
  rownames(data)<-data$Health.State
  data <-data[,-ncol(data)]
  df <- as.data.frame(t(data))
  train <- as.data.frame(t(df))
  train$group <- factor(pull(groups[group_name]))


  fitControl <-
    trainControl(
      method = "repeatedcv",
      number = 5, repeats = 5,
      returnResamp="final",
      classProbs = T,
      savePredictions = T,
      indexFinal=NULL,
      summaryFunction = twoClassSummary
    )

  rf <- train(
    group ~ .,
    data = train,
    method = "rf",
    trControl = fitControl,
    metric = "ROC", verbose = FALSE
  )

  # train
  rocs_train <- roc(response = ifelse(train$group == reference_level, 0, 1), predictor = rf$finalModel$votes[,2])
  ci_auc_train <- suppressWarnings(round(as.numeric(ci.auc(rocs_train)), 3))
  ci_tb_train <- suppressWarnings(as.data.frame(ci.se(rocs_train)))
  ci_tb_train <- suppressWarnings(rownames_to_column(ci_tb_train, var = 'x'))
  ci_tb_train <- as.data.frame(sapply(ci_tb_train, as.numeric))
  names(ci_tb_train) <- c('x', 'low', 'mid', 'high')

  metrics_train <- confusionMatrix(rf$finalModel$predicted, train$group, positive = reference_level, mode = 'everything')

  options(warn = -1)
  g1_out <- ggroc(rocs_train, legacy.axes = TRUE) +
    coord_equal() +
    geom_ribbon(aes(x = 1-x, ymin = low, ymax = high), data = ci_tb_train, alpha = 0.5,
                fill = 'lightblue') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', alpha = 0.7) +
    geom_text(aes(0.5, 0.25, hjust = 0,
                  label = paste0('AUC: ', round(rocs_train$auc, 3), ' 95%CI: ', ci_auc_train[1], ' ~ ', ci_auc_train[3])))  +
    geom_text(aes(0.5, 0.2, hjust = 0,
                  label = paste0('Sensitivity: ', round(as.numeric(metrics_train$byClass[1]), 3)))) +
    geom_text(aes(0.5, 0.15, hjust = 0,
                  label = paste0('Specificity: ', round(as.numeric(metrics_train$byClass[2]), 3)))) +
    geom_text(aes(0.5, 0.1, hjust = 0,
                  label = paste0('F1: ', round(as.numeric(metrics_train$byClass[7]), 3)))) +
    theme_classic() +
    labs(x = '1 - Specificity (% false postivie)',
         y = 'Sensitivity (% true positive')

  pdf(paste0(out_path, 'roc_confusion.pdf'), width = width, height = height)
  print(g1_out)
  dev.off()

  print(g1_out)
  cat("\tDONE.\n")
}

#' @export
X_conf <- function(data = NULL, meta_data = NULL,
                   group_name = 'Group', out_path = './',
                   reference_level = 'A',
                   colors = c('#00BFC4', '#F8766D'),
                   width = 5, height = 5) {

  rownames(data) <- NULL
  groups <- as.data.frame(cbind(sample=paste('sample',rownames(data),sep = ''), group=meta_data[group_name]))
  write.table(groups,paste0(out_path, 'ml_group.txt'),row.names = F,sep = '\t',quote = F)
  data <- cbind(data, Health.State = groups$sample)
  rownames(data)<-data$Health.State
  data <-data[,-ncol(data)]
  df <- as.data.frame(t(data))
  train <- as.data.frame(t(df))
  train$group <- factor(pull(groups[group_name]))


  fitControl <-
    trainControl(
      method = "repeatedcv",
      number = 5, repeats = 5,
      returnResamp="final",
      classProbs = T,
      savePredictions = T,
      indexFinal=NULL,
      summaryFunction = twoClassSummary
    )

  rf <- train(
    group ~ .,
    data = train,
    method = "rf",
    trControl = fitControl,
    metric = "ROC", verbose = FALSE
  )

    ## Confusion matrix plot
  ### train
  cf_mx <- rf$finalModel$confusion[,1:2]
  g2_out <- cf_mx %>%
    reshape2::melt() %>%
    ggplot(aes(x = Var1, y = rev(Var2))) + geom_tile(aes(fill = log(value + 1))) +
    coord_equal() +
    scale_fill_gradient2(low = colors[1], high =colors[2],
                         midpoint = (log(cf_mx[1,1] + 1) +log(cf_mx[1,2] + 1))/2) +
    geom_text(aes(label = value)) +
    theme_minimal() +
    labs(x = 'Predicted', y = 'Reference') +
    theme(legend.position = '')

  pdf(paste0(out_path, 'roc_confusion.pdf'), width = width, height = height)
  print(g2_out)
  dev.off()

  print(g2_out)
  cat("\tDONE.\n")
}

## Code for creating figures from source data

library("tidyverse")
library("readxl")
library("igraph")
library("RColorBrewer")



###############################################################################
## Figure 1
###############################################################################

fig_1a <- read_xlsx("source_data.xlsx", sheet=1)

ggplot(fig_1a, aes(dist)) +
  geom_histogram(binwidth=0.03, colour="white", size=0.003, fill="royalblue3") + 
  theme_bw(base_size=7.5) + 
  labs(x="Squared distance") + 
  scale_y_continuous(breaks=c(0,25,50)) + 
  facet_wrap(~excluded, scales="fixed", ncol=4) + 
  theme(axis.text=element_text(size=5)) 
ggsave("figures_source/fig1a.pdf", width=10, height=8, scale=1, units="cm")


fig_1b <- read_xlsx("source_data.xlsx", sheet=2)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
ggplot(tibble(x=fig_1b$excluded, y=fig_1b$dist), aes(x, y)) +
  geom_line(size=0.4) + 
  geom_point(size=1.75, shape=21, fill="tomato", colour="white", stroke=0.75) +
  theme_bw(base_size=8) + 
  scale_y_continuous(label=scientific_10) + 
  labs(x="Number of excluded lipids", y="Average squared distance") + 
  geom_vline(aes(xintercept=9), linetype="dotted", size=0.4) + 
  geom_vline(aes(xintercept=13), linetype="dotted", size=0.4) + 
  theme(axis.text=element_text(size=6))   
ggsave("figures_source/fig1b.pdf", width=8, height=6, scale=1, units="cm")



###############################################################################
## Figure 2
###############################################################################

fig_2_ausdb <- read_xlsx("source_data.xlsx", sheet=3)
fig_2_lipid <- read_xlsx("source_data.xlsx", sheet=4)


## AusDiab panel ##############################################################

ausdb.pcor <- fig_2_ausdb %>% column_to_rownames("rownames") %>% as.matrix()

g <- graph_from_adjacency_matrix(
  ausdb.pcor,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE,
  add.colnames = NULL,
  add.rownames = NA
)

E(g)[E(g)$weight==1.0]$color = "dodgerblue1"
E(g)[!E(g)$weight==1.0]$color = "grey60"

lipid_family <- str_replace(V(g)$name, pattern="\\(.+\\)$", replacement="")
lipid_family <- str_replace(lipid_family, pattern="\\(.+\\)\\ \\[.+\\]$", replacement="")
length(unique(lipid_family)) # 19

V(g)$label <- NA
V(g)$lipid.family <- lipid_family

lipid.family.index <- factor(lipid_family, labels=as.character(1:19))
V(g)$lipid.family.index <- lipid.family.index

mycolours <- c(brewer.pal(9, "Set1")[1], 
               brewer.pal(9, "Set1")[3], 
               brewer.pal(3, "Greys")[1], 
               brewer.pal(9, "Set1")[2], 
               brewer.pal(8, "Dark2")[1], 
               brewer.pal(9, "Set1")[4], 
               brewer.pal(12, "Set3")[1], brewer.pal(12, "Set3")[7], brewer.pal(12, "Set3")[11], 
               brewer.pal(9, "Set1")[5], brewer.pal(12, "Set3")[6], brewer.pal(8, "Set2")[7], 
               brewer.pal(9, "Set1")[8], brewer.pal(12, "Set3")[12], brewer.pal(12, "Set3")[8], brewer.pal(12, "Set3")[2], 
               brewer.pal(9, "Set1")[7], brewer.pal(12, "Set3")[5], brewer.pal(12, "Set3")[9]) 

V(g)$color <- mycolours[V(g)$lipid.family.index]

modul.louvian <- cluster_louvain(g, weights=abs(E(g)$weight))

pdf(file = "figures_source/fig2ausdb.pdf", width = 7.25, height = 7.25)
plot(g, layout=layout_in_circle, 
     vertex.color=V(g)$color, vertex.frame.color="grey25", vertex.size=2.8, 
     edge.color=E(g)$color, edge.width=(E(g)$weight*1.25), main="AusDiab", 
     mark.groups=modul.louvian, mark.col=adjustcolor(c("lightsteelblue1"), alpha=.001), mark.border=NA, mark.shape=1, 
     margin=rep(-0.04, 4), asp=0)
legend(x=-1.35, y=+0.25, unique(lipid_family), pch=21,
       col="#777777", pt.bg=mycolours, pt.cex=2, cex=1.1, bty="n", ncol=1)
dev.off()


## LIPID panel ################################################################

lipid.pcor <- fig_2_lipid %>% column_to_rownames("rownames") %>% as.matrix()

g <- graph_from_adjacency_matrix(
  lipid.pcor,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE,
  add.colnames = NULL,
  add.rownames = NA
)

E(g)[E(g)$weight==1.0]$color = "dodgerblue1"
E(g)[!E(g)$weight==1.0]$color = "grey60"

lipid_family <- str_replace(V(g)$name, pattern="\\(.+\\)$", replacement="")
lipid_family <- str_replace(lipid_family, pattern="\\(.+\\)\\ \\[.+\\]$", replacement="")
length(unique(lipid_family)) # 19

V(g)$label <- NA
V(g)$lipid.family <- lipid_family

lipid.family.index <- factor(lipid_family, labels=as.character(1:19))
V(g)$lipid.family.index <- lipid.family.index

mycolours <- c(brewer.pal(9, "Set1")[1], 
               brewer.pal(9, "Set1")[3], 
               brewer.pal(3, "Greys")[1], 
               brewer.pal(9, "Set1")[2], 
               brewer.pal(8, "Dark2")[1], 
               brewer.pal(9, "Set1")[4], 
               brewer.pal(12, "Set3")[1], brewer.pal(12, "Set3")[7], brewer.pal(12, "Set3")[11], 
               brewer.pal(9, "Set1")[5], brewer.pal(12, "Set3")[6], brewer.pal(8, "Set2")[7], 
               brewer.pal(9, "Set1")[8], brewer.pal(12, "Set3")[12], brewer.pal(12, "Set3")[8], brewer.pal(12, "Set3")[2], 
               brewer.pal(9, "Set1")[7], brewer.pal(12, "Set3")[5], brewer.pal(12, "Set3")[9]) 

V(g)$color <- mycolours[V(g)$lipid.family.index]

modul.louvian <- cluster_louvain(g, weights=abs(E(g)$weight))

pdf(file = "figures_source/fig2lipid.pdf", width = 7.25, height = 7.25)
plot(g, layout=layout_in_circle, 
     vertex.color=V(g)$color, vertex.frame.color="grey25", vertex.size=2.8, 
     edge.color=E(g)$color, edge.width=(E(g)$weight*1.25), main="LIPID", 
     mark.groups=modul.louvian, mark.col=adjustcolor(c("lightsteelblue1"), alpha=.001), mark.border=NA, mark.shape=1, 
     margin=rep(-0.04, 4), asp=0)
legend(x=-1.35, y=+0.25, unique(lipid_family), pch=21,
       col="#777777", pt.bg=mycolours, pt.cex=2, cex=1.1, bty="n", ncol=1)
dev.off()



###############################################################################
## Figure 3
###############################################################################

fig_3a <- read_xlsx("source_data.xlsx", sheet=5)

my.colors <- c("royalblue", "firebrick1")
subset_label <- c(All="All predictors", PCorr13="Discordant predictors removed")
ggplot(data=fig_3a, 
       aes(x=ausdb, y=lipid)) + 
  geom_point(size=1.5, alpha=0.5, aes(colour=composite)) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  geom_abline(slope=1, intercept=-0.3, colour="grey50", size=0.25) + 
  geom_abline(slope=1, intercept=0.3, colour="grey50", size=0.25) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  scale_y_continuous(breaks=seq(0,1,0.1)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="AusDiab", y="LIPID") + 
  facet_wrap(~subset, scales="fixed", ncol=2, labeller=labeller(subset=subset_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=8), axis.text=element_text(size=6))
ggsave("figures_source/fig3a.pdf", width=12, height=6.6, scale=1, units="cm")


fig_3b <- read_xlsx("source_data.xlsx", sheet=6)

ggplot(fig_3b, aes(delta_improvement)) +
  geom_histogram(binwidth=0.005, colour="white", size=0.1, fill="firebrick2") + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(-0.05,0.15,0.05)) + 
  labs(title="Improvement in prediction accuracy", subtitle="(relative to AusDiab reference)", x="Improvement in correlation") + 
  theme(plot.title=element_text(size=8, hjust = 0.5), plot.subtitle=element_text(size=8, hjust = 0.5), axis.text=element_text(size=6))
ggsave("figures_source/fig3b.pdf", width=6, height=6, scale=1, units="cm")



###############################################################################
## Figure 4
###############################################################################

fig_4a <- read_xlsx("source_data.xlsx", sheet=7)

my.colors <- c("royalblue4", "seagreen4")
subset_label <- c(baseline="Baseline", followup="Follow up")
ggplot(data=fig_4a, 
       aes(x=nontreat, y=treat)) + 
  geom_point(size=1.5, alpha=0.5, aes(colour=timepoint)) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  scale_y_continuous(breaks=seq(0,1,0.1)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="Placebo randomised", y="Pravastatin randomised") + 
  facet_wrap(~timepoint, scales="fixed", ncol=2, labeller=labeller(timepoint=subset_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=8), axis.text=element_text(size=6), legend.position="none")
ggsave("figures_source/fig4a.pdf", width=12, height=7, scale=1, units="cm")


fig_4b <- read_xlsx("source_data.xlsx", sheet=8)

my.colors <- c("royalblue4", "seagreen4")
subset_label <- c(nontreat="Placebo randomised", treat="Pravastatin randomised")
ggplot(data=fig_4b, 
       aes(x=baseline, y=followup)) + 
  geom_point(size=1.5, alpha=0.5, aes(colour=randomisation)) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  scale_y_continuous(breaks=seq(0,1,0.1)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="Baseline", y="Follow up") + 
  facet_wrap(~randomisation, scales="fixed", ncol=3, labeller=labeller(randomisation=subset_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=8), axis.text=element_text(size=6), legend.position="none")
ggsave("figures_source/fig4b.pdf", width=12, height=7, scale=1, units="cm")


###############################################################################
## Figure 5
###############################################################################

fig_5a <- read_xlsx("source_data.xlsx", sheet=9)

target_pct_label <- c('10'="90% of full predictor set", '25'="75% of full predictor set", '50'="50% of full predictor set", '75'="25% of full predictor set")
ggplot(fig_5a, 
       aes(reorder(lipid, rep(r[1:294], 8)), r.1se, colour=data)) + 
  geom_point(size=1, alpha=0.5) + 
  geom_errorbar(aes(ymax=r.1se_max, ymin=r.1se_min), size=0.2) + 
  theme_bw(base_size=8) + 
  scale_colour_manual(values=c("dodgerblue", "tomato")) + 
  scale_x_discrete(breaks=NULL) + 
  labs(x="lipid species", y="Correlation (Range)") + 
  facet_wrap(~target_pct, scales="fixed", ncol=1, labeller=labeller(target_pct=target_pct_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=8), axis.text=element_text(size=6), legend.position="bottom")
ggsave("figures_source/fig5a.pdf", width=12, height=14, scale=1, units="cm")


fig_5b <- read_xlsx("source_data.xlsx", sheet=10)

my.colors <- c("royalblue", "firebrick1")
ggplot(data=fig_5b, 
       aes(x=r.1se_AusDiab, y=r.1se_LIPID)) + 
  geom_point(aes(colour=composite), size=1.25, alpha=0.5) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  geom_abline(slope=1, intercept=-0.3, colour="grey50", size=0.25) + 
  geom_abline(slope=1, intercept=0.3, colour="grey50", size=0.25) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  scale_y_continuous(breaks=seq(0,1,0.1)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="AusDiab", y="LIPID") + 
  facet_wrap(~target_pct, scales="fixed", ncol=2, labeller=labeller(target_pct=target_pct_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=8), axis.text=element_text(size=6))
ggsave("figures_source/fig5b.pdf", width=12, height=10.5, scale=1, units="cm")



###############################################################################
## Figure 6
###############################################################################

fig_6 <- read_xlsx("source_data.xlsx", sheet=11)

ggplot(fig_6, aes(max_corr)) +
  geom_histogram(binwidth=0.0125, colour="white", size=0.1, fill="firebrick2", alpha=0.9) + 
  geom_vline(xintercept=0.6, colour="dodgerblue", size=0.25) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  labs(x="Correlation")
ggsave("figures_source/fig6.pdf", width=9, height=6, scale=1, units="cm")



###############################################################################
## Figure 7
###############################################################################

# Figure 7: generated using Excel template (provided) and source data in sheet 12



###############################################################################
## Figure 8
###############################################################################

# Figure 8: generated using Excel template (provided) and source data in sheet 13



###############################################################################
## Figure 9
###############################################################################

fig_9 <- read_xlsx("source_data.xlsx", sheet=14)

ggplot(data=fig_9, 
       aes(x=AusDiab_test, y=SAFHS)) + 
  geom_point(size=1.75, alpha=0.6, colour="indianred2", stroke=0.4) + 
  geom_label_repel(aes(label=ifelse(abs(AusDiab_test-SAFHS)>0.2, lipid,'')), 
                   box.padding=0.2, label.padding=0.075, point.padding=0.25, label.r=0.1, label.size=0.1, max.time=5, force=2, force_pull=1, 
                   segment.color='grey50', segment.size=0.2, min.segment.length=0.2, size=1.75, nudge_x=-0.01, nudge_y=-0.01, max.overlaps=30) + 
  geom_abline(slope=1, intercept=0, colour="grey30", linewidth=0.4) + 
  geom_abline(slope=1, intercept=-0.2, colour="grey50", linewidth=0.15) + 
  geom_abline(slope=1, intercept=0.2, colour="grey50", linewidth=0.15) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(limits=c(0.08,1), breaks=seq(0.1,1,0.1)) + 
  scale_y_continuous(limits=c(0.08,1), breaks=seq(0.1,1,0.1)) + 
  labs(x="AusDiab", y="SAFHS") + 
  theme(axis.text=element_text(size=6))
ggsave("figures_source/fig9.pdf", width=9, height=9, scale=1, units="cm")







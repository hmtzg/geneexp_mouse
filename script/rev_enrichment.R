library(ggplot2)
library(tidyverse)
library(enrichplot)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(magick)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
ud_gse = readRDS("../../ageing/proc_data/reversal/ud.gse.tissues.rds")
du_gse = readRDS("../../ageing/proc_data/reversal/du.gse.tissues.rds")
ud_gse_shared = readRDS("../../ageing/proc_data/reversal/ud.gse.shared.rds")
du_gse_shared = readRDS("../../ageing/proc_data/reversal/du.gse.shared.rds")

ud_cortex = ud_gse$cortex[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Up-Down/Up-Up Genes in Cortex') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)
  
ggsave("./results/gsea_updown_cortex.pdf",ud_cortex, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_updown_cortex.png",ud_cortex, units = "cm", width = 10, height = 6)

ud_lung = ud_gse$lung[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Up-Down/Up-Up Genes in Lung') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.15,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_updown_lung.pdf",ud_lung, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_updown_lung.png",ud_lung, units = "cm", width = 10, height = 6)

ud_liver = ud_gse$liver[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Up-Down/Up-Up Genes in Liver') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.15,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_updown_liver.pdf",ud_liver, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_updown_liver.png",ud_liver, units = "cm", width = 10, height = 6)

ud_muscle = ud_gse$muscle[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Up-Down/Up-Up Genes in Muscle') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.15,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_updown_muscle.pdf",ud_muscle, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_updown_muscle.png",ud_muscle, units = "cm", width = 10, height = 6)

###
ud_shared = ud_gse_shared[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Shared Up-Down/Up-Up Genes Across Tissues') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_updown_shared.pdf",ud_shared, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_updown_shared.png",ud_shared, units = "cm", width = 10, height = 6)


### du gsea:

du_cortex = du_gse$cortex[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Down-Up/Down-Down Genes in Cortex') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_downup_cortex.pdf",du_cortex, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_downup_cortex.png",du_cortex, units = "cm", width = 10, height = 6)

du_lung = du_gse$lung[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Down-Up/Down-Down Genes in Lung') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_downup_lung.pdf",du_lung, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_downup_lung.png",du_lung, units = "cm", width = 10, height = 6)

du_liver = du_gse$liver[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Down-Up/Down-Down Genes in liver') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_downup_liver.pdf", du_liver, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_downup_liver.png", du_liver, units = "cm", width = 10, height = 6)

# no significant enrichment in muscle:
# du_muscle = du_gse$muscle[,1:10] %>%
#   filter(qvalues < 0.1) %>%
#   arrange(abs(NES)) %>%
#   group_by(sign(NES)) %>% 
#   slice(tail(row_number(),10)) %>%
#   ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
#   geom_bar(stat="identity", width = 0.7) +
#   scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
#   ggtitle('GSEA of Down-Up/Down-Down Genes in muscle') +
#   theme(legend.position = "right",
#         axis.text.y = element_text(size = 12/pntnorm),
#         axis.text.x = element_text(size = 12/pntnorm),
#         plot.title = element_text(size = 12/pntnorm),
#         legend.key.size = unit(.1,"cm")) +
#   #theme_minimal() + 
#   ylab(NULL)
# 
# ggsave("./results/gsea_downup_muscle.pdf",du_muscle, units = "cm", width = 10, height = 6, useDingbats=F)
# ggsave("./results/gsea_downup_muscle.png",du_muscle, units = "cm", width = 10, height = 6)

###
du_shared = du_gse_shared[,1:10] %>%
  filter(qvalues < 0.1) %>%
  arrange(abs(NES)) %>%
  group_by(sign(NES)) %>% 
  slice(tail(row_number(),10)) %>%
  ggplot( aes(NES,fct_reorder(Description, NES),fill = qvalues)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_continuous(low = "#2F2F2F", high = "#CF5C36", guide = guide_colorbar(reverse=TRUE)) +
  ggtitle('GSEA of Shared Down-Up/Down-Down Genes Across Tissues') +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12/pntnorm),
        axis.text.x = element_text(size = 12/pntnorm),
        plot.title = element_text(size = 12/pntnorm),
        legend.key.size = unit(.1,"cm")) +
  #theme_minimal() + 
  ylab(NULL)

ggsave("./results/gsea_downup_shared.pdf",du_shared, units = "cm", width = 10, height = 6, useDingbats=F)
ggsave("./results/gsea_downup_shared.png",du_shared, units = "cm", width = 10, height = 6)


# gsea_ud = ggarrange(ud_cortex, ud_lung, ncol = 2, nrow = 1, common.legend = F,align="v")
# ggsave("./results/gsea_ud_all.pdf", gsea_ud, units = "cm", width = 16, height = 6, useDingbats=F)
# ggsave("./results/gsea_ud_all.png", gsea_ud, units = "cm", width = 16, height = 6)









library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(ggrepel)

aero <- c("Subfossil Vampirovibrionales", "Modern Vampirovibrionales", "Fossil Sericytochromatia_2", "Modern Sericytochromatia", 
          "Modern Bacteriodetes", "Modern Chloroflexi", "Modern Thaumarchaeota")
no_aero <- c("Fossil Cyanobacteriales", "Modern Cyanobacteriales", "Modern Limnotrichales", "Fossil Gloeomargaritales", "Modern Gloeomargaritales")

quasar_normalized_data <- readr::read_csv("processed_data/spectral_quasar_normalized.csv") %>%
  dplyr::filter(sample != "string", sample != "meta") %>%
  dplyr::mutate(sample = ifelse(sample == "Fossil Sericytochromatia...28", "Fossil Sericytochromatia_1",
                                ifelse(sample == "Fossil Sericytochromatia...30", "Fossil Sericytochromatia_2", sample))) %>%
  dplyr::mutate(aerobic_resp = ifelse(sample %in% aero, "YES", 
                                      ifelse(sample %in% no_aero, "NO", "UNKNOWN"))) %>%
  dplyr::select(sample, aerobic_resp, everything())

quasar_normalized_data_filled <- quasar_normalized_data %>%
  tidyr::pivot_longer(cols = 3:last_col(), names_to = "wavelength", values_to = "value") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(wavelength = as.numeric(wavelength)) %>%
  dplyr::arrange(wavelength) %>%
  tidyr::fill(value, .direction = "down") %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = wavelength, values_from = value) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  dplyr::select(sample, aerobic_resp, seq(3, ncol(.), by = 2))


# ====================================================================================================================================================================== #
# Visualizing spectra:
# ====================================================================================================================================================================== #
original <- readr::read_csv("processed_data/midterm_raw_data_transposed.csv") %>%
  dplyr::mutate(sample = ifelse(sample == "Fossil Sericytochromatia...28", "Fossil Sericytochromatia (I)",
                                ifelse(sample == "Fossil Sericytochromatia...30", "Fossil Sericytochromatia (II)", sample))) %>%
  dplyr::mutate(aerobic_resp = ifelse(sample %in% aero, "YES", 
                                      ifelse(sample %in% no_aero, "NO", "UNKNOWN"))) %>%
  dplyr::select(sample, aerobic_resp, everything())

original_filled <- original %>%
  tidyr::pivot_longer(cols = 3:last_col(), names_to = "wavelength", values_to = "value") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(wavelength = as.numeric(wavelength)) %>%
  dplyr::arrange(wavelength) %>%
  tidyr::fill(value, .direction = "down") %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = wavelength, values_from = value) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  dplyr::select(sample, aerobic_resp, seq(3, ncol(.), by = 2)) %>%
  dplyr::mutate(aerobic_resp = factor(aerobic_resp, levels = c("YES","NO","UNKNOWN")))


og_spec <- original_filled %>%
  tidyr::pivot_longer(cols = 3:last_col(), names_to = "wavelength", values_to = "value") %>%
  mutate(
    wavelength = as.numeric(wavelength),
    value = as.numeric(value)
  ) %>%
  ggplot(aes(x = wavelength, y = value, color = aerobic_resp, group = sample)) +
  scale_color_manual(values = c("YES" = "blue", "NO" = "orange", "UNKNOWN" = "grey")) +
  geom_line(linewidth = 0.3) +
  # geom_point() +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.text = element_text(size = 7, color ='black'),
        legend.title = element_text(size = 8, color = 'black'),
        legend.position = "inside",
        legend.position.inside = c(0.2,0.83),
        legend.key.height = unit(0.3, "cm"), 
        # legend.key.width = unit(0.8, "cm"), 
        legend.spacing.y = unit(0.1, "cm"), 
        legend.key = element_blank()) +
  labs(x = "Raman shift", y = "Intensity", color = "Aerobic respiration") +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2))) 
# og_spec


norm_spec <- quasar_normalized_data_filled %>%
  tidyr::pivot_longer(cols = 3:last_col(), names_to = "wavelength", values_to = "value") %>%
  mutate(
    wavelength = as.numeric(wavelength),
    value = as.numeric(value)
  ) %>%
  dplyr::mutate(aerobic_resp = factor(aerobic_resp, levels = c("YES","NO","UNKNOWN"))) %>%
  ggplot(aes(x = wavelength, y = value, color = aerobic_resp, group = sample)) +
  scale_color_manual(values = c("YES" = "blue", "NO" = "orange", "UNKNOWN" = "grey")) +
  geom_line(linewidth = 0.3) +
  # geom_point() +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = "none") +
  labs(x = "Raman shift", y = "Normalized intensity", color = "Aerobic respiration")
# norm_spec


spec_plt <- cowplot::plot_grid(
  og_spec, norm_spec,
  nrow = 1,
  rel_widths = c(1,1),
  labels = c("a","b"),
  align = "h")
spec_plt

ggsave("plots/spectra.png", spec_plt, width = 7.5, height = 3.5, dpi = 500)



# ====================================================================================================================================================================== #
# PCA
# ====================================================================================================================================================================== # 
num_cols <- 3:ncol(quasar_normalized_data_filled)
quasar_normalized_data_filled[num_cols] <- lapply(quasar_normalized_data_filled[num_cols], as.numeric)

pca_result <- prcomp(quasar_normalized_data_filled[num_cols], center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  sample = quasar_normalized_data_filled$sample,
  aerobic_resp = quasar_normalized_data_filled$aerobic_resp,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4]) %>%
  dplyr::mutate(aerobic_resp = factor(aerobic_resp, levels = c("YES","NO","UNKNOWN")))

set.seed(1)
pca_df <- pca_df %>%
  dplyr::mutate(PC1_j = jitter(PC1, amount = 0.5),
                PC2_j = jitter(PC2, amount = 0.5))

var_explained <- summary(pca_result)$importance[2, ] * 100

pca <- ggplot(pca_df, aes(x = PC1_j, y = PC2_j, color = aerobic_resp)) +
  geom_point(size = 5, alpha = 0.7) +
  geom_text_repel(data = dplyr::filter(pca_df, aerobic_resp == "UNKNOWN"), aes(label = sample), color = "black", size = 4, max.overlaps = Inf) +
  scale_color_manual(values = c("YES" = "blue", "NO" = "orange", "UNKNOWN" = "grey")) +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 11, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.text = element_text(size = 10, color ='black'),
        legend.title = element_text(size = 10, color = 'black'),
        legend.position = "inside",
        legend.position.inside = c(0.86,0.9)) +
  labs(x = paste0("PC1 (", round(var_explained[1], 1), "%)"), y = paste0("PC2 (", round(var_explained[2], 1), "%)"), color = "Aerobic respiration") 
pca

ggsave("plots/PCA.png", pca,  width = 7.5, height = 7.5, dpi = 500)

# ====================================================================================================================================================================== #
# Final classification
# ====================================================================================================================================================================== #

classification <- pca_df %>% dplyr::select(sample, aerobic_resp, PC2) %>%
  dplyr::arrange(desc(PC2)) %>%
  dplyr::mutate(inferred = ifelse(aerobic_resp == "UNKNOWN", "YES", "NO")) %>%
  dplyr::mutate(aerobic_resp = ifelse(PC2 > 0, "aerobic", "anaerobic")) %>%
  dplyr::select(-PC2)

write.table(classification, "processed_data/aerobic_species_consensus.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

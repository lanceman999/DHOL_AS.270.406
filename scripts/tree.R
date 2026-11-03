library(ape)
library(phytools)

# Tip names must match your spectral sample names (no spaces)
FSC1 <- "Fossil_Stem_Cyanobacteria_I*";   FSC2 <- "Fossil_Stem_Cyanobacteria_II*"
FCy3 <- "Fossil_Cyanobacteria_III*";      FCya <- "Fossil_Cyanobacteriales"
MCya <- "Modern_Cyanobacteriales"
FLim <- "Fossil_Limnotrichales*";         MLim <- "Modern_Limnotrichales"
FSG  <- "Fossil_Stem_Gloeomargaritales*"; FGlo <- "Fossil_Gloeomargaritales"
MGlo <- "Modern_Gloeomargaritales"
FMel <- "Fossil_Melainabacteria*";        SubV <- "Subfossil_Vampirovibrionales"
MVam <- "Modern_Vampirovibrionales"
FSer <- "Fossil_Sericytochromatia_1*";      FSerS <- "Fossil_Sericytochromatia_2"
MSer <- "Modern_Sericytochromatia"
MBac <- "Modern_Bacteroidetes";          MChl <- "Modern_Chloroflexi"
FDei <- "Fossil_Stem_Deinococcus*"
FTha <- "Fossil_Stem_Thaumarcheota*";     MTha <- "Modern_Thaumarcheota"

cyano_r <- sprintf("((%s,(%s,%s)),(%s,%s))", FCy3, MCya, FCya, FSC2, FSC1)
limno_r <- sprintf("((%s,%s),%s)", MLim, FLim, cyano_r)
gloeo_r <- sprintf("((%s,%s),%s)", MGlo, FGlo, FSG)
upper_r <- sprintf("(((%s,%s),%s),(%s,%s))", MVam, SubV, FMel, gloeo_r, limno_r)
lower_r <- sprintf("(%s,%s,%s,((%s,%s),%s),%s)", FDei, MChl, MBac, MSer, FSerS, FSer, upper_r)
newick_flip <- sprintf("((%s,%s),%s);", MTha, FTha, lower_r)

# Making the blueprint tree!
tree <- read.tree(text = newick_flip)
plot(tree, cex = 0.6, no.margin = TRUE)
nodelabels(cex = 0.5)


# Encoding binary traits of aerobic versus non-aerobic respiration
aerobic_tips <- c(FSC1, FLim, FSer, FSerS, FSG, MVam, FMel, SubV,
                  MBac, MSer, MChl, FDei, FTha, MTha)

stopifnot(all(aerobic_tips %in% tree$tip.label))   # checking that all species exist in the tree

states <- setNames(as.integer(tree$tip.label %in% aerobic_tips), tree$tip.label)
table(states)   # expect 14 ones (aerobic) and 7 zeros (anaerobic)

# Topology only, agnostic to branch length - all branch lengths are set to 1
tree_eq <- compute.brlen(tree, 1)

fit_er  <- fitMk(tree_eq, states, model = "ER")
fit_ard <- fitMk(tree_eq, states, model = "ARD")
AIC(fit_er); AIC(fit_ard)

anc   <- ancr(fit_er)
names(anc)
probs <- anc$ace

# Confirm the matrix has one row per internal node and columns "0" and "1"
dim(probs)          
colnames(probs)       
round(head(probs), 2)

# Save the plot
H <- max(node.depth.edgelength(tree_eq))   # distance from root to deepest tip
H

png("plots/tree.png", width = 7.5, height = 7.5, units = "in", res = 500)

par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = NA)

plot(tree_eq, cex = 0.6, align.tip.label = TRUE, font = 1, label.offset = 0.2,
     x.lim = c(0, H * 1.25))   # increase 1.35 if the longest label is still clipped

ape::nodelabels(node = as.integer(rownames(probs)), pie = probs,
                piecol = c("orange", "blue"), cex = 0.5)
legend("topleft", legend = c("Anaerobic", "Aerobic"),
       fill = c("orange", "blue"), bty = "n", cex = 0.8)
tip_cols <- ifelse(states[tree_eq$tip.label] == 1, "blue", "orange")
tiplabels(pch = 22, bg = tip_cols, col = "black", cex = 1.5)

dev.off()



#step one 
# the essentials ----
# this chunk contains the minimal essential code from this script. Simply uncomment the lines below and run the code.
library(rhdf5)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(biomaRt)
# library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
packageVersion("tidyverse") 
packageVersion("tidyverse") 

targets <- read_tsv("studydesign515.txt")# read in your study design
targets
path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data
path
all(file.exists(path))
lyc.anno <-useMart(biomart="plants_mart", dataset ="slycopersicum_eg_gene", host="plants.ensembl.org")
Tx.lyc <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = lyc.anno) %>%
  as_tibble() %>%
  dplyr::rename(target_id = ensembl_transcript_id, gene_name=ensembl_gene_id) %>% 
  dplyr::select("target_id", "gene_name")

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.lyc, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = FALSE)
#path
#dir()
Txi_gene
sessionInfo()


#######2#222222222222222222222222222222222222222222222#
# the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = CG01:TG02, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=4 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = CG01:TG02, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = CG01:TG02, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

###################33#333333333333333333333333333333333#
# Introduction to this script -----------
# this script walks thorough techniques for data exploration and expands on last week's data wrangling theme
# we'll also continue to create publication-quality graphics
# This script starts with your filtered and normalized abundance data from the Step 2 script.

# Lets set a project-specific library
# Sys.unsetenv("R_LIBS_USER")
# dir.create("RLibrary")
# .libPaths()
# .libPaths(paste(getwd(), "RLibrary", sep="/"))
# setRepositories()

# install.packages('DT')
# install.packages('plotly')
# install.packages('gt')


# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables


# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)

# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
log2.cpm.filtered.norm.df

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) + # color = group
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

# Let's discuss and iteratively refine the PCA code and plot from above
# First, take note of the fact that we can use information from our PCA analysis to label our axes
# Remember that PCA is unsupervised, so knows nothing about group assignment (healthy vs disease)
# But *we* know, and so we can use this knowledge to enhance the plot.  Add a 'color=group' mapping to the aes of the plot above
# Can we figure out the identity of the outlier?  We have already provided samplelabel mapping in aes, so just uncomment the 'geom_label()'
# Uncomment 'coord_fixed()' to apply the correct aspect ratio
# Uncomment 'stat_ellipse()' to see how you can circle clusters on the PCA
# How would this PCA look if you used raw counts (myCounts) instead of log2 CPM?
# What are the disadvantages of looking at a PCA result using such a simple XY plot?

# head back to slides

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

# head back to the slides

# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data
mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(Control.AVG = (CG01+CG02)/2,
         treated.AVG = (TG01+CG02)/2,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = ( treated.AVG- Control.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

# Use dplyr "filter" and "select" functions to pick out genes of interest 
# ways to tweak the 'select' function:
# use ':' between two column names to select all columns between
# use 'contains', 'starts_with' or 'ends_with' to modify how you select
# can refer to columns using exact name or numerical indicator
# use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="MMP1" | geneID=="GZMB" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG"
                | geneID=="CCL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A" ) %>%
  dplyr::select(geneID, Control.AVG, treated.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, Control.AVG, treated.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))

# head back to the slides

# Produce publication-quality tables using the gt package ----
gt(mydata.filter)
# now with a few more options
mydata.filter %>%
  gt() %>%
  fmt_number(columns=1:4, decimals = 1) %>%
  tab_header(title = md("**Regulators of skin pathogenesis**"),
             subtitle = md("*during cutaneous leishmaniasis*")) %>%
  tab_footnote(
    footnote = "Deletion or blockaid ameliorates disease in mice",
    locations = cells_body(
      columns = geneID,
      rows = c(6, 7))) %>% 
  tab_footnote(
    footnote = "Associated with treatment failure in multiple studies",
    locations = cells_body(
      columns = geneID,
      rows = c(2:9))) %>%
  tab_footnote(
    footnote = "Implicated in parasite control",
    locations = cells_body(
      columns = geneID,
      rows = c(2))) %>%
  tab_source_note(
    source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))

mydata.df
# Make an interactive table using the DT package ----
datatable(mydata.df[,c(1,2:8)],  ### CHECK YOUR DATAFRAME COLUMSN
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

# Make an interactive scatter plot with plotly -----
# begin by storing your ggplot object
myplot <- ggplot(mydata.df) + 
  aes(x=Control.AVG, y=treated.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("treated vs. Control") +
  theme_bw()

#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive plot
ggplotly(myplot)

#let's customize this graphic by adding a more informative mouseover tooltip
myplot <- ggplot(mydata.df) +
  aes(x=Control.AVG, y=treated.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("treated vs. Control") +
  theme_bw()

ggplotly(myplot)

# the essentials ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

group <- targets$group
group <- factor(group)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    Control.AVG = (CG01 + CG02)/2, 
                    treated.AVG = (TG01+TG02)/2,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (treated.AVG - Control.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,2:8)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))


#######################4###################################
# Introduction to this script -----------
# the goal of this script is to identify differentially expressed genes (DEGs) and differential transcript usage (DTU)
# you should already know which pairwise comparisons are most important to you
# whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
# if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
# instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 and 4 to identify genes based on log fold-change

# Lets set a project-specific library
# Sys.unsetenv("R_LIBS_USER")
# dir.create("RLibrary")
# .libPaths()
# .libPaths(paste(getwd(), "RLibrary", sep="/"))
# setRepositories()

# Install new packages -----
# install.packages("reshape2")
# install.packages('heatmaply')


# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(reshape2)
library(heatmaply)

# Set up your design matrix ----
group <- factor(targets$group)
design <- model.matrix(~  group)

colnames(design) <- levels(group)

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment)
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)
# Design matrix and linear model: https://youtu.be/R7xd624pR1A

# Contrast matrix ----
contrast.matrix <- makeContrasts(treted =treatment - control,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)
# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
ggplot(myTopHits.df) + # vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  # geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  # geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  # geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  # annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  # annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "slycopersicum",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=4)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="both")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# i missed one Column 
# create interactive tables to display your DEGs ----
#datatable(diffGenes.df,
#  extensions = c('KeyTable', "FixedHeader"),
#   caption = 'Table 1: DEGs in slycopersicum',
#  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
#formatRound(columns=c(2:11), digits=2)
#### fix 
# create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in slycopersicum',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = 2:ncol(diffGenes.df), digits = 2)


####last fix 


#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)


# Create a heatmap of differentially expressed genes ----


##
heatmaply(diffGenes.df[2:5], 
          #dendrogram = "row",
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in slycopersicum",
          scale = "column",
          margins = c(60, 100, 40, 20),
          grid_color = "blue",
          grid_width = 0.0000001,
          titleX = TRUE,
          titleY = TRUE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffGenes.df)[2:5],
          labRow = diffGenes.df$geneID,
          heatmap_layers = theme(axis.line = element_blank()),
          col = colorRampPalette(c("blue", "white", "red"))(50)  # Adjust the colors and number of shades if needed
)




# 1.4  the essentials ----

library(tidyverse)
library(limma)

library(edgeR)
library(gt)
library(DT)
library(plotly)

group <- factor(targets$group)
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
#contrast.matrix <- makeContrasts(infection = disease - healthy,
#  levels=design)#
#Create the design matrix
design <- model.matrix(~0 + group)

# Assign column names to the design matrix
colnames(design) <- levels(group)

# Create a contrast matrix
contrast.matrix <- makeContrasts(treatment - control, levels = design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")


vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels


diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

###4. 2  fix all 
# the essentials ----

library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(heatmaply)

# ... (rest of your code)

# Assume that you have already created diffGenes.df

# Create a heatmap of differentially expressed genes ----
heatmaply(as.matrix(diffGenes.df[, 2:ncol(diffGenes.df)]), 
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in ",
          scale = "column",
          margins = c(60, 100, 40, 20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = TRUE,
          titleY = TRUE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffGenes.df)[2:ncol(diffGenes.df)],
          labRow = diffGenes.df$geneID,
          heatmap_layers = theme(axis.line = element_blank())
)



# The rest of my code for the Volcano plot and interactive table

# The Volcano plot
vplot <- ggplot(myTopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "grey", size = 1) +
  geom_vline(xintercept = 1, linetype = "longdash", colour = "#BE684D", size = 1) +
  geom_vline(xintercept = -1, linetype = "longdash", colour = "#2C467A", size = 1) +
  labs(title = "Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw()

# Convert the ggplot object to a plotly object
vplotly <- ggplotly(vplot)

# Display the interactive Volcano plot
vplotly

# The interactive table
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = 2:11, digits = 2)
# try again heatmap
heatmaply(as.matrix(diffGenes.df[, 2:ncol(diffGenes.df)]), 
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in ",
          scale = "column",
          grid_color = "white",
          titleX = TRUE,
          titleY = TRUE,
          hide_colorbar = TRUE,
          labCol = colnames(diffGenes.df)[2:ncol(diffGenes.df)],
          labRow = diffGenes.df$geneID,
          fontsize_row = 8,
          plot_method = "plotly",
          colors = c("blue", "white", "red"),  # Adjust the colors if needed
          show_hist = FALSE,
          heatmap_layers = theme(axis.line = element_blank())
)
####
#2323#
# the 3.4 essentials ----

library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(heatmaply)

# ... (rest of your code)

# Assume that you have already created diffGenes.df

# Create a heatmap of differentially expressed genes ----
heatmaply(as.matrix(diffGenes.df[, 2:ncol(diffGenes.df)]), 
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in slycopersicum",  # Updated title
          scale = "column",
          grid_color = "white",
          titleX = TRUE,
          titleY = TRUE,
          hide_colorbar = TRUE,
          labCol = colnames(diffGenes.df)[2:ncol(diffGenes.df)],
          labRow = diffGenes.df$geneID,
          fontsize_row = 8,
          plot_method = "plotly",
          colors = c("blue", "white", "red"),  # Adjust the colors if needed
          show_hist = FALSE,
          heatmap_layers = theme(axis.line = element_blank())
)

# ... (rest of your code)

# The Volcano plot
vplot <- ggplot(myTopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "grey", size = 1) +
  geom_vline(xintercept = 1, linetype = "longdash", colour = "#BE684D", size = 1) +
  geom_vline(xintercept = -1, linetype = "longdash", colour = "#2C467A", size = 1) +
  labs(title = "Volcano plot - DEGs in slycopersicum",  # Updated title
       subtitle = "Cutaneous leishmaniasis",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw()

# Convert the ggplot object to a plotly object
vplotly <- ggplotly(vplot)

# Display the interactive Volcano plot
vplotly

# The interactive table
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in slycopersicum',  # Updated title
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = 2:11, digits = 2)


###################################5555555555555555555555555555###
# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(qusage) # Quantitative Set Analysis for Gene Expression
library(heatmaply)


# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis, the list of available species for the "organism" parameters is here: https://biit.cs.ut.ee/gprofiler/page/organism-list 
gost.res <- gost(rownames(myTopHits), organism = "slycopersicum", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = F, capped = T) #set interactive=FALSE to get plot for publications
# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'

###
# List of packages
packages <- c(
  "rhdf5", "tidyverse", "tximport", "ensembldb", "biomaRt",
  "edgeR", "matrixStats", "cowplot", "DT", "plotly", "gt", "limma", 
  "gplots", "heatmaply", "GSEABase", "Biobase", "GSVA", "gprofiler2",
  "clusterProfiler", "msigdbr", "enrichplot", "qusage"
)

# Function to get package version
get_package_version <- function(package_name) {
  package_version <- tryCatch(packageVersion(package_name), error = function(e) NA)
  return(package_version)
}

# Loop through packages and print name along with version
for (package in packages) {
  version <- get_package_version(package)
  cat(sprintf("Package: %s, Version: %s\n", package, version))
}

# List of packages
packages <- c(
  "rhdf5", "tidyverse", "tximport", "ensembldb", "biomaRt",
  "edgeR", "matrixStats", "cowplot", "DT", "plotly", "gt", "limma", 
  "gplots", "heatmaply", "GSEABase", "Biobase", "GSVA", "gprofiler2",
  "clusterProfiler", "msigdbr", "enrichplot", "qusage"
)

# Function to get package version
get_package_version <- function(package_name) {
  package_version <- tryCatch(packageVersion(package_name), error = function(e) NA)
  return(package_version)
}

# Loop through packages and print name along with version
for (package in packages) {
  version <- get_package_version(package)
  if (is.na(version)) {
    cat(sprintf("Package: %s, Version: Not Available\n", package))
  } else {
    cat(sprintf("Package: %s, Version: %s\n", package, version))
  }
}

# Create a data frame with package information
package_data <- data.frame(
  Package = c("rhdf5", "tidyverse", "tximport", "ensembldb", "biomaRt", "edgeR", 
              "matrixStats", "cowplot", "DT", "plotly", "gt", "limma", 
              "gplots", "heatmaply", "GSEABase", "Biobase", "GSVA", 
              "gprofiler2", "clusterProfiler", "msigdbr", "enrichplot", "qusage"),
  Version = c("2.44.0", "2.0.0", "1.28.0", "2.24.1", "2.56.1", "3.42.4", 
              "1.0.0", "1.1.1", "0.29", "4.10.2", "0.9.0", "3.56.2", 
              "3.1.3", "1.4.2", "1.62.0", "2.60.0", "1.48.3", 
              "0.2.2", "4.8.3", "7.5.1", "1.20.3", "2.34.0")
)

# Print the data frame
print(package_data)

write.csv(package_data, "package_versions.csv", row.names = FALSE)
# Save the data frame to a CSV file on the desktop
write.csv(package_data, file.path(Sys.getenv("HOME"), "Desktop", "package_versions.csv"), row.names = FALSE)


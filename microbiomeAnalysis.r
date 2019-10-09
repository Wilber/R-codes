#Use this for for subsequent analyses
#With R version 3.5 or greater, \
 #install Bioconductor packages using BiocManager;\
#see https://bioconductor.org/install 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("phyloseq")

setwd("/Users/ouma.2/qiime/")

#Adapted from: https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html
#Adapted from https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html
#git gub page: https://github.com/Wilber/Microbial-bioinformatics-introductory-course-Material-2018 

# Create Folders as following
# Tables
dir.create("tables")
# Figures
dir.create("figures")
# Phyloseq objects
dir.create("phyobjects")
# Custom codes/notes
dir.create("codes_notes")

#BiocManager::install("microbiome")
#BiocManager::install(c("microbiomeutilities","RColorBrewer","ggpubr"))

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling


#####1: Read input to phyloseq object
#import the biom table (has both OTU table and taxonomy). Make sure it's csv not tsv. 
ps.ng.tax <- read_phyloseq(otu.file = "table-with-taxonomy.biom", 
                           taxonomy.file = NULL, 
                           metadata.file = "sample-metadata.csv", 
                           type = "biom") 
###2: Read tree file:
# Load tree file
library(ape)
treefile_p1 <- read.tree("tree.nwk")

###3: Merge into phyloseq object.
ps.ng.tax <- merge_phyloseq(ps.ng.tax, treefile_p1)
# ps.ng.tax is the first phyloseq object.

print(ps.ng.tax)
rank_names(ps.ng.tax) # we check the taxonomic rank information

otu_table(ps.ng.tax)
tax_table(ps.ng.tax)
sample_data(ps.ng.tax)
datatable(tax_table(ps.ng.tax)) # the table is interactive you can scrol and search thorugh it for details.

#Clean the taxonomy table.
tax_table(ps.ng.tax)[, colnames(tax_table(ps.ng.tax))] <- gsub(tax_table(ps.ng.tax)[, colnames(tax_table(ps.ng.tax))], pattern = "[a-z]__", replacement = "")
tax_table(ps.ng.tax)[tax_table(ps.ng.tax)[, "Phylum"] == "", "Phylum"] <- "Unidentified"

p1 <- plot_taxa_cv(ps.ng.tax, plot.type = "scatter")
p1 + scale_x_log10()##4: #distribution of taxa
# We make a data table with information on the OTUs
ps0_df_taxa <- data.table(tax_table(ps.ng.tax),
                          ASVabundance = taxa_sums(ps.ng.tax),
                          ASV = taxa_names(ps.ng.tax))

ps0_tax_plot <- ggplot(ps0_df_taxa, aes(ASVabundance)) +
  geom_histogram() + ggtitle("Histogram of ASVs (unique sequence) counts (NG-Tax)") +
  theme_bw() + scale_x_log10() + ylab("Frequency of ASVs") + xlab("Abundance (raw counts)")

print(ps0_tax_plot) ##slightly normal dist.?Will affect choise of statistic test later on 

#check for prevalance abundance distribtuion of ASVs and OTUs.
p <- plot_taxa_prevalence(ps.ng.tax, "Phylum")
p

###sequence of mitochondrial and chloroplast origin should be present is low abundance.
ps1.1.ngtax <- subset_taxa(ps.ng.tax, Family != "Mitochondria")
# also check how many lost 
ntaxa(ps.ng.tax)-ntaxa(ps1.1.ngtax)

# save the pseq object
saveRDS(ps.ng.tax, "./phyobjects/ps.ng.tax.rds")
#############################
######Alpha diversities:
#One has to consider the sequencing depth (how much of the taxa have been sampled) for each sample. 
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")
# use print option to see the data saved as phyloseq object.
print(ps1)
#get the dist. of number of reads for all the samples
summary(sample_sums(ps1)) 

#rarefaction curve for the observed ASVs in the entire data set.
otu_tab <- t(abundances(ps1))
p <- vegan::rarecurve(otu_tab, 
                      step = 50, label = FALSE, 
                      sample = min(rowSums(otu_tab), 
                                   col = "blue", cex = 0.6))
#Almost all samples are reaching a plateau and few samples have high number of reads and high number of ASVs.

#We will normalize to the lowest depth of at least 917 reads to keep maximum samples for our anlaysis. 
#Analogous to the sampling depth of sub-sampling???

### Equal sample sums
set.seed(9242)  # This will help in reproducing the filtering and nomalisation. 
ps0.rar <- rarefy_even_depth(ps1, sample.size = 917) #PLAY AROUND WITH SAMPLE SIZE!
saveRDS(ps0.rar, "./phyobjects/ps0.rar.rds")
#Check how much data you have now
ps0.rar <- readRDS("./phyobjects/ps0.rar.rds")
print(ps0.rar)

# quick check for sampling depth
barplot(sample_sums(ps0.rar), las =2)

# quick check taxa prevalence
p.rar <- plot_taxa_prevalence(ps0.rar, "Phylum")
p.rar

####Diversities
##Non-phylogenetic diversities
#For more diversity indices please refer to Microbiome Package
#Goal: calculate alpha diversities, add them to the metadata 
hmp.div <- diversities(ps0.rar, index = "all")
datatable(hmp.div)
datatable(hmp.div)

##plot
# get the metadata out as seprate object
hmp.meta <- meta(ps0.rar)
# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)
# Add the rownames to diversity table
hmp.div$sam_name <- rownames(hmp.div)
# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
# check the tables
colnames(div.df)
head(div.df)
# Now use this data frame to plot 
p <- ggboxplot(div.df, 
               x = "BodySite", 
               y = "diversity_shannon",
               fill = "BodySite", 
               palette = "jco")
print(p)
#OR rotate the x-axis labels/text
p <- p + rotate_x_text()
print(p)

##Alternatively, plot dist of several diversity metrices:
# convert phyloseq object into a long data format.
head(div.df)
div.df2 <- div.df[, c("BodySite", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "diversity_fisher", "diversity_coverage")]
# the names are not pretty. we can replace them
colnames(div.df2) <- c("BodySite", "Inverse Simpson", "Gini-Simpson", "Shannon", "Fisher", "Coverage")
# check
colnames(div.df2)
head(div.df2)

div_df_melt <- reshape2::melt(div.df2)
## Using Location as id variables
head(div_df_melt) #The diversity indices are stored under column named variable.

# Now use this data frame to plot 
p <- ggboxplot(div_df_melt, x = "BodySite", y = "value",
               fill = "BodySite", 
               palette = "jco", 
               legend= "right",
               facet.by = "variable", 
               scales = "free")
p
p <- p + rotate_x_text() 
p
# we will remove the x axis lables
p <- p + rremove("x.text")
p
#save the plot 
ggsave("./figures/Diversities.pdf", height = 4, width = 10)

##plot the dist. of alpha diversity & richness metrics:
#Create 2x2 plot environment so that we can see all 4 metrics at once. 
par(mfrow = c(2, 2))

#Then plot each metric.
hist(div.df$diversity_shannon, main="Shannon diversity", xlab="", breaks=10)
hist(div.df$evenness_simpson, main="Simpson diversity", xlab="", breaks=10)
hist(div.df$chao1, main="Chao richness", xlab="", breaks=15)
hist(div.df$diversity_fisher, main="Fisher", xlab="", breaks=15)

###Plot with pair-wise comparison stats/p-value for different levels of Body site:
lev <- levels(div_df_melt$BodySite) # get the variables

# make a pairwise list that we want to compare.
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

pval <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***", "**", "*", "n.s")
)

###pair-wise comparison of means:
#do we need to check for normality of data first? Looks like it's using wilcox.test: warnings()
#To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.
shapiro.test(div.df$diversity_inverse_simpson) #not normal
shapiro.test(div.df$diversity_gini_simpson) #normal??
shapiro.test(div.df$diversity_shannon) #normal?
shapiro.test(div.df$diversity_fisher) #not normal
shapiro.test(div.df$diversity_coverage) #not normal

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

print(p2)

####Phylogenetic diversity
##calculated using the 'picante' package:
library(picante)
ps0.rar.asvtab <- as.data.frame(ps0.rar@otu_table)
ps0.rar.tree <- ps0.rar@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 

ps0.rar@phy_tree
# it is a rooted tree
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
datatable(df.pd)

# now we need to plot PD

# We will add the results of PD to this file and then plot.
hmp.meta$Phylogenetic_Diversity <- df.pd$PD
hmp.meta

pd.plot <- ggboxplot(hmp.meta,
                     x = "BodySite",
                     y = "Phylogenetic_Diversity",
                     fill = "BodySite",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Body site",
                     legend = "right"
)
pd.plot <- pd.plot + rotate_x_text()

pd.plot + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

#NOTE:
#There are arguments both for and against the use of rarefying to equal library size.
#The application of normalization method will depend on the type of research question. It is always good \
#to check if there is a correlation between increasing library sizes and richness. Observed ASVs and \
#Phylogenetic diversity can be affected by library sizes. It is always good to check for this before making a choice.

lib.div <- diversities(ps1, index = "all")
lib.div2 <- richness(ps1) #i.e chao1 


# let us add number of total reads/samples
lib.div$ReadsPerSample <- sample_sums(ps1)
lib.div$Observed <- lib.div2$observed

colnames(lib.div)

#use ggscatter function from ggpubr package to visualze.
p1 <- ggscatter(lib.div, "diversity_shannon", "ReadsPerSample") +
  stat_cor(method = "pearson")

p2 <- ggscatter(lib.div, "diversity_inverse_simpson", "ReadsPerSample",
                add = "loess"
) +
  stat_cor(method = "pearson")

p3 <- ggscatter(lib.div, "Observed", "ReadsPerSample",
                add = "loess") +
  stat_cor(
    method = "pearson",
    label.x = 100,
    label.y = 50000
  )

ggarrange(p1, p2, p3, ncol = 2, nrow = 2)

##Statistics!!!
##Does BodySite impact the Shannon diversity of the microbiota??
#We choose Shannon because it is slightly normal, see stat test above. 
#Run the ANOVA and save it as an object
aov.shannon.BodySite = aov(diversity_shannon ~ BodySite, data=div.df)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.BodySite)
#To do all the pairwise comparisons between groups and \
#correct for multiple comparisons, we run Tukey’s honest significance test of our ANOVA
TukeyHSD(aov.shannon.BodySite) ##These results are similar to the ones in the box plots above!

#######For non-normal data (eg Gini-Simpson, see above)
#Since BodySite is categorical, we use Kruskal-Wallis \
#(non-parametric equivalent of ANOVA). \
#If we have only two levels, we would run Wilcoxon rank sum test (non-parametric equivalent of t-test)
kruskal.test(diversity_gini_simpson ~ BodySite, data=div.df)

#Test pairwise within the age groups with Wilcoxon Rank Sum Tests.
pairwise.wilcox.test(div.df$diversity_gini_simpson, div.df$BodySite, p.adjust.method="fdr")
#Above results similar to the boxplots stats test 
###########
###ASSIGNMENT: TEST EFFECT OF ANTI-BIOTIC USE ON MICROBIAL DIVERSITY (ACCOUNT FOR BODY SITE?????)

#######For normal continous data (eg Month??)
shapiro.test(div.df$Month)
glm.shannon.Month = glm(diversity_shannon ~ Month, data=div.df)
summary(glm.shannon.Month)
plot(diversity_shannon ~ Month, data=div.df)
#Add the glm best fit line
abline(glm.shannon.Month)


###############
#####COMPOSITION PLOTS:
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")
# use print option to see the data saved as phyloseq object.
print(ps1)

ps1.com <- ps1

# if you have dada2/deblur output and sequences as taxa names, then you can change them as follows
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

# colourCount = length(unique(taxic$Family))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object


# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL

# Second merge at family level

ps1.com.fam <- microbiome::aggregate_top_taxa(ps1.com, "Family", top = 10)

plot.composition.COuntAbun <- plot_composition(ps1.com.fam) + theme(legend.position = "bottom") +
  scale_fill_brewer("Family", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

plot.composition.COuntAbun 

#The above plot is based on the reads per sample. In the next step, we plot the relative abundance.
# the previous pseq object ps1.com.fam is only counts.
# Use traqnsform function of microbiome to convert it to rel abun.
ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")

plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "BodySite",
                                             x.label = "env_material") 
plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Family", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)

# Barplot customize
data.com <- plot.composition.relAbun$data
colnames(data.com)
head(data.com)
p.com <- ggplot(data.com, aes(x = Sample, y = Abundance, fill = Tax))
p.com <- p.com + geom_bar(position = "stack", stat = "identity")
p.com <- p.com + scale_x_discrete(labels = data.com$xlabel, breaks = data.com$Sample)
p.com <- p.com + facet_grid(~xlabel, scales = "free") + theme_bw()
p.com <- p.com + scale_fill_brewer("Family", palette = "Paired")
p.com <- p.com + rremove("x.text")
p.com
ggsave("./figures/Composition plots.pdf", height = 4, width = 6)









##Heatmaps:
# base plot
p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance)) 
# Change color
p.heat <- p.heat + scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw() 

# Make bacterial names italics
p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black', 
                                                    size = 10, 
                                                    face = 'italic')) 
# Make seperate samples based on main varaible
p.heat <- p.heat + facet_grid(~xlabel, 
                              scales = "free") + rremove("x.text") 

p.heat <- p.heat + ylab("Family")

#Clean the x-axis
p.heat <- p.heat + theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()) 

# Clean the facet label box
p.heat <- p.heat + theme(legend.key = element_blank(), 
                         strip.background = element_rect(colour="black", fill="white"))

print(p.heat)

ggsave("./figures/Heatmap.pdf", height = 4, width = 6)

##extra:
# we use count data at family level from the barplot for counts
ps_df <- microbiomeutilities::phy_to_ldf(ps1.com.fam, 
                                         transform.counts = "compositional")
colnames(ps_df)
# this data.frame can be used to customize several plots.  

# example boxplot at family level

p.box <- ggstripchart(ps_df, "BodySite", "Abundance", 
                      facet.by = "Family", color = "ReportedAntibioticUsage",
                      palette = "jco"
)

p.box + rremove("x.text")

#############################
###########################

##BETA DIVERSITY METRICS:
#Beta-diversity: Measures for differences between samples from different groups to identify if \
#there are differences in the overall community composition and structure.

# read non rarefied data
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")

# read rarefied data
ps0.rar.rds <- readRDS("./phyobjects/ps0.rar.rds")
# use print option to see the data saved as phyloseq object.

###Phylogenetic beta-diversity metrics
###1: Unweighted Unifrac

#based on presence/absence of different taxa and abundance is not important. 
#However, it is sensitive to the sequencing depth.
#f a sample is sequenced more than the others then it may have many OTUs \
#(most of them unique) consequently affecting the unifrac dissimilarity estimation.

# if we remove OTUs that are detected atleast 10 times in 5% of the samples
ps0.rar.filtered <- core(ps0.rar.rds, detection = 10, prevalence = 0.05)

summarize_phyloseq(ps0.rar.filtered)
summarize_phyloseq(ps0.rar.rds)
summarize_phyloseq(ps1)
# we reduce the sparsity considerably.

ordu.unwt.uni <- ordinate(ps0.rar.filtered, "PCoA", "unifrac", weighted=F)
# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

unwt.unifrac <- plot_ordination(ps0.rar.filtered, 
                                ordu.unwt.uni, color="BodySite") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("BodySite", palette = "Set2")
print(unwt.unifrac)
#NOTE: Try repeating the above ordination using non-filtered phyloseq object (ps0.rar.rds).

###2: Unweighted Unifrac
#Weighted Unifrac will consider the abundances of different taxa.
ps1.rel <- microbiome::transform(ps1, "compositional")

ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)

# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="BodySite") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("BodySite", palette = "Set2")
print(wt.unifrac)

print(wt.unifrac + stat_ellipse())

###Population-level Density landscapes
p <- plot_landscape(ps1.rel, 
                    "NMDS", 
                    "bray", 
                    col = "BodySite") +
  labs(title = paste("NMDS / Bray-Curtis"))   

p <- p + scale_color_brewer(palette = "Dark2")+ scale_fill_gradient(low = "#e0ecf4", high = "#6e016b") 
p

############################
#############################
###PERMANOVA
##A PERMANOVA lets you statistically determine if the centers (centroids) of the cluster of samples (in an ordination such as PCA) \
#for one groups (e.g treatment) differs from the center of the samples for another group (e.g control). 
library(vegan)
metadf <- data.frame(sample_data(ps1.rel))

unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)


permanova <- adonis(unifrac.dist ~ BodySite, data = metadf)
#permanova <- adonis(unifrac.dist ~ ReportedAntibioticUsage, data = metadf)
permanova


#Checking the homogeneity condition
#Type ?betadisper in R console for more information.
ps.disper <- betadisper(unifrac.dist, metadf$BodySite)
permutest(ps.disper, pairwise = TRUE)

#######################################
##CORE MICROBIOTA!
##Important for differential abundance testing below
# read non rarefied filtered data
ps1 <- readRDS("./phyobjects/ps.ng.tax.rds")

#Subset the data to keep only stool samples.
ps1.gut <- subset_samples(ps1, BodySite == "gut")

# convert to relative abundance  
ps1.gut.rel <- microbiome::transform(ps1.gut, "compositional")
print(ps1.gut.rel)

ps1.gut.rel2 <- prune_taxa(taxa_sums(ps1.gut.rel) > 0, ps1.gut.rel)
print(ps1.gut.rel2)

#Check for the core ASVs/OTUs

core.taxa.standard <- core_members(ps1.gut.rel2, detection = 0.001, prevalence = 50/100)

print(core.taxa.standard)

# Extract the taxonomy table
taxonomy <- as.data.frame(tax_table(ps1.gut.rel2))

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)

DT::datatable(core_taxa_id)

#Total core abundance in each sample (sum of abundances of the core members):
core.abundance <- sample_sums(core(ps1.gut.rel2, detection = 0.001, prevalence = 50/100))
DT::datatable(as.data.frame(core.abundance))

#####Core heatmaps
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(ps1.gut.rel2, 
                    plot.type = "heatmap", 
                    colours = gray,
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p.core)    

# Same with the viridis color palette
library(viridis)
print(p.core + scale_fill_viridis())

##Color change
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette

p.core <- plot_core(ps1.gut.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

print(p.core) 

#Use the format_to_besthit function from microbiomeutilities to get the best classification of the ASVs.
ps1.gut.rel2.f <- microbiomeutilities::format_to_besthit(ps1.gut.rel2)
p.core <- plot_core(ps1.gut.rel2.f, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

p.core + theme(axis.text.y = element_text(face="italic"))
print(p.core)

############################################
#########DIFFERENTIAL ABUNDANCE TESTING FOR UNIVARIATE DATA 
#compares the abundance of a selected bug between two conditions.
#We will use a different set of data, in the dietswap dataset
#NOTE: MAKE SURE DATA ARE NORMALIZED!!!!
theme_set(theme_bw(20))

#Load data
##The diet swap data set represents a study with African and African American \
#groups undergoing a two-week diet swap. Data repository: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
data(dietswap)
d <- dietswap

# Pick microbial abundances for a given taxonomic group
taxa <- "Dialister"

# Construct a data.frame with the selected
# taxonomic group (Dialister) and grouping
df <- data.frame(Abundance = abundances(d)[taxa,],
                 Group = meta(d)$nationality)

library(knitr)
kable(head(df))

#Visualize two groups:
p1 <- ggplot(df, aes(x = Group, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Absolute abundances", y = "Abundance (read count)")

# Let us add the log10(1+x) version:
df$Log10_Abundance <- log10(1 + df$Abundance)
p2 <- ggplot(df, aes(x = Group, y = Log10_Abundance)) +
  geom_boxplot() +
  labs(title = "Log10 abundances", y = "Abundance (log10(1+x) read count)")       


#devtools::install_github("thomasp85/patchwork")
library(patchwork) #combines two ggplots into the same graphics
p1 + p2
## Test whether abundance differences are statistically significant between the two groups
##t-test:

#Significance p-value with t-test:
t.test(Log10_Abundance ~ Group, data = df)
print(t.test(Log10_Abundance ~ Group, data = df)$p.value) #Abundance is not significantly different between the groups (p < 0.05 level)

# Assess whether the abundance data is Gaussian or log-normal within each group 
p <- ggplot(df, aes(fill = Group, x = Log10_Abundance)) +
  geom_density(alpha = 0.5)
print(p) #not Gaussian!!
#Wilcoxon test on log10 data:
print(wilcox.test(Log10_Abundance ~ Group, data = df)$p.value)

#Wilcoxon test on original data:
print(wilcox.test(Abundance ~ Group, data = df)$p.value)

###Fitting Linear Models:
#use Log10 abundances since this is closer to the Gaussian
res <- glm(Log10_Abundance ~ Group, data = df, family = "gaussian")
res
knitr::kable(summary(res)$coefficients, digits = 5)

#The intercept equals to the mean in the first group:
print(mean(subset(df, Group == "AAM")$Log10_Abundance))

#The group term equals to the difference between group means:
print(mean(subset(df, Group == "AFR")$Log10_Abundance) -
          mean(subset(df, Group == "AAM")$Log10_Abundance))

#Note that the linear model (default) significance equals to t-test assuming equal variances.
print(t.test(Log10_Abundance ~ Group, data = df, var.equal=TRUE)$p.value)

##Covariate testing
#Investigate how sex and bmi affect the results
#An important advantage of linear and generalized linear models, \
#compared to plain t-test is that they allow incorporating additional variables, \
#such as potential confounders (age, BMI, gender..):

# Add a covariate:
df$sex <- meta(d)$sex
df$bmi_group <- meta(d)$bmi_group
head(df)

# Fit the model:
res <- glm(Log10_Abundance ~ Group + sex + bmi_group, data = df, family = "gaussian")
knitr::kable(summary(res)$coefficients, digits = 5)
#We can even include interaction terms:
res <- glm(Log10_Abundance ~ Group * sex * bmi_group, data = df, family = "gaussian")
knitr::kable(summary(res)$coefficients, digits = 5)


###########ADVANCED MODELS FOR DIFFERENTIAL ABUNDANCE:
##GLMs are the basis for advanced testing of differential abundance in sequencing data. \
#This is necessary, as the sequencing data sets deviate from symmetric, continuous, Gaussian assumptions in many ways.

#Sequencing data consists of discrete counts:
print(abundances(d)[1:5,1:3])

#The data is sparse:
hist(log10(1 + abundances(d)), 100)

#Rarity
#Long tails of rare taxa:
library(reshape2)
medians <- apply(abundances(d),1,median)/1e3
A <- melt(abundances(d))
A$Var1 <- factor(A$Var1, levels = rev(names(sort(medians))))
p <- ggplot(A, aes(x = Var1, y = value)) +
  geom_boxplot() +
  labs(y = "Abundance (reads)", x = "Taxonomic Group") +
  scale_y_log10()

print(p)

#10.1.4 Overdispersion
#Variance exceeds the mean:
  
means <- apply(abundances(d),1,mean)
variances <- apply(abundances(d),1,var)

# Calculate mean and variance over samples for each taxon
library(reshape2)
library(dplyr)
df <- melt(abundances(d))
names(df) <- c("Taxon", "Sample", "Reads")
df <- df %>% group_by(Taxon) %>%
  summarise(mean = mean(Reads),
            variance = var(Reads))

# Illustrate overdispersion
library(scales)
p <- ggplot(df, aes(x = mean, y = variance)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  labs(title = "Overdispersion (variance > mean)")
print(p)

##########
###(GLM) allows a richer family of probability distributions to describe the data. \
#Intuitively speaking, GLMs allow the modeling of nonlinear, nonsymmetric, and nongaussian associations.

#We fit the abundance (read counts) assuming that the data is Poisson distributed, \
#and the logarithm of its mean, or expectation, is obtained with a linear model.
# Load again the example data
d <- dietswap
df <- data.frame(Abundance = abundances(d)[taxa,],
                 Group = meta(d)$nationality)
head(df)
res <- glm(Abundance ~ 1, data = df, family = "poisson")
knitr::kable(summary(res)$coefficients, digits = 5)
#Note the link between mean and estimated coefficient 
mean(df$Abundance)
exp(coef(res))

####DESeq2: differential abundance testing for sequencing data
# Start by converting phyloseq object to deseq2 format
library(DESeq2)
d <- dietswap # Phyloseq data
ds2 <- phyloseq_to_deseq2(d, ~ group + nationality)

# Run DESeq2 analysis (all taxa at once!)
dds <- DESeq(ds2)

# Investigate results
deseq.results <- as.data.frame(results(dds))
deseq.results$taxon <- rownames(results(dds))

# Sort (arrange) by pvalue and effect size
library(knitr)
deseq.results <- deseq.results %>%
  arrange(pvalue, log2FoldChange)

# Print the result table
# Let us only show significant hits
knitr::kable(deseq.results %>%
               filter(pvalue < 0.05 & log2FoldChange > 1.5),
             digits = 5)

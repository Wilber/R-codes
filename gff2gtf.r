####Convert gff to gtf, only ONCE!
#install.packages("rtracklayer") #install once!
library("rtracklayer")
gff <- import("../reference/ITAG3.2_gene_models.gff")
export(gff,"../reference/ITAG3.2_gene_models.gtf","gtf")

#NOTE: ADD GTF.attrType="ID" argument in Rsubread::featureCounts function, then you don't need the bash command bellow!!!

##the following code uses bash commands:
#then clean the output gff with:
sed -i -e 's/ID=gene/gene_id/' ITAG3.2_gene_models.gtf

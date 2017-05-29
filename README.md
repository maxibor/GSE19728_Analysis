

# GSE19728_Analysis

#### Author : Maxime Borry

Loading the necessary packages for this analysis
``` r
library(affy)
library(hgu133plus2.db)
library(edgeR)
library(gage)
library(coGSEA)
```

Loading the CEL FILES and Normalizing them usin `rma`
``` r
setwd("~/GitLab/GSE19728/data/")
celfiles = ReadAffy()
celfiles = rma(celfiles)
```

Getting the expression data, and the Probe Names
``` r

intensity = exprs(celfiles)
intensity = cbind(rownames(intensity), intensity)
colnames(intensity)[1] = "PROBEID"
intensity = as.data.frame(intensity)
intensity$PROBEID = as.character(intensity$PROBEID)
```

Getting the annotations of Probes IDS to ENTREZ accession numbers
``` r
annots = select(hgu133plus2.db, intensity$PROBEID, "ENTREZID", "PROBEID")
```

Merge expression matrix and annotation
``` r
res = merge(intensity, annots, by= "PROBEID")

```


Getting rid of PROBE ids and type casting
``` r
resmin = res[,2:ncol(res)]
cname = colnames(resmin)
resmin = apply(resmin, 2, as.numeric)
colnames(resmin)= cname
resmin = as.data.frame(resmin)
```

Aggregating the PROBES matching the same ENTREZ accession number by averaging them and applying a log transformation
``` r
result = aggregate(. ~ ENTREZID, resmin, mean)
result$ENTREZID = levels(as.factor(as.character(resmin$ENTREZID)))
rownames(result) = result$ENTREZID
result = result[,-1]
result = log(result)
```

Visualizing in boxplot to check normalization
``` r
boxplot(result, las = 2)
```

Changing column names to remove the `.CEL` filename extension
``` r
colnames(result) = gsub(".CEL","", colnames(result))
```

``` r
colnames(result)

```

Selecting only the samples we are interested in
``` r
result2 = cbind(result$GSM492649_Astrocytomas_N, result$GSM525014, result$GSM525015, result$GSM525016, result$`GSM492662_Astrocytomas_T4-1`, result$`GSM492663_Astrocytomas_T4-2` , result$`GSM492664_Astrocytomas_T4-3`, result$`GSM492665_Astrocytomas_T4-4`, result$`GSM492666_Astrocytomas_T4-5`)
colnames(result2) = c("GSM492649_Astrocytomas_N", "GSM525014", "GSM525015", "GSM525016","GSM492662_Astrocytomas_T4-1", "GSM492663_Astrocytomas_T4-2", "GSM492664_Astrocytomas_T4-3", "GSM492665_Astrocytomas_T4-4", "GSM492666_Astrocytomas_T4-5" )
rownames(result2) = rownames(result)
```

Preparing the design matrix
``` r
Normal = c(rep(1,4),rep(0,5))
Tumor = c(rep(0,4),rep(1,5))
design = cbind(Normal, Tumor)
rownames(design) = colnames(result2)
```

Preparing the contrast matrix
``` r
contr.matrix = makeContrasts(NormalVSTumor = Normal - Tumor, levels = design)
```


Preparing expression list object
``` r
temp = new("EList")
temp$design = design
temp$E = as.matrix(result2)
rownames(temp$E) = as.numeric(rownames(temp$E))
temp$genes$ENTREZ = rownames(result2)
temp$common.dispersion = estimateDisp(temp$E, design = temp$design)$common.dispersion
temp$samples = colnames(result2)
```

Preparing gene set collection
``` r
gs = gage::kegg.gsets(species = "hsa", id.type = "entrez")
geneset = gs$kg.sets
```

Function for simplifying gene sets names
``` r
nameshorter = function(names){
  namemod = c()
  for (i in seq(1,length(names))){
    namemod[i] = paste(strsplit(names[i], " ")[[1]][-1], sep = "", collapse = " ")
    namemod[i] = gsub("/","", names[i])
    namemod[i] = gsub(" ","_", names[i])
  }
  return(namemod)
}
```

Simplifying gene sets names
``` r
names(geneset) = nameshorter(names(geneset))
names(geneset) = gsub("/","_",names(geneset))
```

Saving necessary objects to RDS files
``` r
saveRDS(contr.matrix, "~/GitLab/GSE19728/contrast.rds")
saveRDS(temp, "~/GitLab/GSE19728/elist.rds")
saveRDS(geneset, "~/GitLab/GSE19728/geneset.rds")
```


Reading necessary objects (generated above) from RDS files

``` r
elist = readRDS("~/GitLab/GSE19728/elist.rds")
contrast = readRDS("~/GitLab/GSE19728/contrast.rds")
geneset = readRDS("~/GitLab/GSE19728/geneset.rds")
```



Running coGSEA analysis
``` r
coGSEA(ElistObject = elist, contrastMatrix = contrast, ENTREZGenesIds = elist$genes$ENTREZ, geneSetCollection = geneset,specie = "Homo sapiens", directoryPath = "~/GitLab/GSE19728/results", alpha = 0.05, pvalAdjMethod = "BH", pvalCombMethod = "sumlog",min.intersection.size = 1, GSEA.Methods = c("camera", "gage","globaltest", "gsva", "ssgsea", "zscore", "ora", "padog", "roast","safe"), num.workers = 4, shinyMode = FALSE)

```

``` r
result = read.csv("~/GitLab/GSE19728/results/result_NormalVSTumor.csv", header = TRUE, row.names = 1)
```

``` r
str(cor(result$combined_p.value, result$Avg_Rank, method = "spearman"))
```


``` r
plot(-log10(result$combined_p.value), result$Avg_Rank)

```




``` r
par(mar=c(10,5,3,5))
res_pval_adj = result[,grep("_adj_p.value",colnames(result))]
colnames(res_pval_adj) = gsub("_adj_p.value","", colnames(res_pval_adj))
# png(file = "/home/maxime/GitLab/report/resources/GSSE_number.png", width = 1024)
dtable = apply(res_pval_adj,2, function(x) length(which(x < 0.05)))
bplt = barplot(dtable,las = 2, ylab = "gene sets", ylim = c(0,300),cex.names = 2)
text(y= dtable, x=bplt , pos = 3, cex=2, col = "red", labels = dtable)
# dev.off()
```

``` r
length(unique(unlist(apply(res_pval_adj,2, function(x) which(x < 0.05)))))
```

``` r
colnames(res_pval_adj)[which(res_pval_adj[grep("hsa05214", rownames(res_pval_adj)),] <= 0.05)]
colnames(res_pval_adj[grep("hsa05214", rownames(res_pval_adj)),])[which(res_pval_adj[grep("hsa05214", rownames(res_pval_adj)),]<0.05)]
library(metap)
sumlog(res_pval_adj[grep("hsa05214", rownames(res_pval_adj)),])
```

``` r
length(which(result$combined_p.value < 0.05))
```

``` r
colnames(result)
ranks_table = result[,grep("_Rank",colnames(result))[1:10]]
first5GS = apply(ranks_table,2,function(x) sort(x, index.return = TRUE)$ix)[1:5,]
View(as.data.frame(apply(first5GS, 2, function(x) rownames(result)[x]))
```


``` r
rownames(result)[sort(result$Avg_Rank, index.return = TRUE)$ix]
```

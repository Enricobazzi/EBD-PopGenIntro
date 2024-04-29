# Identificación de señales de selección positiva y análisis de enriquecimiento funcional

En esta práctica, vamos a utilizar datos WGS simulados para 40 individuos de dos poblaciones (20 cada una), a partir de los cuales estimaremos la presencia o no de señales de barrido selectivo en el genoma a través de diferentes aproximaciones que hemos visto en la parte teórica del curso. Los datos, por tanto, consisten en dos VCFs, uno por cada población, que contienen la información de los genotipos en fase de los individuos en cada una de las variantes genéticas detectadas en el genoma.

La estructura de la práctica será la siguiente: 1. Estimación de pi, thetha y Tajima's D 2. Análisis de iHS 3. Enriquecimiento funcional

## iHS

Como ya hemos visto, iHS es un estadístico que mide la asimetría de la extensión de haplotipos en desequilibrio de ligamiento en dos direcciones opuestas (ancestral-derivado). Por lo tanto, cuanto mayor sea la diferencia entre las dos direcciones, mayor será el valor de iHS.

En este caso, vamos a utilizar el paquete rehh para calcular iHS en las dos poblaciones y ver si hay señales de selección positiva en alguna de ellas.

Primero cargamos las librerias necesarias:

```{r message=FALSE, warning=FALSE}
#Dependencies 
library(vcfR) 
library(rehh)
library(ggplot2)
library(tidyverse)
library(topGO)
library(biomaRt)
#BiocManager::install("biomaRt")
```

Cargamos en R los datos para cada población

```{r}
#Read the vcf file 
vcf_pop1 <- data2haplohh(hap_file= "data3_selection_pop1.vcf", allele_coding="01", vcf_reader= "vcfR")
vcf_pop2 <- data2haplohh(hap_file= "data3_selection_pop2.vcf", allele_coding="01", vcf_reader= "vcfR")
```

Vamos a estudiar la POBLACIÓN 1:

```{r}
#Scan the genome: this function calculates iHH (integrated EHH, area under EHH curve) and iES (integrated EHH per Site, which forms the basis for cross-population comparisons)for each SNP. 
scan_pop1 <- scan_hh(vcf_pop1)

#Calculate iHS values 
ihs_pop1 <- ihh2ihs(scan_pop1, freqbin=0.05 , min_maf=0)

#extract candidate SNPs 
ihs_val_pop1 <- ihs_pop1$ihs %>% na.omit() 
candidates_pop1 <- ihs_val_pop1 %>% filter(IHS > 4 | IHS < -4) 

print(candidates_pop1)

#plot iHS results 
pop1 <- ggplot(data = ihs_val_pop1, aes(x = POSITION, y = LOGPVALUE)) +
          geom_point(alpha = 0.5, color = ifelse(ihs_val_pop1$IHS > 4 | ihs_val_pop1$IHS < -4, "red", "grey")) +
          xlab("Position (Mb)") + 
          ylab("log p-value")
          theme_minimal()

print(pop1)
```

Vamos a estudiar la POBLACIÓN 2:

```{r}
#Scan the genome: this function calculates iHH (integrated EHH, area under EHH curve) and iES (integrated EHH per Site, which forms the basis for cross-population comparisons)for each SNP.
scan_pop2 <- scan_hh(vcf_pop2)

#Calculate iHS values
ihs_pop2 <- ihh2ihs(scan_pop2, freqbin=0.05 , min_maf=0)

#extract candidate SNPs
ihs_val_pop2 <- ihs_pop2$ihs %>% na.omit()
candidates_pop2 <- ihs_val_pop2 %>% filter(IHS > 4 | IHS < -4)

print(candidates_pop2)

#plot iHS results
pop2 <- ggplot(data = ihs_val_pop2, aes(x = POSITION, y = LOGPVALUE)) + 
        geom_point(alpha = 0.5, color = ifelse(ihs_val_pop2$IHS > 4 | ihs_val_pop2$IHS < -4, "red", "grey")) + 
        xlab("Position (Mb)") + 
        ylab("log p-value")
        theme_minimal()

print(pop2)
```

## Enriquecimiento funcional

```{r}
#####Get annotation from ensembl.org#####

#read from ensembl.org every felcat ensembl ID:
ensembl <- useMart("ensembl", dataset = "fcatus_gene_ensembl")

#extract the GOterms for every ensembl_id
ensembl_to_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"), mart = ensembl)

#make a list with every GO terms per ENSEMBL_ ID. 
go_list <- split(ensembl_to_go$go_id, ensembl_to_go$external_gene_name)

#Comment: By this we get a total of 29550 ensembl_id corresponding to coding genes (19588) + non-coding genes (9468) + pseudogenes (494) according to https://www.ensembl.org/Felis_catus/Info/Annotation
```

```{r eval=FALSE, include=FALSE}

###But I also want to extract some desired GO IDs (for example that one related with inmune system)
# Define the desired GO IDs
desired_GO_IDs <- c("GO:0016045", "GO:0051607", "GO:0045058")

# Initialize a vector to store the selected gene names
selected_genes <- c()

# Loop through each desired GO ID
for (desired_GO_ID in desired_GO_IDs) {
  # Get gene names associated with the current desired GO ID
  genes_with_desired_GO <- names(go_list)[sapply(go_list, function(x) desired_GO_ID %in% x)]
  
  # Append the gene names to the selected_genes vector
  selected_genes <- c(selected_genes, genes_with_desired_GO)
}

# Remove duplicates in case a gene is associated with multiple desired GO IDs
selected_genes <- unique(selected_genes)

# Set the seed for reproducibility
set.seed(123)

# Define the number of gene names to select randomly
num_genes_to_select <- 10  # Adjust as needed

# Randomly sample gene names from the selected gene names
randomly_selected_genes <- sample(selected_genes, num_genes_to_select, replace = TRUE)

###From the total list of genes (every gene in go_list), I want to extract a random subset
randomly_nonselected_genes <- sample(names(go_list), 100)

###now sum the selected and not selected
genes <- c(randomly_selected_genes, randomly_nonselected_genes)

write.csv(genes, "candidate_genes.csv")


```

```{r}
#read candidate genes csv
genes<- read.csv("candidate_genes.csv")
genes <- genes$x

#cross the felcat annotation (ensembl_id with its GOterms) with my set of genes, to get a logical factor of true (if the gene is in the set) and false (if its not)
allgenes = factor(as.integer(names(go_list) %in% genes))
names(allgenes) <- names(go_list) #ensure names of allgenes are names of go_list

#####Creating the topGOdata object#####
godata <- new("topGOdata", ontology = "BP", allGenes = allgenes, 
              annotationFun = annFUN.gene2GO, gene2GO = go_list, nodeSize=10)


#####run the overrepresentation test#####
over_test <- runTest(godata, statistic = "fisher")


result_table <- GenTable(godata, Fisher=over_test, topNodes=over_test@geneData[2]) %>% 
        as_tibble() %>% 
        mutate(p.adj = round(p.adjust(as.numeric(gsub("<", "", Fisher)), method="BH"), 4), ontology= "BP") %>% 
        filter(p.adj<0.05) 

#write.table(results_table, file ="enrichment.csv", sep =  ";")

print(result_table)

```

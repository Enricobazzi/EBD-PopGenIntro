# Analisis Preliminares de Estructura, Diversidad y Historia Demografica

Con el objetivo de estudiar la historia de una especie, hemos muestreado cuatro regiones de una población de ... (A, B, C, D), extrayendo el ADN de 10 individuos de cada población.

**IMAGEN RESUMEN MUESTREO**

## Exploramos el VCF

Hemos alineado nuestras lecturas al genoma de referencia, y llamado las variantes de los 40 individuos en un [VCF](data/data1.vcf). Por semplicidad, en esta practica solo miraremos las variantes que hemos encontrado en el cromosoma 1 que están guardadas en el file `data1.vcf` en la carpeta `data` de esta practica.

```
# vamos a la carpeta de la practica
cd practica2

# declaramos el vcf como una variable
VCF=data/data1.vcf

# echamos un ojo por encima al vcf
less -S ${VCF}
```

Vemos `##contig=<ID=1,length=10000000>` que nos indica que el tamaño del chromosoma 1 es 10'000'000 pares de bases (10Mb), que es relativamente pequeño. Esto es para acelerar los calculos que haremos en esta practica.

Para sacar la lista de individuos podemos correr el siguiente comando:

```
grep -m1 "#CHROM" $VCF | cut -f10- | tr "\t" " "
```

Y para saber cuantos SNPs tenemos en nuestro VCF:

``` 
grep -v "#" $VCF | wc -l
```

## Analisis de componentes principales (PCA)

Si queremos explorar como están relacionados nuestros individuos, podemos empezar por una PCA basada en sus genotipos.

```
mkdir -m 777 pca
```

El programa [Plink](https://www.cog-genomics.org/plink/) nos ofrece una opción muy facil para realizar una PCA a partir de nuestro VCF.

El primer paso para correr una PCA es la selección de variables no correlacionadas entre si, para que los componente principales no sean distorsionados por relaciones de covarianza entre variables. Como nuestras variables serán los genotipos (SNPs), empezamos llamando la función [indep-pairwise](https://www.cog-genomics.org/plink/1.9/ld#indep) de plink. Le diremos de mirar ventanas de 50 SNPs, moviendose de 1 SNP en cada iteración y eliminando cualquier SNP con coeficiente de correlación cuadrado (r^2) > 0.2:

```
plink --vcf $VCF \
    --indep-pairwise 50 1 0.2 \
    --out data/data1.plink
```

Esta función genera varios files `ls data/`:

- data1.plink.log : contiene la información sobre el comando de plink que acabamos de correr
- data/data1.plink.nosex : contiene la lista de individuos en nuestro análisis por los que no hay información de sexo (todos)
- data1.plink.prune.in : lista de SNPs que pasan los criterios de filtrado (no están correlacionados entre si)
- data1.plink.prune.out : lista de SNPs que NO pasan los criterios de filtrado (están correlacionados entre si)

Eliminando estas ultimas variantes con la función `--exclude`, podemos correr un análisis de componentes principales con la función [--pca](https://www.cog-genomics.org/plink/1.9/strat#pca) de plink:

```
plink --vcf $VCF \
    --exclude data/data1.plink.prune.out \
    --pca 'header' \
    --out pca/data1.plink.pca
```

Esta función genera dos nuevos tipos de file:

- data1.plink.pca.eigenval : contiene información sobre la cantidad de varianza explicada por cada componente principal
- data1.plink.pca.eigenvec : contiene los valores de cada muestra para cada componente principal

Para visualizar estos resultados hemos preparados unos scripts de R:

```{r}
library(ggplot2)

# read eigenvalues
eigenvals <- read.table("pca/data1.plink.pca.eigenval",
                        col.names = c("eigenval"))

# read eigenvectors
eigenvecs <- read.table("pca/data1.plink.pca.eigenvec",
                        header = TRUE)[,-1]

# convert eigenvalues to percentage
percents <- eigenvals$eigenval/sum(eigenvals$eigenval)*100

# scree-plot of eigenvalues
scree <- ggplot() +
  geom_line(aes(x = as.numeric(row.names(eigenvals)), y = percents)) +
  geom_point(aes(x = as.numeric(row.names(eigenvals)), y = percents),
    fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Principal Component")) +
  ylab(paste0("Percentage of variance Explained")) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw()

ggsave(filename = "pca/data1.scree_plot.pdf", plot = scree)

# add the locations
loc <- c(rep("A", 10), rep("B", 10), rep("C", 10), rep("D", 10))

# plot!
pc1pc2 <- ggplot() +
  geom_point(data = eigenvecs,
    aes(x = PC1, y = PC2, fill = loc),
    shape = 21, size = 2.5) +
  xlab(paste0("PC1 - ", round(percents[1], 2), "%")) +
  ylab(paste0("PC2 - ", round(percents[2], 2), "%")) +
  theme_bw()
ggsave(filename = "pca/data1.pc1pc2_plot.pdf", plot = pc1pc2)

pc2pc3 <- ggplot() +
  geom_point(
    data = eigenvecs,
    aes(x = PC2, y = PC3, fill = loc),
    shape = 21, size = 2.5
    ) +
  xlab(paste0("PC2 - ", round(percents[2], 2), "%")) +
  ylab(paste0("PC3 - ", round(percents[3], 2), "%")) +
  theme_bw()
ggsave(filename = "pca/data1.pc2pc3_plot.pdf", plot = pc2pc3)

```

## Admixture para identificar unidades poblacionales

Vamos a correr el software [ADMIXTURE](https://dalexander.github.io/admixture/) para identificar poblaciones en nuestros datos de forma supervisada. Esto significa que tenemos que decidir nosotros el numero de poblaciones en los que nos tiene que dividir nuestros individuos.

Admixture requiere los datos en formato [BED de plink](https://www.cog-genomics.org/plink/1.9/data#make_bed). Para esto usaremos la función `--make-bed` de plink, seleccionando de nuevo SNPs no correlacionados entre si:

```
plink --vcf $VCF \
    --exclude data/data1.plink.prune.out \
    --make-bed \
    --out data/data1.plink
```

Ahora podemos correr ADMIXTURE con diferentes valores de K. Utilizamos [--cv]()

```
# creamos una carpeta para guardar los resultados de admixture
mkdir -m 777 admixture 
cd admixture

# corremos admixture de K=1 a K=8
for k in {1..8}; do
    echo "running admixture for K=${k}"
    admixture --cv ../data/data1.plink.bed ${k} \
        > log${k}.out
done

# resumen del cross validation (cv) en un file
paste <(grep -h CV log*.out | cut -d'=' -f2 | cut -d')' -f1) \
    <(grep -h CV log*.out | cut -d' ' -f4) \
    > cv_table.txt

cd ..
```

Plotear los valores de error del cross validation nos ayuda a tomar decisiones sobre el "correcto" numero de K

```{r}
library(ggplot2)

# read cv file
cv_table <- read.table("admixture/cv_table.txt", col.names = c("K", "CVerr"))

# plot cv values
cvplot <- ggplot() +
  geom_line(data = cv_table,
    aes(x = K, y = CVerr)) +
  geom_point(data = cv_table,
    aes(x = K, y = CVerr),
    fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Populations")) +
  ylab(paste0("Cross Validation Error")) +
  scale_x_continuous(n.breaks = 8) +
  theme_bw()
ggsave(filename = "admixture/data1.cvplot.pdf", plot = cvplot)
```

Diferentes criterios se pueden usar para determinar el mejor numero de K. De los más comunes:
- valor de K donde el error de cv es mínimo
- valor de K donde mejora más el error de cv en un paso

*La información sobre la biología de la especie es importante!*

Ploteamos como nuestras muestras se dividirían con diferentes valores de K:

```{r}
library(ggplot2)
library(tidyr)

# sample names
samples <- read.table("data/data1.plink.fam", sep = " ",
                      col.names = c("FID", "IID", "father",
                                    "mother", "sex", "pheno"))[2]

# loop through k 2 to 8
for (k in 2:8){
  # q values
  qvals <- read.table(paste0("admixture/data1.plink.", k, ".Q"),
                    col.names = paste0("Q",seq(1:k)))
  # join dataframe
  admix <- cbind(samples, qvals) %>% gather(Q, value, -IID)
  # plot
  admix_plot <- ggplot(admix, aes(x = IID, y = value, fill = factor(Q))) +
    geom_bar(stat = "identity", position = "stack") +
    xlab("Sample") + ylab("Ancestry") +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(filename = paste0("admixture/data1.plot.K", k,".pdf"),
         plot = admix_plot, width = 20, height = 5, units = "cm")
}
```

## Estadisticos de Diversidad Poblacionales

Ahora que hemos identificado las poblaciones podemos proceder a calcular algunos estadisticos de diversidad para ellas.

```
mkdir -m 777 stats
```

Primero vamos a dividir el VCF por población y luego usamos las funciones `--window-pi` y `--TajimaD` indicando el tamaño del cromosoma entero para mirar los valores de estos estadisticos dentro de cada una.

```
for pop in A B C D; do
    grep -m1 "#CHROM" $VCF | cut -f10- | tr "\t" "\n" | grep "${pop}" > data/population_${pop}.txt
    
    vcftools --vcf $VCF --maf 0.001 --keep data/population_${pop}.txt \
        --recode --out data/data1.population_${pop}

    vcftools --vcf data/data1.population_${pop}.recode.vcf \
        --window-pi 10000000 --out stats/population_${pop}_pi_20Mb
    
    vcftools --vcf data/data1.population_${pop}.recode.vcf \
        --TajimaD 10000000 --out stats/population_${pop}_tajimaD_20Mb
done
```

Podemos calcular la heterocigosidad individual a partir de los VCFs de cada población en R directamente:

```{r}
library(ggplot2)

pops <- c("A", "B", "C", "D")

het_table <- data.frame()

for (pop in pops){
  
  vcf_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                 "INFO", "FORMAT", c(paste0(pop,seq(1:10))))
   
  vcf <- read.table(paste0("data/data1.population_", pop, ".recode.vcf"),
              col.names = vcf_cols)
  
  samples <- colnames(vcf)[10:19]
  
  pop_data <- data.frame()
  
  for (sample in samples){
    het <- sum(vcf[, sample] == "1|0" | vcf[, sample] == "0|1") / nrow(vcf)
    sample_row <- data.frame(population = pop,
                             individual = sample, heterozygosity = het)
    pop_data <- rbind(pop_data, sample_row)
  }
  
  het_table <- rbind(het_table, pop_data)
}

hetplot <- ggplot() +
  geom_point(data = het_table, aes(individual, heterozygosity, fill = population),
             shape = 21, size = 2.5) +
  xlab("Individual") + ylab("Heterozygosity") +
  theme_bw()
ggsave(filename = "stats/data1.hetplot.pdf")
```

Utilizando los output de los comandos de vcftoos `--window-pi` y `--TajimaD` podemos visualizar estos valores en R:

```{r}
library(ggplot2)
library(tidyr)

pops <- c("A", "B", "C", "D")
diversity_table <- data.frame()

for (pop in pops){
  pi_table <- read.table(paste0("stats/population_", pop, "_pi_20Mb.windowed.pi"),
                         header = TRUE)
  pi <- pi_table$PI
  
  Td_table <- read.table(paste0("stats/population_", pop, "_tajimaD_20Mb.Tajima.D"),
                         header = TRUE)
  Td <- Td_table$TajimaD
  
  pop_row <- data.frame(population=pop,
                        nucleotide_diversity=pi, tajima_d=Td)
  diversity_table <- rbind (diversity_table, pop_row)
}

piplot <- ggplot(diversity_table, 
       aes(x = population, y = nucleotide_diversity)) +
    geom_col() +
    xlab("Population") + ylab("Nucleotide diversity") +
    theme_bw()
ggsave(filename = "stats/data1.piplot.pdf", plot = piplot)

tdplot <- ggplot(diversity_table, 
       aes(x = population, y = tajima_d)) +
    geom_col() +
    xlab("Population") + ylab("Tajima's D") +
    theme_bw()
ggsave(filename = "stats/data1.tdplot.pdf", plot = tdplot)
```

## Análisis de tamaño efectivo

```
chmod +x GONE-Linux/PROGRAMMES/*

VCF=data/data1.big.vcf

vcftools --vcf $VCF --maf 0.001 --keep data/population_A.txt \
    --recode --out data/data1.big.population_A

VCF=data/data1.big.population_A.recode.vcf

plink --vcf $VCF \
    --recode \
    --out GONE-Linux/data1.big.population_A

cd GONE-Linux

bash script_GONE.sh data1.big.population_A


################################################

VCF=data/data1.big.vcf

vcftools --vcf $VCF --maf 0.001 --keep data/population_B.txt \
    --recode --out data/data1.big.population_B

VCF=data/data1.big.population_B.recode.vcf

plink --vcf $VCF \
    --recode \
    --out GONE-Linux/data1.big.population_B

cd GONE-Linux

bash script_GONE.sh data1.big.population_B

################################################

VCF=data/data1.big.vcf

vcftools --vcf $VCF --maf 0.001 --keep data/population_C.txt \
    --recode --out data/data1.big.population_C

VCF=data/data1.big.population_C.recode.vcf

plink --vcf $VCF \
    --recode \
    --out GONE-Linux/data1.big.population_C

cd GONE-Linux

bash script_GONE.sh data1.big.population_C

```

```{r}
library(ggplot2)
library(tidyverse)

popA <- read.table("GONE-Linux/Output_Ne_data1.big.population_A",
                   skip = 1, header = TRUE)

ggplot() +
  geom_step(data = popA %>% filter(Generation <= 150),
            aes(x = Generation, y = Geometric_mean))
```

```{r}

popB <- read.table("GONE-Linux/Output_Ne_data1.big.population_B",
                   skip = 1, header = TRUE)

ggplot() +
  geom_step(data = popB %>% filter(Generation <= 150),
            aes(x = Generation, y = Geometric_mean))
```
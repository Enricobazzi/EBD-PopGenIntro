# Obtención de un dataset de Variantes a partir de datos de Secuenciación - WGS

## Ficheros Fastq

Tenemos unos [fastq](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211) que queremos analizar. Están guardados en la carpeta `data` del directorio de esta practica:

```
cd practica1
ls data/
```

Vamos a mirar uno de ellos. *Tendremos que tener en cuenta que están comprimidos en formato gzip (extensión .gz)*

```
# miramos el file
zcat data/LL240.R1.fq.gz | less -S
```

Cuantas secuencias tenemos?

```
# contamos las secuencias contando el numero de lineas y dividiendo por 4
zcat data/LL240.R2.fq.gz | wc -l
# o directamente
zcat data/LL240.R2.fq.gz | awk 'NR%4==1 {count++} END {print count}'
```

## Control de calidad de las lecturas: Fastqc

Vamos a usar Fastqc para obtener un informe de calidad de nuestras lecturas:

```
# creamos una carpeta para guardar los resultados
mkdir -m 777 fastqc

# corremos fastqc para cada file de lecturas
fastqc data/LL240.R1.fq.gz -o fastqc
fastqc data/LL240.R2.fq.gz -o fastqc
fastqc data/LP220.R1.fq.gz -o fastqc
fastqc data/LP220.R2.fq.gz -o fastqc
```

## Fastp

Vamos a correr fastp sobre nuestras lecturas con unos objetivos comunes en los pasos de pre-procesado:
- eliminación de adaptadores 
- corrección de errores
- eliminación de lecturas demasiado cortas
- eliminación de cadenas de poly-g

```
mkdir -m 777 fastp

fastp --help

fastp \
    -i data/LL240.R1.fq.gz -I data/LL240.R2.fq.gz \
    -o data/LL240.R1.fastp.fq.gz -O data/LL240.R2.fastp.fq.gz \
    -h fastp/LL240.fastp.html -j fastp/LL240.fastp.json \
    --unpaired1 data/LL240.R1.unpaired.fq.gz --unpaired2 data/LL240.R2.unpaired.fq.gz \
    --failed_out data/LL240.R1.failed.fq.gz \
    --trim_poly_g \
    --length_required 30 \
    --correction \
    --detect_adapter_for_pe

fastp \
    -i data/LP220.R1.fq.gz -I data/LP220.R2.fq.gz \
    -o data/LP220.R1.fastp.fq.gz -O data/LP220.R2.fastp.fq.gz \
    -h fastp/LP220.fastp.html -j fastp/LP220.fastp.json \
    --unpaired1 data/LP220.R1.unpaired.fq.gz --unpaired2 data/LP220.R2.unpaired.fq.gz \
    --failed_out data/LP220.R1.failed.fq.gz \
    --trim_poly_g \
    --length_required 30 \
    --correction \
    --detect_adapter_for_pe
```

```
fastqc data/LP220.R1.fastp.fq.gz data/LP220.R2.fastp.fq.gz -o fastqc
fastqc data/LL240.R1.fastp.fq.gz data/LL240.R2.fastp.fq.gz -o fastqc
```

## Alineamiento

Creamos una carpeta donde guardar nuestros alineamientos:

```
mkdir -m 777 bams
```

Vamos a echarle un vistazo a nuestro genoma de referencia:

```
less data/reference_sequence.fa
grep ">" data/reference_sequence.fa
```

Como los genomas de referencia son muy grandes, los algoritmos de alineamiento necesitan indices y diccionarios para ser mas eficientes. Vamos a generarlos facilmente con los siguientes comandos:

```
# indice para BWA
bwa index data/reference_sequence.fa
# indice general
samtools faidx data/reference_sequence.fa
# diccionario para gatk
samtools dict data/reference_sequence.fa -o data/reference_sequence.dict
```

- .amb es un file de texto, registra las N y otros caracteres que no sean ATGC
- .ann es un file de texto, registra las secuencias de referencia (contig/scaffolds/cromosomas), sus nombres, longitud etc.
- .bwt es un file binario, contiene la secuencia transformada por el algoritmo Burrows-Wheeler
- .pac es un file binario, contiene la secuencia enpaquetada
- .sa es un file binario, indice de la collección de sufijos (de la matriz de la secuencia transformada)
- .fai es un file de texto, contiene  la información sobre el genoma de referencia en 5 columnas (nombre, longitud, offset, bases por linea y bytes por linea)
- .dict es un file de texto, parecido a un cabecero de un SAM, contiene información sobre el genóma de referencia

Con los indices listos y nuestras secuencias limpias podemos proceder al alineamiento con BWA-MEM

```
bwa mem \
    data/reference_sequence.fa \
    data/LP220.R1.fastp.fq.gz \
    data/LP220.R2.fastp.fq.gz |
    samtools view -hbS - -o bams/LP220.bam

bwa mem \
    data/reference_sequence.fa \
    data/LL240.R1.fastp.fq.gz \
    data/LL240.R2.fastp.fq.gz |
    samtools view -hbS - -o bams/LL240.bam
```

Para observar nuestro BAM podemos usar una herramienta llamada samtools view, la misma que lo ha comprimido lo puede descomprimir:

```
samtools view -h bams/LP220.bam | less -S
samtools view -h bams/LL240.bam | less -S
```

Muchos de los programas que usaremos para seguir tratando nuestros datos asumen dos cosas principales de nuestro BAM:
1. que las lecturas aparezcan en el orden de su posición en el genoma de referencia
2. que las lecturas estén etiquetadas para poder identificar su grupo de pertenencia (e.g. una carrera de secuenciación en concreto)

Para ordenar nuestras lecturas alineadas en el mismo orden de las posiciones del genoma de referencia podemos usar `samtools sort`:

```
samtools sort bams/LP220.bam -o bams/LP220.sorted.bam
samtools sort bams/LL240.bam -o bams/LL240.sorted.bam

samtools view -h bams/LP220.sorted.bam | less -S
samtools view -h bams/LL240.sorted.bam | less -S
```

Para asignar todas las lecturas dentro de un BAM a un grupo especifico podemos usar `picard AddOrReplaceReadGroups`.

Más información sobre este proceso en:
https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard

Las informaciónes que vamos a añadir al grupo de lectura serán:
- RGID: Read-Group ID
- RGLB: Read-Group library
- RGPL: Read-Group platform
- RGPU: Read-Group platform unit (eg. run barcode)
- RGSM: Read-Group sample name

```
picard AddOrReplaceReadGroups \
    -I bams/LP220.sorted.bam \
    -O bams/LP220.sorted.rg.bam \
    -RGID LP220-seq \
    -RGLB protocolo1 \
    -RGPL Illumina \
    -RGPU LP220-seq-1 \
    -RGSM LP220 \
    -VALIDATION_STRINGENCY SILENT

picard AddOrReplaceReadGroups \
    -I bams/LL240.sorted.bam \
    -O bams/LL240.sorted.rg.bam \
    -RGID LL240-seq \
    -RGLB protocolo1 \
    -RGPL Illumina \
    -RGPU LL240-seq-1 \
    -RGSM LL240 \
    -VALIDATION_STRINGENCY SILENT

# para visualizar el grupo de lecturas en el cabecero:
samtools view -H bams/LP220.sorted.rg.bam
samtools view -H bams/LL240.sorted.rg.bam
```

IGV


----
```
mkdir -m 777 qualimap

qualimap bamqc \
  -bam bams/LP220.sorted.rg.bam \
  --java-mem-size=19G \
  -outfile LP220.html \
  -outformat html \
  -outdir qualimap/LP220

qualimap bamqc \
  -bam bams/LL240.sorted.rg.bam \
  -outfile LL240_qualimap.html \
  -outformat html \
  -outdir qualimap/LL240
```
----

## Llamada de variantes

Utilizaremos GATK para generar un file .g.vcf para cada uno de nuestros individuos usando el flag `-ERC GVCF`:

```
mkdir -m 777 vcfs

# SUPER SLOW!
samtools index bams/LP220.sorted.rg.bam
gatk HaplotypeCaller \
    -R data/reference_sequence.fa \
    -I bams/LP220.sorted.rg.bam \
    -O vcfs/LP220.g.vcf \
    --native-pair-hmm-threads 1 \
    -ERC GVCF

# SUPER SLOW!
samtools index bams/LL240.sorted.rg.bam
gatk HaplotypeCaller \
    -R data/reference_sequence.fa \
    -I bams/LL240.sorted.rg.bam \
    -O vcfs/LL240.g.vcf \
    --native-pair-hmm-threads 1 \
    -ERC GVCF

less -S vcfs/LP220.g.vcf
less -S vcfs/LL240.g.vcf
```

Para juntar los GVCF de diferentes muestras y convertirlos al más compacto VCF usamos `gatk CombineGVCFs` y `gatk GenotypeGVCFs`:

```
gatk CombineGVCFs \
    -R data/reference_sequence.fa \
    --variant vcfs/LP220.g.vcf \
    --variant vcfs/LL240.g.vcf \
    -O vcfs/LP220_LL240.g.vcf

less -S vcfs/LP220_LL240.g.vcf

gatk GenotypeGVCFs \
    -R data/reference_sequence.fa \
    -V vcfs/LP220_LL240.g.vcf \
    -O vcfs/LP220_LL240.vcf

less -S vcfs/LP220_LL240.vcf
```

## Filtrado de variantes

Ahora que tenemos el VCF con nuestras variantes, podemos proceder a varios pasos de filtrado para eliminar las que no nos interesan.

Hay diferentes formas de filtrar los datos, que dependen del criterio que queramos usar:

1. variantes en regiones del genoma que no nos interesan (e.g. baja complejidad, paralogos, genes)
2. variantes de un tipo particular que no vamos a analizar (e.g. INDELs, no-bialelicas, fijadas)
3. variantes con algún valor cuantitativo que no queremos (e.g. calidad, missing data, profundidad, heterozigosidad, hardy-weinberg)

### 1. variantes en regiones del genoma que no nos interesan 

Vamos a ver un ejemplo para eliminar variantes que vienen de regiones repetitivas del genoma. Para esto usamos `bedtools subtract`

https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html

```
less -S data/repetitive_regions.bed

bedtools subtract -a vcfs/LP220_LL240.vcf -b data/repetitive_regions.bed -header |
    uniq > vcfs/LP220_LL240.filter1.vcf

less -S vcfs/LP220_LL240.filter1.vcf
```

### 2. variantes de un tipo particular que no vamos a analizar

Vamos a ver un ejemplo para quedarnos con variantes bialelicas y que no sean INDELs utilizando `vcftools`:

```
vcftools --vcf vcfs/LP220_LL240.filter1.vcf \
    --min-alleles 2 --max-alleles 2 --remove-indels \
    --recode --recode-INFO-all \
    --out vcfs/LP220_LL240.filter2

ls -lhtr vcfs/
 ```

Una forma alternativa es usando `gatk`:

```
gatk SelectVariants \
    -R data/reference_sequence.fa \
    -V vcfs/LP220_LL240.filter1.vcf \
    -O vcfs/LP220_LL240.filter2.gatk.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC

grep -v "#" vcfs/LP220_LL240.filter2.gatk.vcf | wc -l
grep -v "#" vcfs/LP220_LL240.filter2.recode.vcf | wc -l
```

### 3. variantes con algún valor cuantitativo que no queremos

Un ejemplo claro de valor de una variante que puede no interesarnos es su frecuencia alelica. Si una variante está a frecuencia 1.00 en nuestra base de datos, quiere decir que no hay variación entre nuestros individuos y esta variante no será informativa. Para ver cuantas variantes a frecuencia 1.00 tenemos podemos utilizar `bcftools view`, con `-i` para incluir, y luego especificar que en la columna INFO queremos las que tengan AF=1.00:

```
# cuenta el numero de SNPs con AF=1
bcftools view -i 'INFO/AF=1.00' vcfs/LP220_LL240.filter2.recode.vcf | grep -v "#" | wc -l
```

Para eliminarlos podemos usar `bcftools view` pero ahora con `-e` para excluir los que tengan AF=1:

```
bcftools view -e 'INFO/AF=1.00' vcfs/LP220_LL240.filter2.recode.vcf > vcfs/LP220_LL240.filter3.vcf

less -S vcfs/LP220_LL240.filter3.vcf
```

Otras veces queremos eliminar los SNPs que tengan un valor por encima o por debajo de un umbral. Para decidir que valores tienen que tener esos umbrales puede ayudar dibujar su distribucción. Vamos a ver un ejemplo extrayendo los valores de [QD](https://gatk.broadinstitute.org/hc/en-us/articles/360037592191-QualByDepth) y [DP](https://gatk.broadinstitute.org/hc/en-us/articles/13832711528603-Coverage) de todos nuestros SNPs y luego vamos a decidir los umbrales a aplicar al filtrado basados en sus distribuciones.

Empezamos con extraer estos dos valores a una tabla:

```
# creamos la tabla con su cabecero
echo -e "QD\tDP" > qd_dp.table.tsv

# pegamos los valores de QD y DP extrayendolos del vcf directamente y los añadimos a la tabla
paste \
    <(grep -v "#" vcfs/LP220_LL240.filter3.vcf | cut -f8 | grep -o 'QD=[0-9.]\+' | cut -d'=' -f2) \
    <(grep -v "#" vcfs/LP220_LL240.filter3.vcf | cut -f8 | grep -o 'DP=[0-9.]\+' | cut -d'=' -f2) \
    >> qd_dp.table.tsv

less -S qd_dp.table.tsv
```

En R podemos dibujar sus distribuciones y extraer los valores limite para el filtrado:

```{r}
library(ggplot2)
qd_dp.table <- read.table("qd_dp.table.tsv", header = TRUE)

# que valor de QD está debajo de su media - 2 veces su desviación estandar?
qd_lim <- mean(qd_dp.table$QD) - 2 * sd(qd_dp.table$QD)

# vamos a visualizarlo
qd_plot <- ggplot() +
  geom_histogram(data = qd_dp.table, aes(x = QD), bins = 100) +
  scale_x_continuous(breaks = 0:42*2) +
  geom_vline(xintercept = qd_lim, linetype = "dashed",
             colour = "red", size = 0.3) +
  theme_bw()
ggsave(filename = "qd_plot.pdf", plot = qd_plot)

# cuantos SNPs filtramos si aplicamos este valor como limite?
qd_filter <- qd_dp.table[which(qd_dp.table$QD <= qd_lim), ]
nrow(qd_filter)

####

# que valor de DP está por encima de 2 veces su mediana?
dp_lim <- 2 * median(qd_dp.table$DP)

# vamos a visualizarlo
dp_plot <- ggplot() +
  geom_histogram(data = qd_dp.table, aes(x = DP), bins = 100) +
  scale_x_continuous(breaks = 0:300*20, limits = c(0, 300)) +
  geom_vline(xintercept=dp_lim, linetype="dashed", colour="red", size=0.3) +
  theme_bw()
ggsave(filename = "dp_plot.pdf", plot = dp_plot)


# cuantos SNPs filtramos si aplicamos este valor como limite?
dp_filter <- qd_dp.table[which(qd_dp.table$DP <= dp_lim), ]
nrow(dp_filter)

####

# escribimos una tabla para guardar estos limites:
limits <- data.frame(QD_limit = qd_lim, DP_limit = dp_lim)
write.table(limits, file = "qd_dp.limits.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
```

Para aplicar estos filtros podemos utilizar de nuevo `bcftools view -e` leyendo los valores limite en la tabla que hemos generado:

```
qd_lim=$(cut -f1 qd_dp.limits.tsv | tail -1)
dp_lim=$(cut -f2 qd_dp.limits.tsv | tail -1)

bcftools view -e "INFO/QD < ${qd_lim} | INFO/DP > ${dp_lim}" vcfs/LP220_LL240.filter3.vcf > vcfs/LP220_LL240.filter4.vcf

grep -v "#" vcfs/LP220_LL240.filter3.vcf | wc -l
grep -v "#" vcfs/LP220_LL240.filter4.vcf | wc -l
```

Bcftools tiene también una opción para eliminar variantes con missing data. Podemos excluir todas las variantes con más del X% de missing data o con missing data en mas de N individuos de la siguiente forma:

```
# eliminamos lo que tenga missing data en el 10% de individuos:
bcftools view -e 'F_MISSING > 0.1' vcfs/LP220_LL240.filter4.vcf > vcfs/LP220_LL240.filter5.f_missing.vcf
# o en 1 o más individuos
bcftools view -e 'N_MISSING >= 1' vcfs/LP220_LL240.filter4.vcf > vcfs/LP220_LL240.filter5.n_missing.vcf

grep -v "#" vcfs/LP220_LL240.filter5.f_missing.vcf | wc -l
grep -v "#" vcfs/LP220_LL240.filter5.n_missing.vcf | wc -l
```

En nuestra base de datos solo tenemos dos individuos, pero en bases de datos más amplias los umbrales de missing data se podrían calcular de forma empirica como hemos hecho para el QD o el DP

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
mkdir fastqc

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
mkdir fastp

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

## Alinemiento

Creamos una carpeta donde guardar nuestros alineamientos:

```
mkdir bams
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
# indice general
samtools faidx data/reference_sequence.fa
# diccionario para gatk
samtools dict data/reference_sequence.fa -o data/reference_sequence.dict
```

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
----
```
mkdir qualimap

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

## Calling

Utilizaremos GATK para generar un file .g.vcf para cada uno de nuestros individuos usando el flag `-ERC GVCF`:

```
mkdir vcfs

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

gatk GenotypeGVCFs \
    -R data/reference_sequence.fa \
    -V vcfs/LP220_LL240.g.vcf \
    -O vcfs/LP220_LL240.vcf
```

## Filtering

Ahora que tenemos el VCF con nuestras variantes, podemos proceder a varios pasos de filtrado para eliminar las que no nos interesan.

Hay diferentes formas de filtrar los datos, que dependen del criterio que queramos usar:

1. variantes en regiones del genoma que no nos interesan (e.g. baja complejidad, paralogos, genes)
2. variantes de un tipo particular que no vamos a analizar (e.g. INDELs, no-bialelicas, fijadas)
3. variantes con algún valor cualitativo que no queremos (e.g. calidad, missing data, profundidad, heterozigosidad, hardy-weinberg)

### 1. variantes en regiones del genoma que no nos interesan 

Vamos a ver un ejemplo para eliminar variantes que vienen de regiones repetitivas del genoma. Para esto usamos `bedtools subtract`

https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html

```
less -S data/repetitive_regions.bed

bedtools subtract -a vcfs/LP220_LL240.vcf -b data/repetitive_regions.bed -header |
    uniq > vcfs/LP220_LL240.filter1.vcf

less -S vcfs/LP220_LL240.filter1.vcf
```

### biallelic noindel

```
gatk SelectVariants \
    -R data/reference_sequence.fa \
    -V vcfs/LP220_LL240.filter1.vcf \
    -O vcfs/LP220_LL240.filter2.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC
```

### variant

```
bcftools view -e 'INFO/AF=1.00' vcfs/LP220_LL240.filter2.vcf > vcfs/LP220_LL240.filter3.vcf
```

### quality scores

```
# look at distributions first?

bcftools view -e 'INFO/MQRankSum < -12.5 | INFO/ReadPosRankSum < -8.0 | INFO/QD < 6.0 | INFO/FS > 60.0 | INFO/MQ < 40.0' \
    vcfs/LP220_LL240.filter3.vcf > vcfs/LP220_LL240.filter4.vcf
```

### read depth

mosdepth?

```

```

### Final summary of calling and filtering

nvariants and missing data

```
bcftools view -e '' # missing data proportion
```

### Other common filtering

hardy weinberg
heterozygosity
genes

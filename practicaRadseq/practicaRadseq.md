practicaRadseq
================
Mari Jose
5/6/2024

# Análisis de RADseq usando STACKS

## Objetivos

El objetivo de este ejercicio es familiarizar a los estudiantes con el
uso de datos secuenciados de secuenciación masiva producidos a partir de
bibliotecas de representación reducida (RRL), como las etiquetas
asociadas a sitios de restricción (RAD). Estas bibliotecas se utilizan a
menudo para el genotipado por secuenciación, y pueden proporcionar un
conjunto denso de SNPs que se distribuyen uniformemente a través de un
genoma. Los estudiantes adquirirán experiencia con el programa
informático llamado Stacks, diseñado para el análisis de este tipo de
datos. Se analizarán datos de un organismo sin genoma de referencia y se
hará un análisis de novo para identifar SNPs.

Los estudiantes aprenderán a: 1. Procesar datos RAD de Illumina para su
análisis eliminando lecturas de baja calidad y demultiplexando un
conjunto de muestras con barcodes. 2. Utilizar Stacks para ensamblar
sitios RAD de novo. 3. Llamar SNPs, genotipos y haplotipos de estos
individuos dentro de Stacks.

# Datasets

1- Archivo fastq que contiene las secuencias: SRR034310.fastq 2- Archivo
txt que contiene los barcodes usados: Barcode_SRR034310.txt 3- Archivo
txt que contiene los datos de información poblacional popmap.txt

\#Software

Stacks v.2 <https://catchenlab.life.illinois.edu/stacks/>

\#Ejercicio 1. Preparación de datos.

El primer paso en el análisis de todos los datos de secuenciación de
nueva generación, incluidos los datos RAD, es eliminar las secuencias de
baja calidad y separar las lecturas de las distintas muestras que tenían
códigos de barras individuales. Este “desmultiplexado” sirve para
asociar las lecturas a los distintos individuos o muestras de población
de los que proceden.

1- Vemos el tipo de archivo y secuencias que tenemos y analizamos con
FastQC los datos

``` bash
cd practicaRadseq
fastqc SRR034310.fastq
```

2- Una vez que hemos visto como es la calidad, limpiamos los datos
usando el comando process_radtags. Al que le tenemos que decir qué datos
hay que procesar y còmo hay que hacerlo.

Process Radtags para demultiplexar las lecturas: “Archivo con secuencias
Single-end o paired-end”: Single-end “Archivo que contiene las
secuencias”: SRR034310.fastq(.gz) “Archivo con los barcodes”:
Barcodes_SRR034310.tabular “Número de enzimas”: Una “Enzima”: sbfI
“Recuperar secuencias de baja calidad”: No “Formato del Output:” fastq
“Dónde quiero poner el output”: carpeta demultiplexed “Quito reads con
N”: Sí “Quito reads de baja calidad”: Sí

``` bash
mkdir demultiplexed

process_radtags -f ./SRR034310.fastq -o ./demultiplexed/ -b ./Barcode_SRR034310.txt -e sbfI -c -q
```

¿Cuántas lecturas había en el conjunto de datos original? ¿Cuántas se
conservan? ¿Puede explicar por qué perdemos tantas lecturas? ¿Qué tipo
de información proporciona este resultado en relación con el próximo
análisis de datos y el diseño de los códigos de barras en general? La
información se puede encontrar en el archivo de registro de resultados
(.log)

Los parámetros se pueden modificar en función de los datos que tengamos
y la calidad de las secuencias.

Ahora vamos a reejecutar el ódulo jugando con los parámetros “Descartar
lecturas con puntuaciones de calidad bajas”: Sí…. q “score limit”: 20
(por ejemplo; juegue con esto) -score-limit por defecto es 10, lo
pondría en 20 “Establecer el tamaño de la ventana deslizante como una
fracción de la longitud de la lectura, entre 0 y 1”: 0.30 (por ejemplo;
juegue con esto)-window-size (por defecto es 0.15)

``` bash
mkdir demultiplexedquality
process_radtags -f ./SRR034310.fastq -o ./demultiplexedquality/ -b ./Barcode_SRR034310.txt -e sbfI -r -c -q --score-limit 20 --window-size 0.30
```

Volvemos a las preguntas de arriba, ¿ha cambiado algo?

\#Ejercicio 2. Ensamblaje de novo de RADseq

Ahora que tenemos los datos limpios vamos a hacer el ensamblaje de esos
datos. Esto en STACKs se hace usando varios módulos de forma secuencial
que permiten caracterizar los loci que tenemos e identificar los SNPs.

Estos módulos son ustacks, cstacks, sstacks, tsv2bam, gstacks y
populations. Pero se pueden ejecutar todos seguidos usando el comando
denovo_map.pl. En esta práctica vamos a correr este comando pero
mientras corre vamos a explicar que hace cada uno de los módulos y cómo
de optimizan los parámetros a usar. Brevemente, denovo_map realiza
varias etapas, incluyendo:

1-Ejecutar ustacks en cada una de las muestras especificadas,
construyendo loci de novo en cada muestra. 2-Ejecutar cstacks para crear
un catálogo de todos los loci de la población (o sólo de los padres si
se procesa un mapa genético). Los loci se emparejan entre muestras según
la similitud de secuencia. 3-Ejecutar sstacks para cotejar cada muestra
con el catálogo. En el caso de un mapa genético, los padres y la
progenie se comparan con el catálogo. 4-Ejecutar programa tsv2bam para
transponer los datos de la organización por muestras a la organización
por locus RAD. Esto permite a los componentes de la cadena de producción
posterior examinar todos los datos de la metapoblación para cada locus.
Si se dispone de lecturas de extremo pareado, tsv2bam las obtendrá y
almacenará con las lecturas de extremo único para su uso posterior. 5-
Ejecutar gstacks. Si se dispone de lecturas paired-end, se ensamblará un
contig a partir de ellas y se solapará con el locus single-end. Se
llamará a los SNPs, y para cada locus, los SNPs se clasificarán en
haplotipos. 6- Ejecutar el programa de poblaciones para generar
estadísticas resumidas a nivel de población y exportar los datos en
diversos formatos.

Para preparar los datos lo primero hay que tener un archivo popmap, con
los datos de las muestras que se van a analizar. Abrid el archivo
popmap.txt, y mirad que muestras están incluidas. Esas son las que vamos
a analizar.

Estos son los parámetros que vamos a usar para el análisis.

“Número mínimo de lecturas en bruto idénticas necesarias para crear una
pila”: 3 “Número de desajustes permitidos entre loci al construir el
catálogo” 3 “Archivo con datos de la población”popmap.txt

``` bash
mkdir SNPs
denovo_map.pl -M 3 -n 3 -T 6  -o ./SNPs --popmap ./popmap.txt --samples ./demultiplexedquality
```

Mientras corre este comando, vamos a ver cómo se optimizan los
diferentes parámetros, y las utilidades de los diferentes módulos.

Al final podemos analizar los diferentes archivos, y vemos que tenemos
diferentes tipos de archivos que contienen diferentes tipos de
información de la pipeline.

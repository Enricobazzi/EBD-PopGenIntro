# EBD - PopGenIntro

```
docker build -t ebd-popgenintro .
docker run -it -v /Users/enrico/Documents/Curso_Gabinete/EBD-PopGenIntro:/root ebd-popgenintro /bin/bash
```

```
dir=data/practica3

plink --vcf ${dir}/data1.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --indep-pairwise 50 1 0.2 \
    --out ${dir}/data1.plink

plink --vcf ${dir}/data1.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --exclude ${dir}/data1.plink.prune.out --recode A \
    --out ${dir}/data1.plink
```
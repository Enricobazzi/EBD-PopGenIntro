vcftools --vcf data/data_selection.pop1.recode.vcf     --window-pi 50000 --out data/data_selection.pop1
cd practica3/
vcftools --vcf data/data_selection.pop1.recode.vcf     --window-pi 50000 --out data/data_selection.pop1
vcftools --vcf data/data_selection.pop2.recode.vcf     --window-pi 50000 --out data/data_selection.pop2
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt     --recode --recode-INFO-all --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop2.txt     --recode --recode-INFO-all --out data/data_selection.pop2
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt     --recode --recode-INFO-all --out data/data_selection.pop1
vcf_file="data/data_selection.vcf"
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt     --recode --recode-INFO-all --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt --maf=0.001    --recode --recode-INFO-all --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt --maf 0.001    --recode --recode-INFO-all --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt --maf 0.001    --recode  --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt --maf 0.001    --recode --recode-INFO-all --out data/data_selection.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop1.txt --maf 0.001    --recode  --out data/data_selection_norecode.pop1
vcftools --vcf $vcf_file --keep data/data_selection.pop2.txt     --maf 0.001 --recode --recode-INFO-all --out data/data_selection.pop2
vcftools --vcf data/data_selection.pop1.recode.vcf     --window-pi 50000 --out data/data_selection.pop1
vcftools --vcf data/data_selection.pop2.recode.vcf     --window-pi 50000 --out data/data_selection.pop2
#Calculate Tajima D per population
vcftools --vcf data/data_selection.pop1.recode.vcf     --TajimaD 50000 --out data/data_selection.pop1
vcftools --vcf data/data_selection.pop2.recode.vcf     --TajimaD 50000 --out data/data_selection.pop2
vcftools --vcf $vcf_file     --weir-fst-pop data/data_selection.pop1.txt     --weir-fst-pop data/data_selection.pop2.txt     --fst-window-size 50000     --out data/data_selection.pop1_vs_pop2    
pi_pop1 <- read.table("data/data_selection.pop1.windowed.pi", header = TRUE)
Td_pop1 <- read.table("data/data_selection.pop1.Tajima.D", header = TRUE)

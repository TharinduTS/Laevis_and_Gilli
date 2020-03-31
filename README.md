# Laevis_and_Gilli
# Mapped and created VCF file usind the same way in


https://github.com/TharinduTS/allofraseri.md/blob/master/README.md

# Filtered VCF file once again using VCF file

```bash
vcftools --gzvcf ../all_sorted_bam/xlaevis_and_xgilli_sorted.bam.vcf.gz --minGQ 20 --minDP 25 --recode --recode-INFO-all --out ../final_vcf/xlaevis_and_xgilli_all_final
```

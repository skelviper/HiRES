#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate charm

mkdir ./stat

find ./ -name "raw.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#"| wc -l' | sort > ./stat/raw.pairs.stat
find -L ./Rawdata/ -name "*R1*.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/raw.fq.stat
find ./processed -name "*.rna.clean.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/rna.fq.stat
find ./processed -name "*.dna.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/dna.fq.stat
find ./processed -name "contacts.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.dedup.stat
find ./result/cleaned_pairs/c1 -name "*.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#" | wc -l' | sort > ./stat/pairs.c1.stat
find ./result/cleaned_pairs/c12 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c12.stat
find ./result/cleaned_pairs/c123 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c123.stat
for i in `find ./result/cleaned_pairs/c123 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c123.stat

for i in `find ./result/cleaned_pairs/c12 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c12.stat

find processed/ -name "*rms.info" | xargs -I {} grep --with-filename "top3 RMS RMSD" {} | sed -e "s/processed\///g" -e "s/\/3d_info\/50k.align.rms.info:\[M::__main__\] top3 RMS RMSD: /\t/g" | sort > ./stat/rmsd.info

find processed -name "*.yperx.txt" | parallel 'paste <(echo {}) <(cat {})'| sort > ./stat/yperx.stat

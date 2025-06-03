#!/bin/bash

start=$(date +%s.%N)

eval "$($HOME/micromamba/bin/micromamba shell hook -s posix)"

INPUT=$HOME/OneDrive/IGM_PVM_MINION_LIG_LIBRARY20250519/ANALYSIS/02demux/IGM_PVM_MINION_LIG_LIBRARY20250519-DEMUX-sup
ANALYSIS_DIR=$HOME/OneDrive/IGM_PVM_MINION_LIG_LIBRARY20250519/ANALYSIS/03assembly/sup-AB273635.1
DEPTH=20
SAMPLE_ID=MT2-FK
PRIMER_SCHEME=HTLV1_V1
MIN=100
MAX=2000

THREADS=$(nproc)
PRIMER_SCHEME_NAME=$(echo "$PRIMER_SCHEME" | awk -F_ '{print $1}')
PRIMER_SCHEME_DIR=$(echo "$PRIMER_SCHEME" | sed -e 's/ARTIC_V1/ARTIC\/V1/g' -e 's/ARTIC_V2/ARTIC\/V2/g' -e 's/ARTIC_V3/ARTIC\/V3/g' -e 's/ARTIC_V4/ARTIC\/V4/g' -e 's/ARTIC\/V4_1/ARTIC\/V4.1/g' -e 's/ARTIC_V5_3_2/ARTIC\/V5.3.2/g' -e 's/ChikAsianECSA_V1/ChikAsianECSA\/V1/g' -e 's/DENGUESEQ1_V1/DENGUESEQ1\/V1/g' -e 's/DENGUESEQ2_V1/DENGUESEQ2\/V1/g' -e 's/DENGUESEQ3_V1/DENGUESEQ3\/V1/g' -e 's/DENGUESEQ4_V1/DENGUESEQ4\/V1/g' -e 's/FIOCRUZ-IOC_V1/FIOCRUZ-IOC\/V1/g' -e 's/FIOCRUZ-IOC_V2/FIOCRUZ-IOC\/V2/g' -e 's/HTLV1_V1/HTLV1\/V1/g' -e 's/LassaL_V1/LassaL\/V1/g' -e 's/LassaS_V1/LassaS\/V1/g' -e 's/MIDNIGHT_V1/MIDNIGHT\/V1/g' -e 's/MIDNIGHT_V2/MIDNIGHT\/V2/g' -e 's/MPXV_V1/MPXV\/V1/g' -e 's/Nipah_V1/Nipah\/V1/g' -e 's/OROV400L_V1/OROV400L\/V1/g' -e 's/OROV400M_V1/OROV400M\/V1/g' -e 's/OROV400S_V1/OROV400S\/V1/g' -e 's/WNV400_V1/WNV400\/V1/g' -e 's/YFV500_V1/YFV500\/V1/g' -e 's/YFV1000_V1/YFV1000\/V1/g' -e 's/ZaireEbola_V1/ZaireEbola\/V1/g' -e 's/ZikaAsian_V1/ZikaAsian\/V1/g' -e 's/ZikaAsian_V2/ZikaAsian\/V2/g')
if [[ $(echo "$INPUT" | awk -F/ '{print $NF}' | awk -F- '{print $NF}') == "fast" ]]; then
    CLAIR3_MODEL=r1041_e82_400bps_fast_g632
elif [[ $(echo "$INPUT" | awk -F/ '{print $NF}' | awk -F- '{print $NF}') == "hac" ]]; then
    CLAIR3_MODEL=r1041_e82_400bps_hac_v520
elif [[ $(echo "$INPUT" | awk -F/ '{print $NF}' | awk -F- '{print $NF}') == "sup" ]]; then
    CLAIR3_MODEL=r1041_e82_400bps_sup_v520
fi

micromamba activate vigeas-ont

mkdir -p "$ANALYSIS_DIR"
cd "$ANALYSIS_DIR"
echo "$SAMPLE_ID"
echo "$PRIMER_SCHEME_DIR"

samtools fastq "$INPUT"/*unclassified.bam >> "$ANALYSIS_DIR"/"$SAMPLE_ID".fastq

[[ ! -f $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta.fai ]] && samtools faidx $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta

bbduk.sh minlen="$MIN" maxlen="$MAX" minavgquality=10 in="$ANALYSIS_DIR"/"$SAMPLE_ID".fastq out="$ANALYSIS_DIR"/"$SAMPLE_ID".filt.fastq
dedupe.sh in="$ANALYSIS_DIR"/"$SAMPLE_ID".filt.fastq out="$ANALYSIS_DIR"/"$SAMPLE_ID".filtDedup.fastq
minimap2 -t "$THREADS" -ax map-ont $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta "$ANALYSIS_DIR"/"$SAMPLE_ID".filtDedup.fastq | ivar trim -e -m 32 -q 8 -b $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".scheme.bed > "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.bam
samtools sort -o "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.bam
samtools index "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam
run_clair3.sh -t "$THREADS" -p ont --chunk_size=10000 --enable_long_indel --haploid_precise --include_all_ctgs --no_phasing_for_fa --min_coverage="$DEPTH" -m $HOME/rerio/clair3_models/"$CLAIR3_MODEL" --sample_name="$SAMPLE_ID" --bam_fn="$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam --ref_fn=$HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta --bed_fn=$HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".scheme.bed -o "$ANALYSIS_DIR"/vcall
bgzip -dc "$ANALYSIS_DIR"/vcall/merge_output.vcf.gz > "$ANALYSIS_DIR"/"$SAMPLE_ID".vcf
bcftools view -i 'QUAL<10 || DP<=1' "$ANALYSIS_DIR"/"$SAMPLE_ID".vcf > "$ANALYSIS_DIR"/"$SAMPLE_ID".fail.vcf
bcftools view -i 'QUAL>=10 && DP>=20' "$ANALYSIS_DIR"/"$SAMPLE_ID".vcf | awk 'BEGIN {OFS="\t"} /^#/ {print; next} {ref_len=length($4); alt_len=length($5); diff=alt_len-ref_len; if (diff == 0 || diff % 3 == 0) print;}' > "$ANALYSIS_DIR"/"$SAMPLE_ID".pass.vcf
bgzip -f "$ANALYSIS_DIR"/"$SAMPLE_ID".pass.vcf
tabix -p vcf "$ANALYSIS_DIR"/"$SAMPLE_ID".pass.vcf.gz
bcftools query -f '%CHROM\t%POS\n' "$ANALYSIS_DIR"/"$SAMPLE_ID".fail.vcf > "$ANALYSIS_DIR"/"$SAMPLE_ID".nmask.txt
samtools depth -a "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '$3 < 10 {print $1 "\t" $2}' >> "$ANALYSIS_DIR"/"$SAMPLE_ID".nmask.txt && sort -k2,2n "$ANALYSIS_DIR"/"$SAMPLE_ID".nmask.txt -o "$ANALYSIS_DIR"/"$SAMPLE_ID".nmask.txt
bcftools consensus -f $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta "$ANALYSIS_DIR"/"$SAMPLE_ID".pass.vcf.gz -m "$ANALYSIS_DIR"/"$SAMPLE_ID".nmask.txt > "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa
sed -i -e 's/>.*/>'${SAMPLE_ID}'/g' "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa
muscle -align <(awk 'NR==FNR{print; next} {print $1}' "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta) -output "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".aligned.fa
micromamba deactivate
micromamba activate vigeas-stats
samtools depth -a "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam > "$ANALYSIS_DIR"/"$SAMPLE_ID".depth.tsv
[[ ! -f "$ANALYSIS_DIR"/.summary ]] && echo "primer_scheme#sample_id#num_total_reads#num_mapp_reads#avg_depth#depth_20x#depth_100x#depth_1000x#ref_cov#ncount#ncount_perc" | tr '#' '\t' > "$ANALYSIS_DIR"/.summary
echo -n "#" | tr '#' '\n' >> "$ANALYSIS_DIR"/.summary
echo -n "$PRIMER_SCHEME""#" | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
echo -n "$SAMPLE_ID""#" | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
seqkit stats "$ANALYSIS_DIR"/"$SAMPLE_ID".fastq | sed -n '2p' | awk '{print $4}' | sed 's/,//g' | awk '{printf $0"#"}' | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
samtools view -F 4 -c "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
AVG_DEPTH=$(samtools depth "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '{sum+=$3} END {print sum/NR}')
if [[ "$AVG_DEPTH" == "" || "$AVG_DEPTH" == 0 ]]; then
    echo "0.00""#" | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
else
    echo "$AVG_DEPTH" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
fi
paste <(samtools depth "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '{if ($3 > '"20"') {print $0}}' | wc -l) <(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
paste <(samtools depth "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '{if ($3 > '"100"') {print $0}}' | wc -l) <(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
paste <(samtools depth "$ANALYSIS_DIR"/"$SAMPLE_ID".trimmed.sorted.bam | awk '{if ($3 > '"1000"') {print $0}}' | wc -l) <(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$ANALYSIS_DIR"/.summary
COVERAGE=$(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta 2> /dev/null | awk '{print $1}')
if [[ "$COVERAGE" == 0 ]]; then
    echo "0.00" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
    fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk -F" " '{print $1}' | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
    echo "100.00" | awk '{printf $0"\n"}' >> "$ANALYSIS_DIR"/.summary
else
    paste <(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk '{print $1}') <(seqtk comp "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa | awk -F"\t" '{print $9}') | awk -F"\t" '{printf("%0.2f\n", ($1-$2)/$1*100)}' | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
    seqtk comp "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa | awk -F"\t" '{print $9}' | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$ANALYSIS_DIR"/.summary
    paste <(seqtk comp "$ANALYSIS_DIR"/"$SAMPLE_ID".depth"$DEPTH".fa | awk -F"\t" '{print $9}') <(fastalength $HOME/vigeas/primer_schemes/"$PRIMER_SCHEME_DIR"/"$PRIMER_SCHEME_NAME".reference.fasta | awk '{print $1}')| awk -F"\t" '{printf("%0.2f\n", ($1/$2)*100)}'| awk '{printf $0"\n"}' >> "$ANALYSIS_DIR"/.summary
fi
Rscript $HOME/vigeas/scripts/depthCoverage.R "$ANALYSIS_DIR"/"$SAMPLE_ID".depth.tsv "$SAMPLE_ID" "$PRIMER_SCHEME_DIR" "$ANALYSIS_DIR"/"$SAMPLE_ID"
micromamba deactivate
sed '/^[[:space:]]*$/d' "$ANALYSIS_DIR"/.summary > "$ANALYSIS_DIR"/"$SAMPLE_ID".summary.$(uname -n).$(date +'%Y-%m-%d').txt
rm -rf "$ANALYSIS_DIR"/.summary*

end=$(date +%s.%N)
runtime=$(python3 -c "print(${end} - ${start})")

echo "" && echo "Done. The runtime was "$runtime" seconds." && echo ""

#!/bin/bash
#SBATCH --job-name=vigeas
#SBATCH --partition=fat
#SBATCH --cpus-per-task=256
#SBATCH --mem=4000G
#SBATCH --output=fat_%j.out
#SBATCH --error=fat_%j.err
#SBATCH --account=thiago.parente
#SBATCH --ntasks=1

# ill-hyb-direct-assembly.slurm <fastq_folder> <sample_sheet> <ref_seq_list> <se>

INPUT_DIR=$1
SAMPLE_SHEET=$2
REF_SEQ_LIST=$3
SE=$4 # 0 = paired-end; 1 = single-end
DEPTH=10
PIDENT=90
THREADS=$SLURM_CPUS_PER_TASK

module purge
module load bbmap
module load bwa-mem2 # missing
module load blast+
module load entrez-direct # missing
module load exonerate
module load fastp
module load ivar
module load mafft
module load samtools
module load seqkit
module load seqtk # missing
module load spades

LIBRARY_NAME="$(basename $SAMPLE_SHEET | awk -F. '{print $1}')"
ASSEMBLY_DIR="$HOME/proj.laise.moraes/vigeas/${LIBRARY_NAME}_ASSEMBLY"

[[ ! -d "$ASSEMBLY_DIR" ]] && mkdir "$ASSEMBLY_DIR" && chmod 700 -R "$ASSEMBLY_DIR"
echo "" && echo "Output path: $ASSEMBLY_DIR"
echo "" && echo "Log Analysis: $ASSEMBLY_DIR/$LIBRARY_NAME.log.$(uname -n).$(date +'%Y-%m-%d').txt"
echo "" && echo "To display the log output in fulscreen: watch tail -n 20 $ASSEMBLY_DIR/$LIBRARY_NAME.log.$(uname -n).$(date +'%Y-%m-%d').txt" && echo ""

if grep -q "Experiment Name" "$SAMPLE_SHEET"; then
    cat "$SAMPLE_SHEET" | tr -dc '[:print:]\n' | sed -n '/^Sample_ID/,$p' | sed '1d' | awk -F, > "$ASSEMBLY_DIR/.samplesheet"
    SAMPLE_SHEET="$ASSEMBLY_DIR/.samplesheet"
fi

bg() {
    start=$(date +%s.%N)
    module list && echo ""
    while IFS= read -r ID; do
        R1=$(find "$INPUT_DIR" -type f -regex ".*/${ID}.*\(_R1\|_1\|\.R1\).*\.fastq\.gz" -print -quit)
        R2=$(find "$INPUT_DIR" -type f -regex ".*/${ID}.*\(_R2\|_2\|\.R2\).*\.fastq\.gz" -print -quit)
        if [[ -n "$R1" ]]; then
            cp -v "$R1" "$ASSEMBLY_DIR/${ID}.R1.fastq.gz"
            [[ -n "$R2" ]] && cp -v "$R2" "$ASSEMBLY_DIR/${ID}.R2.fastq.gz"
            continue
        fi
    done < <(awk -F, '{print $1}' "$SAMPLE_SHEET")

    [[ ! -f "$ASSEMBLY_DIR/.summary" ]] && {
        echo -e "sample_id\tref_seq\ttarget\tsegment\tflu_vaccine\tnum_total_reads\tnum_mapp_reads\tavg_depth\tref_cov" > "$ASSEMBLY_DIR/.summary"
    }
    while IFS= read -r DESCRIPTION || [ -n "$DESCRIPTION" ]; do
        SAMPLE_ID=$(echo "$DESCRIPTION" | awk -F, '{print $1}')
        PANEL=$(echo "$DESCRIPTION" | awk -F, '{print $2}')
        echo "" && echo "$SAMPLE_ID" && echo "$PANEL"
        mkdir "$ASSEMBLY_DIR"/"$SAMPLE_ID"
        if [[ "$SE" == 1 ]]; then
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" ]] && fastp --thread "$THREADS" -f 0 -t 0 -F 0 -T 0 --cut_front --cut_tail --qualified_quality_phred 20 -l 75 -i "$ASSEMBLY_DIR/$SAMPLE_ID.R1.fastq.gz" -o "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" -h "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.html" -j "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.json"
            [[ ! -d "$ASSEMBLY_DIR/$SAMPLE_ID/denovo" ]] && spades.py -t "$THREADS" --metaviral -s "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" -o "$ASSEMBLY_DIR/$SAMPLE_ID/denovo"
        else
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" ]] && fastp --thread "$THREADS" -f 0 -t 0 -F 0 -T 0 --cut_front --cut_tail --qualified_quality_phred 20 -l 75 -i "$ASSEMBLY_DIR/$SAMPLE_ID.R1.fastq.gz" -I "$ASSEMBLY_DIR/$SAMPLE_ID.R2.fastq.gz" -o "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" -O "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R2.fastq.gz" -h "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.html" -j "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.json"
            [[ ! -d "$ASSEMBLY_DIR/$SAMPLE_ID/denovo" ]] && spades.py -t "$THREADS" --metaviral -1 "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" -2 "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R2.fastq.gz" -o "$ASSEMBLY_DIR/$SAMPLE_ID/denovo"
        fi
        if [ -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fasta" ]; then
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" ]] && blastn -num_threads "$THREADS" -word_size 21 -query "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fasta" -db "$HOME/vigeas/blastdb/$PANEL" -out "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" -outfmt "6 qseqid sseqid pident length"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fastq" ]] && seqtk seq -F '!' "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fasta" > "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fastq"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fastq" ]] && gzip "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fastq"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" ]] && cat "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/transcripts.fastq" > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz"
        elif [ -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fasta" ]; then
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" ]] && blastn -num_threads "$THREADS" -word_size 21 -query "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fasta" -db "$HOME/vigeas/blastdb/$PANEL" -out "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" -outfmt "6 qseqid sseqid pident length"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fastq" ]] && seqtk seq -F '!' "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fasta" > "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fastq"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fastq" ]] && gzip "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fastq"
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" ]] && cat "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" "$ASSEMBLY_DIR/$SAMPLE_ID/denovo/before_rr.fastq.gz" > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz"
        else
            [[ ! -f "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" ]] && seqkit fq2fa "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R1.fastq.gz" | blastn -db "$HOME/vigeas/blastdb/$PANEL" -out "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.blastn.txt" -outfmt "6 qseqid sseqid pident length"
        fi
        while IFS= read -r REFSEQ; do
            echo ""
            {
                printf "%s\t" "${SAMPLE_ID}"
                grep "^$REFSEQ" "$HOME/vigeas/blastdb/$PANEL.tsv" | awk -F"\t" '{if ($3 == "" || $3 == " " || gsub(/\t/, "", $3) == 0) {printf "%s\t%s\t%s\t%s\t", $1, $2, $3, $4} else {printf "%s\t", $0}}'
            } >> "$ASSEMBLY_DIR/.summary"
            bwa-mem2 index -p "$HOME/vigeas/refseq/$REFSEQ" "$HOME/vigeas/refseq/$REFSEQ.fasta"
            if [[ "$SE" == 1 ]]; then
                bwa-mem2 mem -t "$THREADS" "$HOME/vigeas/refseq/$REFSEQ" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" | samtools sort - > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.concat.sorted.bam"
            else
                bwa-mem2 mem -t "$THREADS" "$HOME/vigeas/refseq/$REFSEQ" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R2.fastq.gz" | samtools sort - > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.concat.sorted.bam"
            fi
            samtools index "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.concat.sorted.bam"
            samtools consensus -m simple -a -aa -A -d 1 -@ "$THREADS" --show-del no --show-ins yes -f fasta "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.concat.sorted.bam" -o "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.concat.fasta"
            bwa-mem2 index -p "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.concat" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.concat.fasta"
            if [[ "$SE" == 1 ]]; then
                bwa-mem2 mem -t "$THREADS" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.concat" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" | samtools sort - > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam"
            else
                bwa-mem2 mem -t "$THREADS" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.concat" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.concat.R1.fastq.gz" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.trimmed.R2.fastq.gz" | samtools sort - > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam"
            fi
            samtools index "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam"
            {
                seqkit stats "$ASSEMBLY_DIR/$SAMPLE_ID.R1.fastq.gz" | awk 'NR==2{gsub(/,/, "", $4); printf "%s\t", $4}'
                samtools view -F 4 -c "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam" | awk '{printf "%s\t", $1}'
            } >> "$ASSEMBLY_DIR/.summary"
            samtools depth -a "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam" > "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.depth.tsv"
            samtools mpileup -a -A -d 0 --reference "$HOME/vigeas/refseq/$REFSEQ.fasta" -B "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam" | ivar variants -q 20 -t 0.03 -m "$DEPTH" -p "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ" -r "$HOME/vigeas/refseq/$REFSEQ.fasta"
            samtools mpileup -a -A -d 0 --reference "$HOME/vigeas/refseq/$REFSEQ.fasta" -B "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam" | ivar consensus -t 0 -m "$DEPTH" -p "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.depth${DEPTH}"
            sed -i -e "s/>.*/>${SAMPLE_ID}/g" "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.depth${DEPTH}.fa"
            {
                AVG_DEPTH=$(samtools depth "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.trimmed.sorted.bam" | awk '{sum+=$3} END {print (NR>0) ? sum/NR : 0}')
                printf "%.2f\t" "$AVG_DEPTH"
                REF_LEN=$(fastalength "$HOME/vigeas/refseq/$REFSEQ.fasta" | awk '{print $1}')
                N_COUNT=$(seqtk comp "$ASSEMBLY_DIR/$SAMPLE_ID/$SAMPLE_ID.$REFSEQ.depth${DEPTH}.fa" | awk '{print $9}')
                if [[ "$REF_LEN" -gt 0 ]]; then
                    REF_COV=$(echo "$N_COUNT $REF_LEN" | awk '{printf "%.2f", (($2 - $1) / $2) * 100}')
                else
                    REF_COV="0.00"
                fi
                printf "%s\n" "$REF_COV"
            } >> "$ASSEMBLY_DIR/.summary"
        done < "$REF_SEQ_LIST"
    done < "$SAMPLE_SHEET"
    sed '/^[[:space:]]*$/d' "$ASSEMBLY_DIR/.summary" > "$ASSEMBLY_DIR/$LIBRARY_NAME.summary.$(uname -n).$(date +'%Y-%m-%d').tsv"
    rm -rf "$ASSEMBLY_DIR/.samplesheet" "$ASSEMBLY_DIR/.summary" "$ASSEMBLY_DIR"/*.fastq.gz
    end=$(date +%s.%N)
    runtime=$(python3 -c "print(${end} - ${start})")
    echo "" && echo "Done. The runtime was $runtime seconds." && echo ""
}

bg &>>"$ASSEMBLY_DIR/$LIBRARY_NAME.log.$(uname -n).$(date +'%Y-%m-%d').txt" &

wait
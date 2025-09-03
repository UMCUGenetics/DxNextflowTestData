#!/usr/bin/env bash
# generate_bams.sh — mini RNA-seq BAMs (PAIRED) met 2 junctions binnen één GRCh38-gen
# - 3 samples: sample1/2/3 (+ .bai)
# - AUTO-fallback kiest zelf een geschikt gen
# - BSD-awk compatible; GRCh38 seqlengths via chr_len() (geen associative arrays)
# Usage:
#   bash generate_bams.sh gencode.v44.annotation.gtf[.gz] <GENE|AUTO> [OUTDIR] [PAIRS_PER_SAMPLE]
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <gencode.gtf[.gz]> <GENE|AUTO> [OUTDIR] [PAIRS]" >&2
  exit 1
fi

GTF="$1"
REQ_GENE="$2"              # specifieke gennaam of 'AUTO'
OUT="${3:-.}"
PAIRS="${4:-30}"           # read-pairs per sample
mkdir -p "$OUT"

# -- portable GRCh38 chrom lengths (met en zonder 'chr' prefix) --
chr_len() {
  case "$1" in
    chr1|1) echo 248956422 ;;
    chr2|2) echo 242193529 ;;
    chr3|3) echo 198295559 ;;
    chr4|4) echo 190214555 ;;
    chr5|5) echo 181538259 ;;
    chr6|6) echo 170805979 ;;
    chr7|7) echo 159345973 ;;
    chr8|8) echo 145138636 ;;
    chr9|9) echo 138394717 ;;
    chr10|10) echo 133797422 ;;
    chr11|11) echo 135086622 ;;
    chr12|12) echo 133275309 ;;
    chr13|13) echo 114364328 ;;
    chr14|14) echo 107043718 ;;
    chr15|15) echo 101991189 ;;
    chr16|16) echo 90338345 ;;
    chr17|17) echo 83257441 ;;
    chr18|18) echo 80373285 ;;
    chr19|19) echo 58617616 ;;
    chr20|20) echo 64444167 ;;
    chr21|21) echo 46709983 ;;
    chr22|22) echo 50818468 ;;
    chrX|X)  echo 156040895 ;;
    chrY|Y)  echo 57227415 ;;
    chrM|M|chrMT|MT) echo 16569 ;;
    *) echo "" ;;
  esac
}

# ---- temp files ----
TMP_GTF="$(mktemp)"; EXONS_FILE="$(mktemp)"; EXONS_FILTERED="$(mktemp)"
cleanup() { rm -f "$TMP_GTF" "$EXONS_FILE" "$EXONS_FILTERED"; }
trap cleanup EXIT

# Inlezen (gz of plain)
if [[ "$GTF" =~ \.gz$ ]]; then zcat -- "$GTF" > "$TMP_GTF"; else cat -- "$GTF" > "$TMP_GTF"; fi

# 1) Globale EXONS tabel: gene \t chr \t start \t end  (BSD-awk compatible)
awk 'BEGIN{FS=OFS="\t"}
  $3=="exon" {
    p=index($9,"gene_name \""); if (p==0) next
    rest=substr($9, p+length("gene_name \""))
    q=index(rest,"\""); if (q==0) next
    gene=substr(rest,1,q-1)
    print gene, $1, $4, $5
  }' "$TMP_GTF" \
| LC_ALL=C sort -k1,1 -k2,2 -k3,3n > "$EXONS_FILE"

# 2) Vind eerste bruikbaar exon-triplet (strikt → relaxed)
find_triplet_in_exons() {
  local file="$1" minEx="$2" minIntron="$3"
  awk -v me="$minEx" -v mi="$minIntron" 'BEGIN{FS=OFS="\t"}
    {
      gene=$1; chr=$2; st=$3; en=$4
      if (gene==p2_gene && gene==p1_gene && chr==p2_chr && chr==p1_chr) {
        len1=p1_en-p1_st+1; len2=p2_en-p2_st+1; len3=en-st+1
        i12=p2_st-p1_en-1;  i23=st-p2_en-1
        if (len1>=me && len2>=me && len3>=me && i12>=mi && i23>=mi) {
          print gene, chr, p1_st, p1_en, p2_st, p2_en, st, en, i12, i23
          exit
        }
      }
      p1_gene=p2_gene; p1_chr=p2_chr; p1_st=p2_st; p1_en=p2_en
      p2_gene=gene;    p2_chr=chr;    p2_st=st;    p2_en=en
    }' "$file"
}

CHOSEN_GENE=""
TRIPLET_LINE=""

# 3) Eerst requested gene (tenzij AUTO)
if [[ "$REQ_GENE" != "AUTO" ]]; then
  awk -v g="$REQ_GENE" 'BEGIN{FS=OFS="\t"} $1==g {print}' "$EXONS_FILE" > "$EXONS_FILTERED"
  if [[ -s "$EXONS_FILTERED" ]]; then
    TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILTERED" 30 20 || true)"
    [[ -z "$TRIPLET_LINE" ]] && TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILTERED" 20 8 || true)"
    [[ -z "$TRIPLET_LINE" ]] && TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILTERED" 15 1 || true)"
    [[ -n "$TRIPLET_LINE" ]] && CHOSEN_GENE="$REQ_GENE" || echo "Geen bruikbaar exon-triplet voor '$REQ_GENE'. AUTO-modus…" >&2
  else
    echo "Gen '$REQ_GENE' niet gevonden in GTF. AUTO-modus…" >&2
  fi
fi

# 4) AUTO: doorloop hele tabel
if [[ -z "$TRIPLET_LINE" ]]; then
  TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILE" 30 20 || true)"
  [[ -z "$TRIPLET_LINE" ]] && TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILE" 20 8 || true)"
  [[ -z "$TRIPLET_LINE" ]] && TRIPLET_LINE="$(find_triplet_in_exons "$EXONS_FILE" 15 1 || true)"
  [[ -z "$TRIPLET_LINE" ]] && { echo "Geen bruikbaar exon-triplet gevonden (ook niet relaxed)." >&2; exit 1; }
  CHOSEN_GENE="$(echo "$TRIPLET_LINE" | awk '{print $1}')"
fi

# 5) Parse waarden
read -r CHR EX1_ST EX1_EN EX2_ST EX2_EN EX3_ST EX3_EN INTRON12 INTRON23 \
  <<<"$(echo "$TRIPLET_LINE" | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}')"

LEN1=$((EX1_EN-EX1_ST+1)); LEN2=$((EX2_EN-EX2_ST+1)); LEN3=$((EX3_EN-EX3_ST+1))
SEG=50; ((SEG>LEN1))&&SEG=$LEN1; ((SEG>LEN2))&&SEG=$LEN2; ((SEG>LEN3))&&SEG=$LEN3; ((SEG<15))&&SEG=15

# Junction anchor starts
POS12=$((EX2_ST-SEG))   # ~ EX1_EN-SEG+1
POS23=$((EX3_ST-SEG))   # ~ EX2_EN-SEG+1

# Mate (R2) in exon2
LEN_R2=$((2*SEG)); ((LEN_R2>LEN2))&&LEN_R2=$LEN2
POS_R2=$(( EX2_ST + (LEN2>20 ? 5 : 1) )); (( POS_R2 + LEN_R2 - 1 > EX2_EN )) && POS_R2=$((EX2_EN-LEN_R2+1))

# @SQ LN: echte GRCh38 lengte waar mogelijk
LN="$(chr_len "$CHR")"; [[ -z "$LN" ]] && LN=$((EX3_EN + 1000))

echo "Gekozen gen: ${CHOSEN_GENE}"
echo "Exon-triplet: ${CHR}:${EX1_ST}-${EX1_EN}  ->  ${CHR}:${EX2_ST}-${EX2_EN}  ->  ${CHR}:${EX3_ST}-${EX3_EN}"
echo "Intron12=${INTRON12} bp, Intron23=${INTRON23} bp, SEG=${SEG}, PAIRS/sample=${PAIRS}"

mk_sample() {
  local sample="$1"; local sam="${OUT}/${sample}.sam"; local bam="${OUT}/${sample}.bam"
  {
    printf "@HD\tVN:1.6\tSO:coordinate\n"
    printf "@SQ\tSN:%s\tLN:%d\n" "$CHR" "$LN"

    local n12=$((PAIRS/2)); local n23=$((PAIRS - n12))

    gen_pairs() { # $1=prefix $2=POS_R1 $3=INTRON $4=nPairs
      local pref="$1"; local posr1="$2"; local intron="$3"; local count="$4"
      local i shift p1 max_p1 p2 r1_right r2_right t_start t_end tlen qn seq_r1 qual_r1 seq_r2 qual_r2
      for i in $(seq 0 $((count-1))); do
        shift=$(( (i % 31) - 15 ))
        p1=$(( posr1 + shift ))
        if (( posr1 == POS12 )); then
          ((p1<EX1_ST)) && p1=$EX1_ST
          max_p1=$((EX1_EN-SEG+1))
        else
          ((p1<EX2_ST)) && p1=$EX2_ST
          max_p1=$((EX2_EN-SEG+1))
        fi
        ((p1>max_p1)) && p1=$max_p1

        p2=$(( POS_R2 + (i % 7) ))
        (( p2 + LEN_R2 - 1 > EX2_EN )) && p2=$((EX2_EN-LEN_R2+1))
        (( p2 < EX2_ST )) && p2=$EX2_ST

        if (( posr1 == POS12 )); then r1_right=$(( EX2_ST + SEG - 1 )); else r1_right=$(( EX3_ST + SEG - 1 )); fi
        r2_right=$(( p2 + LEN_R2 - 1 ))
        t_start=$(( p1<p2 ? p1 : p2 )); t_end=$(( r1_right>r2_right ? r1_right : r2_right ))
        tlen=$(( t_end - t_start + 1 ))

        qn="${pref}$(printf '%04d' "$i")"
        seq_r1="$(printf 'A%.0s' $(seq 1 $((2*SEG))))"; qual_r1="$(printf 'I%.0s' $(seq 1 $((2*SEG))))"
        seq_r2="$(printf 'T%.0s' $(seq 1 ${LEN_R2}))";  qual_r2="$(printf 'I%.0s' $(seq 1 ${LEN_R2}))"

        # R1 first-in-pair forward (mate reverse): 99 ; R2 second-in-pair reverse: 147
        printf "%s\t99\t%s\t%d\t60\t%uM%uN%uM\t=\t%d\t%d\t%s\t%s\n" \
          "$qn" "$CHR" "$p1" "$SEG" "$intron" "$SEG" "$p2" "$tlen" "$seq_r1" "$qual_r1"
        printf "%s\t147\t%s\t%d\t60\t%uM\t=\t%d\t-%d\t%s\t%s\n" \
          "$qn" "$CHR" "$p2" "$LEN_R2" "$p1" "$tlen" "$seq_r2" "$qual_r2"
      done
    }

    gen_pairs "a" "$POS12" "$INTRON12" "$n12"
    gen_pairs "b" "$POS23" "$INTRON23" "$n23"

    # klein sampleverschil
    local extra=0; [[ "$sample" == "sample1" ]] && extra=3; [[ "$sample" == "sample2" ]] && extra=6
    local k qn p1 p2 r1_right r2_right t_start t_end tlen seq_r1 qual_r1 seq_r2 qual_r2
    for k in $(seq 1 $extra); do
      qn="x$(printf '%03d' "$k")"
      p1=$(( POS12 + 5 + k )); ((p1 > EX1_EN-SEG+1 )) && p1=$((EX1_EN-SEG+1))
      p2=$(( POS_R2 + 10 + k )); (( p2 + LEN_R2 - 1 > EX2_EN )) && p2=$((EX2_EN-LEN_R2+1))
      r1_right=$(( EX2_ST + SEG - 1 )); r2_right=$(( p2 + LEN_R2 - 1 ))
      t_start=$(( p1<p2 ? p1 : p2 )); t_end=$(( r1_right>r2_right ? r1_right : r2_right ))
      tlen=$(( t_end - t_start + 1 ))
      seq_r1="$(printf 'G%.0s' $(seq 1 $((2*SEG))))"; qual_r1="$(printf 'I%.0s' $(seq 1 $((2*SEG))))"
      seq_r2="$(printf 'C%.0s' $(seq 1 ${LEN_R2}))";  qual_r2="$(printf 'I%.0s' $(seq 1 ${LEN_R2}))"
      printf "%s\t99\t%s\t%d\t60\t%uM%uN%uM\t=\t%d\t%d\t%s\t%s\n" \
        "$qn" "$CHR" "$p1" "$SEG" "$INTRON12" "$SEG" "$p2" "$tlen" "$seq_r1" "$qual_r1"
      printf "%s\t147\t%s\t%d\t60\t%uM\t=\t%d\t-%d\t%s\t%s\n" \
        "$qn" "$CHR" "$p2" "$LEN_R2" "$p1" "$tlen" "$seq_r2" "$qual_r2"
    done
  } > "$sam"

  samtools view -bS "$sam" | samtools sort -o "$bam"
  samtools index "$bam"
  rm -f "$sam"
}

mk_sample "sample1"
mk_sample "sample2"
mk_sample "sample3"

# Manifest
{
  echo -e "gene\tchrom\tex1_start\tex1_end\tex2_start\tex2_end\tex3_start\tex3_end\tintron12\tintron23\tseg\tpairs"
  echo -e "${CHOSEN_GENE}\t${CHR}\t${EX1_ST}\t${EX1_EN}\t${EX2_ST}\t${EX2_EN}\t${EX3_ST}\t${EX3_EN}\t${INTRON12}\t${INTRON23}\t${SEG}\t${PAIRS}"
} > "${OUT}/mini_bams_manifest.tsv"

echo "OK: BAMs + BAI in ${OUT}  |  Gen: ${CHOSEN_GENE}"

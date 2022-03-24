# Raw Data

```bash
cd /tmp
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar xzf sratoolkit.current-ubuntu64.tar.gz

mkdir -p ~/tmp
cd ~/tmp
/tmp/sratoolkit.2.10.5-ubuntu64/bin/sratoolkit.2.10.5-ubuntu64/bin/vdb-config --interactive
```

## Week 2 Part 1: Read  Quality Control

```bash

ILLUMINA_PE_SRA_ACCESSIONS=(
  ERR1949188
)

# Download SRA data
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]}; do
  echo ${ACC}

  /tmp/sratoolkit.2.10.5-ubuntu64/bin/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump --split-files --split-e --gzip ${ACC}
done
```

## Week 3 Part 1 and 2: SARS-CoV-2 Resequencing and Short Read Assembly

```bash

# Download SARS-Cov-2 RefSeq
#####
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512&rettype=fasta&retmode=text" \
  | bgzip \
  > COVID-19.fasta.gz

# Info for Subsampling
#####
GENOME_SIZE=29903
ILLUMINA_TARGET_COVERAGES=(
  10
  100
)
NANOPORE_TARGET_COVERAGES=(
  100
)

# SRA Accessions for Download
#####
ILLUMINA_PE_SRA_ACCESSIONS=(
  SRR11140748
)
NANOPORE_SRA_ACCESSIONS=(
  SRR11140749
)

# Download SRA data
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]} ${NANOPORE_SRA_ACCESSIONS[@]}; do
  echo ${ACC}

  /tmp/sratoolkit.2.10.5-ubuntu64/bin/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump --split-files --split-e --gzip ${ACC}
done

# Illumina PE Subsampling
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]}; do
  echo "Processing: ${ACC}"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${ACC}.tsv ]; then
    pigz -dcp1 ${ACC}_?.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${ACC}.tsv
  fi
  INPUT_READS=$(cut -f2 ${ACC}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${ACC}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${ILLUMINA_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      paste \
        <(pigz -dcp1 ${ACC}_1.fastq.gz | paste - - - -) \
        <(pigz -dcp1 ${ACC}_2.fastq.gz | paste - - - -) \
        | shuf -n $((INPUT_READS/2*TARGET_COVERAGE/INPUT_COVERAGE)) \
        | tee \
          >(cut -f 1-4 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz) \
          | cut -f 5-8 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${ACC}_2_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done

# Nanopore Subsampling
#####
for ACC in ${NANOPORE_SRA_ACCESSIONS[@]}; do
  echo "Processing: $ACC"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${ACC}.tsv ]; then
    pigz -dcp1 ${ACC}.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${ACC}.tsv
  fi
  INPUT_READS=$(cut -f2 ${ACC}.tsv | head -n1)
  INPUT_BP=$(cut -f1 /tmp/${ACC}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${NANOPORE_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      filtlong \
        --assembly COVID-19.fasta.gz \
        --target_bases $((GENOME_SIZE*TARGET_COVERAGE)) \
        ${ACC}.fastq.gz \
      | gzip \
      > ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done
```

TODO: Script (`plot_delta.R`) download

## Week 4 Part 1 and 2: Short and Long Read Alignment, Hybrid Genome Assembly

NC_000913.3.fasta.gz
36_ACGCACCT-GGTGAAGG_L002_R?_001_40x.fastq.gz
lima.bc1106--bc1106_40x.subreadset.fastq.gz

```bash
# Info for Subsampling
#####
GENOME_SIZE=4641652
ILLUMINA_TARGET_COVERAGES=(
  2
  4
  5
  8
  10
  20
  40
)
PACBIO_TARGET_COVERAGES=(
  2
  4
  5
  8
  10
  20
  40
)

# Data for download
#####
declare -A ILLUMINA_PE_SRA_ACCESSIONS
ILLUMINA_PE_SRA_ACCESSIONS=(
  [SRR11075452]=36_ACGCACCT-GGTGAAGG_L002
)
PACBIO_PREFIXES=(
  lima.bc1106--bc1106
)


# Illumina data download
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]}; do
  echo ${ACC}

  /tmp/sratoolkit.2.10.5-ubuntu64/bin/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump --split-files --split-e --gzip ${ACC}
done
for ACC in "${!ILLUMINA_PE_SRA_ACCESSIONS[@]}"; do
  mv ${ACC}_1.fastq.gz ${ILLUMINA_PE_SRA_ACCESSIONS[${ACC}]}_R1_001.fastq.gz
  mv ${ACC}_2.fastq.gz ${ILLUMINA_PE_SRA_ACCESSIONS[${ACC}]}_R2_001.fastq.gz
done

# PacBio data download
#####
URLS=(
  https://downloads.pacbcloud.com/public/dataset/MicrobialMultiplexing_48plex/48-plex_sequences/lima.bc1106--bc1106.subreadset.fastq.gz
)

for URL in ${URLS[@]}; do
  ACC=${URL##*/}
  echo ${ACC}

  wget --continue --no-clobber ${URL}
done



# Illumina PE Subsampling
#####
for ACC in ${!ILLUMINA_PE_SRA_ACCESSIONS[@]}; do
  PREFIX=${ILLUMINA_PE_SRA_ACCESSIONS[${ACC}]}
  echo "Processing: ${PREFIX}"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${PREFIX}.tsv ]; then
    pigz -dcp1 ${PREFIX}_R?_001.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${PREFIX}.tsv
  fi
  INPUT_READS=$(cut -f2 ${PREFIX}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${PREFIX}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${ILLUMINA_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${PREFIX}_R1_001_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      paste \
        <(pigz -dcp1 ${PREFIX}_R1_001.fastq.gz | paste - - - -) \
        <(pigz -dcp1 ${PREFIX}_R2_001.fastq.gz | paste - - - -) \
        | shuf -n $((INPUT_READS/2*TARGET_COVERAGE/INPUT_COVERAGE)) \
        | tee \
          >(cut -f 1-4 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${PREFIX}_R1_001_${TARGET_COVERAGE}x.fastq.gz) \
          | cut -f 5-8 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${PREFIX}_R2_001_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done

# PacBio Subsampling
#####
for PREFIX in ${PACBIO_PREFIXES[@]}; do
  echo "Processing: ${PREFIX}"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${PREFIX}.tsv ]; then
    pigz -dcp1 ${PREFIX}.subreadset.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${PREFIX}.tsv
  fi
  INPUT_READS=$(cut -f2 ${PREFIX}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${PREFIX}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${PACBIO_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${PREFIX}_${TARGET_COVERAGE}x.subreadset.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      pigz -dcp1 ${PREFIX}.subreadset.fastq.gz \
        | paste - - - - \
        | shuf -n $((INPUT_READS/2*TARGET_COVERAGE/INPUT_COVERAGE)) \
        | tr '\t' '\n' \
        | pigz --best --processes 1 \
        > ${PREFIX}_${TARGET_COVERAGE}x.subreadset.fastq.gz
      echo "COMPLETE"
    fi
  done
done
```

# Packaging

```bash
rm *.tsv *.fastq.1 36_ACGCACCT-GGTGAAGG_L002_R?_001.fastq.gz SRR11140749.fastq.gz SRR11140748.fastq.gz lima.bc1106--bc1106.subreadset.fastq.gz

tar -c -f genomics_applications.tar ./*.gz
```

# Download from CloudStor

```bash
mkdir --parents ~/data/nathan
cd ~/data/nathan

curl https://cloudstor.aarnet.edu.au/plus/s/K52Ut5IS4gEuC6U/download \
  > genomics_applications.tar
tar xf genomics_applications.tar && rm genomics_applications.tar
```


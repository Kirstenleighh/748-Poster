# 748-Poster

# Metagenomic Analysis Pipeline

All analysis ran on Linux (PuTTY, CGR server) and R

# DOWNLOADING ZOE DATA – fastq files
   
Bash script for downloading Zoe data using PREDICT accession numbers from https://www.nature.com/articles/s41586-025-09854-7#Sec10  

```
nano Download_ena_fastq.sh
```

```
#!/bin/bash

set -euo pipefail
#Run this directly on the login node on the CBF server

#Project name
PROJECTS=("PRJEB75464") # enter all accession numbers that you want to download

BASE_OUTDIR="$HOME/Projects"

for PROJECT in "${PROJECTS[@]}"; do

  echo "=================================="
  echo "Processing project: $PROJECT"
  echo "=================================="


  OUTDIR="${BASE_OUTDIR}/${PROJECT}"
  mkdir -p "$OUTDIR"


  TSV="${OUTDIR}/${PROJECT}_runs.tsv"

  ENA_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,fastq_ftp,sample_accession,scientific_name&format=tsv&download=true&limit=0"


  echo "Downloading the TSV for ${PROJECT}..."
  wget -qO "${TSV}" "${ENA_URL}"

  if ! head -n 1 "${TSV}" | grep -q "run_accession"; then
    echo "ERROR: Something up with the TSV"
    echo "First few lines:"
    head -n 5 "${TSV}"
    echo
    echo "URL used:"
    echo "${ENA_URL}"
    exit 1
  fi

  echo "TSV saved to: ${TSV}"
  echo "Rows of data: $(wc -l < "${TSV}")"

  echo "Downloading FASTQs listed in fastq_ftp..."

awk -F '\t' -v OUTDIR="${OUTDIR}" '
NR==1 { next }  #skip header
{
    #Columns (per fields order):
    #$1 run_accession
    #$2 fastq_ftp (semicolon-separated), important one 
    #$3 sample_accession
    #$4 scientific_name

    sp = $4
    gsub(/ /, "_", sp)

    sample = $3

    n = split($2, urls, ";")
    for (i = 1; i <= n; i++) {
        url = urls[i]
        if (url == "" || url ~ /^[[:space:]]*$/) continue

        #ENA gives ftp:// 
        #Change it to https://
        gsub(/^ftp:\/\//, "", url)
        gsub(/^https:\/\//, "", url)
        url = "https://" url

        dir = OUTDIR "/" sp "/" sample
        system("mkdir -p \"" dir "\"")

        #-c continue; -nv makes it less verbose which is nicer; -P output dir
        cmd = "wget -c -nv -nc -P \"" dir "\" \"" url "\""
        system(cmd)
    }
}
' "${TSV}"

echo "Done. Output in: ${OUTDIR}"
done
```
Control O, 
enter, and  
control X to exit

Make sure script is executable after editing:
```
chmod +x download.script.file.sh
```
Then use screen session so it will still run even if logged out:
```
screen -S ./downloading file name here 
```

# QUALITY CONTROL
FASTQC performed on raw reads (paired) 

```
nano Run_fastqc_200.sh
```

```
#!/bin/bash

BASE_DIR=~/Projects/PRJEB75460/human_gut_metagenome
OUT_DIR=~/Projects/PRJEB75460/fastqc_results

mkdir -p "$OUT_DIR"

ls -d $BASE_DIR/SAMEA* | head -n 200 | while read sample
do
    echo "Processing $sample"
    fastqc -t 8 $sample/*.fastq.gz -o $OUT_DIR
done

echo "FastQC complete"

Multiqc 
```

Control O, 
enter, and  
control X to exit

Make sure script is executable after editing:
```
chmod +x Run_fastqc_200.sh
```
Then use screen session so it will still run even if logged out:
```
screen -S ./fastqc file name here 
```

# TRIMMING 

```
Bash script for trimming and fastqc:

#!/bin/bash

BASE_DIR=~/Projects/PRJEB75460/human_gut_metagenome
TRIM_DIR=~/Projects/PRJEB75460/Trimmed_reads
TRIM_FASTQC_DIR=~/Projects/PRJEB75460/fastqc_trimmed

mkdir -p "$TRIM_DIR" "$TRIM_FASTQC_DIR"

cd "$BASE_DIR"

for DIR in $(ls -d SAMEA* | head -n 200); do

    echo "Processing $DIR"

    for R1 in "$DIR"/*_1.fastq.gz; do
        [ -e "$R1" ] || continue

        R2=${R1/_1.fastq.gz/_2.fastq.gz}

        echo "Trimming $R1 and $R2"

        trim_galore \
            --paired \
            --quality 20 \
            --output_dir "$TRIM_DIR" \
            "$R1" "$R2"

    done

done
```
# FastQC & MultiQC on Trimmed Reads
```
fastqc "$TRIM_DIR"/*.fq.gz -o "$TRIM_FASTQC_DIR"

multiqc	
```


# HOST REMOVAL on trimmed reads
```
#!/bin/bash

# Set human genome reference
HUMAN_REF=/pub54/edwardco/INTEGRATE/pipeline/env/BMTAGGER_INDEX/Human_genome_2023.fa

THREADS=8

TRIMMED_READS_DIR=~/Projects/PRJEB75460/Trimmed_reads
HOST_REMOVED_DIR=~/Projects/PRJEB75460/Host_removed

mkdir -p "$HOST_REMOVED_DIR"

echo "Unzipping all Trimmed FASTQ .gz files..."
pigz -d -k -p "$THREADS" "$TRIMMED_READS_DIR"/*_1_val_1.fq.gz
pigz -d -k -p "$THREADS" "$TRIMMED_READS_DIR"/*_2_val_2.fq.gz

# Process each sample
for fwd in "$TRIMMED_READS_DIR"/*_1_val_1.fq; do
    base=$(basename "$fwd" _1_val_1.fq)
    rev="$TRIMMED_READS_DIR/${base}_2_val_2.fq"

    # Output folder for this sample
    outdir="$HOST_REMOVED_DIR/${base}_host_removed"

    echo "Processing sample $base"
    mkdir -p "$outdir"

 # Run MetaWRAP read_qc for human read removal (short reads)
    metawrap read_qc -1 "$fwd" -2 "$rev" -o "$outdir" -t "$THREADS" -x "$HUMAN_REF"
    echo "Done: $base"
done

echo "All samples processed."
```

Make sure script is executable after editing:
```
chmod +x run_host_removal.sh
nohup ./run_host_removal.sh > host_removal.log 2>&1 &
```
To check progress:
```
tail -f host_removal.log
```


# Alternative script for trimming and host removal: 
```
for R1_GZ in *_1.fastq.gz; do
    R2_GZ="${R1_GZ%_1.fastq.gz}_2.fastq.gz"
    SAMPLE="${R1_GZ%_1.fastq.gz}"

    TMPDIR=$(mktemp -d)
    R1_FIFO="$TMPDIR/${SAMPLE}_1.fastq"
    R2_FIFO="$TMPDIR/${SAMPLE}_2.fastq"
    mkfifo "$R1_FIFO" "$R2_FIFO"

    pigz -dc "$R1_GZ" > "$R1_FIFO" &
    PID1=$!
    pigz -dc "$R2_GZ" > "$R2_FIFO" &
    PID2=$!

    metawrap read_qc \
      -1 "$R1_FIFO" \
      -2 "$R2_FIFO" \
      -o "${SAMPLE}_HUMAN_REMOVED" \
      -t 16 \
      -x Human_genome_2023

    wait "$PID1" "$PID2"
    rm -rf "$TMPDIR"
done
```

# KRAKEN2 on trimmed reads
```
#!/bin/bash

COUNT=0
mkdir -p kraken2_results_200

for R1 in Trimmed_reads/*_1_val_1.fq.gz
do
    [ -e "$R1" ] || continue
    [ "$COUNT" -ge 200 ] && break
    R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"

    SAMPLE=$(basename "${R1/_1_val_1.fq.gz/}")

    echo "Running Kraken2 on $R1 and $R2"

    kraken2 \
        --paired \
        --db ~/Kraken2-nt_core-251015 \
        --threads 32 \
        --output "kraken2_results_200/${SAMPLE}.kraken" \
        --report "kraken2_results_200/${SAMPLE}.kreport2" \
        "$R1" "$R2"

    COUNT=$((COUNT+1))

    #if [ $COUNT -eq 200 ]; then
     #   break
    #fi

done
```

```
bash run_kraken2_200.sh
```
```
screen -x
```
to leave the screen
crtl a, then d 

```
nohup bash run_kraken2_200.sh > kraken2_200.log 2>&1 &
```

# Krona
```
conda install -c bioconda krona
ktImportTaxonomy *.kraken -o krona.html
```

```
conda install -c bioconda krona
ktUpdateTaxonomy.sh
ktImportTaxonomy -o ~/Projects/PRJEB75460/krona_output.html ~/Projects/PRJEB75460/Bracken_results_200/*.breport2
```

# Pavian plot r code 
```
if (!require(remotes)) install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")

options(shiny.maxRequestSize = 500*1024^2)

pavian::runApp(
  port = 5000,
  launch.browser = TRUE,
  reportDir = "C:/Users/Kirsten/OneDrive/Desktop/MSc_Bioinformatics/MSc Research Project/200 .kreport2s for pavian"
)
```

# Bracken
```
Conda create -n bracken -c bioconda bracken
Conda activate bracken 
```
```
#!/bin/bash

KREPORT_DIR=~/Projects/PRJEB75460/kraken2_results_200
OUTPUT_DIR=~/Projects/PRJEB75460/Bracken_results_200

mkdir -p "$OUTPUT_DIR"

for kreport in "$KREPORT_DIR"/*.kreport2; do
    SAMPLE=$(basename "$kreport" .kreport2)
    bracken \
        -d ~/Kraken2-nt_core-251015 \
        -i "$kreport" \
        -o "${OUTPUT_DIR}/${SAMPLE}.bracken" \
        -w "${OUTPUT_DIR}/${SAMPLE}.breport2" \
        -r 151 \
        -l S \
        -t 5
done
```
```
nohup bash ~/Projects/PRJEB75460/run_bracken_200.sh > ~/Projects/PRJEB75460/bracken.log 2>&1 &
```
```
tail -f ~/Projects/PRJEB75460/bracken.log
```



# Code for subsetting samples for each bacteria into directories 
```
#!/bin/bash

REPORT_DIR="/pub60/kirsten/Projects/PRJEB75460/kraken2_results_200"
OUTPUT_DIR="/pub60/kirsten/Projects/PRJEB75460/pavian_subsets"

mkdir -p "$OUTPUT_DIR"/{ecoli,salmonella,campylobacter,cdiff}

for report in "$REPORT_DIR"/*.kreport2; do
    awk -F'\t' '$3 > 0 && tolower($6) ~ /escherichia coli/'       "$report" | grep -q . && cp "$report" "$OUTPUT_DIR/ecoli/"
    awk -F'\t' '$3 > 0 && tolower($6) ~ /salmonella/'             "$report" | grep -q . && cp "$report" "$OUTPUT_DIR/salmonella/"
    awk -F'\t' '$3 > 0 && tolower($6) ~ /campylobacter/'          "$report" | grep -q . && cp "$report" "$OUTPUT_DIR/campylobacter/"
    awk -F'\t' '$3 > 0 && tolower($6) ~ /clostridioides difficile|clostridium difficile/' "$report" | grep -q . && cp "$report" "$OUTPUT_DIR/cdiff/"
done

for org in ecoli salmonella campylobacter cdiff; do
    echo "$org: $(ls $OUTPUT_DIR/$org | wc -l) samples"
done
```

# Alpha diversity
Shannon index – compare across samples 
The more microbes in a sample there are, the higher the value for the Shannon index is; and the less the inequality of relative abundances is, the higher the Shannon index is too.

Download all bracken reports (.breport2)

R code:
```
# Checking all 200 .breport2 files are in the wd
files <- list.files(pattern="*.breport2")
print(files)
length(files)

# Read in .kreport2 files 
library(dplyr)
library(tidyr)

files <- list.files(pattern="*.breport2")

read_bracken <- function(file) {
  df <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
  colnames(df) <- c("percent", "reads", "assigned", "rank", "taxid", "name")
  
  # keep only species level
  df <- df %>% filter(rank == "S")
  
  df <- df %>% select(name, reads)
  
  df$Sample <- file
  return(df)
}

list_data <- lapply(files, read_bracken)
combined <- bind_rows(list_data)


# Create a feature table
feature_table <- combined %>%
  pivot_wider(names_from = name, values_from = reads, values_fill = 0)

# set rownames
rownames(feature_table) <- feature_table$Sample
feature_table$Sample <- NULL


# Alpha diversity 
#install.packages("vegan")
library(vegan)

shannon <- diversity(feature_table, index = "shannon")

# Plot 
library(ggplot2)

df <- data.frame(
  Sample = rownames(feature_table),
  Shannon = shannon
)

ggplot(df, aes(x = 1, y = Shannon)) +
  geom_violin(fill = "palevioletred", alpha = 0.4, trim = FALSE) +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "Alpha Diversity (Shannon)", x = "", y = "Shannon Index")
```

# Stacked bar chart R code
```
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Read files
files <- list.files(pattern = "\\.breport2$")

read_bracken <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "",
                   fill = TRUE, flush = TRUE)
  df <- df[rowSums(is.na(df)) == 0, ]
  df <- df[, 1:6]
  colnames(df) <- c("percent", "reads_clade", "reads_direct", "rank", "taxid", "name")
  df$reads_clade <- as.numeric(df$reads_clade)
  df$name        <- trimws(df$name)   # remove any leading whitespace from names
  df$Sample      <- gsub("\\.breport2$", "", basename(file))
  df
}

combined_all <- bind_rows(lapply(files, read_bracken))

# Quick sanity check
cat("Rows:", nrow(combined_all), "\n")
cat("Samples:", n_distinct(combined_all$Sample), "\n")
cat("E. coli rows:", sum(combined_all$name == "Escherichia coli"), "\n")

# 2. Assign taxa
target_data <- combined_all %>%
  filter(rank == "S") %>%
  filter(!grepl("phage|virus|uncultured|unclassified|bacterium", name, ignore.case = TRUE)) %>%
  mutate(Taxon = case_when(
    name == "Escherichia coli"                               ~ "Escherichia coli",
    grepl("^Campylobacter ", name)                           ~ "Campylobacter",
    grepl("^Cryptosporidium ", name)                         ~ "Cryptosporidium",
    grepl("^Yersinia ", name)                                ~ "Yersinia",
    name %in% c("Salmonella enterica", "Salmonella bongori") ~ "Salmonella",
    grepl("^Vibrio ", name)                                  ~ "Vibrio",
    grepl("^Shigella ", name)                                ~ "Shigella",
    name == "Clostridioides difficile"                       ~ "Clostridioides difficile",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Taxon))

cat("Target rows found:", nrow(target_data), "\n")

# 3. Summarise
plot_data <- target_data %>%
  group_by(Sample, Taxon) %>%
  summarise(reads = sum(reads_clade), .groups = "drop") %>%
  complete(Sample, Taxon, fill = list(reads = 0)) %>%
  group_by(Sample) %>%
  mutate(rel_abund = reads / sum(reads) * 100) %>%
  ungroup()

cat("Samples with detections:", n_distinct(plot_data$Sample[plot_data$reads > 0]), "\n")
print(plot_data %>% group_by(Taxon) %>% summarise(total = sum(reads)) %>% arrange(desc(total)))

# 4. Plot
sample_order <- plot_data %>%
  group_by(Sample) %>%
  summarise(total = sum(reads)) %>%
  arrange(desc(total)) %>%
  pull(Sample)

plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)


pal <- c(
  "Escherichia coli" = "yellowgreen",
  "Campylobacter"    = "brown4",
  "Cryptosporidium"  = "plum1",
  "Yersinia"         = "indianred3",
  "Salmonella"       = "darkslateblue",
  "Vibrio"           = "orange",
  "Shigella"         = "lightblue",
  "Clostridioides difficile" = "lightpink"
)

tiff("test.tiff", units="in", width=5, height=5, res=300)
p <- ggplot(plot_data, aes(x = Sample, y = rel_abund, fill = Taxon)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal, name = "Pathogen") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = "right",
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Sample", y = "Relative Abundance (%)")
dev.off()

print(p)

library(tidyr)
genus_table <- plot_data %>%
  group_by(Taxon) %>%
  summarise(
    Total_Reads        = sum(reads),
    Samples_Detected   = sum(reads > 0),
    Mean_RelAbund      = round(mean(rel_abund), 4),
    SD_RelAbund        = round(sd(rel_abund), 4),
    Min_RelAbund       = round(min(rel_abund), 4),
    Max_RelAbund       = round(max(rel_abund), 4)
  ) %>%
  arrange(desc(Mean_RelAbund))

print(genus_table)



simple_table <- plot_data %>%
  group_by(Taxon) %>%
  summarise(`Mean % Abundance` = round(mean(rel_abund), 4)) %>%
  arrange(desc(`Mean % Abundance`))

print(simple_table)
```

# Useful commands

To remove a directory
```
rm -rf (name of dir)
```

To check how many of a specific file you have
```
 ls ~/Projects/PRJEB75460/Bracken_results_200/*.bracken | wc -l
```
(change the file ending to search)

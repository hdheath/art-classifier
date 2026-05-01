# Building a Splice Junction Database for flAIr

## TL;DR - Recommended Approach

**You already have excellent junction data!** Your RJunBase file is comprehensive and ready to use.

**Quick Start:**
```bash
# Use your existing RJunBase data (recommended)
python src/extract/standardize_priors.py \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/
```

This will give you ~800,000 high-quality junctions ready for CNN training.

---

## Option 1: Use Your Existing RJunBase Data (RECOMMENDED) ✓

### Why RJunBase is Excellent:

Your `detail_Linear_splice_annotation.txt` contains:
- **~800,000 linear splice junctions**
- Multi-tissue support (Tumor + Normal samples)
- Expression metrics (median TPM)
- Gene annotations
- Already filtered for quality

### What You Have:
```bash
$ head splice_db/detail_Linear_splice_annotation.txt
JunctionID  Junction location  Gene symbol  ...  Tumor median  Normal median
UN_SMIM8_LS0017  chr6:87344607|87345167:+  SMIM8  ...  2.190  0.000
RNF213_LS024  chr17:80328477|80332006:+  RNF213  ...  9.189  5.687
```

### Processing:
The `standardize_priors.py` script already handles this format:
1. Parses `chr:donor|acceptor:strand` format
2. Extracts expression weights (Normal median preferred)
3. Removes duplicates
4. Sorts and indexes for fast lookup

**This is production-ready and requires no additional downloads!**

---

## Option 2: GTEx Direct Download (Alternative)

If you want additional GTEx-specific data:

### Download GTEx v8 Junctions

GTEx provides comprehensive junction data across all tissues.

```bash
# Create download directory
mkdir -p data/downloads/gtex

cd data/downloads/gtex

# Download GTEx junction file (WARNING: Large file ~3GB)
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/\
GTEx_Analysis_2017-06-05_v8_STARjunctions.gct.gz

# Or use smaller per-tissue files
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/\
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
```

### Parse GTEx Junction Format

GTEx uses GCT format (Gene Cluster Text):
```
#1.2
[dimensions]
junction_id  sample1  sample2  sample3  ...
chr1:12345:12567:1  0  5  2  ...
chr1:23456:24789:2  10  8  12  ...
```

Format: `chr:start:end:strand` where strand is 1(+) or 2(-)

---

## Option 3: Use recount3 (Advanced)

### What is recount3?

[recount3](https://rna.recount.bio/) provides uniformly processed RNA-seq data:
- **GTEx**: All tissues (~17,000 samples)
- **TCGA**: Cancer datasets (~11,000 samples)
- **SRA**: Public data

### Access Methods:

#### Method A: Use R package (Recommended)

```r
# Install recount3
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("recount3")

library(recount3)

# Get available projects
human_projects <- available_projects()

# Get GTEx data
gtex_info <- subset(human_projects,
                    file_source == "gtex" &
                    organism == "human")

# Download junctions for a specific tissue
rse_gtex_brain <- create_rse(
    project = "BRAIN",
    project_home = "data_sources/gtex"
)

# Extract junction counts
jxn <- assays(rse_gtex_brain)$counts

# Export to BED format
# ... (R code to parse and export)
```

#### Method B: Direct API Access

```bash
# Get project list
curl https://recount.bio/api/v1/projects > projects.json

# Download specific project junctions
# (API endpoints may vary - check documentation)
```

#### Method C: Pre-computed Files

Check if recount3 provides pre-computed junction files:
```bash
# Example (verify URL)
wget https://recount.bio/data/gtex/junctions.bed.gz
```

---

## Option 4: Extract from GENCODE Annotations

For **annotated junctions only** (high confidence, but incomplete):

```bash
cd data/ref

# Download GENCODE GTF (if not already done)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/\
release_46/gencode.v46.annotation.gtf.gz

# Extract junctions from GTF
zcat gencode.v46.annotation.gtf.gz | \
  awk '$3=="exon"' | \
  # ... parse exon boundaries to get junctions
  # (requires custom script)
```

This gives you ~300,000-400,000 **annotated** junctions.

---

## Comparison of Approaches

| Source | Junction Count | Quality | Coverage | Effort |
|--------|---------------|---------|----------|--------|
| **RJunBase** (yours) | ~800K | High | Multi-tissue | ✓ Ready |
| GTEx Direct | ~500K | High | Tissue-specific | Medium |
| recount3 | ~1M+ | High | Comprehensive | High |
| GENCODE | ~350K | Very High | Annotated only | Low |

---

## Recommended Workflow

### For Most Users (Including You):

**Use RJunBase + GENCODE combo**

```bash
# 1. Process your RJunBase data (primary source)
python src/extract/standardize_priors.py \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/ \
    --prefix rjunbase_

# 2. Extract GENCODE annotated junctions (optional, for validation)
python src/extract/extract_gencode_junctions.py \
    --gtf data/ref/gencode.v46.annotation.gtf.gz \
    --output data/ref/gencode_junctions.tsv

# 3. Merge (optional - only if you want both)
python src/extract/merge_junctions.py \
    --inputs data/ref/rjunbase_junction_prior.tsv.gz \
             data/ref/gencode_junctions.tsv \
    --output data/ref/junction_prior_merged.tsv
```

---

## CNN Training Considerations

### Junction Selection for Training:

**Positives** (real junctions):
- Use junctions with support from multiple samples (min_samples >= 2)
- Weight by expression level (log-transformed read counts)
- Include both annotated and novel junctions

**Negatives** (non-junctions):
- Random genomic positions (>1kb from real junctions)
- "Near-miss" positions (±3-20bp from real junctions)
- Positions with non-canonical motifs

### Expected Counts:

From RJunBase (~800K junctions):
- **Positives for CNN**: ~800K junction sites × 2 (donor + acceptor) = 1.6M training examples
- **Negatives**: Generate 1-5x as many (1.6M - 8M)
- **Total training set**: 3M - 10M examples

### Motif Validation:

Before using junctions for training, validate:
```bash
# Check GT-AG enrichment
python src/extract/validate_junction_motifs.py \
    --junctions data/ref/junction_prior.tsv.gz \
    --genome data/ref/GRCh38.primary_assembly.genome.fa \
    --output qc/junction_motifs.html
```

Expected: >95% of junctions should have canonical splice sites (GT-AG, GC-AG, AT-AC)

---

## Next Steps After Building Database

1. **Validate genome build** (hg38 vs hg19):
   ```bash
   # Compare junction coordinates to GENCODE
   python src/extract/validate_genome_build.py \
       --junctions data/ref/junction_prior.tsv.gz \
       --gencode data/ref/gencode.v46.annotation.gtf.gz
   ```

2. **Generate QC report**:
   ```bash
   python src/extract/junction_qc.py \
       --junctions data/ref/junction_prior.tsv.gz \
       --output qc/junction_summary.html
   ```

3. **Build CNN training windows** (Phase 3):
   ```bash
   python src/cnn/build_training_windows.py \
       --junction-prior data/ref/junction_prior.tsv.gz \
       --genome data/ref/GRCh38.primary_assembly.genome.fa \
       --output data/derived/cnn_train/
   ```

---

## FAQ

### Q: Should I use recount3 or stick with RJunBase?

**A:** RJunBase is excellent and sufficient. Only use recount3 if you need:
- Tissue-specific junction sets
- Cancer-specific junctions (TCGA)
- Extremely rare junctions from SRA

### Q: How many junctions do I need for CNN training?

**A:**
- **Minimum**: 100K junctions (from GENCODE)
- **Good**: 500K junctions (GTEx or RJunBase)
- **Excellent**: 800K+ junctions (your RJunBase data)

More junctions = better CNN performance, especially for novel junction discovery.

### Q: Should I filter for canonical splice sites?

**A:**
- For **positives**: Keep all junctions (CNN learns from real data)
- For **validation**: Check that >95% are canonical
- For **negatives**: Include non-canonical positions

### Q: Genome build verification?

Check coordinates against GENCODE:
```bash
# Get RNF213 gene from GENCODE
zcat data/ref/gencode.v46.annotation.gtf.gz | \
  grep 'gene_name "RNF213"' | head -1

# Should match your junction:
# chr17:80328477|80332006:+ (from RJunBase)
```

If coordinates match GENCODE v46 → it's hg38 ✓

---

## Summary

**For your project:**

1. ✓ **Use RJunBase** (you already have it, it's comprehensive)
2. ✓ **Run standardize_priors.py** (ready to use)
3. ✓ **Validate with GENCODE** (verify genome build)
4. → **Proceed to Phase 3** (CNN training data)

You don't need recount3 unless you have specific requirements for tissue-specific or cancer-specific junctions.

**Your existing data is production-ready!**

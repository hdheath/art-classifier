# Data Validation Report
## flAIr Project Database Assessment

**Date:** 2025-12-26
**Location:** `/private/groups/brookslab/hdheath/projects/flAIr/`

---

## Executive Summary

You have assembled comprehensive prior databases for TSS, TES, and splice junctions. However, **several critical data preparation steps are required** before these can be used for CNN training:

### Critical Issues:
1. **Genome Build Mismatch**: FANTOM CAGE data is hg19, not hg38
2. **Format Heterogeneity**: Multiple file formats need standardization
3. **Missing Tabix Indices**: Most files lack `.tbi` indices for fast querying
4. **Coordinate Extraction Needed**: FANTOM files contain TPM data, need BED coordinates extracted

---

## TSS Database Assessment

### Files Present:
1. **FANTOM_TSS_human.bed** (97MB)
   - Format: BED9 with track header
   - Genome: hg38 (appears correct based on chr10 coordinates)
   - Contains: TSS predictions with confidence scores in name field
   - Status: ✓ Good format, needs tabix indexing

2. **hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz** (339MB)
   - Format: Tab-separated with expression matrix (1800+ samples)
   - Genome: **hg19** (requires liftover!)
   - Contains: CAGE peak IDs in first column (e.g., `chr1:564599:+`)
   - Status: ⚠️ **Needs coordinate extraction and hg19→hg38 liftover**

3. **refTSS_v4.1_human_coordinate.hg38.bed.txt** (9.6MB)
   - Format: Tab-separated BED6 with header
   - Genome: hg38 ✓
   - Contains: ~60,000 refTSS regions
   - Status: ✓ Good, needs header removal and bgzip/tabix

### Recommendations:
- **Use:** refTSS (hg38, ready) + FANTOM BED (after indexing)
- **Transform:** Extract coordinates from FANTOM CAGE file and liftover to hg38
- **Expected TSS count:** ~100,000-200,000 unique peaks after merging

---

## TES/PolyA Database Assessment

### Files Present:
1. **polyadb.hg38.weighted.bed** (14MB)
   - Format: BED6, already lifted to hg38 ✓
   - Contains: Weighted polyA sites (weight in score column)
   - Status: ✓ **Ready to use**, just needs bgzip/tabix

2. **polyadb.hg38.bed** (15MB)
   - Format: BED6, unweighted version
   - Status: ✓ Alternative to weighted version

3. **polyAdb_human.PAS.txt** (43MB)
   - Format: Original polyA_DB format (likely hg19)
   - Status: Already converted above

### Liftover Evidence:
- You have `hg19ToHg38.over.chain` and `liftOver` binary present
- `polyadb.unmapped.bed` (5.9KB) shows minimal liftover failures ✓

### Recommendations:
- **Use:** `polyadb.hg38.weighted.bed` (preferred for CNN training weights)
- **Expected TES count:** ~50,000-100,000 sites

---

## Splice Junction Database Assessment

### Files Present:
1. **detail_Linear_splice_annotation.txt** (137MB)
   - Format: Tab-separated with header
   - Contains: Linear splice junctions with rich annotations
   - Columns: JunctionID, Junction location (chr:donor|acceptor:strand), Gene symbol, Expression metrics
   - Example: `chr6:87344607|87345167:+`
   - Genome: **Likely hg38** (based on RNF213 coordinates match GENCODE)
   - Status: ⚠️ **Needs coordinate parsing**

2. **Alternative_Splice_Junctions.csv** (570MB)
   - Format: CSV, very large
   - Status: Alternative source if needed

3. **detail_Bakc_splice_annotation.txt** (29MB)
   - Format: Back-splice (circular RNA) junctions
   - Status: Not needed for this project

4. **detail_Fusion_splice_annotation.txt** (6.1MB)
   - Format: Fusion gene junctions
   - Status: Not needed for this project

### Coordinate Parsing Required:
The linear splice file has coordinates in format:
`chr6:87344607|87345167:+`

Needs conversion to:
```
chr  donor       acceptor    strand  weight
chr6 87344607    87345167    +       2.190
```

### Recommendations:
- **Use:** `detail_Linear_splice_annotation.txt`
- **Parse:** Extract chr, donor, acceptor, strand from "Junction location" column
- **Weight:** Use "Tumor median" or "Normal median" TPM as weight
- **Expected junction count:** ~500,000-1,000,000 junctions

---

## Data Preparation Pipeline Required

### Phase 1: TSS Prior Standardization

```bash
# 1. FANTOM BED (already hg38)
grep -v "^track" FANTOM_TSS_human.bed | \
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' | \
  sort -k1,1 -k2,2n | \
  bgzip > tss_fantom.bed.gz
tabix -p bed tss_fantom.bed.gz

# 2. refTSS
tail -n +2 refTSS_v4.1_human_coordinate.hg38.bed.txt | \
  sort -k1,1 -k2,2n | \
  bgzip > tss_reftss.bed.gz
tabix -p bed tss_reftss.bed.gz

# 3. FANTOM CAGE peaks (extract coordinates, already hg19→needsliftover)
# Will need custom script to parse annotation column
```

### Phase 2: TES Prior Standardization

```bash
# Already in good shape!
sort -k1,1 -k2,2n polyadb.hg38.weighted.bed | \
  bgzip > tes_prior.bed.gz
tabix -p bed tes_prior.bed.gz
```

### Phase 3: Splice Junction Standardization

```python
# Parse RJunBase format
import pandas as pd

df = pd.read_csv('detail_Linear_splice_annotation.txt', sep='\t')

# Extract coordinates
# Format: chr6:87344607|87345167:+
def parse_junction(junc_str):
    chrom, coords = junc_str.rsplit(':', 1)
    strand = coords[-1]
    donor, acceptor = coords[:-2].split('|')
    return chrom, int(donor), int(acceptor), strand

junctions = df['Junction location'].apply(parse_junction)
# Use Tumor median or Normal median as weight
```

---

## Genome Build Verification Checklist

- [ ] Verify FANTOM BED coordinates match hg38 GENCODE
- [ ] Verify RJunBase coordinates match hg38 GENCODE
- [ ] Download GENCODE v46 GTF for hg38
- [ ] Download GRCh38 primary assembly FASTA

---

## Expected File Outputs

After standardization, you should have:

```
data/ref/
├── genome.fa                          # GRCh38 primary assembly
├── genome.fa.fai                      # samtools faidx index
├── gencode.v46.annotation.gtf.gz      # Gene annotations
├── genes.bed.gz                       # Per-gene boundaries
├── genes.bed.gz.tbi
├── tss_prior.bed.gz                   # Merged TSS peaks
├── tss_prior.bed.gz.tbi
├── tes_prior.bed.gz                   # PolyA sites (weighted)
├── tes_prior.bed.gz.tbi
├── junction_prior.tsv.gz              # Splice junctions
└── junction_prior.tsv.gz.tbi
```

---

## Data Statistics (Estimated)

| Resource | Source | Expected Count | Genome Build |
|----------|--------|---------------|--------------|
| TSS peaks | FANTOM + refTSS | 150,000 | hg38 |
| TES sites | PolyA_DB | 70,000 | hg38 |
| Splice junctions | RJunBase | 800,000 | hg38 (verify) |

---

## Next Steps

1. **Immediate:**
   - Verify genome build for RJunBase (check against GENCODE coordinates)
   - Download GRCh38 reference genome
   - Download GENCODE v46 GTF

2. **Data Processing:**
   - Write standardization scripts (see Phase 1 deliverable)
   - Run coordinate parsing for splice junctions
   - Bgzip and tabix all files

3. **Validation:**
   - Check overlap with GENCODE annotations
   - Verify chromosome naming consistency
   - Run QC metrics (coverage, motif enrichment)

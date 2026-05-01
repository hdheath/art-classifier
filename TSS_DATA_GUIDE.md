# TSS Data Preparation Guide

## What You Have

```
TSS_db/
├── FANTOM_TSS_human.bed (97MB) - Ready to use ✓
├── refTSS_v4.1_human_coordinate.hg38.bed.txt (9.6MB) - Ready to use ✓
└── hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz (339MB) - Needs extraction
```

---

## Quick Decision: Do You Need the CAGE File?

### ✓ **NO - Use simpler files (RECOMMENDED)**

**Why:**
- You already have 2 excellent TSS files (FANTOM BED + refTSS)
- Both are hg38 and ready to use
- Combined, they give you ~100K-150K TSS peaks
- **This is MORE than sufficient for CNN training**

**Quick command:**
```bash
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/
```

**Time:** 5-10 minutes
**Output:** ~150K TSS peaks, ready for CNN training

---

### ✗ **YES - Extract CAGE file (ADVANCED)**

**Only if you need:**
- Maximum TSS coverage (200K+ peaks)
- Cell-type-specific TSS usage
- Expression-weighted TSS peaks

**Problem: It's hg19!**

The CAGE file is in **hg19 coordinates**, which needs liftover to hg38.

**Steps:**

#### 1. Extract TSS coordinates from CAGE file:

```bash
python src/extract/parse_fantom_cage.py \
    --input TSS_db/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz \
    --output data/ref/fantom_cage_tss_hg19.bed \
    --genome-build hg19 \
    --aggregation mean \
    --min-tpm 0.5
```

**Time:** 10-15 minutes (large file!)
**Output:** ~200K TSS peaks in hg19

#### 2. Liftover hg19 → hg38:

**Option A: Using liftOver binary (you have it!)**

```bash
# You already have liftOver in TTS_db/
TTS_db/liftOver \
    data/ref/fantom_cage_tss_hg19.bed \
    TTS_db/hg19ToHg38.over.chain \
    data/ref/fantom_cage_tss_hg38.bed \
    data/ref/fantom_cage_tss_unmapped.bed
```

**Option B: Run parse script with liftover:**

```bash
python src/extract/parse_fantom_cage.py \
    --input TSS_db/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz \
    --output data/ref/fantom_cage_tss.bed \
    --genome-build hg19 \
    --liftover TTS_db/hg19ToHg38.over.chain \
    --aggregation mean \
    --min-tpm 0.5
```

**Expected liftover success:** ~98%
**Output:** ~196K TSS peaks in hg38

#### 3. Merge all TSS files:

```bash
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tss-cage data/ref/fantom_cage_tss.bed \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/
```

**Output:** ~200K-250K unique TSS peaks (after merging nearby peaks)

---

## File Format Breakdown

### 1. FANTOM_TSS_human.bed (Simple, hg38 ✓)

```
chr10  100008587  100008589  p1@CU680531,0.1352  -89  +
chr10  100015362  100015397  p2@LOXL4,0.1291     55   -
```

**Format:** BED9
**Genome:** hg38 ✓
**TSS count:** ~100K
**Quality:** High (DPI-predicted TSS)
**Status:** Ready to use!

---

### 2. refTSS_v4.1_human_coordinate.hg38.bed.txt (Simple, hg38 ✓)

```
chr1  16013  16020  rfhg_1.1  1  -
chr1  29346  29366  rfhg_2.1  1  -
```

**Format:** TSV with header (chr, start, end, refTSS_ID, score, strand)
**Genome:** hg38 ✓
**TSS count:** ~60K
**Quality:** Very high (curated integrative resource)
**Status:** Ready to use!

---

### 3. hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz (Complex, hg19 ✗)

```
00Annotation  ...  tpm.sample1  tpm.sample2  ...
chr10:100013403..100013414,-  ...  2.5  3.1  ...
chr10:100027943..100027958,-  ...  5.2  4.8  ...
```

**Format:** Custom FANTOM format with 1,800+ sample columns
**Genome:** hg19 ✗ (needs liftover!)
**TSS count:** ~200K
**Quality:** Excellent (CAGE-seq validated)
**Status:** Needs extraction + liftover

**Why it's complex:**
- First column has coordinates (`chr:start..end,strand`)
- Next 6 columns are annotations
- Remaining 1,800+ columns are TPM values per sample
- Each row is a CAGE peak (TSS cluster)
- Values are RLE-normalized TPM

---

## Comparison

| File | TSS Count | Genome | Effort | Recommendation |
|------|-----------|--------|--------|----------------|
| FANTOM BED | ~100K | hg38 ✓ | None | ✓ Use |
| refTSS | ~60K | hg38 ✓ | None | ✓ Use |
| CAGE file | ~200K | hg19 ✗ | High | Optional |
| **Combined (simple)** | **~150K** | **hg38 ✓** | **Low** | **✓ Recommended** |
| **Combined (all)** | **~250K** | **hg38** | **High** | Optional |

---

## CNN Training Requirements

### Minimum for good CNN performance:
- **50K TSS peaks** (minimum)
- **100K TSS peaks** (good)
- **150K TSS peaks** (excellent) ← **You have this with simple files!**
- **200K+ TSS peaks** (overkill, diminishing returns)

### Recommendation:

**For most users (including you):** Use FANTOM BED + refTSS
- Quick to process (5 mins)
- Sufficient coverage (150K peaks)
- High quality (both curated)
- Already hg38 (no liftover needed)

**For maximum coverage:** Add CAGE file
- Extra effort (30 mins)
- Slightly more peaks (+50K)
- Requires liftover
- Marginal improvement for CNN

---

## Commands Summary

### Option 1: Simple (RECOMMENDED) ✓

```bash
# Just use the two simple BED files
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/

# Output: data/ref/tss_prior.bed.gz (~150K peaks)
# Time: 5-10 minutes
# Ready for Phase 3: CNN training
```

### Option 2: Maximum coverage (ADVANCED)

```bash
# Step 1: Extract and liftover CAGE file
python src/extract/parse_fantom_cage.py \
    --input TSS_db/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz \
    --output data/ref/fantom_cage_tss.bed \
    --liftover TTS_db/hg19ToHg38.over.chain \
    --min-tpm 0.5

# Step 2: Merge all TSS sources
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tss-cage data/ref/fantom_cage_tss.bed \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/

# Output: data/ref/tss_prior.bed.gz (~250K peaks)
# Time: 30-40 minutes
# Marginal improvement for CNN
```

---

## My Recommendation

**Skip the CAGE file** for now. Use the simple approach:

1. ✓ Fast (5-10 mins vs 30-40 mins)
2. ✓ Sufficient quality (150K peaks is excellent)
3. ✓ No liftover headaches
4. ✓ Can always add CAGE file later if needed

**When would you add CAGE file?**
- If CNN performance is poor on TSS task (unlikely)
- If you need tissue-specific TSS analysis
- If you want maximum possible coverage

**For initial development:** Simple approach is perfect!

---

## Next Steps (Week 1)

```bash
# 1. Process your existing simple files (5-10 mins)
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/

# 2. Verify outputs
ls -lh data/ref/*.gz data/ref/*.tbi

# 3. Test tabix indexing
tabix data/ref/tss_prior.bed.gz chr1:1000000-2000000
tabix data/ref/tes_prior.bed.gz chr1:1000000-2000000
tabix data/ref/junction_prior.tsv.gz chr1:1000000-2000000

# 4. Move to Phase 2: Download GRCh38 genome and GENCODE
```

**You're ready to start Phase 1!**

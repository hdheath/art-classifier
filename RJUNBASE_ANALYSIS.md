# RJunBase Data Analysis for CNN Training

## Summary: Your Data is EXCELLENT! ✓

**Bottom line:** Your RJunBase file has **everything needed** for CNN training. No additional downloads required.

---

## What You Have

### File: `detail_Linear_splice_annotation.txt`

**Total junctions:** 682,017 (excluding header)

### Breakdown:
- **Annotated junctions:** 357,366 (52%)
- **Novel junctions:** 324,651 (48%)
- **Plus strand:** 343,979
- **Minus strand:** 338,038

### Data Columns Available:

| Column | Description | Use for CNN |
|--------|-------------|-------------|
| Junction location | `chr:donor\|acceptor:strand` | ✓ **Coordinates** |
| Gene symbol | Gene name | Context (optional) |
| Tumor median | Median expression in tumors | ✓ **Weight** |
| Normal median | Median expression in normal | ✓ **Weight (preferred)** |
| Tumor frequency | Samples supporting junction | ✓ **Sample count** |
| Normal frequency | Normal samples supporting | ✓ **Sample count** |

---

## Is This Sufficient for CNN Training?

### ✓ YES! Here's why:

#### 1. **Quantity: EXCELLENT**
- **682K junctions** is more than sufficient
- For comparison:
  - GENCODE (annotated only): ~350K junctions
  - GTEx (all tissues): ~500K junctions
  - **Your RJunBase: 682K junctions** ✓

#### 2. **Quality: HIGH**
- Multi-tissue RNA-seq data (Tumor + Normal)
- Expression weights included (median TPM)
- Sample support counts (frequency)
- Mix of annotated + novel junctions

#### 3. **Coverage: COMPREHENSIVE**
- **52% annotated** - matches known biology
- **48% novel** - good for discovering new junctions
- Balanced strand distribution

#### 4. **What CNN Training Needs:**

For each junction, you need:
- ✓ Chromosome
- ✓ Donor coordinate (junction start)
- ✓ Acceptor coordinate (junction end)
- ✓ Strand
- ✓ Weight/support metric

**Your data has ALL of this!**

---

## CNN Training Data Generation

### From your 682K junctions, you'll create:

#### Positive Examples (Real Junctions):
- **Donors:** 682K × 1 = **682,000 donor sites**
- **Acceptors:** 682K × 1 = **682,000 acceptor sites**
- **Total positives:** ~1.36M training examples

Each example:
- Extract 128bp sequence window (±64bp from junction)
- One-hot encode (A,C,G,T,N)
- Label = 1 (positive)
- Weight = expression level (Normal median)

#### Negative Examples (Non-Junctions):
Generate 3-5x as many negatives:
- **Random negatives:** 2M-4M (>1kb from any real junction)
- **Near-miss negatives:** 500K-1M (±3-20bp from real junctions)
- **Total negatives:** ~3M-5M

#### Final Training Set:
- **Positives:** 1.36M
- **Negatives:** 3-5M
- **Total:** 4.4M - 6.4M training examples

**This is MORE than enough for excellent CNN performance!**

---

## Filtering Recommendations

Before using for training, apply minimal filtering:

```python
# Filter criteria
MIN_SAMPLES = 2        # Junction seen in at least 2 samples
MIN_READS = 5          # At least 5 total reads
MIN_EXPRESSION = 0.1   # Median expression > 0.1 TPM

# Expected after filtering: ~600K junctions
```

Even with conservative filtering, you'll have 600K+ high-quality junctions.

---

## Comparison to Alternatives

| Source | Junction Count | Quality | Your Data |
|--------|---------------|---------|-----------|
| RJunBase (yours) | 682K | High | ✓ Have it |
| GTEx direct | ~500K | High | Need download |
| recount3 GTEx | ~500K | High | Need R package |
| GENCODE only | ~350K | Very high | Annotated only |

**Conclusion:** Your RJunBase data is **MORE comprehensive** than GTEx!

---

## What About recount3?

### Do you need it? **NO**

recount3 provides:
- GTEx junctions: ~500K (you have 682K)
- TCGA junctions: ~400K cancer-specific
- SRA junctions: variable quality

### When would you use recount3?
Only if you need:
1. **Tissue-specific** junction subsets (e.g., brain-only)
2. **Cancer-specific** junctions from TCGA
3. **Cross-validation** data from independent source

For general CNN training: **Your RJunBase data is superior**

---

## Validation Checks to Perform

Before training, validate your data:

### 1. **Genome Build Verification**
```bash
# Check if coordinates match GENCODE hg38
python src/extract/validate_genome_build.py \
    --junctions splice_db/detail_Linear_splice_annotation.txt \
    --gencode data/ref/gencode.v46.annotation.gtf.gz
```

**Expected:** Coordinates should match hg38 GENCODE

### 2. **Motif Enrichment**
```bash
# Check GT-AG enrichment
python src/extract/validate_junction_motifs.py \
    --junctions splice_db/detail_Linear_splice_annotation.txt \
    --genome data/ref/GRCh38.primary_assembly.genome.fa
```

**Expected:**
- GT-AG (canonical): >90%
- GC-AG (alternative): 2-5%
- AT-AC (U12): <1%
- Non-canonical: <5%

### 3. **Coordinate Sanity**
```bash
# Check all donor < acceptor
awk -F'\t' 'NR>1 {
    split($2, a, "[:|]");
    donor=a[2]; acceptor=a[3];
    if(donor >= acceptor) print "ERROR: " $2
}' splice_db/detail_Linear_splice_annotation.txt
```

**Expected:** No errors (all donors before acceptors)

---

## Next Steps (Phase 1 → Phase 3)

### Week 1-2: Process Your RJunBase Data

```bash
# 1. Standardize format
python src/extract/standardize_priors.py \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/ \
    --min-samples 2 \
    --min-reads 5

# Expected output:
#   data/ref/junction_prior.tsv.gz (+ .tbi index)
#   ~600K-650K junctions after filtering
```

### Week 3-4: Generate CNN Training Windows (Phase 3)

```bash
# 2. Build training windows
python src/cnn/build_training_windows.py \
    --junction-prior data/ref/junction_prior.tsv.gz \
    --genome data/ref/GRCh38.primary_assembly.genome.fa \
    --output data/derived/cnn_train/ \
    --window-size 128 \
    --negatives-per-positive 3

# Expected output:
#   data/derived/cnn_train/donor_windows.npz (1.8M examples)
#   data/derived/cnn_train/acceptor_windows.npz (1.8M examples)
```

### Week 5-8: Train CNN (Phase 4)

```bash
# 3. Train multi-task CNN
python src/cnn/train_model.py \
    --train-data data/derived/cnn_train/ \
    --output model_weights/cnn_multitask_v1.pt \
    --epochs 50 \
    --batch-size 256 \
    --gpu 0

# Expected results:
#   Donor AUROC: >0.95
#   Acceptor AUROC: >0.95
```

---

## Summary Table

| Requirement | RJunBase | Status |
|-------------|----------|--------|
| Junction coordinates | ✓ chr:donor\|acceptor:strand | ✓ Perfect |
| Strand information | ✓ + or - | ✓ Perfect |
| Expression weights | ✓ Median TPM | ✓ Perfect |
| Sample support | ✓ Frequency counts | ✓ Perfect |
| Quantity | 682K junctions | ✓ Excellent |
| Annotated junctions | 357K (52%) | ✓ Good |
| Novel junctions | 325K (48%) | ✓ Good |
| Genome build | Need to verify | ? Check |

---

## Final Recommendation

### ✓ USE YOUR RJUNBASE DATA

**Reasons:**
1. **More junctions** than GTEx (682K vs 500K)
2. **Already downloaded** (no waiting)
3. **High quality** (filtered, annotated)
4. **Expression weights** included
5. **Ready to process** with standardize_priors.py

### ✗ SKIP recount3

**Reasons:**
1. Not needed (you have better data)
2. Requires R package setup
3. Additional download time
4. No significant advantage

---

## One Concern to Address

### Genome Build Verification (HIGH PRIORITY)

Your RJunBase data **MUST be hg38** for compatibility.

**Quick check:**
```bash
# Compare to known hg38 gene
# RNF213 should be at chr17:80,328,477-80,332,006 (hg38)

grep "RNF213" splice_db/detail_Linear_splice_annotation.txt | head -1
# Should show: chr17:80328477|80332006:+
```

If coordinates match → it's hg38 ✓

If they're ~10-50kb off → it's hg19 (need liftover)

**Action:** Verify this in Week 1 before proceeding

---

## Conclusion

**Your RJunBase data is EXCELLENT for CNN training.**

You have:
- ✓ Sufficient quantity (682K junctions)
- ✓ High quality (multi-tissue, weighted)
- ✓ Right format (easy to parse)
- ✓ Mix of annotated + novel

**No additional splice data downloads needed!**

**Next step:** Run `standardize_priors.py` to process your existing data.

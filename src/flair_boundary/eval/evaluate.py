#!/usr/bin/env python3
"""
Evaluate trained boundary classifiers and compare model performance.

Generates:
  - Precision-Recall curves for each model
  - ROC curves
  - Feature importance plots (XGBoost)
  - Score distribution plots (positives vs negatives)
  - Per-chromosome performance breakdown
  - Calibration plots

Usage:
    python evaluate.py \\
        --features /path/to/data/derived/tss_features.npz \\
        --model-dir /path/to/model_weights/tss/ \\
        --output /path/to/eval_results/tss/ \\
        --type tss \\
        --holdout-chrom chr22
"""

import argparse
import json
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
    f1_score,
    confusion_matrix,
)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Evaluate boundary classifier models",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--features", required=True,
                   help="Feature matrix .npz from prepare_data.py")
    p.add_argument("--model-dir", required=True,
                   help="Directory containing trained model .pkl files")
    p.add_argument("--output", required=True,
                   help="Output directory for evaluation results")
    p.add_argument("--type", required=True, choices=["tss", "tts"],
                   help="Boundary type")
    p.add_argument("--holdout-chrom", default="chr22",
                   help="Chromosome to use as test set")
    p.add_argument("--models", nargs="+", default=None,
                   help="Specific model names to evaluate (default: all in model-dir)")
    p.add_argument("--threshold", type=float, default=0.0,
                   help="Fixed classification threshold (0 = use F1-optimal per model)")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

# Match train.py: exclude database-derived context features to prevent
# soft label leakage during evaluation.
CTX_FEATURE_PREFIX = "ctx_"


def load_features(npz_file: str) -> Tuple:
    data = np.load(npz_file, allow_pickle=True)
    X = data["X"].astype(np.float32)
    y = data["y"].astype(int)
    chroms = data["chroms"]
    positions = data["positions"]
    feature_names = list(data["feature_names"])
    # widths: per-row promoter-region width in bp; sentinel value -1 for
    # rows without a FANTOM5-fair overlap. Used only by the width
    # regressor evaluation; classifier metrics ignore it.
    widths = data["widths"].astype(float) if "widths" in data.files else None

    # Drop ctx_* features (same as train.py)
    keep_mask = [not name.startswith(CTX_FEATURE_PREFIX) for name in feature_names]
    dropped = [name for name in feature_names if name.startswith(CTX_FEATURE_PREFIX)]
    if dropped:
        X = X[:, keep_mask]
        feature_names = [n for n in feature_names if not n.startswith(CTX_FEATURE_PREFIX)]
        print(f"  Dropped {len(dropped)} ctx features: {dropped}")

    return X, y, chroms, positions, feature_names, widths


def load_models(model_dir: Path, model_names: Optional[List[str]] = None) -> Dict:
    models = {}
    for pkl_path in sorted(model_dir.glob("*.pkl")):
        name = pkl_path.stem
        if model_names and name not in model_names:
            continue
        try:
            with open(pkl_path, "rb") as f:
                models[name] = pickle.load(f)
            print(f"  Loaded: {name}")
        except Exception as e:
            print(f"  WARNING: Could not load {pkl_path}: {e}")
    return models


# ---------------------------------------------------------------------------
# Metrics computation
# ---------------------------------------------------------------------------

def compute_full_metrics(
    clf,
    X: np.ndarray,
    y: np.ndarray,
    threshold: float = 0.0,
) -> Dict:
    """Comprehensive metrics for a model on a test set."""
    proba = clf.predict_proba(X)[:, 1]

    # PR curve
    prec, rec, pr_thresholds = precision_recall_curve(y, proba)
    f1_vals = 2 * prec * rec / np.where(prec + rec == 0, 1e-9, prec + rec)

    if threshold > 0:
        opt_thr = threshold
    else:
        opt_idx = np.argmax(f1_vals[:-1])
        opt_thr = float(pr_thresholds[opt_idx])

    pred = (proba >= opt_thr).astype(int)

    # ROC
    fpr, tpr, roc_thresholds = roc_curve(y, proba)

    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y, pred).ravel()

    return {
        "roc_auc": float(roc_auc_score(y, proba)),
        "avg_precision": float(average_precision_score(y, proba)),
        "f1_optimal": float(f1_score(y, pred)),
        "precision_at_opt": float(prec[np.argmax(f1_vals[:-1])]),
        "recall_at_opt": float(rec[np.argmax(f1_vals[:-1])]),
        "optimal_threshold": opt_thr,
        "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        "n_pos": int((y == 1).sum()),
        "n_neg": int((y == 0).sum()),
        # Store curves for plotting
        "_prec": prec,
        "_rec": rec,
        "_fpr": fpr,
        "_tpr": tpr,
        "_proba": proba,
    }


def per_chromosome_metrics(
    clf,
    X: np.ndarray,
    y: np.ndarray,
    chroms: np.ndarray,
) -> pd.DataFrame:
    """Compute AUROC and AUPRC per chromosome."""
    proba = clf.predict_proba(X)[:, 1]
    rows = []
    for chrom in np.unique(chroms):
        mask = chroms == chrom
        y_c = y[mask]
        p_c = proba[mask]
        if (y_c == 1).sum() < 2 or (y_c == 0).sum() < 2:
            continue
        rows.append({
            "chrom": chrom,
            "n_pos": int((y_c == 1).sum()),
            "n_neg": int((y_c == 0).sum()),
            "roc_auc": float(roc_auc_score(y_c, p_c)),
            "avg_precision": float(average_precision_score(y_c, p_c)),
        })
    return pd.DataFrame(rows).sort_values("chrom")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pr_roc_curves(
    model_metrics: Dict[str, Dict],
    boundary_type: str,
    output_dir: Path,
):
    """Plot PR and ROC curves for all models on one figure."""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    for i, (name, metrics) in enumerate(model_metrics.items()):
        color = colors[i % len(colors)]
        prec = metrics["_prec"]
        rec = metrics["_rec"]
        fpr = metrics["_fpr"]
        tpr = metrics["_tpr"]
        auprc = metrics["avg_precision"]
        auroc = metrics["roc_auc"]

        axes[0].plot(rec, prec, color=color, lw=2,
                     label=f"{name} (AUPRC={auprc:.3f})")
        axes[1].plot(fpr, tpr, color=color, lw=2,
                     label=f"{name} (AUROC={auroc:.3f})")

    axes[0].set_xlabel("Recall", fontsize=12)
    axes[0].set_ylabel("Precision", fontsize=12)
    axes[0].set_title(f"{boundary_type.upper()} Precision-Recall Curves", fontsize=13, fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=9)
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].grid(alpha=0.3)

    axes[1].plot([0, 1], [0, 1], "k--", alpha=0.4, label="Random")
    axes[1].set_xlabel("False Positive Rate", fontsize=12)
    axes[1].set_ylabel("True Positive Rate", fontsize=12)
    axes[1].set_title(f"{boundary_type.upper()} ROC Curves", fontsize=13, fontweight="bold")
    axes[1].legend(loc="lower right", fontsize=9)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    out = output_dir / f"{boundary_type}_pr_roc_curves.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def plot_score_distributions(
    model_metrics: Dict[str, Dict],
    y_test: np.ndarray,
    boundary_type: str,
    output_dir: Path,
):
    """Plot score distributions for positives vs negatives."""
    if not HAS_MATPLOTLIB:
        return

    n_models = len(model_metrics)
    fig, axes = plt.subplots(1, n_models, figsize=(5 * n_models, 5))
    if n_models == 1:
        axes = [axes]

    for ax, (name, metrics) in zip(axes, model_metrics.items()):
        proba = metrics["_proba"]
        pos_scores = proba[y_test == 1]
        neg_scores = proba[y_test == 0]

        ax.hist(neg_scores, bins=50, alpha=0.6, color="#d62728",
                density=True, label=f"Non-boundary (n={len(neg_scores):,})")
        ax.hist(pos_scores, bins=50, alpha=0.6, color="#2ca02c",
                density=True, label=f"Boundary (n={len(pos_scores):,})")

        opt_thr = metrics["optimal_threshold"]
        ax.axvline(opt_thr, color="black", linestyle="--", lw=1.5,
                   label=f"Threshold={opt_thr:.3f}")

        ax.set_title(f"{name}\nAUPRC={metrics['avg_precision']:.3f}", fontsize=11)
        ax.set_xlabel("Predicted probability", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.suptitle(f"{boundary_type.upper()} Score Distributions", fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = output_dir / f"{boundary_type}_score_distributions.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def plot_feature_importance(
    clf,
    feature_names: List[str],
    boundary_type: str,
    model_name: str,
    output_dir: Path,
    top_n: int = 40,
):
    """Plot top-N feature importances."""
    if not HAS_MATPLOTLIB:
        return

    importances = None
    if hasattr(clf, "feature_importances_"):
        importances = clf.feature_importances_
    elif hasattr(clf, "named_steps"):
        inner = clf.named_steps.get("clf", None)
        if inner and hasattr(inner, "coef_"):
            importances = np.abs(inner.coef_[0])

    if importances is None:
        return

    df = pd.DataFrame({"feature": feature_names, "importance": importances})
    df = df.sort_values("importance", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(10, top_n * 0.25 + 2))
    y_pos = np.arange(len(df))
    ax.barh(y_pos, df["importance"].values, color="#1f77b4", alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["feature"].values, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Feature Importance", fontsize=11)
    ax.set_title(f"{boundary_type.upper()} Feature Importance\n{model_name} (Top {top_n})",
                 fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    out = output_dir / f"{boundary_type}_{model_name}_feature_importance.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def plot_per_chrom_performance(
    chrom_df: pd.DataFrame,
    boundary_type: str,
    model_name: str,
    output_dir: Path,
):
    """Bar chart of AUROC/AUPRC per chromosome."""
    if not HAS_MATPLOTLIB or chrom_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, metric, label, color in [
        (axes[0], "roc_auc", "ROC AUC", "#1f77b4"),
        (axes[1], "avg_precision", "Average Precision", "#ff7f0e"),
    ]:
        chroms = chrom_df["chrom"].values
        vals = chrom_df[metric].values
        x = np.arange(len(chroms))

        bars = ax.bar(x, vals, color=color, alpha=0.8)
        ax.axhline(vals.mean(), color="black", linestyle="--", lw=1.2,
                   label=f"Mean: {vals.mean():.3f}")
        ax.set_xticks(x)
        ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(f"{boundary_type.upper()} Per-Chromosome {label}\n{model_name}",
                     fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=9)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out = output_dir / f"{boundary_type}_{model_name}_per_chrom.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print("=" * 70)
    print(f"flAIr: Evaluating {args.type.upper()} Boundary Classifiers")
    print("=" * 70)
    print()

    # Load features
    X, y, chroms, positions, feature_names, widths = load_features(args.features)
    X = np.nan_to_num(X, nan=0.0)

    # Split holdout
    if args.holdout_chrom:
        mask = chroms == args.holdout_chrom
        if mask.sum() == 0:
            print(f"WARNING: Chromosome {args.holdout_chrom} not found in data. Using all.")
            X_test, y_test, chroms_test = X, y, chroms
            widths_test = widths
        else:
            X_test = X[mask]
            y_test = y[mask]
            chroms_test = chroms[mask]
            widths_test = widths[mask] if widths is not None else None
            print(f"Using holdout chromosome {args.holdout_chrom}: "
                  f"{mask.sum():,} sites ({(y_test==1).sum():,} pos, {(y_test==0).sum():,} neg)")
    else:
        X_test, y_test, chroms_test = X, y, chroms
        widths_test = widths
        print(f"Evaluating on all data: {len(X_test):,} sites")

    # Load models
    model_dir = Path(args.model_dir)
    print(f"\nLoading models from {model_dir}...")
    models = load_models(model_dir, args.models)

    if not models:
        print("ERROR: No models found")
        sys.exit(1)

    # Create output directory
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Split classifiers vs the (optional) width regressor. The width
    # regressor predicts continuous bp widths and has no .predict_proba();
    # it gets evaluated separately below with regression metrics.
    width_models = {n: m for n, m in models.items() if n.endswith("_width")}
    classifier_models = {n: m for n, m in models.items() if not n.endswith("_width")}

    # Evaluate classifiers
    print(f"\nEvaluating {len(classifier_models)} classifier(s)...")
    model_metrics = {}
    all_results = {}

    for name, clf in classifier_models.items():
        print(f"\n  {name}:")
        metrics = compute_full_metrics(clf, X_test, y_test, args.threshold)
        model_metrics[name] = metrics

        printable = {k: v for k, v in metrics.items() if not k.startswith("_")}
        all_results[name] = printable

        print(f"    AUROC:      {metrics['roc_auc']:.4f}")
        print(f"    AUPRC:      {metrics['avg_precision']:.4f}")
        print(f"    F1 optimal: {metrics['f1_optimal']:.4f}  @ threshold={metrics['optimal_threshold']:.3f}")
        print(f"    Precision:  {metrics['precision_at_opt']:.4f}")
        print(f"    Recall:     {metrics['recall_at_opt']:.4f}")
        print(f"    TP/FP/TN/FN: {metrics['tp']}/{metrics['fp']}/{metrics['tn']}/{metrics['fn']}")

        # Per-chromosome breakdown
        chrom_df = per_chromosome_metrics(clf, X_test, y_test, chroms_test)
        chrom_out = out_dir / f"{args.type}_{name}_per_chrom_metrics.tsv"
        chrom_df.to_csv(chrom_out, sep="\t", index=False)

        # Feature importance
        plot_feature_importance(clf, feature_names, args.type, name, out_dir)
        plot_per_chrom_performance(chrom_df, args.type, name, out_dir)

    # Width regressor evaluation (TSS only; trained against the FANTOM5
    # fair-peak overlap widths). Uses regression metrics — RMSE in log1p
    # space (training scale), MAE in bp, Spearman r — and ignores rows
    # with no width label (sentinel value -1).
    if width_models and widths_test is not None:
        print(f"\nEvaluating {len(width_models)} width regressor(s)...")
        labelled = (widths_test > 0) & (y_test == 1)
        n_labelled = int(labelled.sum())
        if n_labelled < 10:
            print(f"  Skipping: only {n_labelled} holdout positives have a width label")
        else:
            X_w = X_test[labelled]
            w_true = widths_test[labelled]
            log_w_true = np.log1p(w_true)
            for name, reg in width_models.items():
                log_w_pred = reg.predict(X_w)
                w_pred = np.expm1(log_w_pred)
                rmse_log = float(np.sqrt(np.mean((log_w_pred - log_w_true) ** 2)))
                mae_bp = float(np.mean(np.abs(w_pred - w_true)))
                med_ae_bp = float(np.median(np.abs(w_pred - w_true)))
                # Spearman correlation
                from scipy.stats import spearmanr
                rho, _ = spearmanr(w_pred, w_true)
                print(f"  {name}:  n={n_labelled}  RMSE(log1p)={rmse_log:.4f}  "
                      f"MAE={mae_bp:.1f}bp  median|err|={med_ae_bp:.1f}bp  "
                      f"Spearman r={rho:.4f}")
                all_results[name] = {
                    "n_labelled_holdout": n_labelled,
                    "rmse_log1p": rmse_log,
                    "mae_bp": mae_bp,
                    "median_ae_bp": med_ae_bp,
                    "spearman_r": float(rho),
                }

    # Combined plots
    print("\nGenerating combined plots...")
    plot_pr_roc_curves(model_metrics, args.type, out_dir)
    plot_score_distributions(model_metrics, y_test, args.type, out_dir)

    # Save JSON summary
    results_path = out_dir / f"{args.type}_evaluation_summary.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved evaluation summary: {results_path}")

    # Summary table
    print("\n" + "=" * 70)
    print("EVALUATION SUMMARY")
    print("=" * 70)
    print(f"{'Model':<28} {'AUROC':>8} {'AUPRC':>8} {'F1':>8} {'Threshold':>10}")
    print("-" * 70)
    for name, metrics in sorted(model_metrics.items(), key=lambda x: -x[1]["avg_precision"]):
        print(f"{name:<28} {metrics['roc_auc']:>8.4f} {metrics['avg_precision']:>8.4f} "
              f"{metrics['f1_optimal']:>8.4f} {metrics['optimal_threshold']:>10.3f}")

    print(f"\nAll results saved to: {out_dir}/")
    print("Done!")


if __name__ == "__main__":
    main()

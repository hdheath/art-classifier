#!/usr/bin/env python3
"""
Train boundary classifiers (TSS and TTS) using multiple model approaches.

Models trained:
  1. XGBoost (primary) — gradient-boosted trees on k-mer + motif features
  2. Logistic Regression (baseline) — fast, interpretable
  3. Random Forest — ensemble, good OOB calibration
  4. XGBoost + Platt scaling — calibrated probability outputs

All models output probability scores in [0, 1] representing the likelihood
that a given genomic position is a true TSS or TTS boundary.

Usage:
    python train.py \\
        --features /path/to/data/derived/tss_features.npz \\
        --output /path/to/model_weights/tss/ \\
        --type tss \\
        --cv-folds 5 \\
        --model xgboost logistic rf

    # Training on chr1-21 only, evaluating on chr22
    python train.py \\
        --features /path/to/data/derived/tss_features.npz \\
        --output /path/to/model_weights/tss/ \\
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
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    log_loss,
    f1_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import xgboost as xgb

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MODEL_KEYS = ["xgboost", "logistic", "rf", "xgboost_calibrated"]

# Width regression defaults (separate model, regressor not classifier)
XGBOOST_WIDTH_DEFAULTS = {
    "n_estimators": 500,
    "max_depth": 5,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    # log1p target → squared error in log-space ≈ multiplicative error in bp
    "objective": "reg:squarederror",
    "eval_metric": "rmse",
    "tree_method": "hist",
    "random_state": 42,
    "n_jobs": -1,
}

XGBOOST_DEFAULTS = {
    "n_estimators": 500,
    "max_depth": 6,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "gamma": 0.1,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "use_label_encoder": False,
    "eval_metric": "aucpr",
    "tree_method": "hist",
    "random_state": 42,
    "n_jobs": -1,
}

LOGISTIC_DEFAULTS = {
    "C": 0.1,
    "max_iter": 1000,
    "solver": "lbfgs",
    "random_state": 42,
}

RF_DEFAULTS = {
    "n_estimators": 300,
    "max_depth": 10,
    "min_samples_leaf": 5,
    "max_features": "sqrt",
    "oob_score": True,
    "random_state": 42,
    "n_jobs": -1,
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Train TSS/TTS boundary classifiers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--features", required=True,
                   help="Feature matrix .npz from prepare_data.py")
    p.add_argument("--output", required=True,
                   help="Output directory for model weights and metrics")
    p.add_argument("--type", required=True, choices=["tss", "tts"],
                   help="Boundary type (must match feature file)")

    p.add_argument("--model", nargs="+", default=["xgboost"],
                   choices=MODEL_KEYS,
                   help="Which model(s) to train")
    p.add_argument("--cv-folds", type=int, default=5,
                   help="Stratified cross-validation folds")
    p.add_argument("--holdout-chrom", default="",
                   help="Hold out a chromosome for final evaluation (e.g. chr22)")
    p.add_argument("--scale-features", action="store_true",
                   help="Standardize features before training (required for logistic/rf)")
    p.add_argument("--pos-weight", type=float, default=0.0,
                   help="Positive class weight for XGBoost scale_pos_weight "
                        "(0 = auto-compute from class ratio)")
    p.add_argument("--early-stopping", type=int, default=50,
                   help="XGBoost early stopping rounds (0 = disabled)")
    p.add_argument("--feature-importance", action="store_true",
                   help="Save feature importance rankings after training")
    p.add_argument("--width-model", action="store_true",
                   help="Also train a width regressor (positives only) for "
                        "predicting initiation/termination zone width in bp. "
                        "Saved as xgboost_width.pkl alongside the classifier.")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

# Features with this prefix are database-derived context features that
# cause soft label leakage during training (e.g. ctx_nn_dist_log1p encodes
# proximity to the reference boundary set, which trivially separates
# positives from negatives by construction).  They are kept in the .npz
# for potential inference-time use but excluded from training/evaluation.
CTX_FEATURE_PREFIX = "ctx_"


def load_features(npz_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[str], np.ndarray]:
    """
    Load feature matrix from .npz file.

    Database-derived context features (ctx_*) are excluded from training
    to avoid soft label leakage — they remain stored in the .npz for
    inference-time use.

    Returns
    -------
    X, y, chroms, positions, feature_names, widths
    `widths` is the regression target for width prediction (peak width in
    bp for positives, -1 for negatives). Older .npz files without a
    `widths` array return an all-(-1) vector.
    """
    print(f"Loading features from {npz_file}...")
    data = np.load(npz_file, allow_pickle=True)

    X = data["X"].astype(np.float32)
    y = data["y"].astype(int)
    chroms = data["chroms"]
    positions = data["positions"]
    feature_names = list(data["feature_names"])
    if "widths" in data.files:
        widths = data["widths"].astype(float)
    else:
        widths = np.full(len(y), -1.0, dtype=float)
        print("  NOTE: .npz has no 'widths' field — width regression disabled "
              "(re-run prepare_data.py to enable)")

    stored_type = str(data["boundary_type"][0])
    print(f"  Boundary type: {stored_type}")
    print(f"  Samples: {X.shape[0]:,}  Features (raw): {X.shape[1]}")
    print(f"  Positives: {(y==1).sum():,}  Negatives: {(y==0).sum():,}")

    # Drop ctx_* features to prevent leakage during training
    keep_mask = [not name.startswith(CTX_FEATURE_PREFIX) for name in feature_names]
    dropped = [name for name in feature_names if name.startswith(CTX_FEATURE_PREFIX)]
    if dropped:
        X = X[:, keep_mask]
        feature_names = [n for n in feature_names if not n.startswith(CTX_FEATURE_PREFIX)]
        print(f"  Dropped {len(dropped)} ctx features (leakage prevention): {dropped}")
        print(f"  Training features: {X.shape[1]}")

    return X, y, chroms, positions, feature_names, widths


# ---------------------------------------------------------------------------
# Model builders
# ---------------------------------------------------------------------------

def build_xgboost(
    scale_pos_weight: float = 1.0,
    early_stopping: int = 50,
    **kwargs,
) -> xgb.XGBClassifier:
    params = {**XGBOOST_DEFAULTS, "scale_pos_weight": scale_pos_weight}
    params.update(kwargs)
    if early_stopping > 0:
        params["early_stopping_rounds"] = early_stopping
    return xgb.XGBClassifier(**params)


def build_logistic(**kwargs) -> Pipeline:
    params = {**LOGISTIC_DEFAULTS, **kwargs}
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(**params)),
    ])


def build_rf(**kwargs) -> RandomForestClassifier:
    params = {**RF_DEFAULTS, **kwargs}
    return RandomForestClassifier(**params)


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def evaluate_model(clf, X_test: np.ndarray, y_test: np.ndarray) -> Dict:
    """Compute standard classification metrics on a test set."""
    proba = clf.predict_proba(X_test)[:, 1]
    pred = (proba >= 0.5).astype(int)

    prec, rec, thresholds = precision_recall_curve(y_test, proba)

    # F1 at optimal threshold
    f1_scores = 2 * prec * rec / np.where(prec + rec == 0, 1, prec + rec)
    best_thr_idx = np.argmax(f1_scores[:-1])
    best_thr = thresholds[best_thr_idx]
    pred_opt = (proba >= best_thr).astype(int)

    return {
        "roc_auc": roc_auc_score(y_test, proba),
        "avg_precision": average_precision_score(y_test, proba),
        "log_loss": log_loss(y_test, proba),
        "f1_at_0.5": f1_score(y_test, pred),
        "f1_optimal": f1_score(y_test, pred_opt),
        "optimal_threshold": best_thr,
        "n_pos": int((y_test == 1).sum()),
        "n_neg": int((y_test == 0).sum()),
    }


# ---------------------------------------------------------------------------
# Cross-validation
# ---------------------------------------------------------------------------

def cross_validate(
    clf,
    X: np.ndarray,
    y: np.ndarray,
    n_folds: int = 5,
    model_name: str = "",
    is_xgboost: bool = False,
    early_stopping: int = 0,
) -> Dict:
    """Run stratified K-fold cross-validation, return mean metrics."""
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    fold_metrics = []

    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X, y)):
        X_tr, X_val = X[train_idx], X[val_idx]
        y_tr, y_val = y[train_idx], y[val_idx]

        if is_xgboost and early_stopping > 0:
            # Use internal XGBoost eval set
            clf.fit(
                X_tr, y_tr,
                eval_set=[(X_val, y_val)],
                verbose=False,
            )
        else:
            clf.fit(X_tr, y_tr)

        metrics = evaluate_model(clf, X_val, y_val)
        metrics["fold"] = fold_idx
        fold_metrics.append(metrics)

        print(f"    Fold {fold_idx+1}/{n_folds}: "
              f"AUROC={metrics['roc_auc']:.4f}  "
              f"AUPRC={metrics['avg_precision']:.4f}  "
              f"F1={metrics['f1_optimal']:.4f}")

    # Aggregate
    mean_metrics = {}
    for key in fold_metrics[0]:
        if key == "fold":
            continue
        vals = [m[key] for m in fold_metrics]
        mean_metrics[f"{key}_mean"] = float(np.mean(vals))
        mean_metrics[f"{key}_std"] = float(np.std(vals))

    return mean_metrics


# ---------------------------------------------------------------------------
# Feature importance
# ---------------------------------------------------------------------------

def save_feature_importance(
    clf,
    feature_names: List[str],
    output_path: str,
    model_name: str,
):
    """Save feature importance rankings."""
    if hasattr(clf, "feature_importances_"):
        importances = clf.feature_importances_
    elif hasattr(clf, "named_steps"):
        # Pipeline
        inner = clf.named_steps.get("clf", None)
        if inner and hasattr(inner, "coef_"):
            importances = np.abs(inner.coef_[0])
        else:
            return
    else:
        return

    df = pd.DataFrame({
        "feature": feature_names,
        "importance": importances,
    }).sort_values("importance", ascending=False)

    out = Path(output_path) / f"{model_name}_feature_importance.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"    Saved feature importance: {out}")

    # Top 20
    print(f"    Top 20 features:")
    for _, row in df.head(20).iterrows():
        print(f"      {row['feature']:50s} {row['importance']:.6f}")


# ---------------------------------------------------------------------------
# Main training function
# ---------------------------------------------------------------------------

def train_model(
    model_name: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: Optional[np.ndarray],
    y_test: Optional[np.ndarray],
    feature_names: List[str],
    output_dir: Path,
    args,
) -> Dict:
    """Train a single model, cross-validate, save weights."""
    print(f"\n{'='*60}")
    print(f"Training: {model_name.upper()}")
    print(f"{'='*60}")

    # Class weight for imbalanced data
    pos_count = (y_train == 1).sum()
    neg_count = (y_train == 0).sum()
    auto_pos_weight = neg_count / max(pos_count, 1)

    is_xgboost = "xgboost" in model_name

    # Build model
    if model_name == "xgboost":
        spw = args.pos_weight if args.pos_weight > 0 else auto_pos_weight
        clf = build_xgboost(
            scale_pos_weight=spw,
            early_stopping=args.early_stopping,
        )
    elif model_name == "logistic":
        clf = build_logistic()
    elif model_name == "rf":
        clf = build_rf()
    elif model_name == "xgboost_calibrated":
        spw = args.pos_weight if args.pos_weight > 0 else auto_pos_weight
        base = build_xgboost(scale_pos_weight=spw, early_stopping=0)
        clf = CalibratedClassifierCV(base, method="sigmoid", cv=3)
        is_xgboost = False

    print(f"  Class balance: {pos_count:,} pos / {neg_count:,} neg "
          f"(auto weight={auto_pos_weight:.2f})")

    # --- Cross-validation ---
    if args.cv_folds > 1:
        print(f"\n  Cross-validation ({args.cv_folds} folds)...")
        cv_metrics = cross_validate(
            clf, X_train, y_train,
            n_folds=args.cv_folds,
            model_name=model_name,
            is_xgboost=is_xgboost,
            early_stopping=args.early_stopping if is_xgboost else 0,
        )
        print(f"\n  CV summary: AUROC={cv_metrics['roc_auc_mean']:.4f}±{cv_metrics['roc_auc_std']:.4f}  "
              f"AUPRC={cv_metrics['avg_precision_mean']:.4f}±{cv_metrics['avg_precision_std']:.4f}")
    else:
        cv_metrics = {}

    # --- Final fit on full training data ---
    # NOTE: We deliberately do NOT pass (X_test, y_test) as the early-stopping
    # eval_set, because X_test is the held-out chromosome reserved for honest
    # evaluation. Using it for early-stopping would let the model select its
    # n_estimators against the metric we report, which inflates holdout
    # numbers. Instead, carve a small validation slice out of the training
    # set just for early-stopping.
    print("\n  Fitting final model on full training set...")
    if is_xgboost and args.early_stopping > 0:
        from sklearn.model_selection import train_test_split as _tts
        X_tr, X_es, y_tr, y_es = _tts(
            X_train, y_train, test_size=0.1, stratify=y_train, random_state=42
        )
        clf.fit(X_tr, y_tr, eval_set=[(X_es, y_es)], verbose=False)
    else:
        clf.fit(X_train, y_train)

    # --- Holdout evaluation ---
    holdout_metrics = {}
    if X_test is not None and y_test is not None:
        print("\n  Holdout evaluation...")
        holdout_metrics = evaluate_model(clf, X_test, y_test)
        print(f"    AUROC={holdout_metrics['roc_auc']:.4f}  "
              f"AUPRC={holdout_metrics['avg_precision']:.4f}  "
              f"F1={holdout_metrics['f1_optimal']:.4f}  "
              f"Threshold={holdout_metrics['optimal_threshold']:.3f}")

    # --- Save model ---
    model_path = output_dir / f"{model_name}.pkl"
    with open(model_path, "wb") as f:
        pickle.dump(clf, f)
    print(f"\n  Saved model: {model_path}")

    # --- Feature importance ---
    if args.feature_importance:
        save_feature_importance(clf, feature_names, str(output_dir), model_name)

    # Compile all results
    results = {
        "model": model_name,
        "boundary_type": args.type,
        "n_train": len(X_train),
        "n_pos_train": int(pos_count),
        "n_neg_train": int(neg_count),
        "cv_metrics": cv_metrics,
        "holdout_metrics": holdout_metrics,
        "model_path": str(model_path),
    }

    return results


# ---------------------------------------------------------------------------
# Width regression (positives only)
# ---------------------------------------------------------------------------

def train_width_regressor(
    X: np.ndarray,
    y: np.ndarray,
    widths: np.ndarray,
    chroms: np.ndarray,
    feature_names: List[str],
    output_dir: Path,
    boundary_type: str,
    holdout_chrom: str = "",
) -> Dict:
    """Train an XGBoost regressor that predicts peak width (bp) from sequence.

    Trained on positives only (negatives have no defined width). Target is
    log1p(width) so the model handles the long tail of broad CpG-island
    promoters (sometimes >500bp) without being dominated by them.

    Used at inference by ted.py to set per-locus HDBSCAN tss_slack: broad
    promoters get a larger clustering radius, narrow TATA-driven TSSs get
    a tighter one.
    """
    pos_mask = (y == 1) & (widths > 0)
    n_pos = int(pos_mask.sum())
    if n_pos < 100:
        print(f"\n  WIDTH REGRESSOR: skipping — only {n_pos} positives with "
              f"valid widths (need ≥100)")
        return {}

    print("\n" + "=" * 60)
    print("Training: WIDTH REGRESSOR (XGBoost)")
    print("=" * 60)

    Xw = X[pos_mask]
    yw = np.log1p(widths[pos_mask].astype(float))
    cw = chroms[pos_mask]

    # Holdout split on chromosome (same as classifier for honest comparison)
    if holdout_chrom and (cw == holdout_chrom).any():
        train_mask = cw != holdout_chrom
        Xw_train, yw_train = Xw[train_mask], yw[train_mask]
        Xw_test, yw_test = Xw[~train_mask], yw[~train_mask]
        print(f"  Holdout {holdout_chrom}: {len(yw_test):,} positives held out")
    else:
        from sklearn.model_selection import train_test_split as _tts
        Xw_train, Xw_test, yw_train, yw_test = _tts(
            Xw, yw, test_size=0.2, random_state=42
        )
        print(f"  No holdout chrom available; 80/20 random split")

    print(f"  Train: {len(yw_train):,}  Test: {len(yw_test):,}")
    print(f"  Width distribution (bp): "
          f"median={int(np.median(np.expm1(yw_train)))}  "
          f"p10={int(np.percentile(np.expm1(yw_train), 10))}  "
          f"p90={int(np.percentile(np.expm1(yw_train), 90))}")

    # Carve early-stopping val out of train (don't leak holdout)
    from sklearn.model_selection import train_test_split as _tts
    Xw_tr, Xw_es, yw_tr, yw_es = _tts(
        Xw_train, yw_train, test_size=0.1, random_state=42
    )

    params = {**XGBOOST_WIDTH_DEFAULTS, "early_stopping_rounds": 50}
    reg = xgb.XGBRegressor(**params)
    reg.fit(Xw_tr, yw_tr, eval_set=[(Xw_es, yw_es)], verbose=False)

    # Evaluate on holdout
    pred_log = reg.predict(Xw_test)
    pred_bp = np.expm1(pred_log)
    true_bp = np.expm1(yw_test)
    rmse_log = float(np.sqrt(np.mean((pred_log - yw_test) ** 2)))
    mae_bp = float(np.mean(np.abs(pred_bp - true_bp)))
    median_ae_bp = float(np.median(np.abs(pred_bp - true_bp)))
    # Spearman correlation — does the model rank-order narrow vs. broad correctly?
    from scipy.stats import spearmanr
    spearman_r, _ = spearmanr(pred_log, yw_test)

    print(f"  Holdout RMSE (log1p):     {rmse_log:.4f}")
    print(f"  Holdout MAE (bp):         {mae_bp:.1f}")
    print(f"  Holdout median |err| (bp): {median_ae_bp:.1f}")
    print(f"  Holdout Spearman r:       {spearman_r:.4f}")

    # Save
    out_path = output_dir / "xgboost_width.pkl"
    with open(out_path, "wb") as f:
        pickle.dump(reg, f)
    print(f"  Saved: {out_path}")

    return {
        "model": "xgboost_width",
        "boundary_type": boundary_type,
        "n_train": len(yw_train),
        "n_test": len(yw_test),
        "holdout_metrics": {
            "rmse_log1p": rmse_log,
            "mae_bp": mae_bp,
            "median_ae_bp": median_ae_bp,
            "spearman_r": float(spearman_r),
        },
        "model_path": str(out_path),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print("=" * 70)
    print(f"flAIr: Training {args.type.upper()} Boundary Classifier")
    print("=" * 70)
    print()

    # Load data
    X, y, chroms, positions, feature_names, widths = load_features(args.features)

    # Replace NaN with 0 (shouldn't be many after prepare_data filtering)
    X = np.nan_to_num(X, nan=0.0)

    # Split into train/test based on holdout chromosome or stratified split
    if args.holdout_chrom:
        hold_mask = chroms == args.holdout_chrom
        train_mask = ~hold_mask
        print(f"\nHoldout chromosome: {args.holdout_chrom} "
              f"({hold_mask.sum():,} sites)")
        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[hold_mask], y[hold_mask]
        if len(X_test) == 0:
            print(f"WARNING: No sites on {args.holdout_chrom}, using 20% split")
            from sklearn.model_selection import train_test_split
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.2, stratify=y, random_state=42
            )
    else:
        from sklearn.model_selection import train_test_split
        print("\nUsing 80/20 stratified train/test split")
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, stratify=y, random_state=42
        )

    print(f"Train: {len(X_train):,}  Test: {len(X_test):,}")

    # Create output directory
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Train width regressor (positives only) before classifier loop, so it
    # uses the full positive set regardless of which classifier variants
    # are requested.
    width_result = {}
    if args.width_model:
        width_result = train_width_regressor(
            X=X, y=y, widths=widths, chroms=chroms,
            feature_names=feature_names,
            output_dir=out_dir,
            boundary_type=args.type,
            holdout_chrom=args.holdout_chrom,
        )

    # Train each requested model
    all_results = []
    for model_name in args.model:
        result = train_model(
            model_name=model_name,
            X_train=X_train,
            y_train=y_train,
            X_test=X_test if len(X_test) > 0 else None,
            y_test=y_test if len(y_test) > 0 else None,
            feature_names=feature_names,
            output_dir=out_dir,
            args=args,
        )
        all_results.append(result)

    # Save metadata
    meta = {
        "boundary_type": args.type,
        "n_features": len(feature_names),
        "feature_names": feature_names,
        "holdout_chrom": args.holdout_chrom,
        "models": all_results,
        "width_model": width_result if width_result else None,
    }
    meta_path = out_dir / "training_metadata.json"

    def _json_default(obj):
        if isinstance(obj, (np.floating, np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} not JSON serializable")

    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2, default=_json_default)
    print(f"\nSaved training metadata: {meta_path}")

    # Summary table
    print("\n" + "=" * 70)
    print("TRAINING SUMMARY")
    print("=" * 70)
    print(f"{'Model':<25} {'AUROC':>8} {'AUPRC':>8} {'F1_opt':>8}")
    print("-" * 70)
    for r in all_results:
        hm = r.get("holdout_metrics", {})
        auc = f"{hm.get('roc_auc', 0):.4f}" if hm else "N/A"
        apr = f"{hm.get('avg_precision', 0):.4f}" if hm else "N/A"
        f1 = f"{hm.get('f1_optimal', 0):.4f}" if hm else "N/A"
        print(f"{r['model']:<25} {auc:>8} {apr:>8} {f1:>8}")

    print(f"\nAll models saved to: {out_dir}/")
    print("Done!")


if __name__ == "__main__":
    main()

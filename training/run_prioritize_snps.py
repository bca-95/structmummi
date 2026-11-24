#!/usr/bin/env python
# coding: utf-8
import subprocess
import os
from pathlib import Path

import xgboost as xgb
import json
import joblib
from joblib import Parallel, delayed

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, recall_score, confusion_matrix, roc_auc_score, matthews_corrcoef

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

import pandas as pd
import numpy as np
import cupy as cp
import itertools
from kneed import KneeLocator

import shap


from datetime import datetime
import time
import argparse


# --------------------------------
# - This read the config file that contain the input_file in genotype matrix file
# and all parameters that need to be set (see below) 
parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()

# --------------------------------
# --- Load the input variables
with open(args.config) as f:
    config = json.load(f)

genotype_path     = config["genotype_path"]
model_file        = config["model_file"]
features_file     = config["features_file"]
best_tuned_params_file = config["best_tuned_params_file"]
device            = config["device"]
nthread           = config.get("n_jobs_prio", 46)
inter_batch_size  = int(config.get("inter_batch_size", 50))
gpu_interaction   = config["gpu_interaction"] == "yes"

vep_dir           = config["vep_dir"]
fasta_file        = config["fasta_file"]
assembly          = config["assembly"]

repeat_ids        = [int(x) for x in config.get("repeat_ids", range(1, config["nb_repeats"] + 1))]
nb_repeats        = len(repeat_ids)

top_nb_snps       = int(config["top_nb_snps"])
ntop_inter1       = int(config["ntop_inter1"])
ntop_inter2       = int(config["ntop_inter2"])
disease_name      = config.get("disease_name", "diseaseX")

param_out     = config["param_out"]                       # folder with repeat_1 ... repeat_N (models are read here)
result_out    = config.get("result_out", param_out)     # consensus folder: all results are written here
test_size     = config["test_size"]
rd_state_test = config["rd_state_test"]

# --------------------------------


# --------------------------------
# - Record time 
in_time = time.time()
# --------------------------------


# --------------------------------
# - the directory snps dir is automatically created in the .sh file 
snps_dir = f"{result_out}/selected_snps/"
# --------------------------------


# --------------------------------
# In the first step of the model construction, the entry data were
# pre-processed with PLINK which has generated variants name as follow
# chr:pos:ref:alt_alt. Remove the underscore _alt
def strip_snp_suffix(name):
    return name.split("_")[0]
# --------------------------------

# --------------------------------
# - load model
def func_load_model( repeat_id) :
    model = xgb.Booster()
    model.load_model(f"{param_out}/{repeat_id}/{model_file}")
 

    # - load selected features
    top_feature_names = np.load(f"{param_out}/{repeat_id}/{features_file}", allow_pickle=True)
    top_feature_names = np.array([strip_snp_suffix(f) for f in top_feature_names])  


    # - load tuned parameters
    with open(f"{param_out}/{repeat_id}/{best_tuned_params_file}") as f:
        tuned_params = json.load(f)

    best_threshold = tuned_params["threshold_mcc"]

    return model, top_feature_names, best_threshold

# --------------------------------


# --------------------------------
def to_np(x):
    return cp.asnumpy(x) if device == "cuda:0" else x
# --------------------------------


# --------------------------------
# -- Remove correlated SNPs within the same genes

def func_remove_corr_SNPs(X_set):
    # X_set here already include only the top features (not all features)

    # - compute correlation
    corr_matrix_spearman = X_set.corr(method='spearman')

    # - Compute absolute correlation, and highlight SNPs with corr > 0.2
    abs_corr = corr_matrix_spearman.abs()

    # Mask diagonal (self-correlation)
    np.fill_diagonal(abs_corr.values, 0)

    # Identify SNP pairs with high correlation
    high_corr_pairs = np.where(abs_corr > 0.2)
    high_corr_df = pd.DataFrame({
        "SNP1": corr_matrix_spearman.index[high_corr_pairs[0]],
        "SNP2": corr_matrix_spearman.columns[high_corr_pairs[1]],
        "correlation": corr_matrix_spearman.values[high_corr_pairs]
    })

    # Remove duplicates (since matrix is symmetric)
    high_corr_df = high_corr_df[high_corr_df["SNP1"] < high_corr_df["SNP2"]]

    return
# --------------------------------


#----------------------------------------
# Function to find the best threshold for mcc
# 0.5 might not adapted to every case
def func_find_mcc_threshold(y_prob, y_true ):
    thresholds = np.linspace(0.4, 0.68, 15)
    best_mcc = -1
    best_t = 0.5
    for t in thresholds:
        y_pred = (y_prob > t).astype(int)
        mcc = matthews_corrcoef(y_true, y_pred)
        if mcc > best_mcc:
            best_mcc = mcc
            best_t = t
    return best_t, best_mcc

#----------------------------------------


#----------------------------------------
# I do not have a lot of data, for some stability in the search of best_t
# I will not use the all in-sample of val dataset prediction (used during the training of the meta-model)
# but split it to cross validation
def func_find_best_threshold_cv(preds_models_scaled, Y_val_np, C, n_splits=5, random_state=rd_state_test):

    skf_thr = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    oof_preds = np.zeros(len(Y_val_np))

    for train_idx, val_idx in skf_thr.split(preds_models_scaled, Y_val_np):
        X_tr, X_va = preds_models_scaled[train_idx], preds_models_scaled[val_idx]
        y_tr       = Y_val_np[train_idx]

        clf = LogisticRegression(C=C, l1_ratio=0, solver="lbfgs", max_iter=1000)
        clf.fit(X_tr, y_tr)

        oof_preds[val_idx] = clf.predict_proba(X_va)[:, 1]

    best_t, best_mcc_oof = func_find_mcc_threshold(oof_preds, Y_val_np)
    return best_t, best_mcc_oof
#----------------------------------------


#----------------------------------------
# During the hyper-parameter tuning, use mcc as evaluation metric
def func_mcc_eval(model, Y_val, dvalid, num_boost_i):
    # Test multiple threshold, because in some model MCC is weaker, but AUC greater
    # so may be the problem is on the MCC threshold
    thresholds = np.linspace(0.4, 0.68, 10)

    # compute the prediction for each iteration in num_boost_i
    preds = model.predict(dvalid, iteration_range=(0, num_boost_i))

    if device=="cuda:0" :
        y_true = cp.asnumpy(Y_val)
        preds  = cp.asnumpy(preds)

    _, best_mcc = func_find_mcc_threshold(preds, y_true)
    return best_mcc

#----------------------------------------



# --------------------------------
def func_compute_shap(feature_names_subset, explainer_shap, xtest) :
    # top_nb_snps : we need to define the number of snps to prioritize
    shap_values = explainer_shap.shap_values(xtest)

    # mean absolute SHAP value per feature
    mean_abs_shap = np.abs(shap_values).mean(axis=0)

    # make a dictionary of snps with their corresponding shap
    dict_snps_shap = dict(zip(feature_names_subset, mean_abs_shap))

    return dict_snps_shap
# --------------------------------


# --------------------------------
# - Combine SNP-level importance scores across repeats into one consensus
# ranking, weighted by each repeat's meta-model coefficient.

#  Re-normalized internally over the repeats where a given SNP was actually selected,
# so a SNP present in fewer repeats isn't automatically penalized
def func_combine_snp_rankings(nb_repeats, repeat_weights, shap_list, top_nb_snps=top_nb_snps):

    
    #top_shap_list = dict( sorted(shap_list.items(), key=lambda kv: kv[1], reverse=True)[:top_nb_snps] )

    # Weights are normalized to sum to 1 first, so the
    # result has an interpretable, run-independent scale
    repeat_weights = np.asarray(repeat_weights, dtype=float)
    repeat_weights = repeat_weights / repeat_weights.sum()  


    # union of every SNP that appeared in at least one repeat's ranking
    all_snps = set()
    for d in shap_list:
        all_snps.update(d.keys())

    records = []
    for snp in all_snps:
        kappas = np.array([d.get(snp, 0.0) for d in shap_list])  # 0 if absent
        n_present = int(np.sum(kappas > 0))  

        weighted_score = np.sum(repeat_weights * kappas)

        records.append({
            "feature": snp,
            "final_importance": weighted_score,
            "stability_fraction": n_present / nb_repeats,
        })
    full_result = pd.DataFrame(records).sort_values("final_importance", ascending=False)


    # i have to define the top x's snps to prioritize (sorted)
    if top_nb_snps is not None:
        selected_result = full_result.head(top_nb_snps)

    return full_result, selected_result

# -------------------------------- 


# --------------------------------
# SHAP interaction values with shap.TreeExplainer, on GPU or CPU.
# Returns the mean |SHAP interaction value| over the samples, for every pair of features.


def _shap_interactions_sum(explainer, X_chunk):
    # sum over the samples of |SHAP interaction values| for one chunk
    # so later I can easily compute the mean
    # ps: computing in chunck does not modify the values of interaction
    # interaction is computed per samples
    shap_interactions = explainer.shap_interaction_values(X_chunk)
    sum_abs_inter     = np.abs(shap_interactions).sum(axis=0)
    return sum_abs_inter
 
 
def _shap_interactions_sum_cpu_thread(model, X_chunk, threads):
    # To enable compute SHAP interaction using multiple thread,
    # each thread should receive 1 chunck (or by default, shap interaction can only be computed one sample per time)
    # For this parallelism, each thread should receive a copy of the model
    model_copy = model.copy()
    model_copy.set_param({"device": "cpu", "nthread": threads})
    explainer  = shap.TreeExplainer(model_copy)
    return _shap_interactions_sum(explainer, X_chunk)
 
 
def func_compute_shap_interactions(model, X_subset_np, gpu_interaction=gpu_interaction, batch_size=inter_batch_size):
    n_samples, p = X_subset_np.shape
 
    if gpu_interaction:
        # --- GPU: the interaction is computed per chunck
        model.set_param({"device": device})
        explainer = shap.TreeExplainer(model)
        sum_abs = np.zeros((p, p), dtype=np.float64)
        for start in range(0, n_samples, batch_size):
            sum_abs += _shap_interactions_sum(explainer, X_subset_np[start:start + batch_size])
 
    else:
        # --- CPU: the interaction is computed by dividing the samples per number of threads
        # Here we devide the samples by the number of threads

        n_workers = max(1, min(nthread, n_samples))
        threads_per_worker = max(1, nthread // n_workers)
        chunks = np.array_split(X_subset_np, n_workers)

        sums = Parallel(n_jobs=n_workers, prefer="threads")(
            delayed(_shap_interactions_sum_cpu_thread)(model, chunk, threads_per_worker)
            for chunk in chunks)
        sum_abs = np.sum(sums, axis=0)
 
    mean_abs_inter = sum_abs / n_samples                   # same as np.abs(values).mean(axis=0)
    np.fill_diagonal(mean_abs_inter, 0)
    return mean_abs_inter


#def func_compute_shap_interactions_cpu(model, X_subset_np, batch_size=1000):
#
##    # By using set_param, the shap interaction is computed on cpu
##    model.set_param({"device": "cpu", "nthread": nthread})   
#
#    model.set_param({"device": "cuda"})  
#    explainer_cpu = shap.TreeExplainer(model)
#
#    n_samples, p = X_subset_np.shape
#
#    shap_interactions = explainer_cpu.shap_interaction_values(X_subset_np)
#    mean_abs_inter = np.abs(shap_interactions).mean(axis=0)
#
#    np.fill_diagonal(mean_abs_inter, 0)
#    return mean_abs_inter


def func_combine_interaction_rankings(nb_repeats, repeat_weights, interaction_df_list):
    dfs = []
    for i, df in enumerate(interaction_df_list):
        df = df.copy()
        df["repeat_weight"] = repeat_weights[i]
        dfs.append(df)
    all_df = pd.concat(dfs, ignore_index=True)
    #print("func_combine_interaction_rankings_fast -- all_df ----- ", all_df.head(), all_df.shape)

    # vectorized weighted average 
    all_df["weighted_value"] = all_df["value"] * all_df["repeat_weight"]

    grouped = all_df.groupby(["snpA", "snpB"])
    #print("func_combine_interaction_rankings_fast -- grouped ----- ", all_df.head(), all_df.shape)

    weighted_sum = grouped["weighted_value"].sum()
    weight_sum   = grouped["repeat_weight"].sum()
    n_present    = grouped.size()

    result = pd.DataFrame({
        "final_interaction": weighted_sum / weight_sum,   # renormalized over repeats where this pair was present
        "n_repeats_present": n_present,
        "stability_fraction": n_present / nb_repeats,
    }).reset_index()

    # remove zero interaction pair
    result = result[result["final_interaction"] != 0]
    

    return result.sort_values("final_interaction", ascending=False)

# --------------------------------


# --------------------------------
# - Thousands of interaction is generated, need to define a cutoff
# - This script aims to select top interaction SNPs 
#   and SNPs that are interacting with the top SHAP value

def func_select_top_interaction(selected_shap_snp_ranking, final_interaction_ranking, ntop_inter1= ntop_inter1, ntop_inter2=ntop_inter2):
    # - ntop_inter1 : select the top interacting SNPs
    # - ntop_inter2 : select the SNPs that are interacting with the top features with SHAP values 

    # ------  Select top interating SNPs
    top_shap_long = final_interaction_ranking.head(ntop_inter1)

    # ------  Select the SNPs that are interacting with the top SHAP values' SNPs
    top_shap_features = selected_shap_snp_ranking["feature"]

    x = list(range(1, len(final_interaction_ranking) + 1))
    y = final_interaction_ranking["final_interaction"].values

    # Detect elbow
    knee = KneeLocator(x, y, curve="convex", direction="decreasing")

    # Interaction value at the elbow
    threshold_value = y[knee.knee - 1]   # -1 for 0-indexing, same as your original
    #print(f"Interaction threshold detected by elbow: {threshold_value:.6f}")

    # - Select SNPs with strong interaction, where at least one of the
    # SNPs are located in the selected top features
    top_shap_by_inter = final_interaction_ranking[
        final_interaction_ranking["snpA"].isin(top_shap_features) |
        final_interaction_ranking["snpB"].isin(top_shap_features)
    ]

    top_shap_by_inter_sel = top_shap_by_inter[top_shap_by_inter["final_interaction"] >= threshold_value]
    # The selection by elbow might give too many interaction
    # only select the top X interaction
    top_shap_by_inter_sel = top_shap_by_inter_sel.head(ntop_inter2)

    # ----- Combine the selected shap_interaction variables
    final_inter_shap = pd.concat([top_shap_long, top_shap_by_inter_sel], ignore_index=True)
    final_inter_shap = final_inter_shap.drop_duplicates(subset=["snpA", "snpB"])
    final_inter_shap = final_inter_shap.sort_values("final_interaction", ascending=False).reset_index(drop=True)

    return final_inter_shap, threshold_value
# --------------------------------


# -------------------------------- 
# - How to trust the weights from the meta-models
def bootstrap_coef_stability(preds_models_scaled, Y_val_np, C, n_boot=100, random_state=42):
    rng = np.random.RandomState(random_state)
    n = len(Y_val_np)
    all_coefs = []

    for b in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        clf = LogisticRegression(C=C, l1_ratio=0, solver="lbfgs", max_iter=1000)
        clf.fit(preds_models_scaled[idx], Y_val_np[idx])
        all_coefs.append(clf.coef_[0])

    all_coefs = np.array(all_coefs)   # (n_boot, nb_repeats)
    return all_coefs
# -------------------------------- 


# --------------------------------
def func_sort_variant_name(uniq_all_top_SNPs) :
    split_name = [] 
    for snp_i in uniq_all_top_SNPs :
        snp_split = snp_i.split(":")
        chr_ = snp_split[0]
        pos_ = snp_split[1]
        ref_ = snp_split[2]
        alt_ = snp_split[3]

        split_name.append([chr_, pos_, ref_, alt_])
    split_name = np.array(split_name)

    # - Sort chromosome, make sure the values are integer
    int_chr_plink = split_name[:,0].astype(int)
    int_pos_plink = split_name[:,1].astype(int)
    sort_vcf_indx = np.lexsort((int_pos_plink, int_chr_plink))

    sorted_name = split_name[sort_vcf_indx]

    # Join back the name
    sorted_name_snps = [":".join(row) for row in sorted_name]

    return sorted_name_snps
# --------------------------------


# -------------------------------- 
def func_tovcf_format(sorted_all_top_SNPs):
    # The input is a list of SNPs in chr:pos:ref:al that need to be converted to vcf format

    # - Create a VCF file to be read by VEP
    # - info : VEP does not accept rsid offline, use chr:pos:ref:alt
    # - Remove the suffix at the end of each snps id
    # - To avoid failure on VEP, create a structured array to sort by chromosome then position

    all_top_snps_tovcf = []
    for snp_i in sorted_all_top_SNPs :
        snp_split = snp_i.split(":")
        chr_ = snp_split[0]
        pos_ = snp_split[1]
        ref_ = snp_split[2]
        alt_ = snp_split[3]
        all_top_snps_tovcf.append([chr_, pos_, snp_i, ref_, alt_, ".", ".", "."])
        # 3rd column is the name of the variants
    # - Store into file, VEP (to annotate to gene) take fil as input
    # - compress the VCF
    pd_top_snps_tovcf = pd.DataFrame(all_top_snps_tovcf)
    return pd_top_snps_tovcf
# -------------------------------- 


# --------------------------------
def func_map_to_ensg(annot_variants,  ensg_col) :
    # This function map the SNPs identification to the ENSG or ENSP
    #len_snps = len(uniq_all_top_SNPs)
    #variants_ensg_ = np.concatenate(( np.array(uniq_all_top_SNPs).reshape(( len_snps,1 ))
    #                                ,annot_variants[:,ensg_col].reshape(( len_snps,1 ))  
    #                                ), axis=1)

    
    #snp_to_ensg = dict(variants_ensg_)

    annot_df = pd.DataFrame(annot_variants)
    snp_to_ensg = dict(zip(annot_df[0], annot_df[ensg_col]))
    return snp_to_ensg
# --------------------------------


# --------------------------------
def func_map_to_ensp(annot_variants, column_extra_annot_split):
    #  -- Map to dictionary -- SNPs to ENSP if it does exist
    # Not all variants are assigned to ENSP
    # Those will be assigned to the Nan
    len_snps = len(uniq_all_top_SNPs)
    
    all_annotated_ensp = []
    for i, extra in enumerate(column_extra_annot_split) :
        extra = np.array(extra)
        
        indx_feature_gene  = np.where(np.char.startswith(extra, "ENSP="))[0]
        if np.shape(indx_feature_gene)[0] > 0 :
            gene_ = np.char.split( extra[indx_feature_gene],"=" )[0][1]
            all_annotated_ensp.append([i, gene_])
            
        else :
            all_annotated_ensp.append([i, np.nan])
    all_annotated_ensp = np.array(all_annotated_ensp)

    # There are few SNPs, no gene ID, no gene name 

    # -- Map to dictionary
    annot_df = pd.DataFrame(annot_variants)
    snp_to_ensp = dict(zip(annot_df[0], all_annotated_ensp[:,1]))
    return snp_to_ensp
# --------------------------------


# --------------------------------
def func_map_to_swissprot(annot_variants, column_extra_annot_split):
    # In the network visualization, SWISSPROT gene name is easier to read
    # Eventhough, ENSG is more stable annotation

    all_annotated_swissprot  = []
    for i, extra in enumerate(column_extra_annot_split) :
        extra = np.array(extra)
        
        indx_feature_gene  = np.where(np.char.startswith(extra, "SWISSPROT="))[0]
        if np.shape(indx_feature_gene)[0] > 0 :
            gene_ = np.char.split( extra[indx_feature_gene],"=" )[0][1]
            
            all_annotated_swissprot.append([i, gene_.split("_")[0]]) #gene_.split("_")[0]
            
        else :
            all_annotated_swissprot.append([i, np.nan])
    all_annotated_swissprot = np.array(all_annotated_swissprot)

    # -- Map to dictionary
    annot_df = pd.DataFrame(annot_variants)
    all_annotated_swissprot = np.array(all_annotated_swissprot)
    snp_to_swissprot = dict(zip(annot_df[0], all_annotated_swissprot[:,1]))

    return snp_to_swissprot
# --------------------------------


def func_map_to_symbol(annot_variants, column_extra_annot_split):
    all_annotated_gene_name  = []
    for i, extra in enumerate(column_extra_annot_split) :
        extra = np.array(extra)
        
        indx_feature_gene  = np.where(np.char.startswith(extra, "SYMBOL="))[0]
        if np.shape(indx_feature_gene)[0] > 0 :
            gene_ = np.char.split( extra[indx_feature_gene],"=" )[0][1]
            all_annotated_gene_name.append([i, gene_])          
        else :              
            gene_ = annot_variants[i,3] # This is the gene ID
            all_annotated_gene_name.append([i, gene_])
    all_annotated_gene_name = np.array(all_annotated_gene_name)

    # -- Map to dictionary
    annot_df = pd.DataFrame(annot_variants)
    all_annotated_gene_name = np.array(all_annotated_gene_name)
    snp_to_symbol = dict(zip(annot_df[0], all_annotated_gene_name[:,1]))

    return snp_to_symbol
# --------------------------------
# Remove into the list correlated SNPs within the same ENSG

# This functiom compute correlation between SNPs
def func_correlation_snps(X_set, Y, feature_subset, corr_threshold=0.2, random_state=42):
    # - Before computing correlation, create a balance data
    idx_case = np.where(Y.values == 1)[0]
    idx_control = np.where(Y.values == 0)[0]
    minlen_group = min(len(idx_case), len(idx_control))
 
    rng = np.random.RandomState(random_state)

    case_sample    = rng.choice(idx_case, size=minlen_group, replace=False)
    control_sample = rng.choice(idx_control, size=minlen_group, replace=False)
    balanced_idx   = np.concatenate([case_sample, control_sample])

    X_balanced = X_set.iloc[balanced_idx][feature_subset]
    print("np.shape(X_balanced)----", np.shape(X_balanced))
    # - Compute correlation
    corr_matrix_spearman =  X_balanced.corr(method='spearman')

    # - Compute absolute correlation, and highlight SNPs with corr > 0.2
    abs_corr = corr_matrix_spearman.abs()

    # Mask diagonal (self-correlation)
    np.fill_diagonal(abs_corr.values, 0)

    # Identify SNP pairs with high correlation
    high_corr_pairs = np.where(abs_corr > corr_threshold )
    high_corr_df = pd.DataFrame({
        "snpA": corr_matrix_spearman.index[high_corr_pairs[0]],
        "snpB": corr_matrix_spearman.columns[high_corr_pairs[1]],
        "correlation": corr_matrix_spearman.values[high_corr_pairs]
    })

    # Remove duplicates (since matrix is symmetric)
    high_corr_df = high_corr_df[high_corr_df["snpA"] < high_corr_df["snpB"]] 

    return high_corr_df



# --- Remove correlated SNPs within the same gene 

def func_remove_correlation(uniq_SNPs_SHAP, select_corr_snps):
    # -1- Keep only pairs with same gene ID
    same_gene_corr = select_corr_snps.dropna(subset=["Gene1", "Gene2"])
    same_gene_corr = same_gene_corr[same_gene_corr["Gene1"] == same_gene_corr["Gene2"]].copy()

    # - 2- Define removal rule
    def choose_snp_to_remove(row, uniq_SNPs_SHAP, n=100):
        if row["snpA"] not in uniq_SNPs_SHAP[:n]:
            return row["snpA"]   # drop SNP1
        else:
            return row["snpB"]   # drop SNP2
        
    same_gene_corr["SNP_to_remove"] = same_gene_corr.apply(
        lambda r: choose_snp_to_remove(r, uniq_SNPs_SHAP, n=100), axis=1 )

    # -4- Collect SNPs to remove
    snps_to_remove = set(same_gene_corr["SNP_to_remove"].tolist())

    return snps_to_remove


def func_map_interaction_togene(threshold_interaction, interaction_ranking, mapped_final_snps, 
                                snp_to_ensg, snp_to_enst, snp_to_ensp) :
    inter_final = interaction_ranking[ interaction_ranking["snpA"].isin(mapped_final_snps["SNP"]) &\
                                                interaction_ranking["snpB"].isin(mapped_final_snps["SNP"]) ]
    inter_final = inter_final[inter_final["final_interaction"] >= threshold_interaction] 

    inter_final = inter_final.copy()  # avoid SettingWithCopyWarning

    inter_final["ENSG_A"] = inter_final["snpA"].map(lambda s: snp_to_ensg.get(s, np.nan))
    inter_final["ENSG_B"] = inter_final["snpB"].map(lambda s: snp_to_ensg.get(s, np.nan))

    inter_final["ENST_A"] = inter_final["snpA"].map(lambda s: snp_to_enst.get(s, np.nan))
    inter_final["ENST_B"] = inter_final["snpB"].map(lambda s: snp_to_enst.get(s, np.nan))

    inter_final["ENSP_A"] = inter_final["snpA"].map(lambda s: snp_to_ensp.get(s, np.nan))
    inter_final["ENSP_B"] = inter_final["snpB"].map(lambda s: snp_to_ensp.get(s, np.nan))

    ranked_inter_final = ( pd.DataFrame(inter_final)
                    .sort_values("final_interaction", ascending=False)
                    .reset_index(drop=True) ) 

    return ranked_inter_final
# --------------------------------


# --------------------------------
# - load and prepare data
genotype_data = pd.read_csv(genotype_path, sep=r"\s+")
X = genotype_data.iloc[:, 6:]
Y = (genotype_data["PHENOTYPE"] - 1).astype(int)

# drop SNPs that is missing more than 0.98 % in one of the population
cases    = X[Y == 1]
controls = X[Y == 0]

# --- Store the columns names, and convert to idx
# The convertion netween different data format will delete the features name
# Important for the step2
feature_names = X.columns.tolist()
feature_names = [strip_snp_suffix(f) for f in feature_names]
feature_to_idx = {f: i for i, f in enumerate(feature_names)}

# - Fix the column name in X
X.columns = feature_names 
# --------------------------------


# --------------------------------
# - split data similarly like during the training
# to take the test set
# - This test set has been left out during all training
X_train_temp, X_test_temp, Y_train_temp, Y_test_temp = train_test_split(X, Y, test_size=test_size, stratify=Y, 
                                                              random_state=rd_state_test)

# - Here I will split the TEST SET
# for training the metadata and testing the metadata
# 15% for the training and 5% for the test
X_val, X_test, Y_val, Y_test = train_test_split(X_test_temp, Y_test_temp, test_size=0.5,
                                                    stratify=Y_test_temp, random_state=rd_state_test)                                                              

if device == "cuda:0":
    X_val, Y_val = cp.array(X_val.values), cp.array(Y_val.values)
    X_test, Y_test   = cp.array(X_test.values), cp.array(Y_test.values)

# convert to numpy
Y_val_np  = Y_val.copy()
X_val_np  = X_val.copy()

Y_test_np = Y_test.copy()
X_test_np = X_test.copy()
if device == "cuda:0":
    Y_val_np = cp.asnumpy(Y_val)
    X_val_np = cp.asnumpy(X_val)

    Y_test_np = cp.asnumpy(Y_test)
    X_test_np = cp.asnumpy(X_test)

# --------------------------------

# --------------------------------
# - Then I would like to train a meta model, using the validation set
# The coefficent of each model within the meta-model will be used as weight for
# final feature selection across all the models
# --------------------------------


# --------------------------------
# - TRAINING OF THE META-MODELS
# - use the validation set

preds_models_meta = []
models_list  = []
features_per_model_list = [] # here contain the index of SNPs
feat_name_permodel_list = [] # here contain the top name of snps
best_threshold_list = []
explainer_list      = []
val_shap_list = []

for model_i in np.arange(nb_repeats):
    model, top_feature_names, best_threshold = func_load_model(f"repeat_{repeat_ids[model_i]}")

    # -- etract the index of features used in the model
    common_features_idx = [feature_to_idx[f] for f in top_feature_names  if f in feature_to_idx]

    if device == "cuda:0" :
        X_val_set = X_val[:, common_features_idx ]

    elif device == "cpu" :
        X_val_set = X_val.values[:, commin_features_idx]

    # Compute the SHAP
    # - generate an explaner used for diagnisticing the test set
    # - and diagnotiscing the validateion set used for sanity check
    explainer = shap.TreeExplainer(model)

    shap_features = func_compute_shap( top_feature_names, explainer, X_val_set)
    val_shap_list.append(shap_features)
    feat_name_permodel_list.append(top_feature_names)

    # -- Store the model and features and shap values of the models
    models_list.append(model)
    features_per_model_list.append(common_features_idx)
    best_threshold_list.append(best_threshold)
    explainer_list.append(explainer)

    dval = xgb.QuantileDMatrix(X_val_set, label=Y_val)

    preds_model_1 = model.predict(dval)
    if device == "cuda:0":
        preds_model_1 = cp.asnumpy(preds_model_1)

    preds_models_meta.append(preds_model_1)

preds_models_meta = np.column_stack(preds_models_meta)
# -----------------------------------------


# --------------------------------
# - NOW USE THE TEST SET 
# - Aim here is to compute the std error average across all the models, on the test set
# - the  5% of the test set
# It will be compared to the meta-model

# Also compute SHAP on the Xtest to capture SNP-SNP interation

# - Store the model once
mcc_list = []
auc_list = []

features_SHAP_list = []  # here contains the featutes and SHAP values 
preds_models_test  = []
shap_list          = []
shap_interaction_df_list = []

for model_i in np.arange( nb_repeats):

    # -- Performance on individual model computed from the 5% test set
    # -  Store the prediction on the test set
    # -  This prediction will be used to evaluate the meta data

    model = models_list[model_i]
    common_features_idx  = features_per_model_list[model_i]
    subset_feature_names = feat_name_permodel_list[model_i]
    best_threshold = best_threshold_list[model_i]
 
    if device == "cuda:0" :
        X_test_set = X_test[:, common_features_idx ]
    elif devive == "cpu" :
        X_test_set = X_test.values[:, common_features_idx]

    dtest = xgb.QuantileDMatrix(X_test_set, label=Y_test)
    pred  = model.predict(dtest)

    # - Compute SHAP of features in the test set
    shap_features = func_compute_shap(subset_feature_names, explainer_list[model_i], X_test_set)
    shap_list.append(shap_features)
    print(f"SHAP computed on repeat_{repeat_ids[model_i]} ...", flush=True)

    in_time_r = time.time()
    ## - Compute shap interaction (this takes times)
    X_test_set_np = to_np(X_test_set).astype(np.float32)
    mean_abs_inter = func_compute_shap_interactions(model, X_test_set_np, gpu_interaction=gpu_interaction, batch_size=inter_batch_size)
    print(f"SHAP INTERACTION computed on repeat_{repeat_ids[model_i]} ...", flush=True)
    # Each SHAP interaction si computed between 13 to 3.80 mins
    
    records = []
    for a, b in itertools.combinations(range(len(subset_feature_names)), 2):
        snpA, snpB = sorted([subset_feature_names[a], subset_feature_names[b]])
        records.append({"snpA": snpA, "snpB": snpB, "value": mean_abs_inter[a, b]})

    df_this_repeat = pd.DataFrame(records)
    df_this_repeat["repeat_id"] = model_i   # tag which repeat this came from, weight attached later
    df_this_repeat = df_this_repeat[df_this_repeat["value"] != 0]
    shap_interaction_df_list.append(df_this_repeat)

    # I can do global contribution of SNPs, by taking the absolute mean
    # then only pick SNPs with high contribution
    
    # -  Compute performance MCC and AUC
    y_pred_label = (pred > best_threshold ).astype(int)
    # convert to numpy
    if device == "cuda:0":
        pred = cp.asnumpy(pred)
        y_pred_label = cp.asnumpy(y_pred_label )

    mcc = matthews_corrcoef(Y_test_np, y_pred_label )
    mcc_list.append(mcc)

    auc = roc_auc_score(Y_test_np, pred)

    preds_models_test.append(pred)
    auc_list.append(auc)

# combine the predictions of the base models
preds_models_test = np.column_stack(preds_models_test)

# - MCC on test set
mean_mcc = np.mean(mcc_list)
std_mcc  = np.std(mcc_list)
sem_mcc  = std_mcc / np.sqrt(nb_repeats)
# standard error of the mean

# - AUC on test set
mean_auc = np.mean(auc_list)
std_auc  = np.std(auc_list)
sem_auc  = std_auc / np.sqrt(nb_repeats)
# --------------------------------



# --------------------------------


#----------------------------------------
# -- TRAIN META-MODELS USING LOGISTIC REGRESSION 

# - training using validation set
# - do a cross validation to chose the best C
val_mean_mcc_list   = []
val_std_mcc_list    = []

# -- Before fitting the meta-models, I need to
# standardize the base models, because each use different hyper-parameters and fraction of SNPs
# if not done, there will be an artifact for each model where logistic regression will "compensate"
# for a more spread-out column by giving it a smaller coefficient
scaler = StandardScaler()
preds_models_scaled      = scaler.fit_transform(preds_models_meta)  # prediction performed on val, used to train the meta
preds_models_test_scaled = scaler.transform(preds_models_test)


#C_list= np.linspace(-0.5, 1.5, 100)
C_list= np.logspace(-5, 0.25, 130)
skf = StratifiedKFold(n_splits=4, shuffle=True, random_state=42)

for C in C_list:
    mcc_scores = []
    # cv
    for train_idx, val_idx in skf.split(preds_models_scaled, Y_val_np):
        X_tr, X_va = preds_models_scaled[train_idx], preds_models_scaled[val_idx]
        y_tr, y_va = Y_val_np[train_idx], Y_val_np[val_idx]

        # - For each cv, tain a regeression model
        rg_model = LogisticRegression(C=C, 
                                l1_ratio=0, # for penalty="l2",
                                solver="lbfgs", 
                                max_iter=1000)
        rg_model.fit(X_tr, y_tr)

        # - Compute the performance
        preds = rg_model.predict_proba(X_va)[:,1]
        _, mcc = func_find_mcc_threshold(preds, y_va)

        mcc_scores.append(mcc)

    mean_mcc_C = np.mean(mcc_scores)
    val_mean_mcc_list.append(mean_mcc_C)
    val_std_mcc_list.append(np.std(mcc_scores))

val_best_mcc = val_mean_mcc_list[np.argmax(val_mean_mcc_list)]
val_best_C   = C_list[np.argmax(val_mean_mcc_list)]

meta_model = LogisticRegression(C=val_best_C, l1_ratio=0, solver="lbfgs", max_iter=1000)
meta_model.fit(preds_models_scaled, Y_val_np)
#----------------------------------------


#----------------------------------------
# --- plot perfornace vs C
# when select C is small, the mcc is lower but the coefficient between the models are similar
# when is is larger, the mcc improves a bit, and some models contibute more to the prediction than the other

plt.figure(figsize=(8, 6))
plt.plot(C_list, val_mean_mcc_list, marker='o', markersize= 3)

plt.xscale("log")
ax = plt.gca()
# sparsing the xticks for better visualization
ax.xaxis.set_major_locator(mticker.LogLocator(base=10, subs=[1,2,3,5]))

y_lo = min(val_mean_mcc_list) - 0.03
y_hi = max(val_mean_mcc_list) + 0.03
plt.ylim(y_lo, y_hi)

plt.xlabel("C (regularization strength)")
plt.ylabel("Mean MCC")
plt.title("Meta-model tuning (Logistic Regression)")
plt.grid(True)

plt.savefig(f"{result_out}/meta_model{nb_repeats}_C_tuning_reg.png", dpi=150)
plt.close()
#----------------------------------------

#----------------------------------------
# Evaluate the meta data on new dataset --> here the 5% test set

# Make prediction using meta-model
y_meta_pred = meta_model.predict_proba(preds_models_test_scaled)[:, 1]
# -----------------------


# -----------------------
# - Performance of the meta-model on the untouched 5% test set
# - Chose the threshold for mcc for the meta models

# threshold chosen out-of-fold, not in-sample
best_t_meta, best_mcc_oof_meta = func_find_best_threshold_cv(preds_models_scaled, Y_val_np, val_best_C)

# - Perform prediction
y_pred_label = (y_meta_pred > best_t_meta).astype(int)

# convert to numpy
if device == "cuda:0":
    y_pred_label = cp.asnumpy(y_pred_label )

meta_mcc = matthews_corrcoef(Y_test_np, y_pred_label)
meta_auc = roc_auc_score(Y_test_np, y_meta_pred)

meta_mcc_prediction_balance = np.mean(y_pred_label) 

# -----------------------


# -----------------------
# - I - Combine per-repeat SNP rankings, weighted by |meta-model coefficient|
# Meta-model coefficients (interpretability check) meta_model.coef_
repeat_weights = np.abs(meta_model.coef_[0])

full_shap_snp_ranking, selected_shap_snp_ranking = func_combine_snp_rankings( nb_repeats, repeat_weights, shap_list)
print(f"SHAP concensus across models are computed ...", flush=True)
# do i need to save it to file, no for now

# - II - Compute the averaging of SHAP interaction values across SNPs
interaction_ranking = func_combine_interaction_rankings(nb_repeats, repeat_weights, shap_interaction_df_list)
#final_interaction_ranking.to_csv(f"{param_out}/final_consensus_interaction_ranking.tsv", sep="\t", index=False)
#print(f"SHAP INTERACTION concensus across models are computed ...", flush=True)

# - III - Select the top  interactions
top_shap_by_inter_sel, threshold_interaction = func_select_top_interaction(selected_shap_snp_ranking, interaction_ranking)

# - IV - Merge all the SNPs collected from top SHAP values and SHAP interaction values
# - Map the SNPs to genes, using VEP
uniq_SNPs_SHAP = pd.unique(selected_shap_snp_ranking["feature"].values.ravel())
uniq_SNPs_interaction =  pd.unique(top_shap_by_inter_sel[["snpA", "snpB"]].values.ravel())
uniq_all_top_SNPs = list( set(uniq_SNPs_SHAP)| set(uniq_SNPs_interaction) )

# Here, the snsps name is sorted (by chr and pos) to avoid crash when given to VEP
sorted_all_top_SNPs = func_sort_variant_name(uniq_all_top_SNPs)

# - Here, the SNPs are prepared to be read by VEP
# - VEP is used to assign SNPs to genes name
pd_top_snps_tovcf = func_tovcf_format(sorted_all_top_SNPs)
pd_top_snps_tovcf.to_csv(f"{snps_dir}/top_snps_feature_select.vcf",  sep='\t',index=False, header=False)

cmd = f"""
cd {snps_dir} && \
bgzip -f top_snps_feature_select.vcf
"""

vep_input = os.path.join(snps_dir, "top_snps_feature_select.vcf.gz")
vep_output = os.path.join(snps_dir, "top_snps_feature_select_vep.txt")

cmd_vep = [
    "vep",
    "--offline",
    "--cache",
    "--dir", vep_dir,
    "--species", "homo_sapiens",
    "--assembly", assembly,   
    "--fasta", fasta_file,
    "--input_file", vep_input,
    "--output_file", vep_output,
    "--fork", "12",
    "--everything",
    "--pick",
    "--force_overwrite"
]

subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
subprocess.run(cmd_vep, shell=False, check=True)

# - Upload the file generated by VEP
annot_variants = pd.read_csv(vep_output, header = None, sep = "\t", comment="#")
annot_variants = annot_variants.to_numpy()


# Annotation are located in the 13th column
# Gene name is unstable and not all variant are assigned to gene name
# Take ENSG, it is more stable
column_extra_annot = annot_variants[:,13].astype(str)
column_extra_annot_split = np.char.split(column_extra_annot, ";")

#  -- Map to dictionary -- SNPs to GENE/ENSG
ensg_col = 3
enst_col = 4 

#snp_to_swissprot = func_map_to_swissprot(annot_variants, column_extra_annot_split)
snp_to_symbol = func_map_to_symbol(annot_variants, column_extra_annot_split)
snp_to_ensg = func_map_to_ensg(annot_variants,  ensg_col)
snp_to_enst = func_map_to_ensg(annot_variants,  enst_col)
snp_to_ensp = func_map_to_ensp(annot_variants, column_extra_annot_split)
snp_to_shap = dict( zip(full_shap_snp_ranking["feature"], full_shap_snp_ranking["final_importance"]) )
snp_to_stability = dict( zip(full_shap_snp_ranking["feature"], full_shap_snp_ranking["stability_fraction"]) )

# -----------------------------------------
# - Compute correlation between SNPs
# Correlated SNPs within the same genes are discarded
# Here is the best occasion to compute correlation as the feature sets as been already schrinked
##concat_feature_subset = np.unique(np.concat(feat_name_permodel_list))
##correlation_snps = func_correlation_snps(X, Y, concat_feature_subset, corr_threshold=0.2, random_state=42)

corr_step = time.time()
correlation_snps = func_correlation_snps(X, Y, sorted(uniq_all_top_SNPs), 
                                         corr_threshold=0.2, random_state=42)
print(f"Correlation between the {len(uniq_all_top_SNPs)} \
      top SNPs computed in {time.time() - corr_step:.1f}s", flush=True)
# -----------------------------------------


# -- Remove into the list correlated SNPs within the same ENSG
# Annotate correlation SNPs to their gene
select_corr_idx  = correlation_snps["snpA"].isin(uniq_all_top_SNPs) | correlation_snps["snpB"].isin(uniq_all_top_SNPs)
select_corr_snps = correlation_snps[select_corr_idx].copy()
select_corr_snps["Gene1"] = select_corr_snps["snpA"].map(snp_to_ensg)
select_corr_snps["Gene2"] = select_corr_snps["snpB"].map(snp_to_ensg)

# uniq_SNPs_SHAP is used as the SNPs in the top SHAP is priorizited 
# when picking only one SNPs in the correlated interaction
corr_snps_to_remove = func_remove_correlation(uniq_SNPs_SHAP, select_corr_snps)

corr_indices_to_remove = [uniq_all_top_SNPs.index(snp) 
                     for snp in corr_snps_to_remove 
                     if snp in uniq_all_top_SNPs]

top_snps_nocorrelation = np.delete(uniq_all_top_SNPs, corr_indices_to_remove)

# - Build a unique mapping of snps, genes name and Shap value
mapped_final_snps = [ (snp,
                        snp_to_symbol.get(snp, np.nan),
                        snp_to_ensg.get(snp, np.nan), 
                        snp_to_enst.get(snp, np.nan),
                        snp_to_ensp.get(snp, np.nan),
                        snp_to_shap.get(snp, np.nan), 
                        snp_to_stability.get(snp, np.nan) ) 
                        for snp in top_snps_nocorrelation]
columns_mapped = ["SNP", "SYMBOL", "ENSG", "ENST", "ENSP", "SHAP_value", "stability_fraction"]
mapped_final_snps = ( pd.DataFrame(mapped_final_snps, columns=columns_mapped)
                    .sort_values("SHAP_value", ascending=False)
                    .reset_index(drop=True) )                      

# - Write the top SNPs -- correlated within the same genes are remove
# Not alll SNPs might involve in the interaction

mapped_final_snps.to_csv(f"{snps_dir}/final_consensus_shapvalue.tsv", 
                                        sep="\t", index=False)

# - Write the top interaction
mapped_inter_final = func_map_interaction_togene(threshold_interaction, interaction_ranking, mapped_final_snps, 
                                snp_to_ensg, snp_to_enst, snp_to_ensp) 

mapped_inter_final.to_csv(f"{snps_dir}/final_consensus_interaction.tsv", 
                                        sep="\t", index=False)

## Compute the SHAP performed on the val set, for sanity checking
#val_snp_ranking_diagnostic = func_combine_snp_rankings(nb_repeats, repeat_weights, kappa_val_per_model_list)
#val_snp_ranking_diagnostic.to_csv(f"{param_out}/diagnostic_val_snp_ranking.tsv", sep="\t", index=False)

# - How much can I trust the weight from the meta-model to be applied on the SNPs importance
all_coefs = bootstrap_coef_stability(preds_models_scaled, Y_val_np, val_best_C, n_boot=200, random_state=rd_state_test)

# -----------------------

# --------------------------------
# EXPORT: self-contained folder for prediction on new data
# SHAP and interaction values are NOT needed for prediction

meta_bundle = {
    "disease_name": disease_name,
    "date_trained": str(datetime.now()),
    "trained_on": str(genotype_path),
 
    "nb_repeats": int(nb_repeats),
    "repeat_ids": [int(x) for x in repeat_ids],
 
    # relative to the export folder
    #"model_files": [f"repeat_{i+1}/{model_file}" for i in range(nb_repeats)],
    "model_files": [f"repeat_{r}/{model_file}" for r in repeat_ids],

    # ORDERED feature names per repeat, exactly as each model expects them.
    # Saved rather than recomputed: rebuilding this later from the .npy plus
    # a membership filter against a DIFFERENT cohort would silently change
    # the order or the length.
    "features_per_model": [[str(x) for x in f] for f in feat_name_permodel_list],
 
    # per-repeat thresholds, if a single base model is applied on its own
    "thresholds_per_model": [float(t) for t in best_threshold_list],
 
    # standardisation of base-model predictions, fitted on the validation set
    "scaler_mean":  scaler.mean_.tolist(),
    "scaler_scale": scaler.scale_.tolist(),
 
    # fitted aggregation model
    "meta_coef":      meta_model.coef_[0].tolist(),
    "meta_intercept": float(meta_model.intercept_[0]),
    "meta_C":         float(val_best_C),
    "meta_threshold": float(best_t_meta),
 
    # held-out performance, for the user's reference
    "performance": {
        "base_mean_mcc": float(mean_mcc),
        "base_mean_auc": float(mean_auc),
        "meta_mcc": float(meta_mcc),
        "meta_auc": float(meta_auc)},
 
    # how reliable the per-repeat weights are
    "coef_bootstrap_mean": all_coefs.mean(axis=0).tolist(),
    "coef_bootstrap_std":  all_coefs.std(axis=0).tolist()
}

def json_np(o):
    if isinstance(o, (np.integer,)):   return int(o)
    if isinstance(o, (np.floating,)):  return float(o)
    if isinstance(o, (np.ndarray,)):   return o.tolist()
    if isinstance(o, (np.str_,)):      return str(o)
    raise TypeError(f"Not serialisable: {type(o)}")

with open(f"{result_out}/meta_model_bundle.json", "w") as f:
    json.dump(meta_bundle, f, indent=2, default=json_np)

joblib.dump({"scaler": scaler, "meta_model": meta_model},
            f"{result_out}/meta_model_sklearn.joblib")

all_required = sorted({str(f) for feats in feat_name_permodel_list for f in feats})
pd.Series(all_required, name="SNP").to_csv(
    f"{result_out}/required_variants.txt", index=False, header=True)

# --------------------------------


# --------------------------------
# - Record time 
out_time = time.time()
elapsed = out_time - in_time
hours = elapsed // 3600
minutes = (elapsed % 3600) / 60
# --------------------------------

# --------------------------------
# --- Write log file
def write_log():
    with open(f"{result_out}/outlog_{nb_repeats}mcc_mean.log", "w") as f:

        f.write("==== RUN INFO ====\n")
        f.write(f"Date: {datetime.now()}\n\n")

        f.write("\n---- DATA INFO ----\n")
        f.write(f"Nb validation set : S{len(Y_val)}\n")
        f.write(f"Nb test set : S{len(Y_test)}\n")

        f.write("---- INPUT PARAMETERS ----\n")
        f.write(f"These following models were used :\n")

        for r_id in repeat_ids:
            f.write(f"{param_out}/repeat_{r_id}/{model_file} \n")
        f.write("\n\n")

        f.write(f"Number of SNPs (features) per model: {[len(idx) for idx in features_per_model_list] }")

        json.dump(config, f, indent=4)
        f.write("\n\n")

        f.write(f" ---- MODEL PERFORMANCE OF THE BASE MODELS ----\n")
        f.write(f" ---- performed on the test set ----\n")
        f.write(f"MCC: {mean_mcc:.5f} +- {sem_mcc:.5f} \n")
        f.write(f"AUC: {mean_auc:.5f} +- {sem_auc:.5f} \n")
        f.write("\n\n")


        f.write("---- META-MODEL ----\n")
        f.write(f"Meta-models weights {np.abs(meta_model.coef_[0])} \n\n")
        f.write(f"How trustwotry is the weights from the meta-models applied to the SNPs importance\n")
        f.write(f"a large std relative to the mean for any repeat = unstable weight\n")

        f.write(f"Coefficient mean per repeat: {all_coefs.mean(axis=0)}\n")
        f.write(f"Coefficient std per repeat: {all_coefs.std(axis=0)}\n\n")

        f.write(f"Best C: {val_best_C}\n")
    
        f.write(f"Meta-model - test MCC: {meta_mcc:.5f}\n")
        f.write(f"Meta-model - test AUC: {meta_auc:.5f}\n\n")

        f.write("---- SNP CONSENSUS RANKING \n")
        f.write(f"Total SNPs: {len(uniq_all_top_SNPs)}\n")
        f.write(f"Saved to: {snps_dir}/final_consensus_snp_ranking.tsv\n")
        #f.write(f"Val-based diagnostic comparison saved to: {param_out}/diagnostic_val_snp_ranking.tsv\n\n")


        f.write(f"---- Total tuning time : {hours:.0f}h {minutes:.2f}m")

write_log()



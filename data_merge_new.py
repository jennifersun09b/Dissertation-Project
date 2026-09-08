"""Build the baseline and complete-case longitudinal dissertation datasets.

Outputs:
  baseline.csv     -> predictive_model.py
  longitudinal.csv -> EDA_longitudinal.py

Input and output paths can be overridden with BASELINE_SOURCE_PATH,
OUTCOME_SOURCE_PATH, LIFESTYLE_SOURCE_PATH, and OUTPUT_DIR.
"""

# Data Preprocessing & Merging
# ----------------------------

# Import library
# --------------

import os
from pathlib import Path

import pandas as pd
import numpy as np

# Loading datasets
# ----------------

BASELINE_SOURCE_PATH = Path(os.environ.get(
    "BASELINE_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/group_1_clean.tsv",
))
OUTCOME_SOURCE_PATH = Path(os.environ.get(
    "OUTCOME_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/group_1_outcomes_df_without_qc_df.tsv",
))
LIFESTYLE_SOURCE_PATH = Path(os.environ.get(
    "LIFESTYLE_SOURCE_PATH",
    "/myriadfs/projects/ICS_UKB/projects/kasia_project/causal_analysis/ukb_tabular_data_causal_analysis.tsv",
))
OUTPUT_DIR = Path(os.environ.get(
    "OUTPUT_DIR",
    "/home/rmhiund/causal_analysis/data_merge",
))

df = pd.read_csv(BASELINE_SOURCE_PATH, sep='\t', low_memory=False)
df1 = pd.read_csv(OUTCOME_SOURCE_PATH, sep='\t', low_memory=False)
df2 = pd.read_csv(LIFESTYLE_SOURCE_PATH, sep='\t', low_memory=False)

for dataset_name, dataset in {
    "baseline source": df,
    "outcome source": df1,
    "lifestyle source": df2,
}.items():
    if "eid" not in dataset.columns:
        raise ValueError(f"{dataset_name} is missing eid")
    if dataset["eid"].duplicated().any():
        raise ValueError(f"{dataset_name} contains duplicated eid values")

# Data Overview
# -------------

print(df.shape)

print(df1.shape)

print(df2.shape)

df.info()

print(df.isnull().sum()/len(df)*100)

df1.info()

print(df1.nunique())

print(df1.isnull().sum())

outcome_columns = ["eid", "def_CVD_AFTER", "def_CVD_AF_HF_AFTER"]
missing_outcomes = sorted(set(outcome_columns) - set(df1.columns))
if missing_outcomes:
    raise ValueError(f"outcome table is missing columns: {missing_outcomes}")
df1 = df1[outcome_columns]

df1.info()

print(df1.head())

print(df1.isnull().sum())

df2.info()

print(df2.columns)

# Data Merging Strategy
# ---------------------

column = []
for col in df2.columns:
    if col.endswith("-0.0"):
        column.append(col)
    if col.endswith("-2.0"):
        column.append(col)
    elif col == 'eid':
        column.append(col)

df2_i = df2[column]

print(df2_i.head())

# Sleeping Variable
# -----------------
# - 1160 (sleep duration), 1180 (chronotype), 1200 (insomnia), 1210 (snoring), 1220 (daytime dozing) included [literature review]

SLEEP_COLS = ['1160', '1170', '1180', '1190', '1200', '1210', '1220']
INSTANCES = ['0.0', '2.0']

def sleep_score_calculation(row, instance):
    """
    HSC (Healthy Sleep Score) based on:
    - Sleep duration (1160), Chronotype (1180), Insomnia (1200),
      Snoring (1210), Daytime sleepiness (1220)
    Each component weighted equally (weighted version via sensitivity analysis later).
    Higher total score means higher sleep score
    """
    score = 0
    valid_components = 0

    def col(field):
        return f'{field}-{instance}'

    # Chronotype
    if not pd.isna(row[col('1180')]):
        valid_components += 1
        if row[col('1180')] in [1, 2]:
            score += 1

    # Sleep duration
    if not pd.isna(row[col('1160')]):
        valid_components += 1
        if 7 <= row[col('1160')] <= 8:
            score += 1

    # Insomnia
    if not pd.isna(row[col('1200')]):
        valid_components += 1
        if row[col('1200')] == 1:
            score += 1

    # Snoring
    if not pd.isna(row[col('1210')]):
        valid_components += 1
        if row[col('1210')] == 2:
            score += 1

    # Daytime sleepiness
    if not pd.isna(row[col('1220')]):
        valid_components += 1
        if row[col('1220')] in [0, 1]:
            score += 1

    if valid_components == 0:
        return np.nan
    return score / valid_components * 5


def categorise_sleep(x):
    if pd.isna(x):
        return np.nan
    if x >= 4:
        return 2.0
    if x <= 2:
        return 0.0
    return 1.0


# --- Build column list across both instances ---
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in SLEEP_COLS]
# Keep only columns that actually exist in df2_i
all_cols = [c for c in all_cols if c in df2_i.columns]

sleep_i = df2_i[all_cols].replace({-1: np.nan, -3: np.nan})

# --- Calculate scores for each instance ---
for inst in INSTANCES:
    inst_cols = [f'{field}-{inst}' for field in SLEEP_COLS]
    # Skip instance entirely if none of its columns exist
    if not any(c in sleep_i.columns for c in inst_cols):
        print(f"Skipping instance {inst}: no columns found.")
        continue

    sleep_i[f'total_hsc_{inst}'] = sleep_i.apply(
        lambda row: sleep_score_calculation(row, inst), axis=1
    )
    sleep_i[f'sleep_category_{inst}'] = sleep_i[f'total_hsc_{inst}'].apply(categorise_sleep)

print(sleep_i.head())

print(sleep_i['total_hsc_2.0'].value_counts().sort_index())

print(sleep_i['sleep_category_0.0'].value_counts())

print(sleep_i['sleep_category_2.0'].value_counts())

print(sleep_i.columns)

sleep_i = sleep_i[['eid','total_hsc_0.0','sleep_category_0.0', 'total_hsc_2.0', 'sleep_category_2.0']]

print(sleep_i['sleep_category_2.0'].isnull().sum()/len(sleep_i))

both_valid = sleep_i['sleep_category_0.0'].notna() & sleep_i['sleep_category_2.0'].notna()
both_valid.sum()

print(sleep_i.head())

# Smoking variable
# ----------------

SMOKING_COLS = ['20116']
INSTANCES = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in SMOKING_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

smoking_i = df2_i[all_cols].replace({-3: np.nan})

# Rename for clarity
rename_map = {f'20116-{inst}': f'smoking_category_{inst}' for inst in INSTANCES}
smoking_i = smoking_i.rename(columns=rename_map)

# Summary
for inst in INSTANCES:
    col = f'smoking_category_{inst}'
    if col in smoking_i.columns:
        print(f"[{inst}] Null count: {smoking_i[col].isnull().sum()} ({smoking_i[col].isnull().sum()/len(smoking_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{smoking_i[col].value_counts()}\n")

print(smoking_i.head())

print(smoking_i['smoking_category_2.0'].isnull().sum()/len(smoking_i))

both_valid = smoking_i['smoking_category_0.0'].notna() & smoking_i['smoking_category_2.0'].notna()
both_valid.sum()

# Alcohol Variables
# -----------------

# - 20117 has 0.0 (never), 1.0 (previous), current (2.0)
# - 1558 has 6.0 (never), 5.0-1.0 (current but different degree), 5.0 is the least, 1.0 is the most

df2_a = df2_i[['1558-0.0','1558-2.0','20117-0.0','20117-2.0']]

print(df2_a.head())

pairs = [
    ('1558-0.0', '20117-0.0'),   # instance 0 (baseline)
    ('1558-2.0', '20117-2.0'),   # instance 2 (imaging visit)
]

for freq_col, status_col in pairs:
    print(f"\n=== {freq_col}  vs  {status_col} ===")
    print(pd.crosstab(df2_a[freq_col], df2_a[status_col], dropna=False))

# - 152 6.0 (1558) and -3.0 (20117) should be identify as Never

def combine_alcohol(freq, status):
    """Turn one frequency (1558) + one status (20117) into a code 0-4.
    Takes two values so it works for ANY timepoint."""

    # NEW rule: frequency says 'never' but status was 'prefer not to answer'.
    # We trust the frequency and call it Never.
    if freq == 6.0 and status == -3.0:
        return 0

    # 0 = never
    if status == 0.0:
        return 0

    # 1 = previous
    if status == 1.0:
        return 1

    # current drinker -> level comes from frequency
    if status == 2.0:
        if freq in [5.0, 4.0]:   # special occasions, 1-3x/month
            return 2             # low
        if freq in [3.0, 2.0]:   # 1-2x/week, 3-4x/week
            return 3             # moderate
        if freq == 1.0:          # daily or almost daily
            return 4             # high

    # everything else (missing, other odd combos) -> blank
    return np.nan


# Apply the SAME function to each timepoint.
# The lambda just pulls the two values out of the row and hands them over.
df2_i['alcohol_score_0.0'] = df2_i.apply(
    lambda row: combine_alcohol(row['1558-0.0'], row['20117-0.0']), axis=1
)
df2_i['alcohol_score_2.0'] = df2_i.apply(
    lambda row: combine_alcohol(row['1558-2.0'], row['20117-2.0']), axis=1
)

alcohol_i = df2_i[['eid','1558-0.0', '20117-0.0', '1558-2.0','20117-2.0', 'alcohol_score_0.0', 'alcohol_score_2.0']]

print(alcohol_i.drop(columns ='eid')[['alcohol_score_0.0', 'alcohol_score_2.0']].value_counts())

print(alcohol_i.head())

# Diet Variable
# -------------

def safe_sum(values):
    if all(pd.isna(v) for v in values):
        return np.nan
    return sum(v for v in values if not pd.isna(v))

DIET_COLS = ['1289', '1299', '1309', '1329', '1339', '1349', '1369', '1379', '1389']
INSTANCES = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in DIET_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

diet_i = df2_i[all_cols].replace({-1: np.nan, -3: np.nan})

def diet_score_calculation(row, instance):
    def col(field):
        return f'{field}-{instance}'

    fruit      = row[col('1309')]
    vegetable  = safe_sum([row[col('1289')], row[col('1299')]]) / 3
    fish       = safe_sum([row[col('1329')], row[col('1339')]])
    pro_meat   = row[col('1349')]
    red_meat   = safe_sum([row[col('1369')], row[col('1379')], row[col('1389')]])

    # Check if all diet components are missing
    if pd.isna(fruit) and pd.isna(vegetable) and pd.isna(fish) and pd.isna(pro_meat) and pd.isna(red_meat):
        return np.nan

    score = 0

    # Fruit + vegetable condition
    fruit_veg = np.nansum([fruit, vegetable])
    if fruit_veg >= 4.5:
        score += 1

    # Fish condition
    if pd.notna(fish) and fish >= 2:
        score += 1

    # Meat condition
    if pd.notna(pro_meat) and pd.notna(red_meat):
        if pro_meat <= 2 and red_meat <= 5:
            score += 1

    return score


for inst in INSTANCES:
    inst_cols = [f'{field}-{inst}' for field in DIET_COLS]
    if not any(c in diet_i.columns for c in inst_cols):
        print(f"Skipping instance {inst}: no columns found.")
        continue

    diet_i[f'diet_score_{inst}'] = diet_i.apply(
        lambda row: diet_score_calculation(row, inst), axis=1
    )

# Summary
for inst in INSTANCES:
    col = f'diet_score_{inst}'
    if col in diet_i.columns:
        print(f"[{inst}] Null count : {diet_i[col].isnull().sum()} ({diet_i[col].isnull().sum()/len(diet_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{diet_i[col].value_counts()}\n")

print(diet_i.head())

print(diet_i.columns)

diet_i = diet_i[['eid','diet_score_0.0', 'diet_score_2.0']]

print(diet_i.head())

both_valid = diet_i['diet_score_0.0'].notna() & diet_i['diet_score_2.0'].notna()
both_valid.sum()

# Sun Exposure Variables
# ----------------------

SUN_COLS = ['2277']
INSTANCES = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in SUN_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

sun_i = df2_i[all_cols].replace({-1: np.nan, -3: np.nan})

for inst in INSTANCES:
    col_raw = f'2277-{inst}'

    if col_raw not in sun_i.columns:
        print(f"Skipping instance {inst}: {col_raw} not found.")
        continue

    sun_i[f'sunlamp_use_{inst}'] = sun_i[col_raw].apply(
        lambda x: 1.0 if x >= 1 else (0.0 if x == 0 else np.nan)
    )

# Summary
for inst in INSTANCES:
    col = f'sunlamp_use_{inst}'
    if col in sun_i.columns:
        print(f"[{inst}] Null count : {sun_i[col].isnull().sum()} ({sun_i[col].isnull().sum()/len(sun_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{sun_i[col].value_counts()}\n")

print(sun_i.head())

both_valid = sun_i['sunlamp_use_0.0'].notna() & sun_i['sunlamp_use_2.0'].notna()
both_valid.sum()

# Electronic Device Use
# ---------------------

ED_COLS = ['2237']
INSTANCES = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in ED_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

ed_i = df2_i[all_cols]

for inst in INSTANCES:
    col_raw = f'2237-{inst}'

    if col_raw not in ed_i.columns:
        print(f"Skipping instance {inst}: {col_raw} not found.")
        continue

    ed_i[f'ed_category_{inst}'] = ed_i[col_raw]

# Summary
for inst in INSTANCES:
    col = f'ed_category_{inst}'
    if col in ed_i.columns:
        print(f"[{inst}] Null count : {ed_i[col].isnull().sum()} ({ed_i[col].isnull().sum()/len(ed_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{ed_i[col].value_counts()}\n")

print(ed_i.head())

# Physical Activity Variable
# --------------------------

PHYSICAL_COLS = ['1100','2634','1021','894','3647','1001','914','874','981','2624','1011',
                 '3637','943','991','971','884','904','864','1090','1080','1070','6164','6162','924']
REPLACE_COLS  = ['874','864','894','884','914','904']
INSTANCES     = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in PHYSICAL_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

physical_i = df2_i[all_cols].copy()

# Replace only the relevant columns
replace_cols_all = [f'{field}-{inst}' for inst in INSTANCES for field in REPLACE_COLS
                    if f'{field}-{inst}' in physical_i.columns]
physical_i[replace_cols_all] = physical_i[replace_cols_all].replace({-1: np.nan, -3: np.nan, -2: 0})


def process_physical_instance(df, inst):

    def col(field):
        return f'{field}-{inst}'

    # Adjusted minutes
    df[f'walk_mins_adj_{inst}'] = np.where(
        df[col('864')].notna() & df[col('874')].isna(),
        10, df[col('874')]
    )
    df[f'mod_mins_adj_{inst}'] = np.where(
        df[col('884')].notna() & df[col('894')].isna(),
        10, df[col('894')]
    )
    df[f'vig_mins_adj_{inst}'] = np.where(
        df[col('904')].notna() & df[col('914')].isna(),
        10, df[col('914')]
    )

    # METs
    df[f'walking_MET_{inst}']  = 3.3 * df[f'walk_mins_adj_{inst}'] * df[col('864')]
    df[f'moderate_MET_{inst}'] = 4.0 * df[f'mod_mins_adj_{inst}']  * df[col('884')]
    df[f'rigorous_MET_{inst}'] = 8.0 * df[f'vig_mins_adj_{inst}']  * df[col('904')]

    df[f'sum_MET_{inst}'] = df[[f'walking_MET_{inst}', f'moderate_MET_{inst}', f'rigorous_MET_{inst}']].sum(
        axis=1, skipna=True, min_count=1
    )

    # Category function
    def physical_category(row):
        walk_days = row[col('864')]
        walk_mins = row[f'walk_mins_adj_{inst}']
        mod_days  = row[col('884')]
        mod_mins  = row[f'mod_mins_adj_{inst}']
        vig_days  = row[col('904')]
        vig_mins  = row[f'vig_mins_adj_{inst}']

        values = [walk_days, walk_mins, mod_days, mod_mins, vig_days, vig_mins]

        # All NaN or 0 → missing
        if all(pd.isna(x) or x == 0 for x in values):
            return np.nan

        total_days = sum(x for x in [walk_days, mod_days, vig_days] if pd.notna(x))

        if (vig_days >= 3 and row[f'sum_MET_{inst}'] >= 1500) or \
           (row[f'sum_MET_{inst}'] >= 3000 and total_days >= 7):
            return 2

        elif (
            (vig_days >= 3 and vig_mins >= 20) or
            (mod_days >= 5 and mod_mins >= 30) or
            (walk_days >= 5 and walk_mins >= 30) or
            (total_days >= 5 and row[f'sum_MET_{inst}'] >= 600)
        ):
            return 1

        else:
            return 0

    df[f'physical_category_{inst}'] = df.apply(physical_category, axis=1)

    return df


for inst in INSTANCES:
    inst_cols = [f'{field}-{inst}' for field in REPLACE_COLS]
    if not any(c in physical_i.columns for c in inst_cols):
        print(f"Skipping instance {inst}: no columns found.")
        continue
    physical_i = process_physical_instance(physical_i, inst)

# Final column selection
keep_cols = ['eid']
for inst in INSTANCES:
    keep_cols += [
        f'864-{inst}', f'874-{inst}', f'884-{inst}', f'894-{inst}', f'904-{inst}', f'914-{inst}',
        f'walking_MET_{inst}', f'moderate_MET_{inst}', f'rigorous_MET_{inst}',
        f'sum_MET_{inst}', f'physical_category_{inst}'
    ]
keep_cols = [c for c in keep_cols if c in physical_i.columns]
physical_i = physical_i[keep_cols]

# Summary
for inst in INSTANCES:
    col = f'physical_category_{inst}'
    if col in physical_i.columns:
        print(f"[{inst}] Null count : {physical_i[col].isnull().sum()} ({physical_i[col].isnull().sum()/len(physical_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{physical_i[col].value_counts()}\n")

print(physical_i.head())

print(physical_i.columns)

physical_i = physical_i[['eid','physical_category_0.0','physical_category_2.0']]

print(physical_i.head())

both_valid = physical_i['physical_category_0.0'].notna() & physical_i['physical_category_2.0'].notna()
both_valid.sum()

# Mental Health Variable
# ----------------------

MH_COLS   = ['2050', '2060', '2070', '2080', '2020']
INSTANCES = ['0.0', '2.0']

# Build column list across both instances
all_cols = ['eid'] + [f'{field}-{inst}' for inst in INSTANCES for field in MH_COLS]
all_cols = [c for c in all_cols if c in df2_i.columns]

mh_i = df2_i[all_cols].copy()

for inst in INSTANCES:
    scale_cols = [f'{field}-{inst}' for field in ['2050','2060','2070','2080','2020']
                  if f'{field}-{inst}' in mh_i.columns]
    recode_cols = [f'{field}-{inst}' for field in ['2050','2060','2070','2080']
                   if f'{field}-{inst}' in mh_i.columns]

    if not scale_cols:
        print(f"Skipping instance {inst}: no columns found.")
        continue

    # Replace missing codes for all 5 cols
    mh_i[scale_cols] = mh_i[scale_cols].replace({-1: np.nan, -3: np.nan})

    # Recode 1-4 → 0-3 for the 4 depressive/anxiety cols
    mh_i[recode_cols] = mh_i[recode_cols].replace({1.0: 0.0, 2.0: 1.0, 3.0: 2.0, 4.0: 3.0})

    def mental_health_score(row, inst=inst):
        def col(field):
            return f'{field}-{inst}'

        depressive = safe_sum([row[col('2050')], row[col('2060')]])
        anxiety    = safe_sum([row[col('2070')], row[col('2080')]])
        loneliness = row[col('2020')]

        # Check if all components are missing
        if pd.isna(depressive) and pd.isna(anxiety) and pd.isna(loneliness):
            return np.nan

        score = 0

        # Depressive condition
        if pd.notna(depressive) and depressive >= 3:
            score += 1

        # Anxiety condition
        if pd.notna(anxiety) and anxiety >= 3:
            score += 1

        # Loneliness condition
        if pd.notna(loneliness) and loneliness == 1:
            score += 1

        return score

    mh_i[f'mental_score_{inst}'] = mh_i.apply(mental_health_score, axis=1)

# Summary
for inst in INSTANCES:
    col = f'mental_score_{inst}'
    if col in mh_i.columns:
        print(f"[{inst}] Null count : {mh_i[col].isnull().sum()} ({mh_i[col].isnull().sum()/len(mh_i)*100:.2f}%)")
        print(f"[{inst}] Value counts:\n{mh_i[col].value_counts()}\n")

print(mh_i.head())

both_valid = mh_i['mental_score_0.0'].notna() & mh_i['mental_score_2.0'].notna()
both_valid.sum()

print(mh_i.columns)

mh_i = mh_i[['eid','mental_score_0.0', 'mental_score_2.0']]

print(mh_i.head())

# Merging all calculators together
# --------------------------------

final = df.merge(df1, on='eid', how='inner', validate='one_to_one')\
          .merge(sleep_i, on='eid', how='inner', validate='one_to_one')\
          .merge(smoking_i, on='eid', how='inner', validate='one_to_one')\
          .merge(alcohol_i, on='eid', how='inner', validate='one_to_one')\
          .merge(diet_i, on='eid', how='inner', validate='one_to_one')\
          .merge(physical_i, on='eid', how='inner', validate='one_to_one')\
          .merge(mh_i, on='eid', how='inner', validate='one_to_one')

final.info()

print(final.isnull().sum().head()/len(final))

print(final.columns)

print(final.head())

final_filtered = final[['eid', 'genetic_sex',
       'age_defined_baseline', 'merged_bmi', 'def_CVD_AFTER', 'def_CVD_AF_HF_AFTER',
       'sleep_category_0.0', 'sleep_category_2.0',
       'smoking_category_0.0', 'smoking_category_2.0', 'alcohol_score_0.0', 'alcohol_score_2.0',
       'diet_score_0.0', 'diet_score_2.0', 'physical_category_0.0',
       'physical_category_2.0', 'mental_score_0.0', 'mental_score_2.0']]

final_filtered.info()

print(final_filtered.columns)

# Keep only rows where ALL 0.0 columns are non-null
cols_0 = ['sleep_category_0.0', 'smoking_category_0.0', 'alcohol_score_0.0',
          'diet_score_0.0', 'physical_category_0.0', 'mental_score_0.0']

cols_2 = ['sleep_category_2.0', 'smoking_category_2.0', 'alcohol_score_2.0',
          'diet_score_2.0', 'physical_category_2.0', 'mental_score_2.0']

df_complete_both = final_filtered[final_filtered[cols_0 + cols_2].notna().all(axis=1)]

print(df_complete_both.tail(50))

df_complete_both.info()

cols_0 = ['sleep_category_0.0', 'smoking_category_0.0', 'alcohol_score_0.0',
          'diet_score_0.0', 'physical_category_0.0', 'mental_score_0.0']

cols_2 = ['sleep_category_2.0', 'smoking_category_2.0', 'alcohol_score_2.0',
          'diet_score_2.0', 'physical_category_2.0', 'mental_score_2.0']

# Optional diagnostic cohort: keep rows where at least one T2 lifestyle value is present.
has_any_2 = final_filtered[cols_2].notna().any(axis=1)
df_with_2 = final_filtered[has_any_2].copy()

df_with_2.info()

print(df_with_2.isnull().sum()/len(df_with_2)*100)

print(df_with_2.head())

# Saving files
# ------------

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
baseline_path = OUTPUT_DIR / "baseline.csv"
longitudinal_path = OUTPUT_DIR / "longitudinal.csv"
final_filtered.to_csv(baseline_path, index=False)
df_complete_both.to_csv(longitudinal_path, index=False)
print(f"baseline dataset: {len(final_filtered):,} rows -> {baseline_path}")
print(f"longitudinal dataset: {len(df_complete_both):,} rows -> {longitudinal_path}")

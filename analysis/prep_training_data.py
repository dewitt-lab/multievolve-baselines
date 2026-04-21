#!/usr/bin/env python
"""
Extract training data from supplementary tables for MULTI-evolve.

Generates:
    output/apex/apex_training_log2.csv       (log2 75-min FI, A134P excluded)
    output/apex/combo_muts.csv               (15 mutations, A134P excluded)

    output/dcasrx/dcasrx_training_log2.csv   (log2 gTSM fold improvement)
    output/dcasrx/combo_muts.csv

    output/huabc2/huabc2_training_{expr,bind}_log2.csv
    output/huabc2/combo_muts.csv
"""
import os
import argparse
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.join(SCRIPT_DIR, '..')
DATA_DIR = os.path.join(ROOT_DIR, 'data')
OUTPUT_DIR = os.path.join(ROOT_DIR, 'output')


def prep_apex():
    """(−)A134P campaign only: 75-min fold improvement vs WT, log2-transformed."""
    full_df = pd.read_excel(os.path.join(DATA_DIR, 'Table S14 APEX variants.xlsx'))
    full_df = full_df[[x for x in full_df.columns.values if 'Unnamed' not in x]]

    # Doubles (A134P-free only)
    doubles = full_df[full_df['Origin'] == 'Double Mutant'].copy()
    doubles = doubles[~doubles['Mutant'].str.contains('A134P')].copy()
    doubles['muts'] = doubles['Mutant'].apply(lambda x: x.split('/'))

    # Mutation pool (15 mutations, A134P excluded)
    single_muts = [x.split('/') for x in doubles['Mutant'].unique()]
    single_muts = sorted(set([item for sublist in single_muts for item in sublist]))

    wt_row = pd.DataFrame([{'mutation': 'WT', 'property_value': 1.0, 'round': 0}])

    train_singles = full_df[(full_df['Origin'] == 'PLM zeroshot') &
                            (full_df['Mutant'] != 'A134P')][
        ['Mutant', 'Fold Improvement (75 min)']].copy()
    train_singles.columns = ['mutation', 'property_value']
    train_singles['round'] = 0

    train_doubles = doubles[['Mutant', 'Fold Improvement (75 min)']].dropna().copy()
    train_doubles.columns = ['mutation', 'property_value']
    train_doubles['round'] = 1

    training_df = pd.concat([wt_row, train_singles, train_doubles], ignore_index=True)
    training_df = training_df[training_df['property_value'] > 0].reset_index(drop=True)

    out_dir = os.path.join(OUTPUT_DIR, 'apex')
    os.makedirs(out_dir, exist_ok=True)

    training_df['property_value'] = np.log2(training_df['property_value'])
    out_path = os.path.join(out_dir, 'apex_training_log2.csv')
    training_df.to_csv(out_path, index=False)
    print(f"Wrote {len(training_df)} rows to {out_path}  (log2 75-min FI, A134P-free)")

    out_combo = os.path.join(out_dir, 'combo_muts.csv')
    with open(out_combo, 'w') as f:
        for m in single_muts:
            f.write(m + '\n')
    print(f"Wrote {len(single_muts)} mutations to {out_combo}")


def prep_dcasrx():
    full_df = pd.read_excel(os.path.join(DATA_DIR, 'Table S15 dCasRx variants.xlsx'))
    full_df = full_df[[x for x in full_df.columns.values if 'Unnamed' not in x]]
    full_df['logfc'] = full_df['Fold Improvement (gTSM only)'].apply(np.log2)

    doubles = full_df[full_df['Origin'] == 'Double Mutant'].copy()
    doubles['muts'] = doubles['Mutations'].apply(lambda x: x.split('/'))

    single_muts = [x.split('/') for x in doubles['Mutations'].unique()]
    single_muts = list(set([item for sublist in single_muts for item in sublist]))

    wt_row = pd.DataFrame([{'mutation': 'WT', 'property_value': 1.0, 'round': 0}])

    train_singles = full_df[full_df['Origin'] == 'DMS'][
        ['Mutations', 'Fold Improvement (gTSM only)']].copy()
    train_singles.columns = ['mutation', 'property_value']
    train_singles['round'] = 0

    train_doubles = doubles[['Mutations', 'Fold Improvement (gTSM only)']].copy()
    train_doubles = train_doubles.dropna(subset=['Fold Improvement (gTSM only)'])
    train_doubles.columns = ['mutation', 'property_value']
    train_doubles['round'] = 1

    training_df = pd.concat([wt_row, train_singles, train_doubles], ignore_index=True)

    out_dir = os.path.join(OUTPUT_DIR, 'dcasrx')
    os.makedirs(out_dir, exist_ok=True)

    training_df['property_value'] = np.log2(training_df['property_value'])
    out_path = os.path.join(out_dir, 'dcasrx_training_log2.csv')
    training_df.to_csv(out_path, index=False)
    print(f"Wrote {len(training_df)} rows to {out_path}  (log2)")
    print(f"  WT: 1, singles: {len(train_singles)}, doubles: {len(train_doubles)}")

    out_combo = os.path.join(out_dir, 'combo_muts.csv')
    with open(out_combo, 'w') as f:
        for m in sorted(single_muts):
            f.write(m + '\n')
    print(f"Wrote {len(single_muts)} mutations to {out_combo}")


def prep_huabc2():
    """HuABC2 antibody: expression and binding fold improvement vs WT.

    Mutations keep the native MULTI-evolve multi-chain format: colon separates
    chains (VH:VL), slash separates mutations within a chain.  Examples:
        A40P:WT         single VH mutation
        WT:F82G         single VL mutation
        L11Q:Y93W       one VH + one VL (cross-chain double)
        Q13P/S51I:L21I  two VH + one VL (triple)
        WT:WT           wild type

    KFoldProteinSplitter._shift_mutation_position() handles the position
    offset for the second chain when converting to full sequences.
    """
    full_df = pd.read_excel(os.path.join(DATA_DIR, 'Table S16 HuABC2 variants.xlsx'))
    full_df = full_df[[x for x in full_df.columns.values if 'Unnamed' not in x]]

    # Use the original Mutations column directly (colon-separated chain format)
    full_df['mutation'] = full_df['Mutations']

    # Doubles
    doubles = full_df[full_df['Origin'] == 'Double Mutant'].copy()

    # Mutation pool: unique single mutations that appear in doubles,
    # stored in chain notation so downstream scripts know the chain.
    # VH mutations → 'L11Q:WT', VL mutations → 'WT:F82G'
    single_muts_chain = set()   # chain notation (for combo_muts.csv)
    for mut_str in doubles['mutation']:
        parts = mut_str.split(':')
        for m in parts[0].split('/'):
            if m != 'WT':
                single_muts_chain.add(f'{m}:WT')
        for m in parts[1].split('/'):
            if m != 'WT':
                single_muts_chain.add(f'WT:{m}')
    single_muts_chain = sorted(single_muts_chain)

    out_dir = os.path.join(OUTPUT_DIR, 'huabc2')
    os.makedirs(out_dir, exist_ok=True)

    # Write training CSVs for both measurements
    measurements = {
        'expr': 'Fold Improvement (Expression)',
        'bind': 'Fold Improvement (Binding)',
    }

    for tag, col in measurements.items():
        wt_row = pd.DataFrame([{'mutation': 'WT:WT', 'property_value': 1.0, 'round': 0}])

        train_singles = full_df[full_df['Origin'] == 'PLM zeroshot'][
            ['mutation', col]].copy()
        train_singles.columns = ['mutation', 'property_value']
        train_singles['property_value'] = pd.to_numeric(
            train_singles['property_value'], errors='coerce')
        train_singles = train_singles.dropna(subset=['property_value'])
        train_singles['round'] = 0

        train_doubles = doubles[['mutation', col]].copy()
        train_doubles.columns = ['mutation', 'property_value']
        train_doubles['property_value'] = pd.to_numeric(
            train_doubles['property_value'], errors='coerce')
        train_doubles = train_doubles.dropna(subset=['property_value'])
        train_doubles['round'] = 1

        training_df = pd.concat([wt_row, train_singles, train_doubles],
                                ignore_index=True)
        # Drop non-positive values before log2 transform
        training_df = training_df[
            training_df['property_value'] > 0].reset_index(drop=True)

        training_df['property_value'] = np.log2(
            training_df['property_value'])
        out_path = os.path.join(out_dir, f'huabc2_training_{tag}_log2.csv')
        training_df.to_csv(out_path, index=False)
        print(f"Wrote {len(training_df)} rows to {out_path}  (log2)")
        print(f"  WT: 1, singles: {len(train_singles)}, "
              f"doubles: {len(train_doubles)}")

    # combo_muts.csv — chain notation so run_multievolve.py can offset positions
    out_combo = os.path.join(out_dir, 'combo_muts.csv')
    with open(out_combo, 'w') as f:
        for m in single_muts_chain:
            f.write(m + '\n')
    print(f"Wrote {len(single_muts_chain)} mutations to {out_combo}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protein', type=str,
                        choices=['apex', 'dcasrx', 'huabc2', 'all'])
    args = parser.parse_args()

    if args.protein in ('apex', 'all'):
        prep_apex()
    if args.protein in ('dcasrx', 'all'):
        prep_dcasrx()
    if args.protein in ('huabc2', 'all'):
        prep_huabc2()

import csv
import os
import re

import numpy as np
import pandas as pd


def r_array_to_float_array(array, skip_cols):
    """
    Convert from R vector output to a numpy array
    e.g.,

    >>> r_array_to_float_array('[1] "simple prox: -0.171, -0.908"')
    array([-0.171, -0.908])
    """

    array = re.split(r"\s+", array.rstrip('"'))[skip_cols:]
    array = [float(x.strip(",")) for x in array]
    return np.array(array)


def convert_text_to_df(text_fn):
    """
    Read from `mtx_analysis.R` output and produce dat files
        used within `mtx_six_timestep_plot.tex`
    """

    with open(text_fn) as inf:
        # Ignore header
        _ = inf.readline().strip()

        # Counterfactual estimates
        ests = inf.readline().strip()
        ests = r_array_to_float_array(ests, 3)

        # Causal effect estimates
        diffs = inf.readline().strip()
        diffs = r_array_to_float_array(diffs, 3)

        # Treatment values
        treat_vals = re.split(r"\s+", inf.readline().strip())
        treat_vals = np.array(list(map(int, treat_vals)))

        # Counterfactual estimate 2.5% percentiles
        p25s = inf.readline().strip()
        p25s = r_array_to_float_array(p25s, 1)

        # Counterfactual estimate 97.5% percentiles
        p975s = inf.readline().strip()
        p975s = r_array_to_float_array(p975s, 1)

        # Skip repeat of treatment values
        _ = inf.readline()

        # Causal effect esimate 2.5% percentiles
        diff_p25s = inf.readline().strip()
        diff_p25s = r_array_to_float_array(diff_p25s, 1)

        # Causal effect estimate 97.5% percentiles
        diff_p975s = inf.readline().strip()
        diff_p975s = r_array_to_float_array(diff_p975s, 1)

    df = pd.DataFrame({
        "n": treat_vals,
        "p2.5": p25s,
        "est": ests,
        "p97.5": p975s,
        "diff": diffs,
        "diff_p2.5": diff_p25s,
        "diff_p97.5": diff_p975s,
    })

    return df


def get_table_5_row(df):
    a0 = df["n"] == 0
    a1 = df["n"] == 1

    ey0 = df[a0]["est"]
    ey1 = df[a1]["est"]

    effect = df[a1]["diff"]

    effect_p25 = df[a1]["diff_p2.5"]
    effect_p975 = df[a1]["diff_p97.5"]

    row = [ey0, ey1, effect, effect_p25, effect_p975]
    return [x.to_numpy()[0] for x in row]


def main():
    indir = "mtx_results/"

    models = ["naive_frontdoor", "simple_proximal", "proximal_frontdoor"]
    n_bootstraps = [256]

    # Create Table 5 results for binary treatment
    treatment = "binary"
    for model in models:
        for n in n_bootstraps:
            text_fn = os.path.join(indir, f"{model}_{treatment}_{n}.txt")
            if not os.path.exists(text_fn):
                print(f"{text_fn} doesn't exist; skipping")
            else:
                print(text_fn)
                df = convert_text_to_df(text_fn)
                row = get_table_5_row(df)
                print("{} {}".format(model, " ".join(map("{:.3f}".format, row))))

    # Create .dat files for mtx_six_timestep_plot.tex
    treatment = "cat"
    for model in models:
        for n in n_bootstraps:
            text_fn = os.path.join(indir, f"{model}_{treatment}_{n}.txt")
            dat_fn = os.path.join(indir, f"{model}_{treatment}_{n}.dat")
            if not os.path.exists(text_fn):
                print(f"{text_fn} doesn't exist; skipping")
            else:
                start = f"{text_fn} ┠─┐"
                print(start)

                df = convert_text_to_df(text_fn)
                df[["n", "p2.5", "est", "p97.5"]].to_csv(dat_fn, index=False, sep=" ")

                print(" " * (len(start) - 2) + f" └─┨ {dat_fn}")


if __name__ == "__main__":
    main()

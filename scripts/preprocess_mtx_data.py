import argparse
import numpy as np
import os
import pandas as pd

MTX_DATA_DIR = "mtx_data"


def main(args, outfn):
    """
    Prepreprocess the MTX data file for use in analysis
    """
    infn = os.path.join(MTX_DATA_DIR, "mtx.csv.gz")
    df = pd.read_csv(infn, compression="gzip")

    pats = {patkey: df[df["patkey"] == patkey] for patkey in df["patkey"]}
    pats = {x: y for x, y in pats.items() if y.shape[0] > args.ytime}

    rows = []
    for patkey, patdf in pats.items():
        t0 = patdf[patdf["cummonth"] == 0]
        t_m = patdf[patdf["cummonth"] == args.mtime]
        t_w = patdf[patdf["cummonth"] == args.wtime]
        t_y = patdf[patdf["cummonth"] == args.ytime]

        confounders = ["age_0", "sex", "year_0", "onprd2_0", "duration_0", "edu_0",
                       "smoke_0", "dmrd_0", "rapos"]
        x = t0[confounders].to_numpy()
        confounders[confounders.index("rapos")] = "rapos_0"  # rename s/rapos/rapos_0

        if args.binary_treat:
            if args.no_mtx_after < 0:
                # treatment is whether you started MTX from the start 
                a = t0[["mtxspan"]].to_numpy()
            else:
                # treatment is whether you started MTX by the end of this span
                end = patdf[patdf["cummonth"] == (args.no_mtx_after - 1)]
                a = end[["mtxspan"]].to_numpy()

        else:
            # or, treatment is how many months of MTX pre-mediator
            mtx_start = t0[["mtx1stcu"]].to_numpy()
            if np.isnan(mtx_start):
                # never started treatment => A=0
                a = np.array([[0]])
            else:
                a = np.clip(args.mtime - mtx_start, 0, None)

        # optionally exclude people who started treatment "early"
        if args.no_mtx_before > -1:
            early = patdf[patdf["cummonth"] == args.no_mtx_before]
            if early["mtxspan"].to_numpy() > 0:
                continue

        # optionally exclude people who started treatment "late"
        if args.no_mtx_after > -1:
            # did you start mtx before the cutoff
            before = patdf[patdf["cummonth"] == args.no_mtx_after - 1]
            before = before["mtxspan"].to_numpy() > 0

            # had you started mtx by the time we measured y
            later = patdf[patdf["cummonth"] == args.ytime]
            later = later["mtxspan"].to_numpy() > 0
            if not before and later:
                continue

        # define mediators
        mediators = ["onprd2", "dmrd"]
        m_columns = [f"{tmp}_m" for tmp in mediators]
        m = t_m[mediators].to_numpy()

        # define proxies
        zw_columns = ["haqc", "esrc", "jc", "gsc"]
        z_columns = [f"{tmp}_z" for tmp in zw_columns]
        z = t0[zw_columns].to_numpy()
        w_columns = [f"{tmp}_w" for tmp in zw_columns]
        w = t_w[zw_columns].to_numpy()

        # define outcome
        if args.delta_outcome:
            y = t_y[["jc"]].to_numpy() - t0[["jc"]].to_numpy()
        else:
            y = t_y[["jc"]].to_numpy()

        patkey = t0[["patkey"]].to_numpy()
        row = np.concatenate([patkey, x, m, z, w, a, y], axis=1)
        rows.append(row)

    columns = ["patkey"] + confounders + m_columns + z_columns + w_columns + ["a", "y"]
    data = np.concatenate(rows, axis=0)
    df = pd.DataFrame(data=data, columns=columns).astype(int)

    for column in columns:
        uniq = np.unique(df[column])
        # make things binary that can be
        if uniq.shape[0] == 2:
            mapping = dict(zip(uniq, (0, 1)))
            df[column] = df[column].map(mapping)

    if args.a_ints:
        for column in confounders + ["a"]:
            new_col = df[column] * df["a"]
            df[f"a*{column}"] = new_col

    df.to_csv(outfn, index=False, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary_treat", action="store_true")
    parser.add_argument("--delta_outcome", action="store_true")
    parser.add_argument("--a_ints", action="store_true")

    parser.add_argument("--ytime", type=int, default=12)
    parser.add_argument("--wtime", type=int, default=7)
    parser.add_argument("--mtime", type=int, default=6)

    parser.add_argument("--no_mtx_after", type=int, default=-1)
    parser.add_argument("--no_mtx_before", type=int, default=-1)

    args = parser.parse_args()

    bin_str = "bin" if args.binary_treat else "cont"
    timesteps_str = f"m{args.mtime}_w{args.wtime}_y{args.ytime}"

    mtx_limits = [""]
    if args.no_mtx_after > -1:
        mtx_limits.append(f"mtx_to{args.no_mtx_after}")
    if args.no_mtx_before > -1:
        mtx_limits.append(f"mtx_from{args.no_mtx_before}")
    mtx_str = ".".join(mtx_limits)

    delta_str = ".ydelta" if args.delta_outcome else ""

    a_int_str = ""
    if args.a_ints:
        a_int_str = ".a_int"

    outfn = "".join([
        f"mtx_{bin_str}.",
        timesteps_str,
        mtx_str,
        delta_str,
        a_int_str,
        ".csv.gz"
    ])
    outfn = os.path.join(MTX_DATA_DIR, outfn)
    if os.path.exists(outfn):
        print(f"{outfn} exists")
    else:
        main(args, outfn)

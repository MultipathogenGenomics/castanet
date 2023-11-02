import pandas as pd


def combine_output_csvs(fnames, out_fname):
    beeg_df = pd.DataFrame()
    for fname in fnames:
        try:
            df = pd.read_csv(fname)
            df["sample"] = fname.split("/")[-1].split("_")[0]
            beeg_df = pd.concat([beeg_df, df])
        except:
            continue
    beeg_df.to_csv(out_fname)
    return f"Combine analytical output from all samples in batch saved to {out_fname}"

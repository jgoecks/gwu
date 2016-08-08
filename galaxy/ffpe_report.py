import sys
import glob
import pandas as pd

def read_results(sample):
    # Read results into dataframe.
    dfs = []
    for results_file in glob.glob("*%s*results.txt" % sample):
        # Read from file.
        df = pd.DataFrame.from_csv(results_file, sep="\t", index_col=False)

        # Add tissue type to dataframe.
        tissue = results_file.split("-")[-1]
        tissue = tissue[len(sample):tissue.find("_")]
        df["tissue"] = tissue
        cols = df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df = df[cols]
        dfs.append(df)

    # Combine all results into single dataframe.
    all_results = pd.concat(dfs, ignore_index=True)
    return all_results


if __name__ == "__main__":
    # Returns true iff there is a plasma variant that matches variants in other tissues.
    def has_shared_plasma(df):
        return len(df) > 1 and len(df[df["tissue"] == "P"]) > 0

    df = read_results(sys.argv[1])
    # print df.sort_values(["chrom", "start"]).to_csv(sep="\t", index=False)
    print df.groupby(["chrom", "start", "ref", "alt"]).filter(has_shared_plasma).sort_values(["chrom", "start"]).to_csv(sep="\t", index=False)

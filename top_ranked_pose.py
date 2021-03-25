import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os


def rmsd_box(df):
    name = []
    for i, row in df.iterrows():
        if str(row.rescor_func) == "nan":
            name.append(row.dock_func)
        else:
            name.append(f"{row.dock_func}_{row.rescor_func}")

    df["name"] = name
    sns.boxplot(y="r0", x="name", data=df)

    plt.xticks(rotation=90)

    plt.show()


def rmsd_barplot(df):
    g = sns.FacetGrid(df, col="dock_func", legend_out=True)
    g.map(sns.barplot, "rescore_func", "median", "origin",
          palette=sns.color_palette("colorblind"))

    for ax in g.axes[0]:
        ax.tick_params(axis='x', labelrotation=90)

    plt.legend()
    plt.show()


def success_barplot(df):
    g = sns.FacetGrid(df, col="dock_func", legend_out=True, col_order=["plp", "asp", "chemscore", "goldscore"])
    g.map(sns.barplot, "rescore_func", "1.5", "origin",
          palette=sns.color_palette("colorblind"), order=["plp", "asp", "chemscore", "goldscore"])

    for ax in g.axes[0]:
        ax.tick_params(axis='x', labelrotation=90)


    a = g.axes[0][0]


    for p, l in zip(a.patches, a.lines):

        ci = l.get_ydata()
        val = p.get_height()



    plt.legend()
    plt.show()


def percentage_success(df, origin):

    data = {"dock_func": [],
            "rescore_func": [],
            "1.5": [],
            "2.0": [],
            "2.5": [],
            "runid": [],
            "origin": []}

    score_funcs = ["plp", "asp", "chemscore", "goldscore"]
    ri = df.runid.values[0]
    num = len(set(df.pdb))

    for s in score_funcs:
        for r in score_funcs:
            n = df.loc[(df.dock_func == s) &
                       (df.rescore_func == r) &
                       (df.pose_rank == 0)]

            arr = np.array(n.rmsd.values)

            data["1.5"].append(round((arr <= 1.5).sum() * 100/num))
            data["2.0"].append(round((arr <= 2.0).sum() * 100/num))
            data["2.5"].append(round((arr <= 2.5).sum() * 100/num))

            data["dock_func"].append(s)
            data["rescore_func"].append(r)
            data["runid"].append(ri)
            data["origin"].append(origin)

    return pd.DataFrame(data)


def get_origin(fname):
    origins = ["gold", "fhm1", "fhm2", "control"]

    origin = "other"
    for origin_name in origins:
        if origin_name in fname:
            origin = origin_name

    return origin


def median_rmsd(df, origin):

    data = {"dock_func": [],
            "rescore_func": [],
            "median": [],
            "runid": [],
            "origin": []}

    score_funcs = ["plp", "asp", "chemscore", "goldscore"]
    ri = df.runid.values[0]

    for s in score_funcs:
        for r in score_funcs:

            # change the conditions for
            n = df.loc[(df.dock_func == s) &
                       (df.rescore_func == r) &
                       (df.pose_rank == 0)]

            median = np.median(n["rmsd"].values)

            data["dock_func"].append(s)
            data["rescore_func"].append(r)
            data["median"].append(median)
            data["runid"].append(ri)
            data["origin"].append(origin)

    return pd.DataFrame(data)


def fetch_dfs():
    base = "results/"
    fnames = os.listdir(base)
    dfs = []
    pdbs = [p for p in os.listdir("/local/pcurran/leads_frag")]

    for x in ["4J9W", "2RC9", "5F6W"]:
        pdbs.remove(x)

    for fname in fnames:
        df = pd.read_csv(os.path.join(base, fname), index_col=0)
        df = df.loc[df.pdb.isin(pdbs)]
        dfs.append(df)

    return dfs


def rmsd_main():
    dfs = fetch_dfs()
    new_dfs = []
    for df in dfs:
        ri = list(set(df.runid.values))[0]
        origin = get_origin(ri)
        new_dfs.append(median_rmsd(df, origin))

    lit = pd.read_csv("top_pose/median_rmsd_top_pose_literature.csv", index_col=0)
    new_dfs.append(lit)

    out = pd.concat(new_dfs)

    rmsd_barplot(out)


    ndfs = []
    for df in dfs:
        ndf = df.loc[df.pose_rank == 0]
        ndfs.append(ndf)

    all_df = pd.concat(ndfs)

    ax = sns.barplot(x="runid", y="rmsd", data=all_df)

    ax.tick_params(axis='x', labelrotation=90)

    plt.legend()
    plt.show()


def success_table(df):

    data = {"dock_func": [],
            "rescore_func": [],
            "val": [],
            "origin": []}

    for s in ["plp", "asp", "chemscore", "goldscore"]:
        for r in ["plp", "asp", "chemscore", "goldscore"]:
            for o in ["gold", "fhm1"]:
                n = df.loc[(df.dock_func == s) &
                           (df.rescore_func == r) &
                           (df.origin == o)]

                val = f"{round(np.mean(np.array(n['1.5'].values)))} [{round(np.std(np.array(n['1.5'].values)))}]"
                data["val"].append(val)
                data["dock_func"].append(s)
                data["rescore_func"].append(r)
                data["origin"].append(o)

    ndf = pd.DataFrame(data)
    return ndf.pivot(index='origin', columns=["dock_func", "rescore_func"], values="val")


def success_main():
    dfs = fetch_dfs()
    new_dfs = []
    for df in dfs:
        ri = list(set(df.runid.values))[0]
        origin = get_origin(ri)
        new_dfs.append(percentage_success(df, origin))

    out = pd.concat(new_dfs)


    succ_tab = success_table(out)
    succ_tab.to_csv("success_table.csv")


    print(out)
    success_barplot(out)

if __name__ == "__main__":
    rmsd_main()
    success_main()




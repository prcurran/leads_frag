


import pandas as pd


df = pd.read_csv("results/gold_a_poses.csv", index_col = 0)

df = df.loc[df.pose_id == 1]


new = df.pivot(index='pdb', columns=["dock_func", "rescore_func"], values='rmsd')


lev2 = ["plp","plp","plp","plp",
        "asp","asp","asp","asp",
        "chemscore","chemscore","chemscore","chemscore",
        "goldscore","goldscore","goldscore","goldscore"]
lev3 = ["plp", "asp", "chemscore", "goldscore",
        "plp", "asp", "chemscore", "goldscore",
        "plp", "asp", "chemscore", "goldscore",
        "plp", "asp", "chemscore", "goldscore"
        ]



new = new.reindex([lev2,lev3], axis=1)
new.to_csv("test.csv")
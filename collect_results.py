import os
import pandas as pd
import numpy as np

from ccdc.descriptors import MolecularDescriptors
from ccdc.io import MoleculeReader, EntryReader
from tqdm import tqdm


def check_dir(d):
    if not os.path.exists(d):
        os.mkdir(d)
    return d


def remove_dummy_atoms(mol):
    rm = []
    for atm in mol.heavy_atoms:
        if atm.label == "****":
            rm.append(atm)
    mol.remove_atoms(rm)
    return mol


def process_ff_label(func):
    ab = func.split("_")
    ff_a = ab[0]

    if len(ab) == 1:
        ff_b = None

    else:
        ff_b = ab[1]

    return ff_a, ff_b


def rank_array(li):
    array = np.array(li)
    temp = array.argsort()
    temp = temp[::-1]
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(array))
    return list(ranks)


def create_dataframe(base, run_id, pdbs):
    format_dic = {"asp": "ASP",
                  "chemscore": "Chemscore",
                  "goldscore": "Goldscore",
                  "plp": "PLP"}

    data = {"pdb": [],
            "runid": [],
            "pose_id": [],
            "pose_rank": [],
            "dock_func": [],
            "dock_fitness": [],
            "rescore_func": [],
            "rescore_fitness": [],
            "gold_score": [],
            "rmsd": [],
            "rmsd_rank": []}

    pdbs = [pdb for pdb in pdbs if os.path.isdir(os.path.join(base, pdb, run_id))]

    for pdb in tqdm(pdbs):

        dpath = os.path.join(base, pdb, run_id)
        funcs = [d for d in os.listdir(dpath) if not os.path.isfile(os.path.join(dpath, d))]

        for func in funcs:

            ff_a, ff_b = process_ff_label(func)

            s = []          # for the ranking
            r = []
            for i in range(1, 31):
                pose = EntryReader(os.path.join(dpath, func, "data", f"ranked_{pdb}_ligand_m1_{i}.mol2"))[0]
                attr = pose.attributes
                score = float(attr["Gold.Score"].split("\n")[1][:5])
                fit_score = {k.split(".")[1]: attr[k] for k in [a for a in attr.keys() if "Fitness" in a]}
                rmsd = attr["Gold.Reference.RMSD"]

                data["pdb"].append(pdb)
                data["runid"].append(run_id)
                data["pose_id"].append(i)
                data["dock_func"].append(ff_a)
                data["dock_fitness"].append(float(fit_score[format_dic[ff_a]]))
                data["gold_score"].append(score)
                r.append(float(rmsd))

                if ff_b is None:
                    data["rescore_func"].append(ff_a)
                    s.append(float(fit_score[format_dic[ff_a]]))
                else:
                    data["rescore_func"].append(ff_b)
                    s.append(float(fit_score[format_dic[ff_b]]))

            data["rescore_fitness"].extend(s)
            data["pose_rank"].extend(rank_array(s))

            data["rmsd"].extend(r)
            data["rmsd_rank"].extend(rank_array(r))

    return pd.DataFrame(data)


if __name__ == "__main__":

    base = "/local/pcurran/leads_frag"
    pdbs = [p for p in os.listdir(base) if os.path.isdir(os.path.join(base, p))]

    # run_ids = ["gold", "gold_a", "gold_b", "gold_c", "gold_d", "fhm1"]
    # ["fhm2", "fhm2_a", "fhm2_b", "fhm2_c", "fhm2_d"]
    run_ids = ["control"]

    for run_id in run_ids:

        df = create_dataframe(base, run_id, pdbs)

        df.to_csv(f"results/{run_id}_poses.csv")


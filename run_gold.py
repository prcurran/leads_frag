import argparse
import os
from multiprocessing import Pool

import pandas as pd
import shutil
from ccdc.descriptors import MolecularDescriptors
from ccdc.io import MoleculeReader
from tqdm import tqdm

from gold_conf_template import template


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


def do_dock(args):
    pdb, data_dir, fit_pts, fit_pts_path, run_id = args

    protein_path = os.path.join(data_dir, pdb, f"{pdb}_receptor.mol2")
    ligand_path = os.path.join(data_dir, pdb, f"{pdb}_ligand.mol2")
    ref_path = os.path.join(data_dir, pdb, f"{pdb}_ref.mol2")

    if fit_pts == "1":
        fp_path = os.path.join(data_dir, pdb, fit_pts_path)
        assert os.path.exists(fp_path)

    else:
        fp_path = "fit_pts.mol2"

    scoring_funcs = ["goldscore", "chemscore", "asp", "plp"]
    auto_scale = 1
    gold_exe = "/local/pcurran/CCDC/Custom_GOLD/discovery-build-developer/bin/gold_auto"

    dock_func = []
    rescor_func = []
    rmsds = []

    for scor in scoring_funcs:
        for rescor in scoring_funcs:

            outdir = check_dir(os.path.join(data_dir, pdb, run_id))

            if scor == rescor:
                rescor = None
                outdir = check_dir(os.path.join(outdir, f"{scor}"))

            else:
                outdir = check_dir(os.path.join(outdir, f"{scor}_{rescor}"))

            dump_dir = check_dir(os.path.join(outdir, "data"))

            conf_file = template(auto_scale,
                                 ref_path,
                                 ligand_path,
                                 protein_path,
                                 scor,
                                 rescor,
                                 dump_dir=dump_dir,
                                 fit_pts=fp_path,
                                 fit=fit_pts)

            with open(os.path.join(outdir, "gold.conf"), "w") as w:
                w.write(conf_file)

            cmd = f"{gold_exe} {outdir}/gold.conf"
            os.system(cmd)

            ref_path = os.path.join(data_dir, pdb, f"{pdb}_ref-ligand.pdb")
            ref = MoleculeReader(ref_path)[0]
            docks = [remove_dummy_atoms(
                MoleculeReader(os.path.join(dump_dir, f"ranked_{pdb}_ligand_m1_{i}.mol2"))[0]
            )
                for i in range(1, 31)]

            r = [MolecularDescriptors.rmsd(ref, dock, exclude_hydrogens=True) for dock in docks]

            dock_func.append(scor)
            rescor_func.append(rescor)
            rmsds.append(r)

    runs = len(scoring_funcs) ** 2
    ranked_rmsds = zip(*rmsds)

    data = {"pdbs": [pdb] * runs,
            "run": [run_id] * runs,
            "dock_func": dock_func,
            "rescor_func": rescor_func}

    data.update({f"r{i}": x for i, x in enumerate(ranked_rmsds)})

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(data_dir, pdb, run_id, "results.csv"))


class Organiser(argparse.ArgumentParser):
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)

        self.add_argument(
            'run_id',
            help='name of results dir'
        )

        self.add_argument(
            'data_dir',
            help='name of data dir'
        )

        self.add_argument(
            '-f', '--fit_pts', default=0,
            help='read fit points?'
        )

        self.add_argument(
            '-p', '--fit_pts_path', default="hotspot/fitting_pts.mol2",
            help='relative fp path'
        )

        self.args = self.parse_args()

    def run(self):

        pdbs = [p for p in os.listdir(self.args.data_dir)
                if os.path.isdir(os.path.join(self.args.data_dir, p))]

        for x in ["4J9W", "2RC9", "5F6W"]:
            pdbs.remove(x)

        for pdb in pdbs:
            fpath = f"/local/pcurran/leads_frag/{pdb}/{self.args.run_id}"
            if os.path.exists(fpath):
                shutil.rmtree(fpath)

        args = zip(pdbs,
                   [self.args.data_dir] * len(pdbs),
                   [self.args.fit_pts] * len(pdbs),
                   [self.args.fit_pts_path] * len(pdbs),
                   [self.args.run_id] * len(pdbs)
                   )

        processes = 7

        with Pool(processes=processes) as pool:
            list(tqdm(pool.map(do_dock, args), total=len(pdbs)))


if __name__ == "__main__":
    o = Organiser()
    o.run()

import os
from concurrent import futures

import pandas as pd
from tqdm import tqdm

from ccdc.protein import Protein
from ccdc.io import MoleculeReader

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter
from hotspots.result import Extractor, Results


def check_dir(d):
    if not os.path.exists(d):
        os.mkdir(d)
    return d


def calc(args):
    prot_file, hotspot_file = args

    prot = Protein.from_file(prot_file)
    #  pre prepared
    runner = Runner()
    settings = Runner.Settings()
    settings.apolar_translation_threshold = 8
    settings.polar_translation_threshold = 10

    # pdb = os.path.basename(prot_file)[0][:4]
    #
    # mol_path = os.path.join(os.path.dirname(prot_file))

    hr = runner.from_protein(prot, nprocesses=3, settings=settings, probe_size=3)

    for p, g in hr.super_grids.items():
        hr.super_grids[p] = g.dilate_by_atom()

    try:
        e = Extractor(hr)
        bv = e.extract_volume(volume=250)

    except:
        bv = Results(protein=hr.protein.copy(),
                     super_grids={p: g.copy() for p, g in hr.super_grids.items()})

    hr.identifier = "hotspot"
    bv.identifier = "bcv"

    with HotspotWriter(hotspot_file) as w:
        w.write([hr, bv])


def main():

    in_base = "/local/pcurran/leads_frag"
    # pdbs = [p for p in os.listdir(in_base)
    #         if not os.path.exists(os.path.join(in_base, p, "hotspot", "out.zip"))]

    # pdbs = ["4J9W", "2RC9", "5F6W"]
    pdbs = ["1Q11"]
    print(pdbs)
    print(len(pdbs))

    prot_files = [os.path.join(in_base, pdb, f"{pdb}_receptor.mol2") for pdb in pdbs]
    hotspot_files = [check_dir(os.path.join(in_base, pdb, "hotspot")) for pdb in pdbs]

    args = zip(prot_files, hotspot_files)

    # with futures.ProcessPoolExecutor(max_workers=7) as executor:
    #     executor.map(calc, args)

    for arg in args:
        calc(arg)


if __name__ =="__main__":
    main()
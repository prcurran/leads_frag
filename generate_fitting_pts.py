from hotspots.grid_extension import Grid
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Results
from ccdc.io import MoleculeReader
import os
import argparse
from tqdm import tqdm


def check_dir(d):
    if not os.path.exists(d):
        os.mkdir(d)
    return d


def fp_scheme(fpath, percentile, low, high, id):
    fpath = os.path.join(fpath, "out.zip")
    pdb = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(fpath))))

    if os.path.exists(fpath):
        with HotspotReader(fpath) as r:
            hr = r.read()

        hr.docking_fitting_pts(fname=os.path.join(os.path.dirname(fpath), f"fitting_pts_{id}.mol2"),
                               percentile=percentile,
                               low=low,
                               high=high)
    else:
        print(f"FILE NOT FOUND: {pdb}")


def masked_hotspot(base, pdb, hotspot_path):

    assert os.path.exists(hotspot_path)

    with HotspotReader(os.path.join(hotspot_path, "out.zip")) as r:
        hr = [h for h in r.read() if h.identifier == "hotspot"][0]

    b = (hr.buriedness > 3) * hr.buriedness

    crystal_lig = MoleculeReader(os.path.join(base, pdb, f"{pdb}_ref.mol2"))[0]

    g = hr.buriedness.copy_and_clear()

    for atm in crystal_lig.heavy_atoms:
        g.set_sphere(point=atm.coordinates, radius=6, value=1, mode="replace", scaling="None")

    mol_buried = (g & b) * b

    common_mol_buried = hr.super_grids["apolar"].common_boundaries(mol_buried)

    apolar = (common_mol_buried & hr.super_grids["apolar"]) * hr.super_grids["apolar"]
    donor = (common_mol_buried & hr.super_grids["donor"]) * hr.super_grids["donor"]
    acceptor = (common_mol_buried & hr.super_grids["acceptor"]) * hr.super_grids["acceptor"]

    return Results(super_grids={"apolar": apolar, "donor": donor, "acceptor": acceptor},
                   protein=hr.protein,
                   buriedness=common_mol_buried
                   )


class Organiser(argparse.ArgumentParser):
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)

        self.add_argument(
            'percentile',
            help=''
        )

        self.add_argument(
            'low',
            help=''
        )

        self.add_argument(
            'high',
            help=''
        )

        self.add_argument(
            'id',
            help=''
        )

        self.args = self.parse_args()

    def run(self):
        base = "/local/pcurran/leads_frag"

        pdbs = [p for p in os.listdir(base)
                if os.path.isdir(os.path.join(base, p))]

        fails = []
        for pdb in tqdm(pdbs):
            try:
                hotspot_path = os.path.join(os.path.join(base, pdb, "hotspot"))
                masked_path = os.path.join(hotspot_path, "masked_hotspot")
                print(pdb)

                if not os.path.exists(masked_path):
                    masked = masked_hotspot(base, pdb, hotspot_path)

                    with HotspotWriter(masked_path) as w:
                        w.write(masked)

                fp_scheme(fpath=masked_path,
                          percentile=float(self.args.percentile),
                          low=float(self.args.low),
                          high=float(self.args.high),
                          id=self.args.id)
            except:
                print(f"{pdb} FAILED")
                fails.append(pdb)


if __name__ == "__main__":
    o = Organiser()
    o.run()
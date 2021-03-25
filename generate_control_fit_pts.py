from ccdc import io
from hotspots.grid_extension import Grid
import numpy as np
from hotspots.hs_io import HotspotReader
from tqdm import tqdm
import os


def docking_fitting_pts(g, fname, low = 0.1, high = 0.9):
    """

    :return:
    """
    count = 1
    grid_low, grid_high = g.extrema
    nx, ny, nz = g.nsteps

    tail = ""
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                val = g.value(i, j, k)
                if val >= 1:
                    x, y, z = g.indices_to_point(i, j, k)
                    # square to amplify top points
                    scaled_val = (high - low) * (val**2 / grid_high**2) + low
                    tail += f"""      {count} ****        {x:.4f}   {y:.4f}  {z:.4f} Du {scaled_val:.3f}  \n"""
                    count += 1

        num_pts = count - 1
        head = f"""
@<TRIPOS>MOLECULE
GA Fitting Points
{num_pts}    0    0
SMALL
NO_CHARGES


@<TRIPOS>ATOM
"""

    out_str = head + tail

    with open(fname, "w") as w:
        w.write(out_str)


def main():
    base = "/local/pcurran/leads_frag"

    pdbs = [p for p in os.listdir(base)
            if os.path.isdir(os.path.join(base, p))]

    for pdb in tqdm(pdbs):
        fpath = os.path.join(base, pdb, f"{pdb}_ref.mol2")
        mol = io.MoleculeReader(fpath)[0]
        g = Grid.from_molecule(mol, mode='replace', padding=10, scaling=0.5)

        out_path = os.path.join(base, pdb, "control.mol2")
        docking_fitting_pts(g, fname=out_path, high=1)


if __name__ == "__main__":
    main()
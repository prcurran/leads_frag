import json
import os
from urllib.parse import quote_plus
from urllib.request import Request, urlopen

from ccdc.cavity import Cavity
from ccdc.io import MoleculeReader, MoleculeWriter

from hotspots.protein_extension import Protein
from hotspots.wrapper_pymol import PyMOLFile, PyMOLCommands

from scipy.spatial import distance
import numpy as np


def centroid(coords):
    x = set()
    y = set()
    z = set()

    for i in coords:
        x.add(i[0])
        y.add(i[1])
        z.add(i[2])
    return [sum(x) / len(x),
            sum(y) / len(y),
            sum(z) / len(z)]


def ftp_download(pdb_code, out_dir="pdb"):
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = urlopen(Request(url))
    f = response.read().decode("utf-8")
    # write out decoded file
    with open(os.path.join(out_dir, f"{pdb_code}.pdb"), "w") as w:
        w.write(f)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def format_cavity_file(residues):
    # 10 residues per line max (actual max 250 chara but for readability)
    # tag line required
    # last line must be blank
    # seperated by one space
    first_line = "> <Gold.Protein.ActiveResidues>\n"
    lines = "\n".join([" ".join(line) for line in chunks(residues, 10)])
    return first_line + lines + "\n "


def create_pymol_file(residues, pdb):
    f = PyMOLFile()

    f.commands += PyMOLCommands.load(f"{pdb.upper()}.pdb", pdb)
    f.commands += PyMOLCommands.hide("sticks")

    for residue in residues:
        f.commands += PyMOLCommands.select("sele", f"resi {residue[3:]}")
        f.commands += PyMOLCommands.show("stick", "sele")

    f.commands += PyMOLCommands.load(f"{pdb}_ref.mol2", "lig")
    return f


def main():

    data_dir = "/local/pcurran/leads_frag"
    # pdbs = [p for p in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, p))]
    pdbs = ["4P7X"]
    for pdb in pdbs:
        print(pdb)
        #  download
        ftp_download(pdb, out_dir=os.path.join(data_dir, pdb))

        # prepare
        fpath = os.path.join(data_dir, pdb, f"{pdb}.pdb")
        prot = Protein.from_file(fpath)
        prot.remove_all_metals()
        prot.remove_all_waters()

        # cavity reader can not handle incomplete residues
        for r in prot.residues:
            if r.is_incomplete():
                prot.remove_residue(r.identifier)

        # cofactors were removed in the original set
        for cof in prot.cofactors:
            prot.remove_cofactor(cof.identifier)

        for lig in prot.ligands:
            prot.remove_ligand(lig.identifier)

        with MoleculeWriter(fpath) as w:
            w.write(prot)

        #  molecule centre of geometry
        mol_file = os.path.join(data_dir, pdb, f"{pdb}_ref.mol2")
        mol = MoleculeReader(mol_file)[0]
        mol_centroid = mol.centre_of_geometry()

        #  detect cavity
        cavities = Cavity.from_pdb_file(fpath)

        cavity_centroids = []
        for i, c in enumerate(cavities):
            coords = [f.coordinates for f in c.features for c in cavities]
            cavity_centroids.append(centroid(coords))

        index = np.argmin(distance.cdist(cavity_centroids, [mol_centroid]))
        print(index)
        cav = cavities[index]

        for f in cav.features:
            print(f), print(f.residue)

        # Cavity to GOLD cavity file
        residues = list({str(f.residue.identifier).split(":")[1] for f in cav.features})
        cav_str = format_cavity_file(residues)

        with open(os.path.join(data_dir, pdb, "cavity.txt"), "w") as w:
            w.write(cav_str)

        f = create_pymol_file(residues, pdb)
        f.write(os.path.join(data_dir, pdb, "cav_vis.py"))


if __name__ == "__main__":
    main()


import tempfile
import os
from urllib.request import Request, urlopen
from ccdc.protein import Protein
from ccdc.io import MoleculeReader, MoleculeWriter


def ftp_download(pdb_code, hetid):
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = urlopen(Request(url))
    f = response.read().decode("utf-8")
    # write out decoded file
    out_dir = tempfile.mkdtemp()

    fpath = os.path.join(out_dir, f"{pdb_code}.pdb")

    with open(fpath, "w") as w:
        w.write(f)

    prot = Protein.from_file(fpath)
    prot.detect_ligand_bonds()

    return [l for l in prot.ligands if l.identifier.split(":")[1][:3] == hetid][0]


def main():
    base = "/local/pcurran/leads_frag"
    pdbs = [p for p in os.listdir(base) if os.path.isdir(os.path.join(base, p))]

    for pdb in pdbs:
        hetid = MoleculeReader(os.path.join(base, pdb, f"{pdb}_ligand.mol2"))[0].identifier

        mol = ftp_download(pdb, hetid)

        with MoleculeWriter(os.path.join(base, pdb, f"{pdb}_ref.mol2")) as w:
            w.write(mol)


if __name__ == "__main__":
    main()
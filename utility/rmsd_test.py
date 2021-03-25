import os
from ccdc import io
from ccdc.descriptors import MolecularDescriptors


if __name__ == "__main__":

    pdb = "4G46"
    base = f"/local/pcurran/leads_frag/{pdb}"
    #  test
    mol1 = io.MoleculeReader(os.path.join(base, f"{pdb}_ligand.mol2"))[0]
    mol2 = io.MoleculeReader(os.path.join(base, f"{pdb}_ref.mol2"))[0]

    mol3 = io.MoleculeReader(os.path.join(base, "gold/goldscore/data/ranked_4G46_ligand_m1_1.mol2"))[0]

    rm = []

    for atm in mol3.heavy_atoms:
        if atm.label == "****":
            rm.append(atm)

    mol3.remove_atoms(rm)

    print([atm.label for atm in mol1.heavy_atoms])
    print([atm.label for atm in mol2.heavy_atoms])
    print([atm.label for atm in mol3.heavy_atoms])

    a = MolecularDescriptors.rmsd(mol1, mol3)


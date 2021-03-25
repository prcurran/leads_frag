import shutil
import os


src = "/local/pcurran/leads_frag"
dst = "/local/pcurran/for_chris"

pdbs = [p for p in os.listdir(src) if os.path.isdir(os.path.join(src, p))]

for pdb in pdbs:
    print(pdb)
    src_pdb = os.path.join(src, pdb)
    dst_pdb = os.path.join(dst, pdb)

    if not os.path.exists(dst_pdb):
        os.mkdir(dst_pdb)

    fnames = [f"{pdb}_ligand.mol2", f"{pdb}_ref-ligand.pdb", f"{pdb}_ref.mol2", f"{pdb}_receptor.mol2"]

    for fname in fnames:
        shutil.copyfile(os.path.join(src_pdb, fname), os.path.join(dst_pdb, fname))

    dst_hotspot = os.path.join(dst_pdb, "hotspot")
    src_hotspot = os.path.join(src_pdb, "hotspot")

    if not os.path.exists(dst_hotspot):
        os.mkdir(dst_hotspot)

    for fname in ["out.zip", "pymol_file.py"]:
        shutil.copyfile(os.path.join(src_hotspot, fname), os.path.join(dst_hotspot, fname))

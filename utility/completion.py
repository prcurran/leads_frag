import os


def runs(base, run_id, pdbs):

    left = [pdb for pdb in pdbs if not os.path.exists(os.path.join(base, pdb, run_id, "results.csv"))]

    print(left)
    print(len(left))


def hrs(base, pdbs):
    left = [pdb for pdb in pdbs if not os.path.exists(os.path.join(base, pdb, "hotspot/out.zip"))]
    print(left)


def fps(base, pdbs):
    left = [pdb for pdb in pdbs if not os.path.exists(os.path.join(base, pdb, "hotspot/fitting_pts_scheme_1.mol2"))]
    print(left)

if __name__ == "__main__":

    base = "/local/pcurran/leads_frag"
    pdbs = [p for p in os.listdir(base) if os.path.isdir(os.path.join(base, p))]

    # fps(base, pdbs)
    # hrs(base, pdbs)

    runs(base, "fhm1_a", pdbs)
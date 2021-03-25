


from ccdc import io

from scipy.spatial import distance




mol = io.MoleculeReader("/local/pcurran/leads_frag/3CHC/gold/goldscore/data/fit_pts.mol2")[0]


pt_1 = [[a.coordinates.x,
      a.coordinates.y,
      a.coordinates.z]
     for i, a in enumerate(mol.atoms) if i <= 0]

pt_all = [[a.coordinates.x,
      a.coordinates.y,
      a.coordinates.z]
     for i, a in enumerate(mol.atoms) if i > 0]


ds = distance.cdist(pt_1, pt_all)

print(min(ds[0]))


def template(autoscale, ref_ligand, dock_ligand, protein, fitfunc, rescore=None,
             fit_pts="fit_pts.mol2", fit=1, dump_dir="data"):

    if rescore == None:
        rescore_line = ""
        run_flag = ""
        gold_func_line = f"gold_fitfunc_path = {fitfunc}"
    else:
        gold_func_line = f"gold_fitfunc_path = consensus_score"
        rescore_line = f"rescore_fitfunc_path = {rescore}"
        run_flag = "run_flag = CONSENSUS"

    return f"""
  GOLD CONFIGURATION FILE

  AUTOMATIC SETTINGS
autoscale = {autoscale}

  POPULATION
popsiz = auto
select_pressure = auto
n_islands = auto
maxops = auto
niche_siz = auto

  GENETIC OPERATORS
pt_crosswt = auto
allele_mutatewt = auto
migratewt = auto

  FLOOD FILL
radius = 6
origin = 0 0 0
do_cavity = 0
floodfill_atom_no = 0
cavity_file = {ref_ligand}
floodfill_center = cavity_from_ligand

  DATA FILES
ligand_data_file {dock_ligand} 30
ligand_reference_file = {ref_ligand}
param_file = DEFAULT
set_ligand_atom_types = 1
set_protein_atom_types = 0
directory = {dump_dir}
tordist_file = DEFAULT
make_subdirs = 0
save_lone_pairs = 1
fit_points_file = {fit_pts}
read_fitpts = {fit}

  FLAGS

internal_ligand_h_bonds = 0
flip_free_corners = 0
match_ring_templates = 0
flip_amide_bonds = 0
flip_planar_n = 1 flip_ring_NRR flip_ring_NHR
flip_pyramidal_n = 0
rotate_carboxylic_oh = flip
use_tordist = 1
postprocess_bonds = 1
rotatable_bond_override_file = DEFAULT
solvate_all = 1
diverse_solutions = 1
divsol_cluster_size = 3
divsol_rmsd = 1.5

  TERMINATION
early_termination = 0
n_top_solutions = 3
rms_tolerance = 1.5

  COVALENT BONDING
covalent = 0

  SAVE OPTIONS
save_score_in_file = 1 comments
save_protein_torsions = 1
# concatenated_output = {dump_dir}/docked_ligands.mol2
output_file_format = MOL2

  FITNESS FUNCTION SETTINGS
initial_virtual_pt_match_max = 3
relative_ligand_energy = 1
docking_param_file = DEFAULT
docking_fitfunc_path = {fitfunc}
{gold_func_line}
score_param_file = DEFAULT
rescore_param_file = DEFAULT
{rescore_line}
{run_flag}

  PROTEIN DATA
protein_datafile = {protein}

"""

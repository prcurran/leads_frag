#!/bin/bash
conda activate hotspots


python run_gold.py "gold_0" "/local/pcurran/leads_frag"
python run_gold.py "gold_1" "/local/pcurran/leads_frag"
python run_gold.py "gold_2" "/local/pcurran/leads_frag"
python run_gold.py "gold_3" "/local/pcurran/leads_frag"
python run_gold.py "gold_4" "/local/pcurran/leads_frag"


python run_gold.py "fhm1_0" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_1.mol2"
python run_gold.py "fhm1_1" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_1.mol2"
python run_gold.py "fhm1_2" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_1.mol2"
python run_gold.py "fhm1_3" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_1.mol2"
python run_gold.py "fhm1_4" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_1.mol2"


python run_gold.py "fhm2_0" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_2.mol2"
python run_gold.py "fhm2_1" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_2.mol2"
python run_gold.py "fhm2_2" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_2.mol2"
python run_gold.py "fhm2_3" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_2.mol2"
python run_gold.py "fhm2_4" "/local/pcurran/leads_frag" -f 1 -p "hotspot/masked_hotspot/fitting_pts_masked_scheme_2.mol2"


python run_gold.py "control" "/local/pcurran/leads_frag" -f 1 -p "control.mol2"
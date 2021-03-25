#!/bin/bash
conda activate hotspots

#python run_gold.py "gold" "/local/pcurran/leads_frag"
#python run_gold.py "gold2" "/local/pcurran/leads_frag"


#python generate_fitting_pts.py 10 0.1 0.9 "masked_scheme_1"

#python generate_fitting_pts.py 50 0.1 0.9 "masked_scheme_2"

python generate_fitting_pts.py 10 0.1 2 "masked_scheme_3"
python generate_fitting_pts.py 10 0.1 5 "masked_scheme_4"
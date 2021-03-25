"""
A module for key PDB data collection tools
"""
import json
import os
from urllib.parse import quote_plus
from urllib.request import Request, urlopen
from pprint import pprint
import pandas as pd
from tqdm import tqdm


def data_query(pdb):
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb}'
    # call the url
    response = urlopen(Request(url))

    if response.status == 200:
        # response received in bytes
        decoded = response.read().decode("utf-8")
        # covert json str to dict
        return json.loads(decoded)
    else:
        print(f"Status:{response.status}\n{response.msg}")


if __name__ == "__main__":
    # Step 1: Run the UniProt Search

    base = "/local/pcurran/leads_frag"
    pdbs = os.listdir(base)

    data = {'pdb': [],
            'comp_id': [],
            'link': [],
            'provenance_code': [],
            'type': [],
            'unit': [],
            'value': []}

    for pdb in tqdm(pdbs):
        results = data_query(pdb.upper())
        try:
            pieces = results['rcsb_binding_affinity']
            for x in pieces:
                data['pdb'].append(pdb)
                for k in data.keys():
                    if k == "pdb":
                        continue
                    elif k not in x:
                        data[k].append(None)
                    else:
                        data[k].append(x[k])



        except KeyError:
            data['pdb'].append(pdb)
            for k in data.keys():
                if k != "pdb":
                    data[k].append(None)

    pprint(data)
    df = pd.DataFrame(data)
    df.to_csv("data/binding_data.csv")
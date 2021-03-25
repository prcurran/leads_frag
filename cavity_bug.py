import os
from urllib.request import Request, urlopen
from ccdc.cavity import Cavity
import tempfile


def ftp_download(pdb_code, out_dir="pdb"):
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = urlopen(Request(url))
    f = response.read().decode("utf-8")
    # write out decoded file
    with open(os.path.join(out_dir, f"{pdb_code}.pdb"), "w") as w:
        w.write(f)


def main():

    print(os.environ["CSDHOME"])
    # from ccdc import io
    # print(io.EntryReader('CSD'))
    # pdb = "4P7X"
    # tmp = tempfile.mkdtemp()
    # ftp_download(pdb, out_dir=tmp)
    #
    # fpath = os.path.join(tmp, f"{pdb}.pdb")
    # cavities = Cavity.from_pdb_file(fpath)
    # cav = cavities[0]
    #
    # for feature in cav.features:
    #     print(feature)
    #     print(feature.residue)




if __name__ == "__main__":
    main()
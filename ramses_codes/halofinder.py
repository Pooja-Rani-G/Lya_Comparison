import yt
import numpy as np
import argparse
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="path of info_*.txt file")
parser.add_argument("-o","--output", help="folder path for outputs")

args = parser.parse_args()
info = args.input
direc= args.output

ds=yt.load(info)

os.makedirs(direc, exist_ok=True)

hc = HaloCatalog(data_ds=ds, finder_method="hop",output_dir=direc,
				finder_kwargs={"threshold": 100, "ptype":"all"})
hc.create()
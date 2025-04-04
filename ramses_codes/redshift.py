import numpy as np
import yt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-len","--boxsize", help="box size",default=100)
parser.add_argument("-lm","--levelmin",help="Base level of simulation",default=8)
parser.add_argument("-o","--output", help="number of outputs")

args = parser.parse_args()
length= int(args.boxsize)
lev=int(args.levelmin)
n_out=int(args.output)

os.makedirs(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/",exist_ok=True)

f=open(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/redshift.txt","x")
for i in range(n_out):
    j=i+1
    path=f"/user1/poojarani/Lya_Comparison/ramses_data/lev{lev:03}_len{length}/output_{j:05}/"
    ds=yt.load(path)
    z = ds.current_redshift
    f.write(str(z)+"\n")
f.close()


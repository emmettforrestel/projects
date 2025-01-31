from podio import root_io
import glob
import hist
import functions
import pickle
import argparse
import math
from itertools import combinations
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plots', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23"
folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"

files = glob.glob(f"{folder}/*.root")

# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

z_step = 2

if args.calculate:
    Z_plus_values = []
    Z_neg_values = []
    
    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            nEvents += 1
            all_hits =  event.get("VertexBarrelCollection")
            good_hits = functions.filter_hits(all_hits, layer_radii)
            for hit in good_hits:
                edep=1e6*hit.getEDep() # convert to keV
                mc = hit.getMCParticle()
                if mc.getGeneratorStatus() != 1:  # only primaries
                    continue
                z = hit.getPosition().z
                if z >= 0:
                    Z_plus_values.append(z)
                else:
                    Z_neg_values.append(abs(z))
                
        if i > args.maxFiles:
            break
if args.plots:
    outdir = "./"
    plt.hist(Z_plus_values, bins=100, color='r', edgecolor='r', histtype='step', label="pos") # Add labels and title 
    plt.hist(Z_neg_values, bins=100, color='b', edgecolor='b', histtype='step', label='neg') # Add labels and title 
    plt.xlabel('Z_comp') 
    plt.ylabel('Events') 
    plt.savefig('Z_comp.png')
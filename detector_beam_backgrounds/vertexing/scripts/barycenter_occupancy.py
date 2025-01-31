from podio import root_io
import glob
import hist
import functions
import pickle
import argparse
import math
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()


##########################################################################################
#  this file is for calculating the average and max occupancies in the first layer
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23"
#folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"

files = glob.glob(f"{folder}/*.root")


# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

if args.calculate:
    cellCounts = {}
    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            nEvents += 1
            hits = event.get("VertexBarrelCollection")
            good_hits = functions.filter_hits(hits, layer_radii)
            barycenters = functions.barycenter_energy(good_hits)
            for b in barycenters:
                ##### ENERGY-CALC #####
                edep = b[1]
                bx = b[2]
                by = b[3]
                bz = b[4]
                microcount = b[5]
                # the cellID corresponds to a cluster of pixels = module
                # the readout is done per module
                cellID = b[6]
                if not cellID in cellCounts:
                    cellCounts[cellID] = 0
                # increment the hit count for this module
                # strictly speaking all the primaries+secondaries within the dR should be treated as 1 hit
                cellCounts[cellID] += 1
        if i > args.maxFiles:
            break

    # normalize the hits over the number of events
    print(cellCounts)
    cellCounts_norm = [cellCounts[c]/nEvents for c in cellCounts]
    
    max_hits = max(cellCounts_norm)
    avg_hits = np.mean(cellCounts_norm)

    safety_factor = 3
    cluster_size = 5

    # this corresponds to the occupancy per module
    print(f"Maximum occupancy: {max_hits*safety_factor*cluster_size}")
    print(f"Average occupancy: {avg_hits*safety_factor*cluster_size}")
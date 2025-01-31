from podio import root_io
import glob
import hist
import functions
import pickle
import argparse
import math
from itertools import combinations

parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plots', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()


##########################################################################################
# this file is for plotting the number of hits in a 2D map of phi and z and purely as a
# function of phi and theta
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23"
folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"

files = glob.glob(f"{folder}/*.root")


# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

z_step = 2

if args.calculate:
    axis_z = hist.axis.Regular(int(max_z/2), -max_z, max_z, name = "z", underflow=False, overflow=False)
    axis_theta = hist.axis.Regular(int(180/5), 0, 180, name = "theta", underflow=False, overflow=False)
    axis_phi = hist.axis.Regular(int(360/5), 0, 360, name = "phi", underflow=False, overflow=False)
    axis_microhits = hist.axis.Regular(15, 0, 15, name = "microhits", underflow=False, overflow=False)
    
    axis_energy = hist.axis.Regular(100, 0, 100, name = "energy")
    axis_distance = hist.axis.Regular(40, 0, 40, name = "distance")

    hist_theta_phi = hist.Hist(axis_theta, axis_phi)
    hist_z_phi = hist.Hist(axis_z, axis_phi)
    hist_energy = hist.Hist(axis_theta, axis_energy)
    hist_distance = hist.Hist(axis_energy, axis_distance)
    hist_microhits = hist.Hist(axis_microhits)
    
    microcounts = []

    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            nEvents += 1
            all_hits =  event.get("VertexBarrelCollection")
            good_hits = functions.filter_hits(all_hits, layer_radii)
            barycenters = functions.barycenter_energy(good_hits)
                     
            for b in barycenters:
                ##### ENERGY-CALC #####
                edep = b[1]
                bx = b[2]
                by = b[3]
                bz = b[4]
                microcount = b[5]
                ph = functions.phi(bx, by) * (180 / math.pi)
                
                hist_z_phi.fill(bz, ph)
                hist_microhits.fill(microcount)
        if i > args.maxFiles:
            break

    # normalize the histograms over the number of events, to get the average number of hits per event
    hist_z_phi /= nEvents
    hist_microhits /=nEvents
    
    hists = {}
    
    hists['hist_z_phi'] = hist_z_phi
    hists["hist_microhits"] = hist_microhits

    with open("output_hitmaps.pkl", "wb") as f:
        pickle.dump(hists, f)

if args.plots:

    outdir = "./"
    with open("output_hitmaps.pkl", "rb") as f:
        hists = pickle.load(f)    
    hist_z_phi = hists['hist_z_phi']
    hist_microhits = hists['hist_microhits']
    
    functions.plot_hist(hist_microhits, f"{outdir}/bary_microhits.png", f"Hit count first layer", xMin=0, xMax=15, xLabel="Number of Microhits", yLabel="Microhit Count Freqency")
    functions.plot_2dhist(hist_z_phi, f"{outdir}/bary_z_phi.png", f"Hit maps first layer", xMin=-max_z, xMax=max_z, xLabel="z (mm)", yLabel="Azimuthal angle (deg)")
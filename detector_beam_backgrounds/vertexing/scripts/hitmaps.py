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
    axis_theta = hist.axis.Regular(180, 0, 180, name = "theta", underflow=False, overflow=False)
    axis_phi = hist.axis.Regular(360, 0, 360, name = "phi", underflow=False, overflow=False)
    
    axis_energy = hist.axis.Regular(100, 0, 100, name = "energy")
    axis_distance = hist.axis.Regular(40, 0, 40, name = "distance")

    hist_theta_phi = hist.Hist(axis_theta, axis_phi)
    hist_z_phi = hist.Hist(axis_z, axis_phi)
    
    hist_energy = hist.Hist(axis_theta, axis_energy)
    hist_distance = hist.Hist(axis_energy, axis_distance)

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
                ##### ENERGY-CALC #####
                edep = 1e6*hit.getEDep() # convert to keV
                mc = hit.getMCParticle()
                if mc.getGeneratorStatus() != 1:  # only primaries
                    continue
                x = hit.getPosition().x
                y = hit.getPosition().y
                z = hit.getPosition().z
                ph = functions.phi(x, y) * (180 / math.pi)
                th = functions.theta(x,y,z) * (180 / math.pi)
                r = functions.r(x, y, z)
                hist_z_phi.fill(z, ph)
                hist_theta_phi.fill(th, ph)
                hist_energy.fill(th, edep)
                hist_distance.fill(edep, r)
                
        if i > args.maxFiles:
            break

    # normalize the histograms over the number of events, to get the average number of hits per event

    hist_z_phi /= nEvents
    hist_theta_phi /= nEvents
    
    hist_energy /= nEvents
    hist_distance /=nEvents
    
    hists = {}
    hists['hist_z_phi'] = hist_z_phi
    hists['hist_theta_phi'] = hist_theta_phi
    
    hists['hist_energy'] = hist_energy 
    hists['hist_distance'] = hist_distance

    with open("output_hitmaps.pkl", "wb") as f:
        pickle.dump(hists, f)

if args.plots:

    outdir = "./"
    with open("output_hitmaps.pkl", "rb") as f:
        hists = pickle.load(f)

    # get the histograms
    hist_energy = hists['hist_energy'] #added
    hist_distance = hists['hist_distance'] #added
    
    hist_z_phi = hists['hist_z_phi']
    hist_theta_phi = hists['hist_theta_phi']

    functions.plot_2dhist(hist_z_phi, f"{outdir}/z_phi.png", f"Hit maps first layer", xMin=-max_z, xMax=max_z, xLabel="z (mm)", yLabel="Azimuthal angle (deg)")
    functions.plot_2dhist(hist_theta_phi, f"{outdir}/theta_phi.png", f"Hit maps first layer", xMin=0, xMax=180, xLabel="Theta (deg)", yLabel="Azimuthal angle (deg)")
    
    functions.plot_2dhist(hist_energy, f"{outdir}/theta_energy_hits.png", f"Hit maps theta-energy", xMin=0, xMax=180, xLabel="Theta (deg)", yLabel="Energy (keV)")
    functions.plot_2dhist(hist_distance, f"{outdir}/distance_energy_hits.png", f"Hit maps distance-energy", xMin=0, xMax=40, yMin=10, yMax=40, xLabel="Energy (keV)", yLabel="Distance")
    # plot the projections on phi and theta
    hist_theta = hist_theta_phi[:, sum]
    hist_phi = hist_theta_phi[sum, :]
    functions.plot_hist(hist_theta, f"{outdir}/theta.png", f"Hit maps first layer", xMin=0, xMax=180, xLabel="Theta (deg)")
    functions.plot_hist(hist_phi, f"{outdir}/phi.png", f"Hit maps first layer", xMin=0, xMax=360, xLabel="Azimuthal angle (deg)")
from podio import root_io
import glob
import hist
import functions
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plots', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default="-1", help="Maximum files to run over")
args = parser.parse_args()

##########################################################################################
# this file is for plotting the number of hits over energy deposited in each layer
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23"
files = glob.glob(f"{folder}/*.root")


#layer_radii = [14, 23, 34.5, 141, 316] # IDEA vertex approximate layer radii
layer_radii = [14, 36, 58] # CLD vertex approximate layer radii
hits = {r: [] for r in layer_radii}


if args.calculate:
    axis_energy = hist.axis.Regular(200, 0, 200, name = "energy")
    axis_dedx = hist.axis.Regular(100, 0, 1000, name = "dedx")
    axis_layers = hist.axis.Regular(len(layer_radii), 0, len(layer_radii), name = "layers", underflow=False, overflow=False)

    hist_energy = hist.Hist(axis_layers, axis_energy)
    hist_dedx = hist.Hist(axis_layers, axis_dedx)
    hist_dedx_noSecondary = hist.Hist(axis_layers, axis_dedx)

    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            for hit in event.get("VertexBarrelCollection"):
                radius_idx = functions.radius_idx(hit, layer_radii)

                edep = 1000000*hit.getEDep() # convert to keV
                path_length = hit.getPathLength() # mm
                mc = hit.getMCParticle()

                if mc.getGeneratorStatus() != 1:
                    continue # particles not input into geant

                hist_energy.fill(radius_idx, edep)
                hist_dedx.fill(radius_idx, edep/path_length)
                if not hit.isProducedBySecondary(): # mc particle not tracked
                    hist_dedx_noSecondary.fill(radius_idx, edep/path_length)

        if i > args.maxFiles:
            break

    hists = {}
    hists['hist_energy'] = hist_energy
    hists['hist_dedx'] = hist_dedx
    hists['hist_dedx_noSecondary'] = hist_dedx_noSecondary

    with open("output.pkl", "wb") as f:
        pickle.dump(hists, f)


if args.plots:

    outdir = "./"
    with open("output.pkl", "rb") as f:
        hists = pickle.load(f)

    # get the histograms
    hist_energy = hists['hist_energy']
    hist_dedx = hists['hist_dedx']
    hist_dedx_noSecondary = hists['hist_dedx_noSecondary']

    # plot them for all layers
    for i,r in enumerate(layer_radii):
        # we slice the histograms for each layer
        hist_energy_layer = hist_energy[i,:]
        hist_dedx_layer = hist_dedx[i,:]
        hist_dedx_noSecondary_layer = hist_dedx_noSecondary[i,:]
    
        functions.plot_hist(hist_energy_layer, f"{outdir}/energy_layer{i}.png", f"Energy layer {i}", xMin=0, xMax=100, xLabel="Energy (keV)")
        functions.plot_hist(hist_dedx_layer, f"{outdir}/dedx_layer{i}.png", f"Energy layer {i}", xMin=0, xMax=1000, xLabel="dE/dx (keV/mm)")
        functions.plot_hist(hist_dedx_noSecondary_layer, f"{outdir}/dedx_noSecondary_layer{i}.png", f"Energy layer {i}", xMin=0, xMax=1000, xLabel="dE/dx (keV/mm)")

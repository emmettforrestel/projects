from podio import root_io
import glob
import hist
import functions
import pickle
import argparse
import math
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plot', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()
folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_wz3p6_ee_qq_ecm91p2"
files = glob.glob(f"{folder}/*.root")

# Example approximate radii for CLD
layer_radii = [14, 36, 58]
max_z = 110

if args.calculate:
    # Initialize histogram
    axis_energy = hist.axis.Regular(200, 0, 200, name="energy")
    axis_dR     = hist.axis.Regular(100, 0, 30, name="dR")
    hist_dR_energy = hist.Hist(axis_dR, axis_energy)
    hist_energies = hist.Hist(axis_energy)
    hist_dR_values = hist.Hist(axis_dR)

    # Array of all dR and energy values in all events
    dR_values = []
    energies  = []

    # Collect nEvents to normalize plot later
    nEvents   = 0

    # Go thtough file and get events in each file
    for i, filename in enumerate(files):
        if i >= args.maxFiles:
            break
        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)
        events = podio_reader.get("events")


        for event in events:
            nEvents += 1

            # All hits in first layer in this event
            all_hits = event.get("VertexBarrelCollection")
            good_hits = functions.filter_hits(all_hits, layer_radii)

            # hits_by_mc: key = MCParticle; list of all hits from that particle
            hits_by_mc = {}

            # Go through all hits in this event and determine what MC Particle they're coming from
            # Also find position of max momentum hit in event

            #print("Length of good_hits",len(good_hits))
            for hit in good_hits:
                mc = hit.getMCParticle()

                # If mc not already a key, add it
                #for hit in hits_by_mc:

                if mc.getEnergy() not in hits_by_mc:
                    hits_by_mc[mc.getEnergy()] = []
                
                # Add hit to the list associated with mc
                hits_by_mc[mc.getEnergy()].append(hit)

            
            #print("Length of hits_by_mc",len(hits_by_mc))
            # Print number of keys in this event
            #print("# of MC Particles in event ", nEvents, ": ", len(hits_by_mc))


            # iterates over the key-value pairs in the dictionary
            # mc_particle is each mc particle in the event
            # hits_for_mc is the array of hits associated with that particular mc particle
            for mc_particle, hits_for_mc in hits_by_mc.items():

                #print(len(hits_for_mc))

                # Max momentum hit of this mc particle
                x_of_max_momentum_hit_in_mc = 0
                y_of_max_momentum_hit_in_mc = 0
                z_of_max_momentum_hit_in_mc = 0
                max_momentum = -1

                for hit in hits_for_mc:

                    momentum = functions.r(hit.getMomentum().x, hit.getMomentum().y, hit.getMomentum().z)

                    if momentum > max_momentum:
                        x_of_max_momentum_hit_in_mc = hit.getPosition().x
                        y_of_max_momentum_hit_in_mc = hit.getPosition().y
                        z_of_max_momentum_hit_in_mc = hit.getPosition().z
                        max_momentum = momentum
                
                #if len(hits_for_mc) != 1:
                    #print(len(hits_for_mc))

                for h in hits_for_mc:
                    dx = h.getPosition().x - x_of_max_momentum_hit_in_mc
                    dy = h.getPosition().y - y_of_max_momentum_hit_in_mc
                    dz = h.getPosition().z - z_of_max_momentum_hit_in_mc

                    dR = functions.r(dx, dy, dz)
                    edep = 1e6 * h.getEDep()  # keV

                    if dR != 0:
                        dR_values.append(dR)
                        
                    energies.append(edep)
                    #print(dR)

                    hist_dR_energy.fill(dR, edep)
                    hist_energies.fill(edep)
                    hist_dR_values.fill(dR)
        # end for event in events
    # end for files
        if i > args.maxFiles:
            break


    hist_dR_energy /= nEvents
    hist_energies /= nEvents
    hist_dR_values /= nEvents

    hists = {}
    hists["hist_dR_energy"] = hist_dR_energy
    hists["hist_energies"] = hist_energies
    hists["hist_dR_values"] = hist_dR_values

    with open("output_hitmaps.pkl", "wb") as f:
        pickle.dump(hists, f)



if args.plot:
    outdir = "./"
    with open("output_hitmaps.pkl", "rb") as f:
        hists = pickle.load(f)

    hist_dR_energy = hists['hist_dR_energy'] #added
    hist_energies = hists["hist_energies"]
    hist_dR_values = hists["hist_dR_values"]

    # Plot dR-energy 2D histogram
    functions.plot_2dhist(hist_dR_energy, f"{outdir}/mc__hist_dR_energy_2d.png", f"Hit maps dR-energy", xMin=0, xMax=4, yMin =0, yMax=200, yLabel="Energy (keV)", xLabel="dR")
    
    #Plot dR 1D histogram
    functions.plot_hist(hist_dR_values, f"{outdir}/mc_hist_dR_1d.png", f"Hit map: dR", xMin=0, xMax=10, yMin =0, yMax=22, yLabel="# of hits", xLabel="dR")

    #Plot energy 1D histogram
    functions.plot_hist(hist_energies, f"{outdir}/mc_hist_energy_1d.png", f"Hit map: energy", xMin=0, xMax=200, yMin =0, yMax=3, yLabel="# of hits", xLabel="Energy (keV)")
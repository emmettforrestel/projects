from podio import root_io
import glob
import hist
import functions
import pickle
import argparse
import math
import pprint
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

def barycenter_energy(hits):
    """
    Returns a 2D array of particles with their associated energy, barycenter.x, barycenter.y, barycenter.z as an array:
    particle_energy_barycenter[0][0] = particle energy
    particle_energy_barycenter[0][1] = x-coord of barycenter of hits
    particle_energy_barycenter[0][2] = y-coord of barycenter of hits
    particle_energy_barycenter[0][3] = z-coord of barycenter of hits
    """ 
    particle_energy_barycenter = []
    # hits_by_mc: key = MCParticle; list of all hits from that particle
    hits_by_mc = {}
    for hit in hits:
        mc = hit.getMCParticle()
            # If mc not already a key, add it
        if mc.getEnergy() not in hits_by_mc:
            hits_by_mc[mc.getEnergy()] = []    
        hits_by_mc[mc.getEnergy()].append(hit)
           
    for mc_particle_energy, hits_for_mc in hits_by_mc.items():
        b = functions.geometric_baricenter(hits_for_mc)
        particle_energy_barycenter.append([mc_particle_energy, b[0], b[1], b[2]])
    print(particle_energy_barycenter)
    return(particle_energy_barycenter)

if args.calculate:
    # Initialize histogram
    axis_energy = hist.axis.Regular(200, 0, 200, name="energy")
    axis_dR = hist.axis.Regular(100, 0, 30, name="dR")
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
            barycenter_energy(good_hits)
            
            for hit in good_hits:
                mc = hit.getMCParticle()

                # If mc not already a key, add it
                #for hit in hits_by_mc:

                if mc.getEnergy() not in hits_by_mc:
                    hits_by_mc[mc.getEnergy()] = []
                
                # Add hit to the list associated with mc
                hits_by_mc[mc.getEnergy()].append(hit)

        
            # iterates over the key-value pairs in the dictionary
            # mc_particle is each mc particle in the event
            # hits_for_mc is the array of hits associated with that particular mc particle
            for mc_particle, hits_for_mc in hits_by_mc.items():

                #b = array of position: b[0] = b.x, b[1] = b.y, b[2] = b.z
                b = functions.geometric_baricenter(hits_for_mc)
                
                for h in hits_for_mc:
                    dx = h.getPosition().x - b[0]
                    dy = h.getPosition().y - b[1]
                    dz = h.getPosition().z - b[2]

                    dR = functions.r(dx, dy, dz)
                    edep = 1e6 * h.getEDep()  # keV

                    if dR != 0:
                        dR_values.append(dR)
                        
                    energies.append(edep)

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
    functions.plot_2dhist(hist_dR_energy, f"{outdir}/barycenter_hist_dR_energy_2d.png", f"Hit maps dR-energy", xMin=0, xMax=2, yMin =0, yMax=50, yLabel="Energy (keV)", xLabel="dR")
    
    #Plot dR 1D histogram
    functions.plot_hist(hist_dR_values, f"{outdir}/barycenter_hist_dR_1d.png", f"Hit map (barycenter): dR", xMin=0, xMax=10, yMin =0, yMax=16, yLabel="# of hits", xLabel="dR")

    #Plot energy 1D histogram
    functions.plot_hist(hist_energies, f"{outdir}/barycenter_hist_energy_1d.png", f"Hit map (barycenter): energy", xMin=0, xMax=200, yMin =0, yMax=3, yLabel="# of hits", xLabel="Energy (keV)")
      
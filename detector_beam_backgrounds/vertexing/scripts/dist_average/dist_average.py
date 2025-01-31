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
    axis_x = hist.axis.Regular(200, -100, 100, name = "x")
    axis_y = hist.axis.Regular(200, -100, 100, name = "y")
    axis_z = hist.axis.Regular(200, -100, 100, name = "z")
    hist_xy = hist.Hist(axis_x, axis_y)
    hist_xz = hist.Hist(axis_x, axis_z)
    X_values =[]
    Y_values = []
    Z_values = []
    
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
                x = hit.getPosition().x
                y = hit.getPosition().y
                z = hit.getPosition().z
                #ph = functions.phi(x, y) #* (180 / math.pi)
                #th = functions.theta(x,y,z)# * (180 / math.pi)
                #pseudorapidity = functions.pseudorapiditiy(th)
                #dR = functions.dR(ph, th)
                hist_xy.fill(x, y)
                hist_xz.fill(x, z)
                X_values.append(x)
                Y_values.append(y)
                Z_values.append(z)
                
        if i > args.maxFiles:
            break
        
    hist_xy /=nEvents
    hist_xz /=nEvents
    hists = {}    
    hists['hist_xy'] = hist_xy
    hists['hist_xz'] = hist_xz

    with open("output_hitmaps.pkl", "wb") as f:
        pickle.dump(hists, f)

if args.plots:
    outdir = "./"
    with open("output_hitmaps.pkl", "rb") as f:
        hists = pickle.load(f)
        
    plt.figure(figsize=(8,6))
    plt.scatter(X_values, Y_values, color='blue', s=1, alpha=0.7, label='Hits')
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.xlabel('X') 
    plt.ylabel('Y') 
    plt.title('Scatter Plot of X vs Y') 
    plt.grid(True) 
    plt.legend() 
    plt.savefig("scatter_XY.png")
    
    plt.figure(figsize=(8,6))
    plt.scatter(X_values, Z_values, color='blue', s=1, alpha=0.7, label='Hits')
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.xlabel('X') 
    plt.ylabel('Z') 
    plt.title('Scatter Plot of X vs Z') 
    plt.grid(True) 
    plt.legend() 
    plt.savefig("scatter_XZ.png") 
    plt.show()

    hist_xy = hists['hist_xy'] #added
    hist_xz = hists['hist_xz'] #added
    functions.plot_2dhist(hist_xy, f"{outdir}/hist_xy.png", f"hist_xy", xMin=-100, xMax=100, yMin=-100, yMax=100, yLabel="Y", xLabel="X")
    functions.plot_2dhist(hist_xz, f"{outdir}/hist_xz.png", f"hist_xz", xMin=-100, xMax=100, yMin=-100, yMax=100, yLabel="Z", xLabel="X")
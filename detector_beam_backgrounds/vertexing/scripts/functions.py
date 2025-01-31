from podio import root_io
import ROOT
import math
import numpy as np
import os
import matplotlib.pyplot as plt
import mplhep as hep
import pickle
import hist
hep.style.use(hep.style.ROOT)
def r(x,y,z):
    """
    Calculates the magnitude of the radius.
    """
    return math.sqrt(x**2 + y**2 + z**2)
def r_comp(x1, y1, z1, x2, y2, z2):
    """
    Computes the distances between two points.
    """
    x_dif = x2 - x1
    y_dif = y2 - y1
    z_dif = z2 - z1
    return math.sqrt(x_dif**2 + y_dif**2 + z_dif**2)
def pseudorapiditiy(theta):
    return -math.log(math.tan((theta/2)))

def phi(x,y):
    """
    Calculates phi of particle.
    Inputs: x,y floats.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    phi = math.atan(y/x)
    if x < 0:
        phi +=  math.pi
    elif y < 0:
        phi += 2*math.pi
    return phi
def theta(x,y,z):
    """
    Calculates theta of particle.
    Inputs: x,y,z floats.
    Output: theta, float representing angle in radians from 0 to pi.
    """
    return math.acos(z/np.sqrt(x**2 + y**2 + z**2))
def radius_idx(hit, layer_radii):
    """
    Calculates polar radius of particle.
    Inputs: hit, SimTrackerHit object.
    Output: r, int representing polar radius in mm.
    """
    true_radius = hit.rho()
    for i,r in enumerate(layer_radii):
        if abs(true_radius-r) < 4:
            return i
    raise ValueError(f"Not close enough to any of the layers {true_radius}")

def filter_hits(hits, layer_radii):
    """
    Filter hits:
    - only on the first layer (radius_idx == 0),
    - only primary particles ( mc.getGeneratorStatus() == 1 ),
    - remove secondaries (hit.isProducedBySecondary() == False),
    etc.
    """
    filtered = []
    for h in hits:
        rad = radius_idx(h, layer_radii)
        if rad != 0:  # only first layer
            continue
        #if h.isProducedBySecondary():
            #continue
        filtered.append(h)
    return filtered


def geometric_baricenter(hits):
    """
    Returns energy barycenter of hits for one mc
    """
    sum_energy = 0.0
    sum_energy_x = 0.0
    sum_energy_y = 0.0
    sum_energy_z = 0.0

    for h in hits:
        edep = 1e6 * h.getEDep()
        x = h.getPosition().x
        y = h.getPosition().y
        z = h.getPosition().z

        sum_energy += edep
        sum_energy_x += edep*x
        sum_energy_y += edep*y
        sum_energy_z += edep*z

    b_x = sum_energy_x/sum_energy
    b_y = sum_energy_y/sum_energy
    b_z = sum_energy_z/sum_energy
   
    return [sum_energy, b_x, b_y, b_z]

def max_momentum_hit_cell(hits_for_mc):
    """
    returns max momentum of hits for an mc and cell id of max momentum hit
    """
    cell_of_max = None
    max_momentum = -1
    for hit in hits_for_mc:
        momentum = r(hit.getMomentum().x, hit.getMomentum().y, hit.getMomentum().z)
        if momentum > max_momentum:
            max_momentum = momentum
            cell_of_max = hit.getCellID()
    return([max_momentum, cell_of_max])

def barycenter_energy(hits):
    """
    Returns a 2D array of particles with their associated energy, barycenter.x, barycenter.y, barycenter.z as an array:
    particle_energy_barycenter[0][0] = particle energy
    particle_energy_barycenter[0][1] = summed energy of barycenter of hits
    particle_energy_barycenter[0][2] = x-coord of barycenter of hits
    particle_energy_barycenter[0][3] = y-coord of barycenter of hits
    particle_energy_barycenter[0][4] = z-coord of barycenter of hits
    particle_energy_barycenter[0][5] = count of microhits
    particle_energy_barycenter[0][6] = cellID of max_momentum hit
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
        b = geometric_baricenter(hits_for_mc)
        max_cellid = max_momentum_hit_cell(hits_for_mc) 
        particle_energy_barycenter.append([mc_particle_energy, b[0], b[1], b[2], b[3], len(hits_for_mc), max_cellid[1]])
    return(particle_energy_barycenter)

def plot_hist(h, outname, title, xMin=-1, xMax=-1, yMin=-1, yMax=-1, xLabel="", yLabel="Events", logY=False):
    fig = plt.figure()
    ax = fig.subplots()
    hep.histplot(h, label="", ax=ax, yerr=False)
    ax.set_title(title)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    #ax.legend(fontsize='x-small')
    if logY:
        ax.set_yscale("log")
    if xMin != -1 and xMax != -1:
        ax.set_xlim([xMin, xMax])
    if yMin != -1 and yMax != -1:
        ax.set_ylim([yMin, yMax])
    fig.savefig(outname, bbox_inches="tight")
    fig.savefig(outname.replace(".png", ".pdf"), bbox_inches="tight")
def plot_2dhist(h, outname, title, xMin=-1, xMax=-1, yMin=-1, yMax=-1, xLabel="", yLabel="Events", logY=False):
    fig = plt.figure()
    ax = fig.subplots()
    hep.hist2dplot(h, label="", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    #ax.legend(fontsize='x-small')
    if logY:
        ax.set_yscale("log")
    if xMin != -1 and xMax != -1:
        ax.set_xlim([xMin, xMax])
    if yMin != -1 and yMax != -1:
        ax.set_ylim([yMin, yMax])
    fig.savefig(outname, bbox_inches="tight")
    fig.savefig(outname.replace(".png", ".pdf"), bbox_inches="tight")

def fill_MC_Particles_collection(hit, mc, MCParticles_collection):

                mc = hit.getMCParticle()
                #if mc.getGeneratorStatus() != 1: #(mc is None) or
                #    continue
                if mc not in hits_by_mc:
                    hits_by_mc[mc] = []
                hits_by_mc[mc].append(hit)

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
        b = geometric_baricenter(hits_for_mc)
        particle_energy_barycenter.append([mc_particle_energy, b[0], b[1], b[2]])
    return(particle_energy_barycenter)

def geometric_baricenter(hits):
    """
    Iterates over an array of hits and calculates their barycenter.
    
    It returns an array where 
    [0] is the x-coord of the barycenter
    [1] is the y-coord of the barycenter
    [2] is the z-coord of the barycenter
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
   
    return [b_x, b_y, b_z]
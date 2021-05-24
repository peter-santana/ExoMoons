# Simulation of an exomoon around a gaint planet in an eccentric orbit
# Function returns orbit data for star, planet, and moon. 
# PHL @ UPR Arecibo
# Uses REBOUND: https://github.com/hannorein/rebound

import rebound
import numpy as np

def exomoon_orbit1(mass_star = 1.0, mass_planet = 1.0, a_planet = 1.0,
                   e_planet = 0.0, mass_moon = 1.0, a_moon = 1.0, e_moon = 0.0,
                   n_orbits = 100, n_points_orbit = 100):

    # constants

    mass_jupiter = 9.55e-4 # solar units
    mass_earth = 3e-6 # solar units
    r_earth_moon = 0.002687 # AU
    
    # initialize simulation variables

    n = n_orbits*n_points_orbit
    
    torb = 2*np.pi*np.sqrt(a_planet**3/mass_star)
    times = np.linspace(0, n_orbits*torb, n)

    class data:
        
        # planet position
        xp = np.zeros(n)
        yp = np.zeros(n)
        a = np.zeros(n)
        e = np.zeros(n)        

        # star position
        xs = np.zeros(n)
        ys = np.zeros(n)

        # moon position
        x = np.zeros(n)
        y = np.zeros(n)
        d = np.zeros(n)
        
        t = times
    
    # initialize simulation

    sim = rebound.Simulation()

    sim.integrator = 'ias15' # orbit integrator
    sim.dt = 0.1 # time step

    m_star = mass_star
    m_planet = mass_planet*mass_jupiter
    m_moon = mass_moon*mass_earth
    ar_moon = a_moon*r_earth_moon

    r_planet = a_planet*(1. + e_planet)
    vy_planet = np.sqrt(m_star*(2./r_planet - 1./a_planet))

    r_moon = ar_moon*(1. + e_moon)
    vy_moon = np.sqrt(m_planet*(2./r_moon - 1./ar_moon))

    sim.add(m = m_star, hash = 'star') # star
    sim.add(m = m_planet, x = r_planet, vy = vy_planet, hash = 'jovian') # planet
    sim.add(m = m_moon, x = r_planet + r_moon, vy = vy_planet + vy_moon, hash = 'moon') # moon

    # sim.status()
    
    sim.move_to_com() # move coordinates to center of mass

    print('-' * (4*15+3))
    print('%15s %15s %15s %15s' % ('Orbit', 'Time (yrs)', 'a (AU)', 'd (LD)' ))
    print('-' * (4*15+3))
    
    for i, time in enumerate(times):
        sim.integrate(time)
        data.a[i] = sim.particles['jovian'].a
        data.e[i] = sim.particles['jovian'].e
        data.x[i] = sim.particles['moon'].x
        data.y[i] = sim.particles['moon'].y
        data.xp[i] = sim.particles['jovian'].x
        data.yp[i] = sim.particles['jovian'].y
        data.xs[i] = sim.particles['star'].x
        data.ys[i] = sim.particles['star'].y
        data.d[i] = np.sqrt((data.x[i]-data.xp[i])**2+(data.y[i]-data.yp[i])**2)/r_earth_moon
        if ((i/(n/10)) % 1) == 0:
            print('%15i %15i %15.2f %15.2f' % (i / n_points_orbit,times[i]/(2*np.pi),data.a[i],data.d[i]))

    print('%15i %15i %15.2f %15.2f' % (n_orbits,times[i]/(2*np.pi),data.a[i],data.d[i]))
    print('-' * (4*15+3))
    
    # converted to eart-moon distance
    # data.d = np.sqrt((data.x-data.xp)**2+(data.y-data.yp)**2)/r_earth_moon 
    
    return data()
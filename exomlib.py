import exomlib
import numpy as np
import matplotlib.pyplot as plt


# main simulation code
# HD 20782 b

# star properties

mass_star = 0.34 # solar masses

# planet properties

mass_planet = 0.3 # jupiter masses
a_planet = 0.166 # AU, astronomical units
e_planet = 0.38 # eccentricity

# moon properties

#a_moon = 1.0 # lunar distances
e_moon = 0.0 # moon eccentricity

rh = (a_planet/0.002695)*(1-e_planet)*(mass_planet/(3.*mass_star*1048.))**(1/3.)
rhp = (a_planet/0.002695)*(1-e_planet)*(mass_planet/((3.+e_planet)*mass_star*1048.))**(1/3.)
rha = (a_planet/0.002695)*(1+e_planet)*(mass_planet/((3.-e_planet)*mass_star*1048.))**(1/3.)

mass_moon_values = [0.1, 0.5, 1.0, 5.0, 10.0]

print('Planet Mass   = ', mass_planet, 'Mj')
print('Planet eccen. = ', e_planet)
print('Star Mass     = ', mass_star, 'Ms')
print('Hill Sphere   = ', round(rh,2), 'LD')
print('Hill Sphere P = ', round(rhp,2), 'LD')
print('Hill Sphere A = ', round(rha,2), 'LD')
print()

for mass_moon in mass_moon_values: # earth masses
    
    print()
    print('Moon Mass = ', mass_moon, ' Me')
    
    rh_moon_values = [1.0, 0.80, 0.60, 0.40, 0.20]
    
    for dist in rh_moon_values: #lunar distances
        
        a_moon = dist*rhp
        
        # integration properties

        n_orbits = 1000000
        n_points_orbit = 10

        print()
        print('m_moon = ', mass_moon, 'Me, a_moon = ', dist, 'RH (', round(a_moon,2), 'LD)')
        data = exomlib.exomoon_orbit1(mass_star = mass_star, 
                                  mass_planet = mass_planet, a_planet = a_planet, e_planet = e_planet, 
                                  mass_moon = mass_moon, a_moon = a_moon, e_moon = e_moon, 
                                  n_orbits = n_orbits, n_points_orbit = n_points_orbit)

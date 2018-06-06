# PCA-for-sattelite
My course project at institute. I should to create a PCA algorhitm in order to estimate an orbit of satellite

based on apriory analysis

__Model__

A sattelite which moves around the Earth orbit at _normal Gravity Field_ 

Initial data are specified in the osculating elements

1. Orbital Inclination is 42 degrees

2. Semimajor axis we can count from specified axes of ellips 

        re = 6371  # Earth radius, km
        h_pi = 21000
        h_alpha = 970
        r_pi = h_pi + re
        r_alpha = h_alpha + re
        self.a = (r_pi + r_alpha) / 2  # Semimajor axis(km)
        
3.


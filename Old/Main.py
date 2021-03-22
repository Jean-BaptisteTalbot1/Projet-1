# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 14:29:29 2021

@author: jbtalbot
"""
import math


N1 = 1750 # Newton
H01 = 5000 # Heures

param_fiabilite = {
    "defaut": [4.48, 0, 1.5],
    "autre" : [4.439, 0.02, 1.483]
        }
p = "defaut"

# Fiabilite
R = 0.99

#%% Charge dynamique

# Billes 
a = 3
# Manufacturier
theta   = param_fiabilite[p][0]
x0      = param_fiabilite[p][1]
b       = param_fiabilite[p][2]

# Durée de vie pour une fiabilité 1%
L01 = 60 * H01 * N1
L10 = L01 / (x0 + theta*(math.log(1/R))**(1/b))

# Force radiale
Fe = 26.7328 #TODO : prendre les valeurs du design 

# Charge dynamique
C10 = Fe * (L10/10**6)**(1/a)


# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 21:22:59 2021

@author: jbtalbot
"""

import math
from array import *
from sympy import solve, symbols, Symbol, Eq


# R : fiabilité du roulement
# L1R : Durée pour une fiabilité R
# L10 : Durée pour une fiabilité 90%
# p : pouvant être 1 ou 2 pour le premier jeu de données ou le second


param_fiabilite = {
    "defaut": [4.48, 0, 1.5],
    "autre" : [4.439, 0.02, 1.483]
        }



#%% Relation charge-durée de vie

def fiabilite(L1R, L10, p):
    # Dépendamment des manufacturiers, les valeurs suivantes sont suggérées
    theta   = param_fiabilite[p][1]
    x0      = param_fiabilite[p][2]
    b       = param_fiabilite[p][3]
    
    try:
        R = math.exp( -pow((L1R/L10 - x0) / (theta - x0), b) )
        return R
    except ValueError:
        print("Erreur dans le calcul de la fiabilité R")
                
        
def chargeDureeVie(L1R,R,p):
    
    # Dépendamment des manufacturiers, les valeurs suivantes sont suggérées
    theta   = param_fiabilite[p][1]
    x0      = param_fiabilite[p][2]
    b       = param_fiabilite[p][3]
        
    try:
        R = math.exp( -pow((L1R/L10 - x0) / (theta - x0), b) )
        return R
    except ValueError:
        print("Erreur dans le calcul de la dureeVie90")
     
        
     
def chargeDyn_C10(Fe, typeRoulement):
    
    if (typeRoulement == "billes"):
        a = 3
    elif (typeRoulement == "rouleaux"):
        a = 10/3
    else:
        print("erreur de type de roulements")
        return 1
    
    C10 = Fe * pow(L10 / pow(10,6), 1/a)
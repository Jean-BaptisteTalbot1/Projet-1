# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:44 2021

@author: Catherine
"""

import math
from math import pi
from sympy import Eq, solve, symbols 


# Nettoyage de la console
print("\033[H\033[J")

##############################################################################
##############################################################################
#%%                                                                        ###
###                     ENTRÉES OPÉRATEUR                                  ###
##############################################################################
##############################################################################

# Diamètre de la poulie 1, soit la poulie pour la courroir du moteur (en mm)
D_poulie1       = 10 # mm
# Distance entre la poulie du moteur et le boîtier
l1              = 5 # mm

# Diamètre de la poulie 2, soit la poulie de sortie du boîtier (en mm)
D_poulie2       = 10 # mm
# Distance entre la poulie de sortie et le boîtier (en mm)
l2              = 5 # mm

# Angle d'orientation des brins des courroies par rapport au boîtier (en rad)
theta_courroie  = math.radians(45) # rad

# Puissance du moteur d'entrée (en hp)
P_hp            = 1 #hp

# Vitesse de rotation du moteur (en tr/min)
n1              = 1750  # tr/min

##############################################################################
##############################################################################
##############################################################################



##############################################################################
#%% CONSTANTES

#Acier AISI 1035 EF - Matériau ductile
Sut     = 550   #MPa
Sy      = 460   #MPa
Kf      = 3     #CC en fatigue
Kfs     = 3     #CC en fatigue
FS      = 3
ratio   = 4
hp2watts = 745.7
##############################################################################


##############################################################################
#%% Vitesse de sortie n2
n2 = n1/ratio
  

##############################################################################
#%% Conversion puissance en watts
P_watt = hp2watts * P_hp


##############################################################################
#%% Couples
# Couple T1
T1 = (60 * P_watt) / (2*pi*n1) * 1000 # Nmm
# Couple T2
T2 = (60 * P_watt) / (2*pi*n2) * 1000 # Nmm

##############################################################################
#%% Entraxe, diamètres des engrenages

# Entraxe a
a = 15.4 * pow(T2/1000,0.337) # a en mm - T2 en Nmm (Nm * 1000)

# Dimensions des engrenages (pignon et roue)
# Sachant que le ratio de vitesse est de 4, le diamètre suit le même ratio
d_p, d_r = symbols('d_p d_r')

#(d_pignon/2) + (d_roue/2) == a
eq1 = Eq(- a + d_p/2 + d_r/2) 
# d_roue == 4*d_pignon
eq2 = Eq(4 * d_p - d_r)

resultDiametreEngrenage = solve((eq2,eq1),(d_r,d_p))
d_pignon    = resultDiametreEngrenage[d_p]
d_roue      = resultDiametreEngrenage[d_r]

l = 0.5 * d_pignon

del d_p, d_r, eq1, eq2
##############################################################################
#%% DCL - Détermination des forces et des moments

# Puisque Rcy = Rdy ET Rcz = Rdz ET Ray = Rby ET Raz = Rbz :
Ry, Rz, Fr, Ft = symbols('Ry Rz Fr Ft')

# Relation de la force transmise et radiale
FrFt    = Eq(Fr/Ft-math.tan(math.radians(20)))
# Somme des forces en Y
Fy      = Eq(2*(Ry) - Fr)                    
# Somme des forces en Z
Fz      = Eq(2*(Rz) - Ft)                
# Somme des moments sur l'axe des X 
Mx      = Eq(-T1 + Ft * (d_pignon / 2))

# Résolution des équations du système
DCL     = solve((Fy,Fz,FrFt,Mx) , (Ry, Rz, Fr,Ft))


print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("     PREMIÈRE PARTIE     ")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("")
print("-------------------------------------------")
print("Données du système de la première partie")
print("")
for F in DCL:
    print(str(F) + " : " + str(DCL[F]) + " N")


# Saut de ligne
print("-------------------------------------------")
print("")

Ry = round(DCL[Ry],5)
Rz = round(DCL[Rz],5)
Ft = round(DCL[Ft],5)
Fr = round(DCL[Fr],5)

My = (l * Ry) #moment au point critique
Mz = (l * Rz) #moment au point critique

M_resultant = math.sqrt(My**2 + Mz**2)



##############################################################################
#Largeur des engrenages w
w = 0.75*d_pignon #

##############################################################################
#%% Courbe SN pour le diamètre requis pour l'arbre CD

Se_CD, Se_AB = symbols('Se_CD Se_AB')
  
# Vérification des conditions pour le Se_Prime
if (Sut < 1400):
    Se_prime = 0.5*Sut
else:
    Se_prime = 700

# Paramètre a,c,d et e selon l'énoncé
ka = 4.51*pow(Sut,-0.265) #(usinage, laminage à froid)
kc = 1 #chargement combiné
kd = 1 #pas d'info donc 20°C
ke = 0.814 #fiabilité de 99%

# Valeur tampon de départ conservateur juste au dessus de la valeur pour entrer
# dans la boucle WHILE
d_buffer_CD = 255
d_ArbreCD = 254

# Tant que la différence entre la valeur précédente (buffer) et la nouvelle
# (d_Arbre) est plus grande que le seuil, une nouvelle itération est lancée
print("-------------------------------------------")
print("Itérations pour le diamètre de l'arbre CD :")
print("")
while (math.fabs(d_buffer_CD - d_ArbreCD >= 10**(-6))):
    
    # Mise en mémoire de la valeur de départ
    d_buffer_CD = d_ArbreCD
    
    # Calcul du paramètre kb selon le diamètre d'entrée
    if(d_ArbreCD <= 51):
        kb = 1.24 * d_ArbreCD**-0.107
    else:
        kb = 1.51 * d_ArbreCD**-0.157 #d_ArbreCD inconnu donc pire cas à d=254
    
    # Calcul du Se
    Se_CD = ka*kb*kc*kd*ke*Se_prime
    
    # Calcul du diamètre requis. Si la différence entre cette nouvelle valeur
    # et le buffer du départ est plus grande que le seuil, l'algorithme est 
    # relancé
    d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)
    print(str(d_ArbreCD) + " mm")

# Saut de ligne
print("-------------------------------------------")
print("")

##############################################################################
#%% Courbe SN pour le diamètre requis pour l'arbre AB

# Valeur tampon de départ conservateur juste au dessus de la valeur pour entrer
# dans la boucle WHILE
d_buffer = 255
d_ArbreAB = 254

# Tant que la différence entre la valeur précédente (buffer) et la nouvelle
# (d_Arbre) est plus grande que le seuil, une nouvelle itération est lancée
print("-------------------------------------------")
print("Itérations pour le diamètre de l'arbre AB :")
print("")
while (math.fabs(d_buffer - d_ArbreAB >= 10**(-6))):
    
    # Mise en mémoire de la valeur de départ
    d_buffer = d_ArbreAB
    
    # Calcul du paramètre kb selon le diamètre d'entrée
    if(d_ArbreAB <= 51):
        kb = 1.24 * d_ArbreAB**-0.107
    else:
        kb = 1.51 * d_ArbreAB**-0.157 #d_ArbreCD inconnu donc pire cas
    
    # Cacul du Se
    Se_AB = ka*kb*kc*kd*ke*Se_prime
    
    # Calcul du diamètre requis. Si la différence entre cette nouvelle valeur
    # et le buffer du départ est plus grande que le seuil, l'algorithme est 
    # relancé
    d_ArbreAB = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_AB) + (math.sqrt(3)*Kfs*T2 / Sut)) )**(1/3)
    print(str(d_ArbreAB) + " mm")
    
# Saut de ligne
print("-------------------------------------------")
print("")
    
##############################################################################
#%%Type roulements

# R : fiabilité du roulement
# L1R : Durée pour une fiabilité R
# L10 : Durée pour une fiabilité 90%
# p : pouvant être 1 ou 2 pour le premier jeu de données ou le second


#%% Relation charge-durée de vie

# Dépendamment des manufacturiers, les valeurs suivantes sont suggérées
theta       = 4.439
x0          = 0.02
b           = 1.483
a_roul      = 3     #billes
Xi          = 1     #billes à gorges profondes
V           = 1     #bague interne du roulement tourne
H01         = 5000  #heures 

# Fiabilite
R = 0.99

#%% Charge dynamique CD

# Durée de vie pour une fiabilité 99%
L01_CD  = 60 * H01 * n1
L10     = L01_CD / (x0 + theta*(math.log(1/R))**(1/b))

# Force radiale
Fr = math.sqrt((Ry**2) + (Rz**2))

#Force radiale éqivalente
Fe = Xi*V*Fr

# Charge dynamique
C10_CD = (Fe * (L10/10**6)**(1/a_roul)) / 1000 #kN

print("-------------------------------------------")
print("Charges dynamiques des roulements :")
print("")
print("C10_CD : " + str(C10_CD) + " kN")

#%% Charge dynamique CD

# Durée de vie pour une fiabilité 99%
L01_AB  = 60 * H01 * n2
L10     = L01_AB / (x0 + theta*(math.log(1/R))**(1/b))

# Force radiale
Fr = math.sqrt((Ry**2) + (Rz**2))

#Force radiale éqivalente
Fe = Xi*V*Fr 

# Charge dynamique
C10_AB = (Fe * (L10/10**6)**(1/a_roul)) / 1000 #kN

print("C10_AB : " + str(C10_AB) + " kN")


print("-------------------------------------------")


#%%











print("")
print("")

print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("     DEUXIÈME PARTIE     ")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")

##############################################################################
##############################################################################
#%%                                 PARTIE 2
##############################################################################
##############################################################################
#%% Forces poulies
F2 = (15 * P_watt)/(n1 * (D_poulie1*(10**-3)) * math.pi)
F1 = 5 * F2

F4 = (15 * P_watt)/(n1 * (D_poulie2*(10**-3)) * math.pi)
F3 = 5 * F4


##############################################################################
#%% Torques poulies

T_entree = P_watt / ((2*math.pi*n1)/60)
T_sortie = T_entree*4

##############################################################################
##%% Calcul équilibre des forces arbre CD

# Arbre CD
Rcy, Rdy, Rcz, Rdz = symbols('Rcy Rdy Rcz Rdz')

# Somme des forces en Y
FyCD = Eq(Rcy + Rdy - Fr + (F1+F2)*math.cos(theta_courroie))
# Somme des forces en Z
FzCD = Eq(Rcz + Rdz + Ft + (F1+F2)*math.sin(theta_courroie))
# Moment par rapport à D
Myd = ((l * Ft) + ((F1+F2)*math.cos(theta_courroie)*(l1 + 2*l)) - Rcz*l*2)
# moment par rapport à D
Mzd = ((l * Fr) + ((F1+F2)*math.sin(theta_courroie)*(l1 + 2*l)) - Rcy*l*2)

# Résolution des équations du système pour l'arbre CD
DCL = solve((FyCD,FzCD,Myd,Mzd), (Rcy,Rdy,Rcz,Rdz))

##############################################################################
##%% Calcul équilibre des forces arbre AB

# Arbre AB
Ray, Rby, Raz, Rbz = symbols('Ray Rby Raz Rbz')

# Somme des forces en Y
FyAB = Eq(Ray + Rby - Fr + (F3+F4)*math.cos(theta_courroie))
# Somme des forces en Z
FzAB = Eq(Raz + Rbz + Ft + (F3+F4)*math.sin(theta_courroie))
# Moment par rapport à A
Mya = (-(l * Ft) - ((F3+F4)*math.cos(theta_courroie)*(l2 + 2*l)) + Rbz*l*2)
# Moment par rapport à A
Mza = (-(l * Fr) - ((F3+F4)*math.sin(theta_courroie)*(l2 + 2*l)) + Rby*l*2)

# Résolution des équations du système pour l'arbre AB
DCL1 = solve((FyAB,FzAB,Mya,Mza), (Ray,Rby,Raz,Rbz))

# Saut de ligne
print("")
print("-------------------------------------------")

print("Valeur du second système : ")
print("")
for F in DCL:
    print(str(F) + " : " + str(DCL[F]) + " N")

Rdz = round(DCL[Rdz],5)
Rdy = round(DCL[Rdy],5)
Rcz = round(DCL[Rcz],5)
Rcy = round(DCL[Rcy],5)

for F in DCL1:
    print(str(F) + " : " + str(DCL1[F]) + " N")

Rbz = round(DCL1[Rbz],5)
Rby = round(DCL1[Rby],5)
Raz = round(DCL1[Raz],5)
Ray = round(DCL1[Ray],5)


##############################################################################
#%% Courbe SN

# Calcul de Se pour arbre CD
# Calcul du paramètre kb selon le diamètre d'entrée
if(d_ArbreCD <= 51):
    kb = 1.24 * d_ArbreCD**-0.107
else:
    kb = 1.51 * d_ArbreCD**-0.157 

Se_CD = ka*kb*kc*kd*ke*Se_prime


# Calcul de Se pour arbre AB
# Calcul du paramètre kb selon le diamètre d'entrée
if(d_ArbreAB <= 51):
    kb = 1.24 * d_ArbreAB**-0.107
else:
    kb = 1.51 * d_ArbreAB**-0.157

Se_AB = ka*kb*kc*kd*ke*Se_prime

##############################################################################
#%% Calcul des moments résultants

# Moment au point critique pour arbre CD
MyCD = ((F1+F2)*math.cos(theta_courroie) * ((2*l) + l1)) # Nmm
MzCD = ((F1+F2)*math.sin(theta_courroie) * ((2*l) + l1)) # Nmm

# Moment au point critique pour arbre CD
MyAB = ((F3+F4)*math.cos(theta_courroie)* ((2*l) + l2)) # Nmm
MzAB = ((F3+F4)*math.sin(theta_courroie)* ((2*l) + l2)) # Nmm


# Moment résultant pour les deux arbres
M_resultant1 = math.sqrt(MyCD**2 + MzCD**2) # Nmm
M_resultant2 = math.sqrt(MyAB**2 + MzAB**2) # Nmm


##############################################################################
#%% Calcul des Facteurs de sécurité

# Arbre CD
FS1 = ((d_ArbreCD**3)*Se_CD*Sut*math.pi)/(16*((Kfs*Se_CD*math.sqrt(3)*(T_entree*1000))+(2*M_resultant1*Kf*Sut)))

# Arbre AB
FS2 = ((d_ArbreAB**3)*Se_AB*Sut*math.pi)/(16*((Kfs*Se_AB*math.sqrt(3)*(T_sortie*1000))+(2*M_resultant2*Kf*Sut)))


# Saut de ligne
print("-------------------------------------------")
print("")
print("-------------------------------------------")

print("Facteurs de sécurité des arbres")
print("")
print("FS1 : " +str(FS1))
print("FS2 : " +str(FS2))
print("-------------------------------------------")

##############################################################################
#%% Print durée de vie roulement

# Saut de ligne
print("")
print("-------------------------------------------")
print("Durée de vie des roulements")
print("")
print("L01_CD : " + str(int(L01_CD)) + " cycles") #cycles
print("L01_AB : " + str(int(L01_AB)) + " cycles") #cycles
print("-------------------------------------------")

   
   







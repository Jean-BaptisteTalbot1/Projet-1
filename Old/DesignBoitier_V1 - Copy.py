# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:44 2021

@author: Catherine
"""

import math
from sympy import Eq, solve, symbols 


##############################################################################
#%% Variables

#Acier AISI 1035 EF - Matériau ductile
Sut     = 550   #MPa
Sy      = 460   #MPa
Kf      = 3     #CC en fatigue
Kfs     = 3     #CC en fatigue
FS      = 3

##############################################################################
#%%Interface utilisateur

D_poulie1       = 10 #mm
D_poulie2       = 10 #mm
l1              = 5 #mm
l2              = 5 #mm
theta_courroie  = math.radians(45) #rad
P_hp            = 1 #hp
n1              = 1750  #tr/min - moteur

##############################################################################
#%% Vitesse de sortie n2
n2 = n1/4
  

##############################################################################
#%% Conversion puissance en watts
P_watt = 745.7*P_hp


##############################################################################
#%% Couples
# Couple T1
T1 = (60*P_watt)/(2*math.pi*n1) * 1000 # Nmm
# Couple T2
T2 = (60*P_watt)/(2*math.pi*n2) * 1000 # Nmm

T1 = round(T1,5)
T2 = round(T2,5)

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
# RcdY, RabY, RcdZ, RabZ

Ry, Rz, Fr, Ft = symbols('Ry Rz Fr Ft')

FrFt    = Eq(Fr/Ft-math.tan(math.radians(20)))
Fy      = Eq(4*(Ry) - Fr)    # Somme des forces en Y
Fz      = Eq(4*(-Rz) + Ft)    # Somme des forces en Z
Mx      = Eq(-T1 + Ft * (d_pignon / 2))


DCL     = solve((Fy,Fz,FrFt,Mx) , (Ry, Rz, Fr,Ft))


for F in DCL:
    print(str(F) + " : " + str(DCL[F]))

Ry = round(DCL[Ry],5)
Rz = round(DCL[Rz],5)
Ft = round(DCL[Ft],5)
Fr = round(DCL[Fr],5)

My = (l * Ry) #moment au point critique
Mz = (l * Rz) #moment au point critique

M_resultant = math.sqrt(My**2 + Mz**2)


#%% Diamètre requis

Se_CD = symbols('Se_CD')
Se_AB = symbols('Se_AB')

d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)
d_ArbreAB = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_AB) + (math.sqrt(3)*Kfs*T2 / Sut)) )**(1/3)

##############################################################################
#Largeur des engrenages w
w = 0.75*d_pignon #

##############################################################################
#%% Courbe SN

d_buffer_CD = 255
  
if (Sut < 1400):
    Se_prime = 0.5*Sut
else:
    Se_prime = 700

ka = 4.51*pow(Sut,-0.265) #(usinage, laminage à froid)
kc = 1 #chargement combiné
kd = 1 #pas d'info donc 20°C
ke = 0.814 #fiabilité de 99%
 
d_ArbreCD = 254

while (math.fabs(d_buffer_CD - d_ArbreCD >= 0.1)):
    d_buffer_CD = d_ArbreCD
    if(d_ArbreCD <= 51):
        kb = 1.24 * d_ArbreCD**-0.107
    else:
        kb = 1.51 * d_ArbreCD**-0.157 #d_ArbreCD inconnu donc pire cas à d=254
    Se_CD = ka*kb*kc*kd*ke*Se_prime
    d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)


##############################################################################
#%% Courbe SN

d_buffer = 255
  
if (Sut < 1400):
    Se_prime = 0.5*Sut
else:
    Se_prime = 700  
 
d_ArbreAB = 254

while (math.fabs(d_buffer - d_ArbreAB >= 0.1)):
    d_buffer = d_ArbreAB
    if(d_ArbreAB <= 51):
        kb = 1.24 * d_ArbreAB**-0.107
    else:
        kb = 1.51 * d_ArbreAB**-0.157 #d_ArbreCD inconnu donc pire cas
    Se_AB = ka*kb*kc*kd*ke*Se_prime
    d_ArbreAB = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_AB) + (math.sqrt(3)*Kfs*T2 / Sut)) )**(1/3)
    
    
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
a_roul      = 3 #billes
Xi          = 1
H01         = 5000 #heures 

# Fiabilite
R = 0.99

#%% Charge dynamique CD

# Durée de vie pour une fiabilité 1%
L01_CD = 60 * H01 * n1
L10 = L01_CD / (x0 + theta*(math.log(1/R))**(1/b))

# Force radiale
Fe = Xi*Fr #Prendre les valeurs du design 

# Charge dynamique
C10_CD = (Fe * (L10/10**6)**(1/a_roul)) / 1000 #kN

print("C10_CD : " + str(C10_CD))

#%% Charge dynamique CD

# Durée de vie pour une fiabilité 99%
L01_AB = 60 * H01 * n2

L10 = L01_AB / (x0 + theta*(math.log(1/R))**(1/b))

# Force radiale
Fe = Xi*Fr #Prendre les valeurs du design 

# Charge dynamique
C10_AB = (Fe * (L10/10**6)**(1/a_roul)) / 1000 #kN

print("C10_AB : " + str(C10_AB))

##############################################################################
##############################################################################
##%%PARTIE 2
##############################################################################
##############################################################################
#%% Force poulie
F2 = (15 * P_watt)/(n1 * (D_poulie1*(10**-3)) * math.pi)
F1 = 5 * F2

F4 = (15 * P_watt)/(n1 * (D_poulie2*(10**-3)) * math.pi)
F3 = 5 * F4


##############################################################################
# #%% Torque Poulies

T_entree = P_watt / ((2*math.pi*n1)/60)
T_sortie = T_entree*4

##############################################################################
##%% Calcul équilibre des forces

#Arbre CD
Rcy, Rdy, Rcz, Rdz = symbols('Rcy Rdy Rcz Rdz')

FyCD = Eq(Rcy + Rdy - Fr + (F1+F2)*math.cos(theta_courroie))    # Somme des forces en Y
FzCD = Eq(Rcz + Rdz + Ft + (F1+F2)*math.sin(theta_courroie))    # Somme des forces en Z
Myd = ((l * Ft) + ((F1+F2)*math.cos(theta_courroie)*(l1 + 2*l)) - Rcz*l*2) #moment par rapport à D
Mzd = ((l * Fr) + ((F1+F2)*math.sin(theta_courroie)*(l1 + 2*l)) - Rcy*l*2) #moment par rapport à D

DCL = solve((FyCD,FzCD,Myd,Mzd), (Rcy,Rdy,Rcz,Rdz))

#Arbre AB
Ray, Rby, Raz, Rbz = symbols('Ray Rby Raz Rbz')

FyAB = Eq(Ray + Rby - Fr + (F3+F4)*math.cos(theta_courroie))    # Somme des forces en Y
FzAB = Eq(Raz + Rbz + Ft + (F3+F4)*math.sin(theta_courroie))    # Somme des forces en Z
Mya = (-(l * Ft) - ((F3+F4)*math.cos(theta_courroie)*(l2 + 2*l)) + Rbz*l*2) #moment par rapport à A
Mza = (-(l * Fr) - ((F3+F4)*math.sin(theta_courroie)*(l2 + 2*l)) + Rby*l*2) #moment par rapport à A

DCL1 = solve((FyAB,FzAB,Mya,Mza), (Ray,Rby,Raz,Rbz))


for F in DCL:
    print(str(F) + " : " + str(DCL[F]))

Rdz = round(DCL[Rdz],5)
Rdy = round(DCL[Rdy],5)
Rcz = round(DCL[Rcz],5)
Rcy = round(DCL[Rcy],5)

for F in DCL1:
    print(str(F) + " : " + str(DCL1[F]))

Rbz = round(DCL1[Rbz],5)
Rby = round(DCL1[Rby],5)
Raz = round(DCL1[Raz],5)
Ray = round(DCL1[Ray],5)


##############################################################################
#%%Calcul SN 

#Calcul Se arbre CD
if (Sut < 1400):
    Se_prime = 0.5*Sut
else:
    Se_prime = 700

ka = 4.51*pow(Sut,-0.265) #(usinage, laminage à froid)
kc = 1 #chargement combiné

if(d_ArbreCD <= 51):
    kb = 1.24 * d_ArbreCD**-0.107
else:
    kb = 1.51 * d_ArbreCD**-0.157 #d_ArbreCD inconnu donc pire cas

kd = 1 #pas d'info donc 20°C
ke = 0.814 #fiabilité de 99%
Se_CD = ka*kb*kc*kd*ke*Se_prime


#Calcul Se arbre AB
if (Sut < 1400):
    Se_prime = 0.5*Sut
else:
    Se_prime = 700

ka = 4.51*pow(Sut,-0.265) #(usinage, laminage à froid)
kc = 1 #chargement combiné

if(d_ArbreCD <= 51):
    kb = 1.24 * d_ArbreAB**-0.107
else:
    kb = 1.51 * d_ArbreAB**-0.157 #d_ArbreCD inconnu donc pire cas

kd = 1 #pas d'info donc 20°C
ke = 0.814 #fiabilité de 99%
Se_AB = ka*kb*kc*kd*ke*Se_prime

##############################################################################
#%%Calcul des moments résultants

#Moment au point critique pour arbre CD
MyCD = ((F1+F2)*math.cos(theta_courroie))
MzCD = ((F1+F2)*math.sin(theta_courroie))

#Moment au point critique pour arbre CD
MyAB = ((F3+F4)*math.cos(theta_courroie))
MzAB = ((F3+F4)*math.sin(theta_courroie))


#Moment résultant pour les deux arbres
M_resultant1 = math.sqrt(MyCD**2 + MzCD**2)
M_resultant2 = math.sqrt(MyAB**2 + MzAB**2)

##############################################################################
#%%Calcul des Facteurs de sécurité

#Arbre CD
FS1 = ((16*math.pi*d_ArbreCD)*(((2*Kf*M_resultant1)/Se_CD)*((math.sqrt(3)*Kfs*T_entree)/Sut)))**-1

#Arbre AB
FS2 = ((16*math.pi*d_ArbreAB)*(((2*Kf*M_resultant2)/Se_AB)*((math.sqrt(3)*Kfs*T_sortie)/Sut)))**-1

print("FS1 : " +str(FS1))
print("FS2 : " +str(FS2))

##############################################################################
#%%Print durée de vie roulement

print("L01_CD : " + str(L01_CD) + " cycles") #cycles
print("L01_AB : " + str(L01_AB) + " cycles") #cycles

   
   







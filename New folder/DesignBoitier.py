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
n1      = 1750  #tr/min - moteur
Kf      = 3     #CC en fatigue
Kfs     = 3     #CC en fatigue
P_hp    = 1
FS      = 3


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


print("l : " + str(l))
print("a : " + str(a))
print("d_pignon : " + str(d_pignon))
print("d_roue : " + str(d_roue))

del d_p, d_r, eq1, eq2
##############################################################################
#%% DCL - Détermination des forces et des moments

# Puisque Rcy = Rdy ET Rcz = Rdz ET Ray = Rby ET Raz = Rbz :
# RcdY, RabY, RcdZ, RabZ

Ry, Rz, Fr, Ft = symbols('Ry Rz Fr Ft')

FrFt = Eq(Fr/Ft-math.tan(math.radians(20)))
Fy = Eq(2*(Ry) - Fr)    # Somme des forces en Y
Fz = Eq(2*(Rz) - Ft)    # Somme des forces en Z
Mx  = Eq(-T1 + Ft * (d_pignon / 2))


DCL = solve((Fy,Fz,FrFt,Mx) , (Ry, Rz, Fr,Ft))


for F in DCL:
    print(str(F) + " : " + str(DCL[F]))

Ry = round(DCL[Ry],5)
Rz = round(DCL[Rz],5)
Ft = round(DCL[Ft],5)
Fr = round(DCL[Fr],5)

My = (l * Ry)
Mz = (l * Rz)

M_resultant = math.sqrt(My**2 + Mz**2)

print("My : " + str(My))
print("Mz : " + str(Mz))



#%% Diamètre requis

Se_CD = symbols('Se_CD')
Se_AB = symbols('Se_AB')

d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)
d_ArbreAB = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_AB) + (math.sqrt(3)*Kfs*T2 / Sut)) )**(1/3)

print("")
print("P hp : " + str(P_hp))
print("")

##############################################################################
#Largeur des engrenages w
w=0.75*d_pignon #

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
        kb = 1.51 * d_ArbreCD**-0.157 #d_ArbreCD inconnu donc pire cas
    Se_CD = ka*kb*kc*kd*ke*Se_prime
    d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)
    print("d_CD : " + str(d_ArbreCD))



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
    print("d_AB : " + str(d_ArbreAB))
    
    
    
    
##############################################################################
#%%Interface utilisateur

D_poulie1 = 10 #mm
D_poulie2 = 10 #mm
l1 = 5 #mm
l2 = 5 #mm
n_entree = 1000 #tr/min
H_entree = 745.7 #watt
theta_courroie = 45 #degrée

    
##############################################################################
#%% Force poulie
F2 = (15 * H_entree)/(n_entree * (D_poulie1) * math.pi)
F1 = 5 * F2


##############################################################################
# #%% Torque Poulie 1

T= H_entree / ((2*math.pi*n_entree)/60)


Rcy, Rdy, Rcz, Rdz, Fr, Ft = symbols('Rcy Rdy Rcz Rdz Fr Ft')

FrFt = Eq(Fr/Ft-math.tan(math.radians(20)))
Fy = Eq(Rcy+ Rdy - Fr + (F1+F2)*math.cos(math.radians(theta_courroie)))    # Somme des forces en Y
Fz = Eq(Rcz+ Rdz - Ft - (F1+F2)*math.sin(math.radians(theta_courroie)))    # Somme des forces en Z
Mx  = Eq(-T + Ft * (d_pignon / 2))

DCL = solve((FrFt,Ft,Fz,Mx), (Rcy,Rdy,Rcz,Rdz))

My = (-(l * Fr) - ((F1+F2)*math.cos(math.radians(theta_courroie))*(l1 + 2*l)) + Rcy*l*2)
Mz = (-(l * Ft) + ((F1+F2)*math.sin(math.radians(theta_courroie))*(l1 + 2*l)) + Rcz*l*2)

DCL = solve((Fy,Fz,FrFt,Mx,My,Mz) , (Rcy, Rdy, Rcz, Rdz, Fr,Ft))


for F in DCL:
    print(str(F) + " : " + str(DCL[F]))

Rdz = round(DCL[Rdz],5)
Rdy = round(DCL[Rdy],5)
Rcz = round(DCL[Rcz],5)
Rcy = round(DCL[Rcy],5)
Ft = round(DCL[Ft],5)
Fr = round(DCL[Fr],5)

My = eval(str(My))
Mz = eval(str(Mz))

M_resultant = math.sqrt(My**2 + Mz**2)

print("My : " + str(My))
print("Mz : " + str(Mz))

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
        kb = 1.51 * d_ArbreCD**-0.157 #d_ArbreCD inconnu donc pire cas
    Se_CD = ka*kb*kc*kd*ke*Se_prime
    d_ArbreCD = ( (16*FS/math.pi) * ((2*Kf*M_resultant / Se_CD) + (math.sqrt(3)*Kfs*T1 / Sut)) )**(1/3)
    print("d_CD : " + str(d_ArbreCD))






    

















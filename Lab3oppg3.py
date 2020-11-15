"TILNÆRMING AV KRITISK VINKEL FOR RULLENDE OBJEKT NED KVARTSIRKEL. RULLING OG SLURING"

import numpy as np
import matplotlib.pyplot as plt


#Startsverdier
G=9.81
R=0.54
r=0.02
c=1/2
mu_s=1000
mu_k=0.204
m=0.5

#Startsbetingelser
t_start=0
t_slutt=2
t_delta=0.0001
t_verdier=np.linspace(t_start,t_slutt,int(1/t_delta))

theta_init=0.2617993878 #her må vi ha theta>0 slik at løkken senere ikke repeterer tomme lister
tol=0.01 #Toleransen er avhengig av t_delta. 
        #Jo mindre t_delta, altså flere iterasjoner, jo mindre kan man sette tol, og man får et mer nøyaktig svar.

#Lager tomme lister for alle variablene
N=np.zeros(len(t_verdier))
f=np.zeros(len(t_verdier))

a_verdier=np.zeros(len(t_verdier))
v_verdier=np.zeros(len(t_verdier))

alpha_verdier=np.zeros(len(t_verdier))
omega_verdier=np.zeros(len(t_verdier))
theta_verdier=np.zeros(len(t_verdier))

alpha_objekt = np.zeros(len(t_verdier))
omega_objekt = np.zeros(len(t_verdier))

EK_trans=np.zeros(len(t_verdier))
EK_rot=np.zeros(len(t_verdier))
EP=np.zeros(len(t_verdier))

#Startsbetingelser for de tomme listene
theta_verdier[0]=theta_init
omega_verdier[0]=0

a_verdier[0]=G*np.sin(theta_init)/(1+c) #Akselerasjonsutrykket for ren rulling.
v_verdier[0]=0

N[0]=m*(G*np.cos(theta_verdier[0])-(v_verdier[0]**2/(R+r)))
f[0]=mu_s*N[0]

alpha_verdier[0] = a_verdier[0]/(r+r)
alpha_objekt[0] = a_verdier[0]/r

EP[0] = m*G*(R+r)*np.cos(theta_verdier[0])

#Løkken fyller de tomme listene for akselerasjon, vinkel, translatorisk fart og vinkelfart.
#alle disse listene er avhengige av startsbetingelsene og verdier fra forrige iterasjon

for i in range(len(t_verdier)-1):
    if mu_s*(N[i+1]) >= c*m*a_verdier[i+1]:

        v_verdier[i+1]=v_verdier[i]+a_verdier[i]*t_delta

        omega_verdier[i+1]=(v_verdier[i+1])/(R+r)

        theta_verdier[i+1]=theta_verdier[i]+omega_verdier[i]*t_delta

        a_verdier[i+1]=G*np.sin(theta_verdier[i+1])/(1+c)
        
        N[i+1]=m*(G*np.cos(theta_verdier[i+1])-(v_verdier[i+1]**2/(R+r)))
        
        f[i] = N[i]*mu_s
        
        alpha_objekt[i] = f[i]/(c*m*r)
        
        omega_objekt[i+1] = omega_objekt[i] + alpha_objekt[i]*t_delta
        
        EK_trans[i+1] = 1/2*m*v_verdier[i+1]**2
        
        EK_rot[i+1] = 1/2*c*m*v_verdier[i+1]**2
        
        EP[i+1] = m*G*(R+r)*np.cos(theta_verdier[i+1])
    
    if mu_s*(N[i+1])<=c*m*a_verdier[i+1]:
        
        a_verdier[i+1]=G*np.sin(theta_verdier[i+1])-(mu_k*N[i+1])/m
        
        v_verdier[i+1]=v_verdier[i]+a_verdier[i]*t_delta

        omega_verdier[i+1]=(v_verdier[i+1])/(R+r)

        theta_verdier[i+1]=theta_verdier[i]+omega_verdier[i]*t_delta
        
        N[i+1]=m*(G*np.cos(theta_verdier[i+1])-(v_verdier[i+1]**2/(R+r)))
        
        f[i] = N[i]*mu_k
        
        alpha_objekt[i] = f[i]/(c*m*r)
        
        omega_objekt[i+1] = omega_objekt[i] + alpha_objekt[i]*t_delta
        
        EK_trans[i+1] = 1/2*m*v_verdier[i+1]**2
        
        EK_rot[i+1] = 1/2*c*m*v_verdier[i+1]**2
        
        EP[i+1] = m*G*(R+r)*np.cos(theta_verdier[i+1])
        
        if N[i+1]<tol:
            #Løkken stopper dersom normalkraften blir mindre enn toleransen
            print("Kritisk verdi: ",round(theta_verdier[i+1]*180/np.pi,4),"grader")
            global b
            b = i+2
            break

#Utregning for den totale energien
#EK_trans = 0.5*m*(v_verdier**2)
#EK_rot = 0.5*c*m*(r**2)*(omega_objekt**2)
#EP = m*G*(R+r)*np.cos(theta_verdier)

#E_tot = EK_trans+EK_rot+EP

#print(E_tot)

plt.figure(0)
plt.plot(t_verdier,theta_verdier*180/np.pi,label="Theta verdier")
plt.xlabel("Tid")
plt.ylabel("Vinkel (grader)")
plt.legend()
plt.show()

plt.figure(0)
plt.plot(t_verdier, EK_trans + EK_rot + EP, color='orange',label="Total energi")
plt.xlabel("Tid")
plt.ylabel("Energi (J)")
#plt.axis([0,2,0,5])
plt.legend()
plt.show()
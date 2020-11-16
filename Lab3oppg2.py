"TILNÆRMING AV KRITISK VINKEL FOR RULLENDE OBJEKT NED KVARTSIRKEL. REN RULLING"

import numpy as np
import matplotlib.pyplot as plt

#Startsverdier
G=9.81
R=0.54
r=0.0305
c=1 #Treghetskonstanten til objektet
m=1 #Massen er kun nødvendig for utregning av energi. Kan være hva som helst.

#Startsbetingelser
t_start=0
t_slutt=2
t_delta=0.0001
t_verdier=np.linspace(t_start,t_slutt,int(t_slutt/t_delta))

theta_init=0.115 #Her må vi ha en vinkel>0 slik at løkken ikke repeterer tomme lister
tol=0.01 #Toleransen er avhengig av t_delta. 
        #Jo mindre t_delta, altså flere iterasjoner, jo mindre kan man sette tol, og man får et mer nøyaktig svar.

#Lager tomme lister for alle verdiene
a_verdier=np.zeros(len(t_verdier))
v_verdier=np.zeros(len(t_verdier))

omega_verdier=np.zeros(len(t_verdier))
alpha_verdier=np.zeros(len(t_verdier))
theta_verdier=np.zeros(len(t_verdier))

EK_trans=np.zeros(len(t_verdier))
EK_rot=np.zeros(len(t_verdier))
EP=np.zeros(len(t_verdier))

#Startsbetingelser for de tomme listene
theta_verdier[0]=theta_init
a_verdier[0]=G*np.sin(theta_init)/(1+c) #Akselerasjon med hensyn på friskjon og rotasjonbevegelse.
#Dette er også oppdatert for løkken nedenfor.
v_verdier[0]=0
omega_verdier[0]=0

#Løkken fyller de tomme listene for akselerasjon, vinkel, translatorisk fart og vinkelfart.
#Alle disse listene er avhengige av startsbetingelsene og verdier fra forrige iterasjon
for i in range(len(t_verdier)-1):

        v_verdier[i+1]=v_verdier[i]+a_verdier[i]*t_delta

        omega_verdier[i+1]=(v_verdier[i+1])/(R+r)

        theta_verdier[i+1]=theta_verdier[i]+omega_verdier[i]*t_delta

        a_verdier[i+1]=G*np.sin(theta_verdier[i+1])/(1+c)

        if abs(G*np.cos(theta_verdier[i])-(v_verdier[i]**2/(R+r)))<tol:
            #Løkken stopper dersom normalkraften blir mindre enn toleransen
            print("Kritisk verdi: ",round(theta_verdier[i]*180/np.pi,4),"grader")
            z=i+1
            break
   
#Utregning for den totale energien
EK_trans=0.5*m*(v_verdier**2)
EK_rot=0.5*c*m*(v_verdier**2)
EP=m*G*(R+r)*np.cos(theta_verdier)

#Grafer vinkelen over tidsintervallet
plt.figure(0)
plt.plot(t_verdier[0:z],theta_verdier[0:z]*180/np.pi,color='black',label="Theta verdier")
plt.xlabel("Tid")
plt.ylabel("Vinkel (grader)")
plt.grid(True)
plt.legend()
plt.show()

#Grafer den totale energien
plt.figure(0)
plt.plot(t_verdier[0:z+1],EK_trans[0:z+1]+EP[0:z+1]+EK_rot[0:z+1],color='orange',label="Total energi")
plt.plot(t_verdier[0:z+1],EK_trans[0:z+1],color='purple',label="Kinetisk translatorisk energi")
plt.plot(t_verdier[0:z+1],EK_rot[0:z+1],color='blue',label="Kinetisk rotasjonsenergi")
plt.plot(t_verdier[0:z+1],EP[0:z+1],color='green',label="Potensiell energi")
plt.grid(True)
plt.xlabel("Tid (s)")
plt.ylabel("Energi (J)")
plt.axis([0,t_verdier[z+1],0,10]) #Her kan man bestemme aksene for energi-grafen
plt.legend()
plt.show()
plt.show()

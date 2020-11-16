"TILNÆRMING AV KRITISK VINKEL FOR RULLENDE OBJEKT NED KVARTSIRKEL. RULLING OG SLURING"

import numpy as np
import matplotlib.pyplot as plt

#Startsverdier
G=9.81
R=0.54
r=0.03
c=1 #Treghetskonstanten til objektet.
mu_s=0.6
mu_k=0.3
m=1 #Massen er kun nødvendig for utregning av energi. Kan være hva som helst.

#Startsbetingelser
t_start=0
t_slutt=2
t_delta=0.0001
t_verdier=np.linspace(t_start,t_slutt,int(t_slutt/t_delta))

theta_init=0.2 #Her må vi ha en vinkel>0 slik at løkken ikke repeterer tomme lister
tol=0.01 #Toleransen er avhengig av t_delta. 
        #Jo mindre t_delta, altså flere iterasjoner, jo mindre kan man sette tol, og man får et mer nøyaktig svar.

#Lager tomme lister for alle variablene
#Her er det lagt til nye lister for friksjon og normalkraft. Dette er oppdatert fra koden i oppg. 2
N=np.zeros(len(t_verdier))
f=np.zeros(len(t_verdier))

a_verdier=np.zeros(len(t_verdier))
v_verdier=np.zeros(len(t_verdier))

alpha_verdier=np.zeros(len(t_verdier))
omega_verdier=np.zeros(len(t_verdier))
theta_verdier=np.zeros(len(t_verdier))

#Nye lister for vinkelfart og vinkelakselerasjon til objektet. Nødvendig for energiutregning.
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

alpha_verdier[0] = a_verdier[0]/(R+r)
alpha_objekt[0] = a_verdier[0]/r


#Løkken fyller de tomme listene for akselerasjon, vinkel, translatorisk fart og vinkelfart.
#alle disse listene er avhengige av startsbetingelsene og verdier fra forrige iterasjon

for i in range(len(t_verdier)-1):
    if mu_s*(N[i+1]) >= c*m*a_verdier[i+1]: #Denne løkken kjører så lenge objektet ruller rent.

        v_verdier[i+1]=v_verdier[i]+a_verdier[i]*t_delta

        omega_verdier[i+1]=(v_verdier[i+1])/(R+r)

        theta_verdier[i+1]=theta_verdier[i]+omega_verdier[i]*t_delta

        a_verdier[i+1]=G*np.sin(theta_verdier[i+1])/(1+c)
        
        N[i+1]=m*(G*np.cos(theta_verdier[i+1])-(v_verdier[i+1]**2/(R+r)))
        
        f[i] = m*(G*np.sin(theta_verdier[i])-a_verdier[i])
        
        alpha_objekt[i] = f[i]/(c*m*r)
        
        omega_objekt[i+1] = omega_objekt[i] + alpha_objekt[i]*t_delta
        
        EK_trans[i+1] = (1/2)*m*v_verdier[i+1]**2
        
        EK_rot[i+1] = (1/2)*c*m*r**2*omega_objekt[i+1]**2
        
        EP[i+1] = m*G*(R+r)*np.cos(theta_verdier[i+1])
    
    if mu_s*(N[i+1])<=c*m*a_verdier[i+1]: #Dette kravet er oppfylt når objektet begynner å slure.
            
        a_verdier[i+1]=G*np.sin(theta_verdier[i+1])-(mu_k*N[i+1])/m #Nytt utrykk for akselerasjonen.
        
        v_verdier[i+1]=v_verdier[i]+a_verdier[i]*t_delta

        omega_verdier[i+1]=(v_verdier[i+1])/(R+r)

        theta_verdier[i+1]=theta_verdier[i]+omega_verdier[i]*t_delta
        
        N[i+1]=m*(G*np.cos(theta_verdier[i+1])-(v_verdier[i+1]**2/(R+r)))
        
        f[i] = N[i]*mu_k #Friksjonen kan utrykkes slik ved sluring.
        
        alpha_objekt[i] = f[i]/(c*m*r)
        
        omega_objekt[i+1] = omega_objekt[i] + alpha_objekt[i]*t_delta
        
        EK_trans[i+1] = (1/2)*m*v_verdier[i+1]**2
        
        EK_rot[i+1] = (1/2)*c*m*r**2*omega_objekt[i+1]**2
        
        EP[i+1] = m*G*(R+r)*np.cos(theta_verdier[i+1])
        
        if N[i+1]<tol:
            #Løkken stopper dersom normalkraften blir mindre enn toleransen
            print("Kritisk verdi: ",round(theta_verdier[i+1]*180/np.pi,4),"grader")
            z=i+1
            break

#Grafer vinkelen over tidsintervallet
plt.figure(0)
plt.plot(t_verdier[0:z+1],theta_verdier[0:z+1]*180/np.pi,color='black',label="Theta verdier")
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

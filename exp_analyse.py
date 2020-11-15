import numpy as np
import matplotlib.pyplot as plt

filnavn = 'HUL_SYL2_f_data.csv' #Her kan man teste de ulike filene. Husk at filene og python
#programmet må være lagret i samme mappe for at koden skal kjøre de rette filene.

fildata = np.loadtxt(filnavn, delimiter=';')

tol=0.496 #Her må man finne riktig verdi for de ulike filene for å avgjøre når objektet faller av kvartsirkel.
#Dette kan man lese av grafen etter du har kjørt programmet. Tips: sett tol=1 for å sjekke for alle verdiene,
#så velger man en toleranse verdi som man mener passer grafen og sjekker hvilken kritisk vinkel programmet printer.

#Disse listene blir fyllt med dataen fra filen man velger
tid = fildata[:,0]
x_verdier = fildata[:,1]
y_verdier = fildata[:,2]
r=np.zeros(len(tid))
dr=np.zeros(len(tid))
theta=np.arctan(x_verdier/y_verdier)*180/np.pi

#Denne regner distansen objektet har fra origo i kvartsirkelen
for z in range(len(tid)):
    r[z]=np.sqrt(x_verdier[z]**2+y_verdier[z]**2)
    if tid[z]>tol:
        print(f'Faktiske: {theta[z]}')
        print(f'Usikkerhet: {(theta[z+1]-theta[z]+theta[z]-theta[z-1])/2}')
        print(f'Startsvinkel: {np.arctan(x_verdier[0]/y_verdier[0])} radianer')
        print(f'Startsvinkel: {np.arctan(x_verdier[0]/y_verdier[0])*180/np.pi} grader')
        break


    
plt.figure(0)
plt.title('Endring av radius over tid', fontsize=12)
plt.plot(tid,r,'.',label="radius")
#plt.axvline(x=1.078, color='r', linestyle='--')
plt.xlabel("Tid (s)")
plt.ylabel("Radius (r)")
plt.legend()
plt.show()




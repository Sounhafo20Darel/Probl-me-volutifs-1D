
import matplotlib.pyplot as plt
import numpy as np

E=200 #prix strike
a, b=-4, 4 # intervalle que parcour x=[-4,4]
sigma=30/100 # volatilite
r=3/100 #taux dinteret fixe
k=(2*r)/sigma**2
alpha=-(k-1)/2
beta=-((k+1)**2)/4
T=3
Z=(sigma**2)*T/2
N, M=100, 100 #les nombres d'intervalles

#discretisation

x=np.linspace(a,b,N+1)
t=np.linspace(0,Z,M+1)
c, d=(b-a)/N, Z/M # c le pas en espace, d le pas en temps
F= d/c**2 # represente le lambda 

#conditions initiale
def I(x):
    if x>0:
        return (np.exp(-alpha*x))*(np.exp(x)-1)
    else:
        return 0
u0=np.zeros(N+1)
for i in range (N+1):
    u0[i]=I(x[i])
print(u0)

#conditions aux limites

def g1(t):
    return 0
def g2(t):
    return np.exp(-alpha*b-beta*t)*(np.exp(b)-np.exp(-k*t))

#matrice solution

Usol=np.zeros((M+1,N+1))
Usol[0,:]=u0

H=np.zeros(N-1)

#matrices A

Aa=-F*np.eye(N-1,k=-1)
Ab=-F*np.eye(N-1,k=1)
Ac=(1+2*F)*np.identity(N-1)
A= Aa+Ab+Ac

for n in range(1, M+1):
    Usol[n,0]=g1(t[n])
    Usol[n,N]=g2(t[n])
    H=Usol[n-1,1:N].copy()
    H[0]+=F*Usol[n,0]
    H[-1]+=F*Usol[n,N]
    Usol[n, 1:N]=np.linalg.solve(A,H)
    
#retour aux variables financieres

S=E*np.exp(x) # prix sous-jacent issus de S a partir de x
U_final=Usol[M,:] # on concidere le prix a la maturite
V=E*np.exp(alpha*x+beta*Z)*U_final 
payoff=np.maximum(S-E,0)

#figure

#premiere figure
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
for i in range (0, M+1, M//5):
    plt.plot(x, Usol[i,:], label=f"t={t[i]:.3f}")
plt.grid(True, alpha=0.5)
plt.xlabel("espace adimentionne (x)")
plt.ylabel("espace adimensionne (u)")
plt.title("Evolution de u(x,t) dans le maillage")
plt.legend()

#figure2 graphique financier
plt.subplot(1,2,2)

plt.plot(S, V, label="prix call", color="blue", lw=2)
plt.plot(S, payoff, "--", label="payoff a l echeance", color="red")
plt.grid(True, alpha=0.5)
plt.xlabel("prix de laction (S)")
plt.ylabel("prix de loption (V)")
plt.title("Restitution complete de Black-Scholes")
plt.legend()
plt.xlim(0,500)
plt.ylim(0,500)


plt.tight_layout()
plt.show

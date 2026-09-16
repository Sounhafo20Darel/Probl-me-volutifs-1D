print('*****APPROXIMATION DE LA SOLUTION DE L EQUATION DE LA CHALEUR EN 1D*****')
print('Methode des differences finies schema implicite')
T=float(input('entrez la borne superieure l axe temps'))
M=int(input('entrez le nombre d instant sur lesquels vous souhaitez visualiser l evolution du phenomene'))
a=float(input('entrez la borne inferieure sur l axe spacial'))
b=float(input('entrez la borne superieure sur l axe spacial'))
N=int(input('entrez le nombre de subdivision de notre intervalle sur l axe spacial'))
µ=float(input('entrez la constante de l equation de la chaleur'))
from numpy import*
import matplotlib.pyplot as plt
Dx=(b-a)/(N+1)
Dt=T/M
t=zeros(M)
for i in range(M):
    t[i]=i*Dt
    
x=zeros(N+1)
for i in range(N+1):
    x[i]=a+i*Dx
    
U=zeros((M,N+1))
for i in range(1,N):
    if -1<=x[i]<=1:
        U[0,i]=1
    else:
        U[0,i]=0
    
for i in range(M):
    U[i,0]=0
    U[i,N]=0
    

Vs=(-µ/(Dx)**2)*ones(N-1)
Vi=(-µ/(Dx)**2)*ones(N-1)
Vm=((1/Dt)+2*µ/(Dx)**2)*ones(N)
A=diag(Vm)+diag(Vi,k=-1)+diag(Vs,k=1)
B=(1/Dt)*eye(N)

F=zeros((M,N))
for n in range(M):
    for i in range(N):
        F[n,i]=0

    
V=zeros(N)
for i in range(N):
    V[i]=U[0,i+1]
    
for n in range(1,M):
    V=linalg.solve(A,B.dot(V)+F[n,:])
    for i in range(N):
        U[n,i+1]=V[i]

for n in range(M):
    plt.plot(x,U[n,:],label=f'instant{n}')

plt.legend()
plt.show()

        
        

    

      
    



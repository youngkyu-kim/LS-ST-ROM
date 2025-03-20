import numpy as np; import scipy as sp
from scipy.sparse.linalg import spsolve,splu

# train parameter space
mu1,mu2=np.meshgrid(np.linspace(-0.9,-0.5,2),np.linspace(-0.9,-0.5,2),indexing='ij')
mu1=mu1.flatten(); mu2=mu2.flatten(); no_para=np.size(mu1,axis=0)

# space and time domain
N=70; Ns=(N-1)**2; Nt=50; h=1/N; k=2/Nt
x1,x2=np.meshgrid(np.linspace(0,1,N+1)[1:-1],np.linspace(0,1,N+1)[1:-1],indexing='ij')
x1=x1.flatten(); x2=x2.flatten(); t=np.linspace(0,2,Nt+1)

# basic operators
e=np.ones(N-1); I=sp.sparse.eye(Ns,format='csc')
A1D=1/h**2*sp.sparse.spdiags(np.vstack((e,-2*e,e)),[-1,0,1],N-1,N-1,format='csc')
A2D=sp.sparse.kron(A1D,sp.sparse.eye(N-1,format='csc'),format='csc')\
    +sp.sparse.kron(sp.sparse.eye(N-1,format='csc'),A1D,format='csc')

# snapshot, free DoFs, No IC included
U=np.zeros((Ns,no_para*Nt))
for p in range(no_para):
    A=A2D-sp.sparse.diags(1/np.sqrt((x1-mu1[p])**2+(x2-mu2[p])**2),format='csc'); invA=splu(I-k*A)
    f=(1/np.sqrt((x1-mu1[p])**2+(x2-mu2[p])**2).reshape(-1,1))@(np.sin(2*np.pi*t).reshape(1,-1))
    u=np.zeros((Ns,Nt+1)) # IC is zero
    for n in range(Nt):
        u[:,n+1]=invA.solve(u[:,n]+k*f[:,n+1])      
    U[:,np.arange(p*Nt,(p+1)*Nt)]=u[:,1:]
    
W,S,VT=np.linalg.svd(U) # POD of solution snapshot

# test parameter space
nparam1=15; nparam2=15; ParamD1,ParamD2=np.meshgrid(np.linspace(-1.7,-0.2,nparam1),np.linspace(-1.7,-0.2,nparam2))
ParamD1=ParamD1.flatten(); ParamD2=ParamD2.flatten(); num_test=np.prod(ParamD1.shape)

# variables to store (also store ParamD1 and ParamD2)
UROM_g2=np.zeros((num_test,Nt+1,N+1,N+1)); UFOM_2=np.zeros((num_test,Nt+1,N+1,N+1))
AvgRelErr_g2=np.zeros((num_test)); STResROM_g2=np.zeros((num_test))

# Construct Phi_s
ns=5; PHIs=W[:,:ns]; PHIsT=PHIs.T

# construct Djk
nt=3; D=np.zeros((ns,nt*Nt))
for i in range(ns):
    Ri=VT[i,:]; Ri=Ri.reshape(-1,no_para,order='F')
    Wi,Si,ViT=np.linalg.svd(Ri); PHIti=Wi[:,:nt]
    for j in range(nt):
        D[i,Nt*j:Nt*(j+1)]=PHIti[:,j]

# construct PHIst for post processing, i.e., reconstruction
PHIst=np.zeros((Ns*Nt,ns*nt))
for i in range(Nt):
    for j in range(nt):
        PHIstij=np.zeros((Ns,ns))
        Dij=sp.sparse.diags(D[:,j*Nt+i],format='csc')
        PHIst[Ns*i:Ns*(i+1),ns*j:ns*(j+1)]=PHIs@Dij  
PHIstT=PHIst.T

for test_count in range(num_test): # predictive cases
    param1,param2=ParamD1[test_count],ParamD2[test_count] # set parameter

    # set model
    A=A2D-sp.sparse.diags(1/np.sqrt((x1-param1)**2+(x2-param2)**2),format='csc')
    f=(1/np.sqrt((x1-param1)**2+(x2-param2)**2).reshape(-1,1))@(np.sin(2*np.pi*t).reshape(1,-1))

    # construct Ast and fst using block structure (ust0 is zero)
    As=PHIsT@A@PHIs; fs=PHIsT@f
    Ast=np.zeros((ns*nt,ns*nt))
    for i in range(nt):
        for j in range(nt):
            Astij=np.zeros((ns,ns))
            for kk in range(Nt):
                Dik=sp.sparse.diags(D[:,i*Nt+kk],format='csc')
                Djk=sp.sparse.diags(D[:,j*Nt+kk],format='csc')
                Astij+=Dik@Djk-k*Dik@As@Djk
                if kk != Nt-1:
                    Dik_next=sp.sparse.diags(D[:,i*Nt+kk+1],format='csc')
                    Astij-=Dik_next@Djk   
            Ast[ns*i:ns*(i+1),ns*j:ns*(j+1)]=Astij  
    fst=np.zeros(ns*nt)
    for j in range(nt):
        for kk in range(Nt):
            Djk=sp.sparse.diags(D[:,j*Nt+kk],format='csc')
            fst[ns*j:ns*(j+1)]+=k*Djk@fs[:,kk+1]

    # Space-Time ROM (online phase)
    UstROM=PHIst@np.linalg.solve(Ast,fst)

    # FOM
    u=np.zeros((Ns,Nt+1)); invA=splu(I-k*A)
    for n in range(Nt):
        u[:,n+1]=invA.solve(u[:,n]+k*f[:,n+1])
    UstFOM=u[:,1:].flatten(order='F')
    
    AvgRelErr_g2[test_count]=np.linalg.norm(UstROM-UstFOM)/np.linalg.norm(UstFOM) # Avg. rel. error
    
    rstROM=np.zeros(Ns*Nt) # ST residual
    for n in range(1,Nt+1):
        if n==1:
            rstROM[(n-1)*Ns:n*Ns]=k*f[:,n]-(I-k*A).dot(UstROM[(n-1)*Ns:n*Ns])
        else:
            rstROM[(n-1)*Ns:n*Ns]=k*f[:,n]+UstROM[(n-2)*Ns:(n-1)*Ns]-(I-k*A).dot(UstROM[(n-1)*Ns:n*Ns])
    STResROM_g2[test_count]=np.linalg.norm(rstROM)
    
    # store solutions for each param
    urom=np.zeros((Nt+1,N+1,N+1)); ufom=np.zeros((Nt+1,N+1,N+1))
    urom[0]=np.zeros((N+1,N+1)); ufom[0]=np.zeros((N+1,N+1))
    for n in range(1,Nt+1):
        urom[n,1:-1,1:-1]=UstROM[(n-1)*Ns:n*Ns].reshape((N-1,N-1))
        ufom[n,1:-1,1:-1]=UstFOM[(n-1)*Ns:n*Ns].reshape((N-1,N-1))
    UROM_g2[test_count]=urom; UFOM_2[test_count]=ufom
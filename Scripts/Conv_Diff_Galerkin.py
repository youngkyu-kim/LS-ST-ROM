import numpy as np; import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, splu

# train parameter space
mu1,mu2=np.meshgrid(np.linspace(0.03,0.05,2),np.linspace(0.33,0.35,2),indexing='ij')
mu1=mu1.flatten(); mu2=mu2.flatten(); no_para=np.size(mu1,axis=0) 

# space and time domain
N=70; Ns=(N-1)**2; Nt=50; h=1./N; Tfinal = 1.; k=Tfinal/Nt
x1,x2=np.meshgrid(np.linspace(0,1,N+1)[1:-1],np.linspace(0,1,N+1)[1:-1],indexing='ij')

# IC
u0=100*np.sin(2*np.pi*x1)**3*np.sin(2*np.pi*x2)**3; u0[np.nonzero(x1>0.5)]=0.0; u0[np.nonzero(x2>0.5)]=0.0; u0=u0.flatten()

# basic operators
e=np.ones(N-1); I=sp.eye(Ns,format='csc')
A1D_diff=1/h**2*sp.spdiags(np.vstack((e,-2*e,e)),[-1,0,1],N-1,N-1,format='csc') 
A2D_diff=sp.kron(A1D_diff,sp.eye(N-1,format='csc'),format='csc')\
    +sp.kron(sp.eye(N-1,format='csc'),A1D_diff,format='csc')
A1D_conv=1/(h)*sp.spdiags(np.vstack((-1*e,1*e,0*e)),[-1,0,1],N-1,N-1,format='csc') 
A2D_conv=sp.kron(A1D_conv,sp.eye(N-1,format='csc'),format='csc')\
    +sp.kron(sp.eye(N-1,format='csc'),A1D_conv,format='csc')

# snapshot, free DoFs, No IC included
U=np.zeros((Ns,no_para*Nt))
for p in range(no_para):
    A=-mu1[p]*A2D_conv + mu2[p]*A2D_diff; invA=splu(I-k*A) 
    u=np.zeros((Ns,Nt+1)); u[:,0]=u0
    for n in range(Nt):
        u[:,n+1]=invA.solve(u[:,n])     
    U[:,np.arange(p*Nt,(p+1)*Nt)]=u[:,1:]

# POD of solution snapshot
W,S,VT=np.linalg.svd(U)

# test parameter space
nparam1=12; nparam2=12; ParamD1,ParamD2=np.meshgrid(np.linspace(0.01,0.07,nparam1),np.linspace(0.31,0.37,nparam2))
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
        Dij=sp.diags(D[:,j*Nt+i],format='csc')
        PHIst[Ns*i:Ns*(i+1),ns*j:ns*(j+1)]=PHIs@Dij  
PHIstT=PHIst.T 

# precompute (parameter seperation)
A1s=PHIsT@A2D_conv@PHIs; A2s=PHIsT@A2D_diff@PHIs; 
DD=np.zeros((ns*nt,ns*nt)); DA1sD=np.zeros((ns*nt,ns*nt)); DA2sD=np.zeros((ns*nt,ns*nt))
for i in range(nt):
    for j in range(nt):
        DDij=np.zeros((ns,ns)); DA1sDij=np.zeros((ns,ns)); DA2sDij=np.zeros((ns,ns))
        for kk in range(Nt):
            Dik=sp.diags(D[:,i*Nt+kk],format='csc'); Djk=sp.diags(D[:,j*Nt+kk],format='csc')
            DDij+=Dik@Djk; DA1sDij-=k*Dik@A1s@Djk; DA2sDij-=k*Dik@A2s@Djk 
            if kk != Nt-1:
                Dik_next=sp.diags(D[:,i*Nt+kk+1],format='csc')
                DDij-=Dik_next@Djk   
        DD[ns*i:ns*(i+1),ns*j:ns*(j+1)]=DDij; DA1sD[ns*i:ns*(i+1),ns*j:ns*(j+1)]=DA1sDij
        DA2sD[ns*i:ns*(i+1),ns*j:ns*(j+1)]=DA2sDij
        
# construct Ust
ust=np.zeros(ns*nt); us0=PHIsT@u0
for j in range(nt):
    Djk=sp.diags(D[:,j*Nt],format='csc')
    ust[ns*j:ns*(j+1)] = Djk@us0   

# predictive cases
for test_count in range(num_test):
    param1,param2=ParamD1[test_count],ParamD2[test_count] # set parameter mu

    # construct Ast (fst is zero)
    Ast=DD-param1*DA1sD+param2*DA2sD  

    A = -param1*A2D_conv + param2*A2D_diff # set model

    # Space-Time ROM (online phase)
    UstROM=PHIst@np.linalg.solve(Ast,ust)

    # FOM
    u=np.zeros((Ns,Nt+1)); u[:,0] = u0; invA=splu(I-k*A)
    for n in range(Nt):
        u[:,n+1]=invA.solve(u[:,n])
    UstFOM=u[:,1:].flatten(order='F')

    AvgRelErr_g2[test_count]=np.linalg.norm(UstROM-UstFOM)/np.linalg.norm(UstFOM) # Avg. rel. error

    rstROM=np.zeros(Ns*Nt) # ST residual
    for n in range(1,Nt+1):
        if n==1:
            rstROM[(n-1)*Ns:n*Ns]=u0-(I-k*A).dot(UstROM[(n-1)*Ns:n*Ns])
        else:
            rstROM[(n-1)*Ns:n*Ns]=UstROM[(n-2)*Ns:(n-1)*Ns]-(I-k*A).dot(UstROM[(n-1)*Ns:n*Ns])
    STResROM_g2[test_count]=np.linalg.norm(rstROM)

    # store solutions for each param
    urom=np.zeros((Nt+1,N+1,N+1)); ufom=np.zeros((Nt+1,N+1,N+1))
    urom[0,1:-1,1:-1]=u0.reshape((N-1,N-1)); ufom[0,1:-1,1:-1]=u0.reshape((N-1,N-1))
    for n in range(1,Nt+1):
        urom[n,1:-1,1:-1]=UstROM[(n-1)*Ns:n*Ns].reshape((N-1,N-1))
        ufom[n,1:-1,1:-1]=UstFOM[(n-1)*Ns:n*Ns].reshape((N-1,N-1))
    UROM_g2[test_count]=urom; UFOM_2[test_count]=ufom
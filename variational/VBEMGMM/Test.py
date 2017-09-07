import numpy as np
import pandas as pandas

K = 6

of_df = pandas.read_csv("/Users/dom/Dropbox/Private/NumMethHaskell/variational/VBEMGMM/faithful.txt", sep = ' ', header=None)

N = of_df.shape[0]

np.random.seed(42)
Z = np.array([np.random.dirichlet(np.ones(K)) for _ in range(N)])

Y = pandas.DataFrame(Z)

Y.to_csv("/Users/dom/Dropbox/Private/NumMethHaskell/variational/VBEMGMM/faithZ.txt")

YHask = pandas.read_csv("/Users/dom/Dropbox/Private/NumMethHaskell/variational/bigR.csv", sep = ',', header=None)

NK = Z.sum(axis=0)

ZHask = pandas.DataFrame.as_matrix(YHask)

NKHask = ZHask.sum(axis=0)

def calcXd(Z,X):
    #weighted means (by component responsibilites)
    (N,XDim) = np.shape(X)
    (N1,K) = np.shape(Z)
    NK = Z.sum(axis=0)
    assert N==N1
    xd = np.zeros((K,XDim))
    for n in range(N):
        for k in range(K):
            xd[k,:] += Z[n,k]*X[n,:]
    #safe divide:
    for k in range(K):
        if NK[k]>0: xd[k,:] = xd[k,:]/NK[k]

    return xd

X = of_df.as_matrix()

xd = calcXd(ZHask,X)

def calcS(Z,X,xd,NK):
    (N,K)=np.shape(Z)
    (N1,XDim)=np.shape(X)
    assert N==N1

    S = [np.zeros((XDim,XDim)) for _ in range(K)]
    for n in range(N):
        for k in range(K):
            B0 = np.reshape(X[n,:]-xd[k,:], (XDim,1))
            L = np.dot(B0,B0.T)
            assert np.shape(L)==np.shape(S[k]),np.shape(L)
            S[k] += Z[n,k]*L
    #safe divide:
    for k in range(K):
        if NK[k]>0: S[k] = S[k]/NK[k]
    return S

S = calcS(Z,X,xd,NK)

# -*- coding: utf-8 -*-
"""
Editor Spyder

Este é um arquivo de script temporário.
"""
import numpy as np
#ITEM A: rotina que resolve o sistema pelo método de Crank Nicolson para os pontos ps dados, retornando u1 que armazena as soluções uk(T,xi) para T=1 
def Decomposicao(N, M, dt, dx, item, nf, P, erro):
    u1 = []
    Mat = []
    for k in range (N-1):
        Mat.append([0]*(N-1))
    Mat[0][0] = 1 + erro
    for i in range (1,N-1):
        Mat[i-1][i] = -erro/2
        Mat[i][i] = 1 + erro
        Mat[i][i-1] = -erro/2
    L = []
    for k in range (N-1):
        L.append([0]*(N-1))
    D = []
    for k in range (N-1):
        D.append([0]*(N-1))
        D[0][0] = 1 + erro
    for j in range (N-1):
        for i in range (j):
            h = Mat[j][i]
            L[j][i] = h/D[i][i]
            for k in range(i+1, j+1):
                Mat[j][k] = Mat[j][k] - h*L[k][i]
        D[j][j] = Mat[j][j]
    for l in range (nf):
        a = P[l]
        
        U = []
        for k in range (M+1):
            U.append([0]*(N+1))
        F = []
        for k in range (M+1):
            F.append([0]*(N+1))
        for k in range (M+1):
            for i in range (N):
                if a - dx/2 <= i*dx and i*dx <= a + dx/2: 
                    F[k][i] = (10 * (1 + np.cos(5*k*dt)))/dx
                else:
                    F[k][i]= 0
        for k in range (M):    
            B = []
            for q in range (N-1):
                B.append(0)
            
            for i in range(N-1):
                if (i==0):
                    B[i] = (erro/2)*U[k+1][i] + (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
                elif (i==N-2):
                    B[i] = (erro/2)*U[k+1][N] + (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
                else:
                    B[i] = (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
        
            Z = []
            for v in range (N-1):
                Z.append(0)
            C = []
            for a in range (N-1):
                C.append(0)
            for j in range (N-1):
                Z[j] = B[j]
                for i in range (j):
                    Z[j] = Z[j] - L[j][i]*Z[i] 
                C[j] = Z[j]/D[j][j]
            X = []
            for o in range (N-1):
                X.append(0)  
            for t in reversed(range(N-1)):
                X[t] = C[t]
                for i in range (t+1, N-1):
                    X[t] = X[t] - L[i][t]*X[i]
            for i in range (1,N):
                U[k+1][i] = X[i-1]
        u1.append(U[N]) 
    return u1
#ITEM B: rotina que monta as matrizes A e B do sistema normal do problema de mínimos quadrados
def matrizAeB(N, M, dt, dx, item, nf, P, graficasso):
    ut= []
    if item == 'a':
        for i in range (N):
            ut.append(7*graficasso[0][i])
    if item == 'b':
        for t in range (N):
            ut.append(2.3*graficasso[0][t] + 3.7*graficasso[1][t] + 0.3*graficasso[2][t] + 4.2*graficasso[3][t])
    A = []
    for i in range (nf):
        t = []
        for j in range (nf):
            p = 0
            for k in range (N-1):
                p = graficasso[i][k]*graficasso[j][k] + p
            t.append(p)
        A.append(t)
    B = []
    for i in range (nf):
        p = 0
        for k in range (N-1):
            p = ut[k]*graficasso[i][k] + p
        B.append(p)
    return (A, B)
#ITEM C: rotina que realiza a decomposição LDLt da matriz simétrica A e resolve o sistema 
def rotinaLDLt(N, M, dt, dx, item, nf, P, graficasso, erro, A, B):
    L = []
    for i in range(nf):
        Prov = []
        for j in range(nf):
            if i==j:
                Prov.append(1)
            else:
                Prov.append(0)
        L.append(Prov)
    D = []
    for t in range(nf):
        D.append(0)
    M = []
    for t in range(nf):
        M.append(0)
    for i in range (nf):
        soma1=0
        for j in range (nf):
            M[j] = L[i][j] * D[j]
            soma1 += L[i][j]*M[j]
        D[i] = A[i][i] - soma1
        for j in range (i+1, nf):
            soma2=0
            for k in range (nf):
                soma2 += L[j][k]*M[k]
            L[j][i] = (A[j][i] - soma2)/D[i]
    Lt = np.transpose(L)
    Y = []
    Z = []
    X = []
    for i in range (nf):
        soma3=0 
        for j in range (i):
            soma3 += L[i][j]*Y[j]
        Y.append(B[i] - soma3)
    for i in range (nf):
        Z.append(Y[i]/D[i])
    for r in range(nf):
        X.append(0)
    for i in reversed(range(nf)):
        soma4=0
        for j in reversed(range(i+1,nf)):
            soma4 += Lt[i][j]*X[j] 
        X[i] = Z[i] - soma4
    return X
def main():
    N = int(input ('N:'))
    item = input('Qual item (a ou b)?: ')
    nf = int(input ('nf: '))
    P = []
    #TESTE A
    if item == 'a':
        p1 = float(input ('p1: '))
        P.append(p1)
    #TESTE B 
    if item == 'b':
        p1 = float(input ('p1: '))
        P.append(p1)
        p2 = float(input ('p2: '))
        P.append(p2)
        p3 = float(input ('p3: '))
        P.append(p3)
        p4 = float(input ('p4: '))
        P.append(p4)
    erro = N
    M = N
    dt = 1/N
    dx = dt
    A = []
    B = []
    X = []
    graficasso = Decomposicao(N, M, dt, dx, item, nf, P, erro)
    A, B = matrizAeB(N, M, dt, dx, item, nf, P, graficasso)
    X = rotinaLDLt(N, M, dt, dx, item, nf, P, graficasso, erro, A, B)
    print(X)
main()
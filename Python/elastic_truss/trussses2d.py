import pandas as pd
import numpy as np
from Constitutrel import *
global S, Breq, nmemb, L, A, funcname, member, deltasupp, ns, multiplier

inpuflname = pd.ExcelFile("compound planar.xlsx")
node = pd.read_excel(inpuflname, sheet_name='Node Data').values
member = np.delete(pd.read_excel(inpuflname, sheet_name='Member Data').values, 0, axis = 1)
funcname = pd.read_excel(inpuflname, sheet_name='Member Data', header=None).applymap(lambda x: str(x) if isinstance(x, str) else np.nan).dropna(axis=1, how='all').values
bc = pd.read_excel(inpuflname, sheet_name='Boundary Conditions').values
Force = pd.read_excel(inpuflname, sheet_name='Force Data').values
Suppsettle = pd.read_excel(inpuflname, sheet_name='Support Settlements').values
nnode = np.shape(node)[0]
nmemb = np.shape(member)[0]
L = np.zeros((nmemb,1))
A = member[:,3].reshape(-1,1)
F = Force[:, 1:3].reshape(-1,1)
Fm = np.max(np.abs(F))
avgforce = Fm
# Calculating the length for each member
for i in range(nmemb):
    lx = node[member[i, 2]-1, 1] - node[member[i, 1]-1, 1]
    mx = node[member[i, 2]-1, 2] - node[member[i, 1]-1, 2]
    L[i] = np.sqrt(lx ** 2 + mx ** 2)
# Setting up the equilibrium equations by arranging the corresponding equations for each node by tension coefficient method and assembling them in final master matrix(B)
B = np.zeros((2 * nnode, nmemb))
for i in range(nnode):
    mem, nodesn = np.where(member[:, 1:3] == i + 1)
    s = mem.size
    for k in mem:
        B[2*i, k] = (-node[member[k, 1]-1, 1] + node[member[k, 2]-1, 1]) / L[k]  # direction cosines lx, mx, and nx
        B[2*i +1, k] = (-node[member[k, 1]-1, 2] + node[member[k, 2]-1, 2]) / L[k]
    mem, nodesn = np.where(member[:, 2].reshape(-1,1) == i+1)
    for k in mem:
        B[2*i, k] = (node[member[k, 1]-1, 1] - node[member[k, 2]-1, 1]) / L[k]  # direction cosines lx, mx, and nx
        B[2*i+1, k] = (node[member[k, 1]-1, 2] - node[member[k, 2]-1, 2]) / L[k]
B = -B
# Finding positions of restraints in movement in different direction for all nodes, storing respective equations in Br and removing them from the main system of equations to be solved.
nodenum, xyz = np.where(bc[:, 1:3].T == 0)
ns = nodenum.shape[0]
Br = np.zeros((ns, nmemb))
Br = B[2 * (xyz) + nodenum, :]
Breq = B[np.setdiff1d(np.arange(1, 2*nnode+1), 2*(xyz)+nodenum+1)-1, :]
qs = Suppsettle.T.flatten()[nnode*(nodenum+1)+xyz]
deltasupp = np.dot(Br.T, qs)
BC = np.hstack((B, np.zeros((2*nnode, 2*nnode))))
nodepos = 2*(xyz) + nodenum 
reactpos = np.arange(nmemb + 2 * nnode - ns, nmemb + 2 * nnode)
BC[(2*nnode*reactpos+nodepos)%BC.shape[0], (2*nnode*reactpos+nodepos)//B.shape[0]] = -1
# Assembling the master matrix(S), the system of non-linear equations to be solved
sb = np.shape(Breq)
S = np.hstack((np.zeros((sb[1], sb[1])), Breq.T, np.zeros((sb[1], ns))))       # The matrix with B and B' as diagonal elements for the final set of equations to be solved
# For initial guess d0
d0 = avgforce*np.ones((sb[0] + sb[1]+ns, 1))
Fm = np.max(np.abs(F))
avgforce = Fm
numsigdigf = 1+ np.floor(np.log10(avgforce))
f = np.ones((nmemb, 1))
f = f*avgforce
deltacheck = Constitutrel(L,A,funcname,member,nmemb,f)[0]
de1 = np.max(np.abs(deltacheck))
nosigdigd = 1+ np.floor(np.log10(de1))
multiplier = 10**(numsigdigf-nosigdigd)
# Minization of error function
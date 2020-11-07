import numpy as np
import matplotlib.pyplot as plt


class Dislocation:
    def __init__(self, burgerVec, lineVec, position):
        self.BurgerVec = burgerVec
        self.LineVec = lineVec
        self.Position = position


class MaterialParameters:
    def __init__(self, nu, G):
        self.nu = nu
        self.G = G  # in GPa

#given a dislocation dis, this function returns the sigmaXX stress at point x0, y0
def sigmaXX(dis, x0, y0, materialParams):
    b = np.sqrt((dis.BurgerVec).dot(dis.BurgerVec))
    G = materialParams.G
    nu = materialParams.nu
    D = (G*b)/((2*np.pi)*(1-nu))
    disX = x0 - dis.Position[0]
    disY = y0 - dis.Position[1]
    sXX = (-1*D*(disY))*(((3*disX*disX)+(disY**2))/(((disX**2)+(disY**2))**2))
    return sXX

#given a mesh and a dislocation this function returns another 2d array (mesh) with values of sigmaXX at each point in mesh
#this calls sigmaXX for each point in mesh
def plotSigmaXX(dis, materialParams, meshX, meshY):
    sigmaXXMesh = np.zeros_like(meshX, dtype = np.float64)
    for i in range(meshY.shape[0]):
        for j in range(meshY.shape[1]):
            if ((meshX[i][j]==dis.Position[0]) and (meshY[i][j]==dis.Position[1])):
                continue
            else:
                a = sigmaXX(dis, meshX[i][j], meshY[i][j], materialParams)
                sigmaXXMesh[i][j] = a
    return sigmaXXMesh

#initialise material constants
materialParams = MaterialParameters(0.3, 30e9)  # (nu, G)
#make burger vec, line vec, and the position of dislocation
b1 = np.array([-2e-10, 0, 0])
l1 = np.array([0, 0, -1])
p1 = np.array([1e-6, 1e-6, 0])

#create a dislocation using above parameters
d1 = Dislocation(b1, l1, p1)

#create a mesh
#refer to a detailed answer about whe meshGrids: https://stackoverflow.com/a/49439331/9239725
xx, yy = np.meshgrid(np.linspace(0, 2e-6, 100), np.linspace(0, 2e-6,120))

#get sigmaXX for each point in mesh
stressXX = plotSigmaXX(d1, materialParams, xx, yy)
print(stressXX)

#plot it
plt.figure()  # 10 inch by 10 inch
#main plot
plt.pcolormesh(xx, yy, stressXX, vmin=-1e7, vmax=1e7, cmap='Spectral_r', shading='gouraud')
#include the side colorbar
plt.colorbar()

#add the edge Dislocation sign at the position of dislocation
Tsign = ((-7, 0), (0, 0), (0, -8), (0, 0), (7, 0))
Trevsign= ((-7, 0), (0, 0), (0, 8), (0, 0), (7, 0))
plt.plot(d1.Position[0],d1.Position[1],marker=Trevsign,markeredgecolor='white',markersize=10,markeredgewidth=2)
#include grids in  the figure
plt.grid()
plt.savefig("xx")

import numpy as np
import matplotlib.pyplot as plt


class Dislocation:
    def __init__(self, burgerVec, sign, position):
        self.BurgerVec = burgerVec
        self.Sign = sign
        self.Position = position


class MaterialParameters:
    def __init__(self, nu, G):
        self.nu = nu
        self.G = G  # in GPa

#given a dislocation dis, this function returns the sigmaXY stress at point x0, y0
def sigmaXY(dis, x0, y0, materialParams):
    b = np.sqrt((dis.BurgerVec).dot(dis.BurgerVec))
    G = materialParams.G
    nu = materialParams.nu
    D = (G*b)/((2*np.pi)*(1-nu))
    delX = x0 - dis.Position[0]
    delY = y0 - dis.Position[1]
    sXY = ((dis.Sign)*D*(delX))*(((delX**2)-(delY**2))/(((delX**2)+(delY**2))**2))
    return sXY

#given a mesh and a dislocation this function returns another 2d array (mesh) with values of sigmaY at each point in mesh
#this calls sigmaXY for each point in mesh
def plotSigmaXY(dis, materialParams, meshX, meshY):
    sigmaXYMesh = np.zeros_like(meshX, dtype = np.float64)
    for i in range(meshY.shape[0]):
        for j in range(meshY.shape[1]):
            if ((meshX[i][j]==dis.Position[0]) and (meshY[i][j]==dis.Position[1])):
                continue
            else:
                a = sigmaXY(dis, meshX[i][j], meshY[i][j], materialParams)
                sigmaXYMesh[i][j] = a
    return sigmaXYMesh

def overallStress(disArr, materialParams,meshX, meshY):
    sigmaXYMesh = np.zeros_like(meshX, dtype = np.float64)
    #for every Point
    for i in range(meshY.shape[0]):
        for j in range(meshY.shape[1]):
            #for every dislocation
            cumulativeForce = 0
            for dis in disArr:
                if ((meshX[i][j]==dis.Position[0]) and (meshY[i][j]==dis.Position[1])):
                    continue
                else:
                    cumulativeForce += sigmaXY(dis, meshX[i][j], meshY[i][j], materialParams)
            sigmaXYMesh[i][j] = cumulativeForce
    return sigmaXYMesh

def overallDisForcesXY(disArr,  materialParams):
    forceArr = np.zeros_like(disArr)
    #for every Dislocation
    for i in range(len(disArr)):
        #we will calculate force on dis due to all other dislocations
        dis = disArr[i]
        modB = np.sqrt((dis.BurgerVec).dot(dis.BurgerVec))
        sign = dis.Sign
        cumulativeForce = 0
        xPos = dis.Position[0]
        yPos = dis.Position[1]
        #pick all other dislocations and calculate the force becasue of them
        for dis2 in disArr:
            if ((xPos==dis2.Position[0]) and (yPos==dis2.Position[1])):
                #ignore the same dislocation because there is no force on a dislocation because of itself
                continue
            else:
                #else add the force due to dis2 on original dis
                #print("stress on dis ", i, " from other is ", sigmaXY(dis2, xPos, yPos, materialParams))
                cumulativeForce += sign*modB*sigmaXY(dis2, xPos, yPos, materialParams)
        forceArr[i] = cumulativeForce
    #now we have filled up the array containing forec on all dislocations due to all other dislocations
    return forceArr

def randomDislocationArr(numDislocation, xLimit, yLimit):
    disArr = []
    rng = np.random.default_rng()
    #we will generate random floats by first generating random integer and adding a random fraction value to it

    #we also maintain a set of already selected random points just so we dont accidentaly choose the same point again
    alreadySelectedPoints = set()

    #generate n random dislocaions
    while(len(alreadySelectedPoints)!= numDislocation):
        randomXposI = rng.integers(int(xLimit/(1e-6)))
        randomYposI = rng.integers(int(yLimit/(1e-6)))
        randomXposF = rng.random()
        randomYposF = rng.random()

        #add the integer and fraction part to get a random decimal value within limits
        xPos = randomXposI + randomXposF
        yPos = randomYposI + randomYposF

        if (xPos, yPos) not in alreadySelectedPoints:
            alreadySelectedPoints.add((xPos, yPos))
            b1 = np.array([-2e-10, 0, 0])
            #assign half negative half positive signs
            if (len(alreadySelectedPoints))%2==0:
                sign = +1
            else:
                sign = -1
            p = np.array([xPos*1e-6, yPos*1e-6, 0])
            dis = Dislocation(b1, sign, p)
            disArr.append(dis)
        else:
            continue
    return disArr

def printDislocationArr(darr):
    print('{:15} {:6}  {:6}   {:6}'.format('', "X-Cordinate", "Y - Cordinate", "Sign"))
    for i in range(len(darr)):
        print('{:15} {:11.4}  {:11.4}  {:11}'.format('Dislocation '+str(i), darr[i].Position[0], darr[i].Position[1], darr[i].Sign))

def plotStress(darr, xx, yy, stress,fileName):
    fig = plt.figure()  # 10 inch by 10 inch
    #main plot
    plt.pcolormesh(xx, yy, stress, vmin=-1e7, vmax=1e7, cmap='Spectral_r', shading='gouraud')
    #include the side colorbar
    plt.colorbar()
    #add the edge Dislocation sign at the position of dislocation
    Tsign = ((-7, 0), (0, 0), (0, -8), (0, 0), (7, 0))
    Trevsign= ((-7, 0), (0, 0), (0, 8), (0, 0), (7, 0))
    for dis in darr:
        if dis.Sign == +1:
            plt.plot(dis.Position[0],dis.Position[1],marker=Trevsign,markeredgecolor='white',markersize=10,markeredgewidth=2)
        if dis.Sign == -1:
            plt.plot(dis.Position[0],dis.Position[1],marker=Tsign,markeredgecolor='white',markersize=10,markeredgewidth=2)
    #include grids in  the figure
    plt.grid()
    plt.savefig(fileName)
    plt.close(fig)
#initialise material constants
materialParams = MaterialParameters(0.3, 30e9)  # (nu, G)
#make burger vec, line vec, and the position of dislocation
#b1 = np.array([-2e-10, 0, 0])
#l1 = np.array([0, 0, -1])
#sign = +1
#p1 = np.array([1e-6, 1e-6, 0])
#p2 = np.array([1e-6, 1.25e-6, 0])
xLimit = 3e-6
yLimit = 3e-6
#create a dislocation using above parameters
#d1 = Dislocation(b1, sign, p1)
#d2 = Dislocation(b1, sign, p2)
numDis = int(input("No of Dilocations?"))
darr = randomDislocationArr(numDis, xLimit, yLimit)
darrinit = darr
printDislocationArr(darr)
forcearr = overallDisForcesXY(darr, materialParams)
print("Force on dislocaions")
print(forcearr)
#create a mesh
#refer to a detailed answer about whe meshGrids: https://stackoverflow.com/a/49439331/9239725
xx, yy = np.meshgrid(np.linspace(0, xLimit, 100), np.linspace(0, yLimit,120))

#get sigmaXY for each point in mesh
stressXY = overallStress(darr, materialParams, xx, yy)
#print(stressXY)

plotStress(darr,xx,yy, stressXY,"start")

maxTimeSteps = int(input("MaxTimeSteps"))
maxDelX = ((np.float64(input("Maximum delta x? in terms of percentage of width")))/100)*xLimit
print(maxDelX)
for itr in range(maxTimeSteps):
    print('{:15} {}  {:4}'.format('', "Iteration", itr))
    forceXY = overallDisForcesXY(darr, materialParams)
    print("Force on dislocations")
    print(forceXY)
    maxForce = np.min(abs(forceXY))
    delT = maxDelX/maxForce
    #now we move each dislocation
    for i in range(len(darr)):
        dis = darr[i]
        #no change in y cordinate as no climb consideration assuming no diffusion at room temp
        delt = dis.Position[0] + (delT*forceXY[i])
        #print(delt)
        dis.Position[0] += delt
        if dis.Position[0] >xLimit:
            dis.Position[0] =  dis.Position[0] - xLimit*np.round(dis.Position[0]/xLimit)
        if dis.Position[0] <0:
            dis.Position[0] = dis.Position[0] + xLimit*np.round(abs(dis.Position[0]/xLimit)) + xLimit
    printDislocationArr(darr)
    #plot the iteration
    stressXY = overallStress(darr, materialParams, xx, yy)
    filename = "xy"+str(itr)
    plotStress(darr,xx,yy, stressXY,filename)

printDislocationArr(darr)
plotStress(darr,xx,yy, stressXY,"final")
maxDisp = np.float64(0)
for i in range(len(darr)):
    disp = darrinit[i].Position[0] - darr[i].Position[0]
    disp = disp - xLimit*np.round(disp/xLimit)
    if abs(disp) > maxDisp:
        maxDisp = disp

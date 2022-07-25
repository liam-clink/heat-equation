import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

## Functions
def Positions(xrange,yrange,resolution):
    '''
    Creates an arbitrary array of positions at a specified resolution.
    xrange&yrange must be a sequence (I used tuples)
    '''
    
    xaxis = np.arange(xrange[0],xrange[1]+resolution,resolution)
    yaxis = np.arange(yrange[0],yrange[1]+resolution,resolution)
    xvals,yvals = np.meshgrid(xaxis,yaxis,sparse=True)
    return xvals,yvals
            
def f(positions):
    '''
    Input array of independent var values and output an array
    '''
    z = np.zeros((positions[0].size,positions[1].size))
    for i in range(positions[0].size):
        for j in range(positions[1].size):
            x,y = positions[0][0][i], positions[1][j][0]
            
            z[i][j] = np.sin(x)+np.cos(y)
            #z[i][j] = -(np.cos(x)**2+np.cos(y)**2)
            #z[i][j] = np.exp(-x**2-y**2)
    return z

def Divergence(z,xy):
    '''
    Takes in 2darray and outputs array of divergence of the gradient
    '''
    
    newz = np.zeros_like(z)
    for i in range(xy[0].size):
        for j in range(xy[1].size):
            xgrad, ygrad = 0, 0
            if (i == 0 or i == (xy[0].size-1) or j == 0 or j == (xy[1].size-1)) == False:
                # Second derivative in the x direction then y direction
                # This weird indexing (not what I started with) is because of what the 3D plotter needs to work afaik, yes, it's annoying
                xgrad = (z[i+1][j]-2*z[i][j]+z[i-1][j])/((xy[0][0][i+1]-xy[0][0][i])*(xy[0][0][i]-xy[0][0][i-1]))
                ygrad = (z[i][j+1]-2*z[i][j]+z[i][j-1])/((xy[1][j+1][0]-xy[1][j][0])*(xy[1][j][0]-xy[1][j-1][0]))
                newz[i][j] = (xgrad+ygrad) # Because divergence is sum of second partial derivatives
    for i in range(xy[0].size):
        for j in range(xy[1].size):
            #Edges (Heat transfer = 0 (insulated), so derivative = 0 at edge)
#            if i == 0 and j not in (0,(xy[1].size-1)):
#                newz[i][j] = newz[i+1][j]
#            elif i == (xy[0].size-1) and j not in (0,(xy[1].size-1)):
#                newz[i][j] = newz[i-1][j]
            if j == 0 and i not in (0,(xy[0].size-1)):
                newz[i][j] = newz[i][j+1]
            elif j == (xy[1].size-1) and i not in (0,(xy[0].size-1)):
                newz[i][j] = newz[i][j-1]
    for i in range(xy[0].size):
        for j in range(xy[1].size):
            #Corners (average of two neighbors, probably a better way to do this)
            if i == 0 and j == 0:
                newz[i][j] = newz[i+1][j+1]
                print(newz[i][j])
            elif i == 0 and j == (xy[1].size-1):
                newz[i][j] = newz[i+1][j-1]
            elif i == (xy[0].size-1) and j == 0:
                newz[i][j] = newz[i-1][j+1]
            elif i == (xy[0].size-1) and j == (xy[1].size-1):
                newz[i][j] = newz[i-1][j-1]
    return newz
    
def Evolve(z,xy,dt,a=1.0):
    '''
    Code to evolve the distribution
    '''
    thermalConductivity = a
    newz = np.zeros_like(z)
    divergence = Divergence(z,xy)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            newz[i][j] = z[i][j] + (thermalConductivity*divergence[i][j])*dt
    return newz

def AxesSetup(xmin,xmax,ymin,ymax,zmin,zmax):
    ax.set_xlim3d([xmin,xmax])
    ax.set_xlabel('X')
    
    ax.set_ylim3d([ymin,ymax])
    ax.set_ylabel('Y')
    
    ax.set_zlim3d([zmin, zmax])
    ax.set_zlabel('Z')
    
    ax.set_title('3D Test')
    
# Functions for different types of plots
def SurfacePlot(frameNumber,xmin,xmax,ymin,ymax,zmin,zmax):    
    if frameNumber <= len(data):
        Z = data[frameNumber]
        #Refreshes the plot in between frames (prevents smear)
        ax.clear()   
        AxesSetup(xmin,xmax,ymin,ymax,zmin,zmax)
        #Drawing the plot
        ax.plot_surface(X, Y, Z, rstride=5,cstride=5,cmap=cm.jet,linewidth=0)

def TriSurfPlot(frameNumber,xmin,xmax,ymin,ymax,zmin,zmax):
    if frameNumber <= len(data):
        x = positions[0][0]
        y = positions[1].flatten()
        z = data[frameNumber].flatten()
        #Refreshes the plot in between frames (prevents smear)
        ax.clear()   
        AxesSetup(xmin,xmax,ymin,ymax,zmin,zmax)
        #Drawing the plot
        ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
        
def ContourPlot(frameNumber,xmin,xmax,ymin,ymax,zmin,zmax):
    if frameNumber <= len(data):
        Z = data[frameNumber]
        #Refreshes the plot in between frames (prevents smear)
        ax.clear()   
        AxesSetup(xmin,xmax,ymin,ymax,zmin,zmax)
        #Draws the plot
        ax.plot_surface(X, Y, Z, rstride=2, cstride=2, alpha=0.3)
        ax.contour(X, Y, Z, zdir='z', offset=zmin, cmap=cm.coolwarm)
        ax.contour(X, Y, Z, zdir='x', offset=xmin, cmap=cm.coolwarm)
        ax.contour(X, Y, Z, zdir='y', offset=ymax, cmap=cm.coolwarm)

## Main code
# Initial conditions
xmin,xmax,ymin,ymax,resolution = -3.14,3.14,-3.14,3.14,0.2
positions = Positions((xmin,xmax),(ymin,ymax),resolution)
temperature = f(positions)
zmin,zmax = np.amin(temperature),np.amax(temperature)
# Constants
a = 1.0
dt = 0.01
ti, tf = 0, 2

data = [] # List to store all the data points for later animation
counter = 0
for t in np.arange(ti,tf+dt,dt):
    temperature = Evolve(temperature,positions,dt,a)
    if counter % 10 == 0: #for settting fps
        data.append(temperature)        

## Code for making plots
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111, projection='3d')
X,Y = np.meshgrid(positions[0],np.transpose(positions[1]))

# Creating the Animation object
tempDistribution = animation.FuncAnimation(fig, SurfacePlot,len(data),
                                           fargs=(xmin,xmax,ymin,ymax,zmin,zmax),
                                           repeat=True,repeat_delay=0,
                                           interval=1000*(tf-ti)/len(data), blit=False)
#interval is in milliseconds
plt.show()

#tempDistribution.save('Heat Equation.mp4',codec='mpeg4',fps=len(data)/(tf-ti))
#print('Ding!')

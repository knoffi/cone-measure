import numpy as np
import matplotlib.pyplot as plt
import coneVol3 as cV
import matplotlib.pyplot as plt
import numpy as np
import matrix as M
import polygonTests as pT
center = [0,0]
lowerBound = -5
upperBound = 5
density = 0.5
X, Y = np.meshgrid(np.arange(center[0]+lowerBound, center[0]+upperBound, density), np.arange(center[1] + lowerBound, center[1] + upperBound,density))
x_shape = X.shape

U = np.zeros(x_shape)
V = np.zeros(x_shape)

#cD_test = cV.getConeVol( [[1,0],[0,1],[-1,0],[0,-1]])
#cD_test = cV.getConeVol( [[1,1],[-1,1],[-1,-1],[1,-1]])
P = [ [ 1,1 ], [-1,1] , [-1,-1 ], [2,-1]]
cD_test = cV.getConeVol( P )
for i in range(x_shape[0]):
    for j in range(x_shape[1]):
        point = [ X[i][j], Y[i][j]]
        try:
            v = M.subVek(cV.getConeVolIterator(cD_test,point),point)
            if(M.norm(v) != 0):
                v=M.scaleVek(0.8/M.norm(v) , v)
            else:
                v=[0,0]
                print(point)
            if v[0]**2+v[1]**2 < 0:
                v=[0,0]
        except ZeroDivisionError:
            #print(point)
            v=[0,0]
        U[i,j] = v[0]
        V[i,j] = v[1]

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, units='xy' ,scale=2, color='red')

ax.set_aspect('equal')

plt.xlim(center[0]+lowerBound,center[0]+upperBound)
plt.ylim(center[1]+lowerBound,center[1]+upperBound)

#plt.title('How to plot a vector field using matplotlib ?',fontsize=10)

#plt.savefig('how_to_plot_a_vector_field_in_matplotlib_fig1.png', bbox_inches='tight')
plt.show()
#plt.close()
#print(cV.getConeVolIterator(cD_test,[1,1]))
#print(cD_test)
#u= [cD_test[0][0],cD_test[0][1]]
#V = cD_test[0][2]
#print(u)
#print(V)
#pT.plotPoly(P,'r')
#print(cV.getNextVertex([cD_test[0][0],cD_test[0][1]],cD_test[0][2],[1,1]))
#print(cV.getConeVolIterator( cD_test , [1,1]))

point = [-0.5, 4.0]
P2 = [[-0.5, 4.0],[-1.0, 4.0], [-1.0, 2.0], [0.5, 2.0]]
print(P2)
pT.plotPoly(P2 , 'b')
pT.plotPoly(P,'r')
print(cV.getConeVol(P2))
print(cV.getConeVol(P))
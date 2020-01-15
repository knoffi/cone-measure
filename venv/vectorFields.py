import numpy as np
import matplotlib.pyplot as plt
import coneVol3 as cV
import matplotlib.pyplot as plt
import numpy as np
import matrix as M
import polygonTests as pT

const = 5
center = [const,const]
lowerBound = -const
upperBound = const
density = const / 10
X, Y = np.meshgrid(np.arange(center[0]+lowerBound, center[0]+upperBound, density), np.arange(center[1] + lowerBound, center[1] + upperBound,density))
x_shape = X.shape

U = np.zeros(x_shape)
V = np.zeros(x_shape)

x_fix = []
y_fix = []
x_sing = []
y_sing = []
#cD_test = cV.getConeVol( [[1,0],[0,1],[-1,0],[0,-1]])
#cD_test = cV.getConeVol( [[1,1],[-1,1],[-1,-1],[1,-1]])
P = [ [ 1,1 ], [-1,1] , [-1,-1 ], [2,-1]]
cD_test = cV.getConeVol( P )
for i in range(x_shape[0]):
    for j in range(x_shape[1]):
        point = [ X[i][j], Y[i][j]]
        try:
            v = M.subVek(cV.getConeVolIterator(cD_test,point),point)
            if v[0]**2+v[1]**2  < 0:
                v=[0,0]
                x_sing.append(point[0])
                y_sing.append(point[1])
            else:
                if(M.norm(v) != 0):
                    v=M.scaleVek(density/M.norm(v) , v)
                else:
                    v=[0,0]
                    x_fix.append(point[0])
                    y_fix.append(point[1])
        except ZeroDivisionError:
            #print(point)
            x_sing.append(point[0])
            y_sing.append(point[1])
            v=[0,0]
        U[i,j] = v[0]
        V[i,j] = v[1]

fig, ax = plt.subplots()
ax.quiver(X, Y, U, V, units='xy' ,scale=2, color='red')

ax.set_aspect('equal')
ax.plot(x_fix, y_fix, 'go')
ax.plot(x_sing, y_sing, 'bo')

plt.xlim(center[0]+lowerBound,center[0]+upperBound)
plt.ylim(center[1]+lowerBound,center[1]+upperBound)

#plt.title('How to plot a vector field using matplotlib ?',fontsize=10)

plt.savefig('trapez05x55.png', bbox_inches='tight')
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

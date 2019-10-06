import math
import matrix as M
import machEps as mE
import matplotlib.pyplot as mp
import randomPolygon3Centered as rP

# only works in R^2
# is this counterclockwise? like my new polygon orderings... ?
#changed

cD_test = [ [ 1 , 0 , 1] ,  [ 0 , 1 , 1] , [ -1 , 0 , 1 ] , [ 0 , -1 , 1 ] ]
point_test = [ 1 , -1 ]
def nextVertex(u, volume, point):
    orth = []
    # gegen uhrzeigersinn um 90 grad drehen
    orth.append( - u[1] )
    orth.append( u[0] )

    # assumes that 0 lies in K (i.e. point is not orthogonal)
    d = M.scaleVek(2 * volume / M.scal(point, u), orth)
    temp = M.addVek(point, d)
    point[0] = temp[0]
    point[1] = temp[1]

def getNextVertex(u, volume, point):
    result = []
    nextVertex(u, volume, point)
    for x in point:
        result.append(x)
    return result


nextPoint_test = getNextVertex( [ 1 , 0 ] , 1 , point_test)
if M.dist( nextPoint_test , [ 1 , 1] ) > 0.00001:
    print( nextPoint_test)
    print( 'Fehler bei NextPoint' )

def domain(start, end, step):
    result = []
    i = 0
    while (start + i * step <= end):
        result.append(start + i * step)
        i += 1
    return result

# IMPORTANT: does nextVertex work with local variable rules of python?
def coneVolIterator(coneData, point):
    result = []
    n = len(coneData)

    for i in range( n ):
        facet = coneData[ ( i + 1 ) % n ]
        nextVertex([facet[0], facet[1]], facet[2], point)
    for x in point:
        result.append(x)
    return result

def getConeVolIterator(coneData, point):
    result = []
    for x in point:
        result.append(x)
    coneVolIterator(coneData, result)
    return result


if M.dist(getConeVolIterator( cD_test , point_test ) , point_test ) > 0.00001 :
    print( 'Fehler bei coneVolIterator' )

# berechnet abstand zum quadrat
def polyFunctional(coneData, point):
    return M.norm(M.subVek(point, getConeVolIterator(coneData, point))) ** 2

if polyFunctional( cD_test , point_test ) > 0.0001:
    print( ' Fehler bei polyFunctional ')

def arrayPolyFunctional1(coneData, point, param):
    result = []
    for i in param:
        result.append(polyFunctional(coneData, M.scaleVek(i, point)))
        if polyFunctional(coneData, M.scaleVek(i, point)) < eps_graph:
            print(i)
    return result

# changed
def diffGetNextVertex(u, volume, point):
    B = M.Matrix([ [ -u[1] ] , [ u[0] ] ])
    A = M.Matrix( [ [ -u[0] , -u[1] ] ] )
    C = B.mult(A)
    c = M.scal(point, u)
    c = 2 * volume / (c**2)
    C.scale(c)
    # Matrix.plus oder Matrix.add verwenden?
    C.add(M.idMatrix(2))
    return C

def arrayPolyFunctional2(coneData, param):
    result = []
    for i in param:
        result.append(polyFunctional(coneData, [1, i]))
    return result

#changed
def gradGetNextVertex(u, volume, point):
    result = diffGetNextVertex( u , volume , point)
    return result.copyTrans()

#gradGetNextVertex([1,1],1,[1,1]).list()



# Does the for-loop start with the correct facet? It depends on coneData...
def diffConeVolumeIterator(coneData, point):
    result = M.idMatrix(2)
    n = len( coneData )

    for i in range( n ):
        facet = coneData[ (i + 1) % n ]
        result = diffGetNextVertex([facet[0], facet[1]], facet[2], point).mult(result)
        nextVertex([facet[0], facet[1]], facet[2], point)
    return result

    diff_test = diffConeVolumeIterator( cD_test , point_test )
    exspectedDiff_test = Matrix ( [ [ 0 , 0 ] , [ 0 , 0 ] ] )
    if M.distMatrix( diff_test , exspectedDiff_test ) > 0.00001:
        print( 'Fehler bei diffConeVolume' )

# changed
def gradConeVolumeIterator(coneData, point):
    result = diffConeVolumeIterator( coneData , point )
    return result.copyTrans()

def gradPolyFunctional(coneData, point):
    B = M.idMatrix(2)
    A = gradConeVolumeIterator(coneData, point)
    A.scale(-1)
    B.add(A)

    b = M.subVek(point, getConeVolIterator(coneData, point))
    b = M.scaleVek(2, b)

    result = B.image(b)
    return result


def gamma(coneData, params):
    # d_0 und d_1 verdrehen?
    d_0 = [-coneData[0][1], coneData[0][0]]
    d_1 = [ coneData[1][1], -coneData[1][0]]
    v = M.scaleVek(params[0], d_0)
    w = M.scaleVek(params[1], d_1)
    return M.addVek( v , w )


def gradGamma(coneData, params):
    d_0 = [-coneData[0][1], coneData[0][0]]
    d_1 = [coneData[1][1], -coneData[1][0]]

    # ist das richtig?
    return M.Matrix([d_0, d_1])


def diffGamma(coneData, params):
    d_0 = [-coneData[0][1], coneData[0][0]]
    d_1 = [ coneData[1][1], -coneData[1][0]]

    return M.Matrix([d_0, d_1]).copyTrans()


def phi(params,cD):
    return 1 * polyFunctional(cD, gamma(cD,params))

# is this the correct gradient?
def gradPhi(params,cD):

    A = gradGamma(cD, params)
    v = gradPolyFunctional(cD, gamma(cD, params))

    return M.scaleVek( 1 , A.image(v))

def gradPhiApprox( params , cD , h ):
    e_1 = [1, 0]
    e_2 = [0, 1]

    point_diff1 = M.addVek(params, M.scaleVek(h, e_1))
    quotient_diff1 = phi(point_diff1 , cD ) / h  - phi(params , cD ) / h

    point_diff2 = M.addVek(params, M.scaleVek(h, e_2))
    quotient_diff2 = phi(point_diff2, cD) / h - phi(params, cD) / h

    return [ quotient_diff1 , quotient_diff2 ]



def g(a, x):
    return a * x[0] ** 2





def gradG(x):
    return [2 * x[0]]

def F(cD,params):
    return M.subVek(getConeVolIterator(cD,gamma(cD,params)),gamma(cD,params))

def psi(cD , params):
    return M.norm( F(cD,params))

def diffF(cD, params):
    A=diffGamma(cD,params)
    B=diffConeVolumeIterator(cD,gamma(cD,params))
    C=M.idMatrix(2)
    C.scale(-1)
    B.add(C)
    return B.mult(A)
# hier einfach mal die Norm hoch Anzahl der nextPolyIteratorSteps multiplizieren
def regulator(cD , params):
    return M.norm(F(cD,params)) #* M.norm(params)

def coneVolumeNewton(cD,params):
    value=float("inf")
    print( ' newton begins ')
    print(params)
    n = 0
    while(regulator( cD , params ) > eps_newton ):
        #print(regulator(cD , params))
        v=F(cD,params)
        C=diffF(cD,params)
        params=M.subVek(params,C.lgsSolve(v))
        #print(value)
        if n % 50 == 0:
            print(params)
        n = n + 1
    return params

# es ist wichtig, dass u_1 in der richtigen reihenfolge ist und als erster rechts von p steht
# p=[1,10] hat eine sehr interessanten Minimizer auf spann([1,10])

#x = domain(0.7, 19, 0.03)
#y = arrayPolyFunctional2(cD, x)
# IMPORTANT: algorithm throws exception if one of the p_i is orthogonal to a u_i.
#print(coneVolumeNewton(cD,p))


def getClockRotation(vector):
    return[ vector[ 1 ] , -vector[ 0 ] ]


# exspect vertices to be ordered counter clockwise
# muss überprüft werden
def getConeVol(vertices):
    result=[]
    n=len(vertices)
    for i in range( n ):
        v = M.subVek( vertices[ i ] , vertices[ (i - 1 + n) % n ])
        u_tilde = getClockRotation( v )
        divisor = M.norm(v)
        if u_tilde[0] == 0 and u_tilde[1] == 0:
            divisor = mE.getMachEps()
        u = M.scaleVek( 1.0 / divisor , u_tilde )

        height= M.scal(u,vertices[i-1])
        if height < 0:
            print('height is negativ at ')
            print(i)
        facetVolume= M.norm(v)
        coneVolume= 0.5 * height * facetVolume

        result.append([u[0],u[1],coneVolume])
    return result

A_test = M.Matrix( cD_test )
cD_resultTest = getConeVol( [[ 1 , 1 ] , [ -1 , 1 ] , [ -1 , -1 ], [ 1 , -1 ] ] )
B_test = M.Matrix( cD_resultTest )
if M.distMatrix( A_test , B_test ) > 0.00001:
    print( ' Fehler bei getConeVol in preSolver2 ')
    print( cD_resultTest )


# ziemlich unsymmetrisches fünfeck, 5. und 3. facette sind parallel

eps_newton=0.0000000000001
eps_graph= 0.0000001
const_step= 0.1
eps_descend = 0.1
#cD_test = [[0, 1, 1], [1, 0, 1], [ 0 , -1, 1],[-1,0,1],[-1, 2,1]]
p=[1,0.5]

# cD2 wird gar nicht bei start verwendet, sondern cD....
# cD2=getConeVol([[1,1],[1,-1],[-1,-1],[-1,1]])
# print(coneVolumeDescend(cD,startPoint()[0]))

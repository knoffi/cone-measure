from cmath import sqrt
import math
import numpy
import random
import coneVol as cV
import matrix as M
import matplotlib.pyplot as mp
import machEps as mE

eps_max = 4096 * mE.getMachEps()
eps_min = 4096 * mE.getMachEps()

generalMaxRadiusBound = 50

def getRandomPolygon(number):
    polarResult = getRandomNoncenteredPolarPolygon(number)
    result = getCartesian( polarResult )
    # maybe this should be done in every step. But then again,
    makeCentered( result )
    #makeUnitVolume( result )

    return result

def getRandomNoncenteredPolarPolygon(m):
    polarResult = getRandomOriginPolarTriangle()
    steps = 3
    while (steps < m):
        # maybe it works fine, if the start triangle is just centered....
        # only works, if polarResult is centered, right? ... this means I have to switch back and forth from polar to cartesian ... xS
        angle = getRandomAngle(polarResult)
        orderedPolarInsertion(polarResult, angle)
        steps = steps + 1
    return polarResult

def makeUnitVolume( orderedPolytop ):
    coneVolume = cV.getConeVol( orderedPolytop )
    volume = 0
    for data in coneVolume:
        volume = volume + data[2]

    for i in range( len(orderedPolytop)):
        orderedPolytop[i-1] = M.scaleVek( 1.0 / math.sqrt( volume ) , orderedPolytop[i-1] )




def angleDist ( alpha, beta ):
    if ( alpha <= math.pi and beta <= math.pi ):
        return abs( alpha - beta )
    if ( alpha >= math.pi and beta >= math.pi ):
        return abs( alpha - beta )
    if abs( alpha - beta ) <= math.pi :
        return abs( alpha - beta )
    return 2 * math.pi - abs( alpha - beta )

def noOrigin( alpha, beta, gamma ):
    count = 0
    if ( angleDist( alpha, beta ) + angleDist( beta , gamma ) > math.pi ):
        count = count+1
    if (angleDist(beta, alpha) + angleDist(alpha, gamma) > math.pi):
            count = count + 1
    if (angleDist(alpha, gamma) + angleDist(beta, gamma) > math.pi):
        count = count + 1

    if count == 3:
        return False
    return True

#alpha = 0.5 * math.pi
#beta = 1 * math.pi
#gamma = 1.5 * math.pi
# print( angleDist( alpha , beta) / math.pi )
#print( noOrigin( alpha , beta , gamma ))

def firstComponent( point ):

    return point[0]

def getRandomOriginPolarTriangle():
    alpha = random.random() * 2 * math.pi
    beta = random.random() * 2 * math.pi
    gamma = random.random() * 2 * math.pi
    #angles = randomAngles( 3 )
    #alpha = angles[ 0 ]
    #beta = angles[ 1 ]
    #gamma = angles[ 2 ]
    r_1 = 1 / random.random()
    r_2 = 1 / random.random()
    r_3 = 1 / random.random()

    while noOrigin( alpha, beta, gamma):
        alpha = random.random() * 2 * math.pi
        beta = random.random() * 2 * math.pi
        gamma = random.random() * 2 * math.pi

    result= [  [ alpha , r_1 ] , [ beta , r_2 ] , [ gamma , r_3 ] ]
    result.sort( key = firstComponent )
    return result



def getCartesian( vertices ):
    result = [ ]
    for v in vertices:
        result.append( polar( v[0] , v[1] ) )
    return result

# seems to work well, was tested in polygonTest...
def makeCentered( vertices ):

    center = [ 0 , 0 ]
    i=0

    for v in vertices:
        center = M.subVek( center , v )
    center = M.scaleVek( 1 / len(vertices) , center )

    while i < len(vertices):
        vertices[ i ] = M.addVek( vertices[ i ] , center )
        i = i + 1



# better distribution: do not set angle of first point as 0
def orderedPolarInsertion( polarVertices , angle ):

    neighbours = getCartesian(getPolarNeighbours( polarVertices , angle ))
    left2 = neighbours[0]
    left1 = neighbours[1]
    right1 = neighbours[2]
    right2 = neighbours[3]

    max = maxRadius( [ left2, left1 , right1 , right2 ] , angle ) - eps_max
    min = minRadius( left1, right1 , angle ) + eps_min

    if min > max:
        print('oh shit')
        print(neighbours)
        print( max )
        print( min )
        print( angle )

    r = randomRadius( min , max )
    newVertex = [ angle , r ]
    # UPDATE : Habe momentan nicht mehr zu wenige Vertices... höchstens zu viele in manchmal Fällen...
    # gibt es in dieser Schleife einen Fehler, der für zu wenige oder zu viele Angles sorgt? (manchmal könnte der gleiche Punkt mehrmals hnzugefügt werden, was für Normalenvektor gleich 0 sorgen könnte ..)
    orderedFirstComponentInsertion( polarVertices , newVertex )

def orderedFirstComponentInsertion( orderedVectorArray, vector ):
    index = 0
    for elements in orderedVectorArray:
        if elements[0] < vector[0]:
            index = index + 1
        else:
            break
    orderedVectorArray.insert( index , vector )

#test_array = [ [ 0 , 1 ] , [ 2 , 0 ] , [ 63 , -1 ] ]
#test_vector = [ 0 , 19 ]
#orderedFirstComponentInsertion( test_array , test_vector)
#print( test_array )

def getAngle(point):
    x=point[0]
    y=point[1]
    if x  == 0 and y > 0:
        return 0.5 * math.pi
    if x == 0 and y < 0:
        return 1.5 * math.pi
    if x > 0 and y > 0:
        return math.atan(y/x)
    if x > 0 and y < 0:
        return math.atan(y/x) + 2* math.pi
    if x < 0:

        return math.atan(y/x) + math.pi
    return 0

#better distribution: select random n random angles between 0 and pi, then order them and return an array with every angle
# berücksichtigt das hier die periodizität, also pi = 3*pi


# can be improved if ' u ' is already in polar coordinates
def supportFunction( K_polar , u ):
    angle = getAngle( u )
    neighbours = getPolarNeighbours( K_polar , angle )
    v = getCartesian([neighbours[1]])[0]
    w = getCartesian([neighbours[2]])[0]
    h_tilde = getMin( v , w , angle)
    divisor = M.norm(u)
    if divisor == 0:
        print('warum ist u als Normalenvektor so klein?')
        print( u )
        return h_tilde / mE.getMachEps()
    return h_tilde / M.norm(u)


def angleIsContained( polarVertices , angle ):
    for vector in polarVertices:
        if vector[0] == angle:
            return True
    return False



def getRandomAngle( polarVertices ):
    angle = random.random() * 2 * math.pi

    while angleIsContained( polarVertices , angle ):
        angle = random.random() * 2 * math.pi
    return angle

def randomAngles( n ):
    result = [ 0 ]
    steps = 1

    while steps < n :
        value = random.random() * 2 * math.pi

        #while( isContained( result , value ) ):
        #    value = random.random() * 2 * math.pi
        result.append( value )
        steps = steps + 1

    result.sort()
    return result



# print( randomAngles( 5 ) )


# has to be improved, so that it will give a good distributed number between 0 and +infty
def randomRadius(min, max):
    #randomRadius will always be below 10 + min . Is that good?
    if( max == math.inf ):
        if( min < generalMaxRadiusBound ):
            return random.random() * ( generalMaxRadiusBound - min ) + min
        else:
            return random.random() *  min + min
    return min + ( max - min ) * random.random()

def polar(angle, radius):
    x= radius * math.cos(angle)
    y= radius * math.sin(angle)
    return[x,y]

# may be better to save points in polar coordinates, because maxRadius would not be forced to calculate so often arctan(y/x)
def maxRadius( neighbours , angle ):
    points = neighbours
    r1 = getMax( points[0] , points[1] , angle )
    r2 = getMax( points[2] , points[3] , angle )
    if r1 < r2:
        return r1
    return r2



def minRadius( left1, right1, angle ):
    return getMin( left1 , right1 , angle)


# may be better, if centeredVertices are represented in polar coordinates
def getPolarNeighbours( centeredPolarVertices , angle ):
    index = 0
    length = len( centeredPolarVertices )

    for v in centeredPolarVertices:
        if v[0] > angle:
            right2 = centeredPolarVertices[ ( index - 2 ) % length  ]
            right1 = centeredPolarVertices[ ( index - 1 ) % length  ]
            left1 = centeredPolarVertices[ ( index ) % length  ]
            left2 = centeredPolarVertices[ ( index + 1 ) % length  ]
            return [ left2 , left1 , right1 , right2 ]
        index = index + 1

    return [ centeredPolarVertices[ 1 ] , centeredPolarVertices[ 0 ] , centeredPolarVertices[ length - 1 ] , centeredPolarVertices[ length -2 ] ]









# test = [ [ 1 , 1 ] , [ 1 , -1 ] , [ -1 , 1 ] , [ -1 , -1 ] ]
#print(getNeighbours( test , 0))

def getMax(v,w,angle):
    u_1 = polar(angle,1)
    u_2 = M.subVek(v,w)
    #print(u_1)
    U= M.Matrix([u_1,u_2]).copyTrans()
    #U.list()
    result=U.lgsSolve(v)
    # my lgsSolver returns a 'solution' even if v is not in im(M) . Therefore the if-check needs to get adjusted
    if(len(result) < 2):
        return math.inf
    if( result[0] <= 0 ):
        return math.inf
    return result[0]

def getMin(v,w,angle):
    u_1 = polar(angle,1)
    u_2 = M.subVek(v,w)
    U = M.Matrix([u_1,u_2]).copyTrans()
    result=U.lgsSolve(v)
    if( len(result) < 2 ):
        return 0
    if ( result[0] <= 0):
        return 0
    return result[0]



# test = [ [ 0 , 2 ] , [ 0.25 * math.pi , 1.42  ] , [ 0.75 * math.pi , 1.42 ] , [ 1 * math.pi , 2 ] ]
# neighbours = getNeighbours( test , 0.5 * math.pi )
# print( maxRadius( neighbours , 0.5 * math.pi ) )

v=[0,1]
w=[1,0]
angle= math.pi/4
vertices=[ [ 1, 0 ],  [ 0 , -1 ] , [ -1 , 0 ] ]
# ein result mit zwei gleichen Punkten: [[1.3506563545534096, 0.0], [1.8164857151592047, 0.0], [1.3506563545534096, 0.0], [1.4946652011162949, 0.0]]
#P = getRandomOriginPolarTriangle()
# - minRadius testen und einbauen in getNextPoint von getRandomPolygon()
# - minRadius ist notwendig, aber vllt muss man bei maxRadius und minRadius nur die Vorgänger überprüfen ... -> mache 2
# - habe einmal [[44.44310714415645, 0.0], [-2.868745831991136, 18.486864345855214], [-14.100910361720999, -23.4225736149311], [-13.38842882794231, -23.35396107847429], [-14.085022122502007, 28.289670347550175]]
#   rausbekommen. Das ist kein 5-Eck, zwei Punkte sind in der konvexen Hülle der anderen enthalten...



# oh shit
# [[1, 0], [0, -1], [-1, 0]]
# 1.4898145505571896
# 1.5915425615297067
# 3.404131410794002
# oh shit
# [[1, 0], [0, -1], [-1, 0]]
# 1.7818867557071751
# 1.8200557773619472
# 3.6352952855406166

# vertices = [ [ 1 , 1 ] , [ -1 , 1 ] , [ 0 , 1 ] ]
# testPoint = [ 1 , 1 ]
# print( containsPoint( vertices , testPoint ) )
# [[2.6280147103737663, -3.9985536987429766], [-1.2515947854797322, 0.22196113640884813], [0.7102207863932186, 0.856232362200955]] may have been angles were the random triangle did not contain the origin! (mixed up checks, I am not sure about the tests...)
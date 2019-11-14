from cmath import sqrt
import math
import numpy
import random
import coneVol as cV
import matrix as M
import matplotlib.pyplot as mp
import machEps as mE
import fractions as f

eps_max = 4096 * mE.getMachEps()
eps_min = 4096 * mE.getMachEps()

generalMaxRadiusBound = 1

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

def getVolume( orderedPolygon ):
    coneVolume = cV.getConeVol(orderedPolygon)
    volume = 0
    for data in coneVolume:
        volume = volume + data[2]
    return volume

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

def getPolar( vertices_polar ):
    result = [ ]
    for v in vertices_polar:
        result.append( reversePolar( v[0] , v[1] ) )
    return result

# seems to work well, was tested in polygonTest...
def makeCentered( vertices ):

    center = getCenter(vertices)

    for point in vertices:
        point[0] = point[0] - center[0]
        point[1] = point[1] - center[1]

def translateRational( verticesRational , translationRational):

    for point in verticesRational:
        point[0] = point[0] - translationRational[0]
        point[1] = point[1] - translationRational[1]

def getCenter( vertices ):

    P = vertices
    x = 0
    y = 0
    n = len(P)
    volume = 0

    for i in range(n):
        p = P[i]
        q = P[(i + 1) % n]

        determinante = (p[0] * q[1] - q[0] * p[1])

        x += (p[0] + q[0]) * determinante
        y += (p[1] + q[1]) * determinante

        volume += 0.5 * determinante

    x = x / ( 6.0 * volume )
    y = y / ( 6.0 * volume )
    center = [ x , y ]
    return center

def getCenterAndVolumeRational( verticesRational ):

    P = verticesRational
    x = f.Fraction( 0 , 1)
    y = f.Fraction( 0 , 1)
    n = len(P)
    volume = f.Fraction( 0 , 1)

    for i in range(n):
        p = P[i]
        q = P[(i + 1) % n]

        determinante = (p[0] * q[1] - q[0] * p[1])

        x += (p[0] + q[0]) * determinante
        y += (p[1] + q[1]) * determinante

        volume += f.Fraction( 1 , 2) * determinante

    x = x / ( f.Fraction( 6 , 1) * volume )
    y = y / ( f.Fraction( 6 , 1) * volume )
    center = [ x , y ]
    return [ center , volume ]


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

def downRoundVertices( vertices, digits ):

    for point in vertices:

        factor = 10**digits

        point[0] = math.floor( point[0] * factor ) / factor
        point[1] = math.floor( point[1] * factor ) / factor

def roundVertices(vertices, digits):
    for point in vertices:
        factor = 10 ** digits
        stretched_x = point[0] * factor
        stretched_y = point[1] * factor
        if math.floor( stretched_x * 10 ) % 10 < 5:
            point[0] = math.floor(stretched_x) / factor
        else:
            point[0] = ( math.floor(stretched_x) + 1 ) / factor

        if math.floor( stretched_y * 10 ) % 10 < 5:
            point[1] = math.floor(stretched_y) / factor
        else:
            point[1] = ( math.floor(stretched_y) + 1 ) / factor


tester_1 = [ [ 1.1 , 1.5 ] , [ 1.6 , 1.1 ] , [ 1.1 , 1.1 ] , [ 1.7 , 1.7 ] ]

tester_2 = [ [ 1.1 , 1.5 ] , [ 1.6 , 1.1 ] , [ 1.1 , 1.1 ] , [ 1.7 , 1.7 ] ]

tester_roundDownResult = [ [ 1 , 1 ] , [ 1 , 1 ] , [ 1 , 1 ] , [ 1 , 1 ] ]
tester_roundResult = [ [ 1 , 2 ] , [ 2 , 1 ] , [ 1 , 1 ] , [ 2 , 2 ] ]
downRoundVertices(tester_1 , 0 )
roundVertices( tester_2 , 0)
if(M.distMatrix( M.Matrix(tester_roundDownResult) , M.Matrix(tester_1)) > 0 ):
    print('Fehler bei roundDownVertices')
if( M.distMatrix( M.Matrix(tester_roundResult) , M.Matrix(tester_2))) > 0:
    print('Fehler bei roundVertices')


def makeBaryCentered( cartesianPolygon ):
    P = cartesianPolygon
    barycenter = [ 0 ,  0 ]
    n = len(P)

    for vertex in P:
        barycenter = M.addScaleVek( barycenter , 1 / n , vertex )

    for vertex in P:
        vertex[0] -= barycenter[0]
        vertex[1] -= barycenter[1]

def getBaryCenter( cartesianPolygon ):
    P = cartesianPolygon
    barycenter = [ 0 ,  0 ]
    n = len(P)

    for vertex in P:
        barycenter = M.addScaleVek( barycenter , 1 / n , vertex )

    return barycenter

# can be improved if ' u ' is already in polar coordinates
def supportFunctionPolar( K_polar , u ):
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

# can be improved if ' u ' is already in polar coordinates
def supportFunctionCartesianCentered( K_cartesianCentered , u ):
    if math.fabs( M.norm( u ) - 1 ) > 0.0001:
        print( ' u is not normed for support function, result needs to be divided by norm ')
        return []
    P = K_cartesianCentered
    n = len(P)
    result = math.inf

    for i in range(n):
        scaling = getMinForSupport( P[ i % n] , P[ i-1 ] , u )
        if scaling > 0 and scaling < result:
            result = scaling

    return result

def supportFunctionCartesianCenteredRational(K_cartesianCenteredRational, u_Rational):

    P = K_cartesianCenteredRational
    n = len(P)
    result = math.inf

    # luckily, norm of u_i cancels in the fraction of the logMin product inequality...
    # hier muss ich noch dran arbeiten...

    for i in range(n):
        # I do not need a special rational support function, do I?
        scaling = getMinForSupport(P[i % n], P[i - 1], u_Rational)
        if scaling > 0 and scaling < result:
            result = scaling

    return result


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

def reversePolar( x_koord , y_koord ):

    angle = getAngle( [ x_koord , y_koord ])
    r = math.sqrt( x_koord *+2  + y_koord ** 2 )

    return [ angle , r ]

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

def getMinForSupport( v , w , u ):
    u_1 = u
    u_2 = M.subVek(v,w)
    #print(u_1)
    U= M.Matrix([u_1,u_2]).copyTrans()
    #U.list()
    result=U.lgsSolve(v)
    # my lgsSolver returns a 'solution' even if v is not in im(M) . Therefore the if-check needs to get adjusted
    if( M.dist( U.image(result) , v ) > 0.0001 ):
        return -math.inf
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


# works like round but is rounding -0.5 to -1 isntead of 0
def roundAwayFromZero( number, digits):
    if( number >= 0):
        x = number * 10**digits
        y = number * 10**(digits + 1) % 10

        if y < 0:
            print( 'error in round?')
        if y < 5:
            return math.floor(x) / ( 10 ** digits)
        else:
            return math.floor(x + 1 ) / ( 10 ** digits)
    else:
        return - roundAwayFromZero( -number, digits )



def getRationalWithDigits(K , digits ):

    result = []
    for v in K:
        a_0 = roundAwayFromZero( v[0] , digits ) * 10**digits
        a_1 = roundAwayFromZero( v[1] , digits ) * 10**digits
        numerator_0 = math.floor( a_0 )
        numerator_1 = math.floor( a_1 )

        if numerator_0 - a_0 > 0.5 :
            print( 'approximated float number is lower than a0, I am trying to fix this')
            numerator_0 = math.floor( a_0  + 0.5 )

        if numerator_1 - a_1 > 0.5 :
            print( 'approximated float number is lower than a0, I am trying to fix this')
            numerator_1 = math.floor( a_1  + 0.5 )


        denominator_0 = 10 ** digits
        denominator_1 = 10 ** digits

        result.append(  [ f.Fraction( numerator_0 , denominator_0) , f.Fraction( numerator_1 , denominator_1 ) ] )

    return result
def printRationalPolygon( rationalPolygon ):
    string_result = ' [ '
    for point in rationalPolygon:
        frac0 = point[0]
        frac1 = point[1]
        if rationalPolygon.index(point) < len(rationalPolygon) - 1 :
            string_result = string_result + ' [ %d/%d , %d/%d ] ,' % (frac0.numerator , frac0.denominator , frac1.numerator , frac1.denominator )
        else:
            string_result = string_result + ' [ %d/%d , %d/%d ] ' % (frac0.numerator, frac0.denominator, frac1.numerator, frac1.denominator)

    print(string_result + ' ] ')

test_getRational = [ [ 1.19 , 1.37 ] , [ -0.901 , 0.72 ] , [ -0.5 , -0.5  ] , [ 1.2075639364836 , -1.17947621917]]
cube_rational = getRationalWithDigits( test_getRational , 0 )




M_test = M.Matrix( [ [ f.Fraction( 0 , 1) , f.Fraction( 2 , 1) ] , [ f.Fraction( 2 , 1) , f.Fraction( 0 , 1) ]])
b_test = [ f.Fraction( 10 , 1) , f.Fraction( 8 , 1) ]
solutuion_test = M_test.lgsSolve( b_test )
exspectedSolution_test =  [ f.Fraction( 4 , 1) , f.Fraction( 5 , 1) ]

if M.dist( exspectedSolution_test , solutuion_test ) > 0.001 :
    print( exspectedSolution_test )
    print( solutuion_test )
    print( ' Fehler in lgs solve für rationale Zahlen')


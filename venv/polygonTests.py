from cmath import sqrt
import math
import numpy
import random
import coneVol as cV
import matrix as M
import matplotlib.pyplot as mp
import randomPolygon3Centered as rP
import machEps as mE
import randomRationalPolygon as rRP
import fractions as f


# eventuell wird ein tester besser, wenn ich die polarkoordinaten nicht komplett umwandle
# HIER SOLLTE ALS NÄCHSTES GEARBEITET WERDEN; UM TESTER ZU VERBESSERN ::: WELCHE SIND DIE SINNVOLLEN EPSILONS FÜR DIE TEST-SCHRANKENß SOLLTE MAN SCHNELL MAL MACHINE EPSILON ZUVOR BERECHNEN? UND KONDITIONIEREUNG BESTIMMEN?
# warum kriege ich so oft  benachbarte vertices mit fast gleichem Winkel (Zufall?) und (unnötig) fast gleichem Radius (Fehler bei minRadius, maxRadius) ?

containsPoint_eps = 0.1 #1000000 * mE.getMachEps()
convexRun_eps = containsPoint_eps

def getPolyDomain(vertices):
    a = vertices[0][0]
    b = vertices[0][0]
    c = vertices[0][1]
    d = vertices[0][1]

    for v in vertices:
        if v[0] < a :
            a = v[0]
        if v[0] > b :
            b = v[0]
        if v[1] < c :
            c = v[1]
        if v[1] > d :
            d = v[1]
    return [ a -1 , b +1  , c -1 , d +1 ]

def plotPoly(cartesianVertices, colorString):
    x=[]
    y=[]

    for v in cartesianVertices:
        x.append(v[0])
        y.append(v[1])
    # make a line between last point and first point
    x.append(cartesianVertices[0][0])
    y.append(cartesianVertices[0][1])
    #x.append( 0 )
    #y.append( 0 )
    mp.plot( x , y , colorString + '--' )
    mp.plot( x , y , colorString + 'o' )
    mp.plot( 0 , 0 )
    mp.axis(getPolyDomain(cartesianVertices))
    #print(vertices)
    mp.show()

def plotPolyInArea( cartesianVertices, colorString , lowBound_x , upBound_x , lowBound_y , upBound_y ):
    x=[]
    y=[]

    for v in cartesianVertices:
        x.append(v[0])
        y.append(v[1])
    # make a line between last point and first point
    x.append(cartesianVertices[0][0])
    y.append(cartesianVertices[0][1])
    #x.append( 0 )
    #y.append( 0 )
    mp.plot( x , y , colorString + '--' )
    mp.plot( x , y , colorString + 'o' )
    mp.plot( 0 , 0 )
    mp.axis([ lowBound_x , upBound_x , lowBound_y , upBound_y])
    #print(vertices)
    mp.show()

# because of machine epsilon this wont work if point is on the boundary
def containsPoint( orderedCartesianVertices , point ):
    if ( len( orderedCartesianVertices )  == 1 ):
        if( M.norm( M.subVek( point , orderedCartesianVertices[1] ) ) <= test_eps ):
            return True
        else:
            return False

    if( len(orderedCartesianVertices) == 2 ):

        d = M.subVek(orderedCartesianVertices[0] , orderedCartesianVertices[1] )
        normal = [ d[1] , - d[0] ]

        if( math.fabs( M.scal( point , normal ) - M.scal( orderedCartesianVertices[1] , normal ) ) <= containsPoint_eps ):
            return True
        else:
            return False
    i=0
    n = len( orderedCartesianVertices )
    for v in orderedCartesianVertices:
        w = orderedCartesianVertices[ (i+1) % n ]
        d = M.subVek( w , v )
        normal = [ d[1] , -d[0] ]
        if ( M.scal( normal , point) - containsPoint_eps > M.scal( normal , v) ):
            return False
        i = i + 1
    return True

#P = getRandomPolygon(8)
#plotPoly(P)

# buggy
def isConvex( orderedCartesianVertices ):
    for v in orderedCartesianVertices:
        test = orderedCartesianVertices.copy()
        test.remove(v)
        if( containsPoint( test , v ) ):
            return False

    return True

# maybe worse, because here we have to calculate arctan ( x / y )...

# TODO: check, if eps at bound checks are making condition stronger (desirable for logMin counterexample) or easier to fulfill
def isConvexRun( orderedCartesianVertices ):
    n = len(orderedCartesianVertices)
    for i in range( n ):
        v= orderedCartesianVertices[i % n]
        w= orderedCartesianVertices[ ( i + 1 ) % n]
        x = orderedCartesianVertices[ ( i + 2 ) % n]
        d= M.subVek( w , v )
        e = M.subVek( x , w )
        if( isNotStrictlyLeft( d , e ) ):
            return False
    return True



def isNotStrictlyLeft( u , v_notLeft ):
    alpha = rP.getAngle( u )

    beta = rP.getAngle( v_notLeft )

    if alpha < math.pi:
        if alpha < beta - convexRun_eps and beta < alpha + math.pi - convexRun_eps:
            return False
        else:
            return True
    else:
        if alpha < beta - convexRun_eps or beta < alpha - math.pi - convexRun_eps:
            return False
        else:
            return True

#ConvexRunTester =  [ [ 1 , 1 ] , [ 0 , 2 ] , [ -1 , - 1] , [ 1 , 0 ] ]
#plotPoly( ConvexRunTester )
#print( isConvexRun( ConvexRunTester ) )


# oh shit
# dreieck = [[0.43463668469106065, 4.0069543735500055], [3.355131889073447, 2.1534620225883176], [6.236445733444501, 11.26009871536809]]# maxRadius = 1.003015893755336
# maxRadius = 3.265324153125075
# minRadius = 5.926254308177569
# angle = 1.9198063967206336
# result = [[0.9256533780571976, 0.17646635275518419], [-4.651411054314339, 3.8275672327520813], [-4.813297192669639, -1.9671473960856636], [8.539054868926781, -2.036886189421601]]


def makeTest( n ):
    rounds = 0
    wrongs = []
    counting = [ 0 , 0 , 0 , 0 ]
    while rounds < test_repeats:
        P = rP.getRandomPolygon( n )
        convexTest = isConvex( P )
        convexRunTest = isConvexRun( P )
        if not ( convexTest) and convexRunTest:
            print( 'okay')
            #plotPoly(P, 'b')
            counting[1] += 1
        if convexTest and convexRunTest:
            print('nice')
            #plotPoly(P, 'g')
            counting[0] += 1
        if convexTest and not ( convexRunTest ):
            print( 'strange' )
            #plotPoly(P, 'y')
            counting[2] += 1
        if not ( convexTest ) and not ( convexRunTest ):
            print( 'this is bad' )
            wrongs.append(P)
            plotPoly( P  , 'r')
            counting[3] += 1
        rounds = rounds + 1
    return counting

def makeEdgeNumberTest( repeats , edgeNumber ):
    fail_less = 0
    fail_more = 0
    for i in range(repeats):
        P = rP.getRandomPolygon( edgeNumber )
        if len(P) < edgeNumber:
            print( 'not enough edges' )
            fail_less = fail_less + 1
        if len(P) > edgeNumber:
            print('too many edges')
            fail_more = fail_more + 1
    #print(fail_less)
    #print(fail_more)

def makeVertexDistTest( polygon, epsilon):
    n = len(polygon)
    for i in range(n):
        if M.dist( polygon[ i % n] ,  polygon[i-1] ) >= epsilon:
            return False
    return True

def getMinVertexDist( polygon):
    result = math.inf
    n = len(polygon)
    for i in range(n):
        value_1 = math.fabs(polygon[ i % n ][0] -  polygon[ i - 1 ][0] )
        value_2 = math.fabs(polygon[i % n][1] - polygon[i - 1][1])
        if value_1 < result:
            result = value_1
        if value_2 < result:
            result = value_2
    return result

def makeNormalsTest( polygon , eps_norm, eps_direction ):
    n = len( polygon )
    cD = cV.getConeVol( polygon )

    for i in range(n):
        u = [ cD[i % n ][0] , cD[ i % n ][1] ]
        if math.fabs( M.norm(u) - 1) > eps_norm:
            return False

        v = polygon[ i - 1 ]
        w = polygon[ i % n ]

        diff = math.fabs( M.scal( v , u ) - M.scal( w , u ) )

        if diff > eps_direction:
            return False

    return True

def getWorsteNormalsDirection( polygon ):
    n = len( polygon )
    result = 0
    cD = cV.getConeVol( polygon )

    for i in range(n):
        u = [ cD[ i % n ][0] , cD[ i % n ][1] ]

        v = polygon[ i - 1 ]
        w = polygon[ i % n ]

        diff = math.fabs( M.scal( v , u ) - M.scal( w , u ) )

        if diff > result:
            result = diff

    return result

dist = 1.0 / 2.0
P_test = [ [ 2 , 5 ] , [ 0 , 5 ] , [ 0 , -3  ] , [ 2 , -3 ] ]

rP.makeCentered( P_test )


repeats = 10000
eps_normTest = 0.00000001
eps_directionTest = 0.00002


# for eps_directionTest = 0.00002
def getWorstNormDirection( repeats ,  eps_direction ):
    success = 0
    fails = 0
    value_max = 0
    worstPoly = []
    for n in range(repeats):
        K = rP.getRandomPolygon( 4 )
        value = getWorsteNormalsDirection(K)
        if value < eps_directionTest:
            success += 1
        else:
            if value > value_max:
                value_max = value
                worstPoly = K
            minVertDist = getMinVertexDist( K )
            if minVertDist > 1:
                print( minVertDist )

    print( [ success , fails ])
    plotPoly( worstPoly , 'r')
    print( value_max )

#getWorstNormDirection( 1000000 , eps_directionTest )




def supportTest( repeats):
    tooSmall = 0
    tooSmall2 = 0
    tooSmall3 = 3
    tooBig = 0
    easyFail = 0
    n = 0
    while( n < repeats):

        P = rP.getRandomPolygon(3)
        Q = rP.getRandomPolygon(4)
        cD_P = cV.getConeVol(P)

        for i in range(len(P)):
            u = [ cD_P[i-1][0] , cD_P[i-1][1] ]
            alpha = rP.supportFunctionCartesianCentered( Q , u)
            beta = rP.supportFunctionCartesianCentered2(Q, u)
            gamma = trySupport( Q , u , 0.0000001)

            if math.fabs( M.scal( u , P[ i - 1 ]) - rP.supportFunctionCartesianCentered2( P , u) ) > 1 :
                easyFail += 1

            if not ( containsPoint( Q , M.scaleVek( alpha - math.exp( -15) , u ) )):
                tooBig += 1

            if ( containsPoint(Q, M.scaleVek(alpha + 0.01, u))):
                #plotPoly( P , 'r' )
                #plotPoly( Q , 'y')

                tooSmall += 1

            if ( containsPoint(Q, M.scaleVek(beta + 0.01, u))):
                #plotPoly( P , 'r' )
                #plotPoly( Q , 'y')

                tooSmall2 += 1

            if ( containsPoint(Q, M.scaleVek(gamma + 0.01, u))):
                #plotPoly( P , 'r' )
                #plotPoly( Q , 'y')

                tooSmall3 += 1
        n += 1

    print( [ easyFail , tooBig , tooSmall , tooSmall2 , tooSmall3 ])

def supportRationalTest ( repeats):
    tooSmall = 0
    tooSmall2 = 0
    tooBig = 0
    easyFail = 0
    n = 0
    while( n < repeats):

        P = rRP.getRandomRationalPolygon(3 , 3)
        Q = rRP.getRandomRationalPolygon(4 , 3)
        #rP.supportFunctionCartesianCenteredRational()
        for i in range(len(P)):
            u = [ P[i % len(P)][1] - P[i - 1][1]  , P[ i % len(P)][0] - P[ i - 1][0] ]
            alpha = rP.supportFunctionCartesianCenteredRational( Q , u)
            beta = rP.supportFunctionCartesianCenteredRational(Q, u)

            #pos_diff = M.scal( u , P[ i - 1 ]) - rP.supportFunctionCartesianCentered2( P , u)

            #if   pos_diff > 0.1  or - pos_diff > 0.1:
            #   easyFail += 1

            if not ( containsPoint( Q , M.scaleVek( alpha - f.Fraction( 1, 10) , u ) )):
                tooBig += 1

            if ( containsPoint(Q, M.scaleVek(alpha + f.Fraction( 1, 1000), u))):
                #plotPoly( P , 'r' )
                #plotPoly( Q , 'y')

                tooSmall += 1

            if ( containsPoint(Q, M.scaleVek(beta + f.Fraction( 1, 1000), u))):
                #plotPoly( P , 'r' )
                #plotPoly( Q , 'y')

                tooSmall2 += 1
        n += 1

    print( [ easyFail , tooBig , tooSmall , tooSmall2 ])

def trySupport( Q , u , eps):
    alternative = 1
    exp = 0
    if containsPoint(Q, u):
        while containsPoint(Q, M.scaleVek(alternative + eps, u)) or not containsPoint(Q,
                                                                                      M.scaleVek(alternative - eps, u)):
            if containsPoint(Q, M.scaleVek(alternative + 10 ** exp, u)):
                alternative += 10 ** exp
            else:
                exp -= 1

    else:
        while containsPoint(Q, M.scaleVek(alternative + eps, u)) or not containsPoint(Q,
                                                                                      M.scaleVek(alternative - eps, u)):

            if not containsPoint(Q, M.scaleVek(alternative + 10 ** exp, u)):
                alternative -= 10 ** exp
            else:
                exp -= 1
    return alternative


supportTest(1)


def compareSupportFunctions( repeats , eps ):
    P = rP.getRandomPolygon(3)
    Q = rP.getRandomPolygon(4)
    cD_P = cV.getConeVol(P)

    for i in range(len(P)):
        u = [cD_P[i - 1][0], cD_P[i - 1][1]]
        alpha = rP.supportFunctionCartesianCentered(Q, u)
        beta = rP.supportFunctionCartesianCentered2(Q, u)
        alternative = trySupport( Q , u )




        print(alternative)
        print( alpha)
        print( beta )










test_repeats= 10000
test_verticeNumber = 7
#makeEdgeNumberTest( test_repeats , test_verticeNumber )
#print(makeTest(5))

test = [ [ 8 , 1 ] , [ 8 , 8 ] , [ 0 , 2 ] , [ -2 , -2 ] ]


dist = math.sqrt(2)
pi = math.pi

Quadrat_polar = [ [ pi / 4 , dist ] , [ 3*pi/4 , dist] , [ 5*pi/4 , dist ] , [  7*pi / 4 , dist ] ]
Diamond_polar = [ [ pi / 2 , 1 ] , [ pi , 1 ] , [ 3 * pi / 2 , 1 ] , [ 0 , 1 ] ]
Quadrat = rP.getCartesian( Quadrat_polar )
Diamond = rP.getCartesian( Diamond_polar )

u = [ math.sqrt(0.5) , math.sqrt(0.5) ]

if ( math.fabs( rP.supportFunctionPolar( Diamond_polar , u ) - math.sqrt(0.5) ) > 0.00000001  ):
    print( ' Fehler bei support function test ')
if ( math.fabs( rP.supportFunctionPolar( Quadrat_polar , u ) - math.sqrt(2) ) > 0.00000001  ):
    print( ' Fehler bei support function test ')
alpha = rP.supportFunctionCartesianCentered( Diamond, u )
if math.fabs( alpha - math.sqrt(0.5) ) > 0.00000001  :
    print( ' Fehler bei supportCartesianCentered function test ')
    print(alpha)
beta = rP.supportFunctionCartesianCentered( Quadrat , u )
if math.fabs( beta - math.sqrt(2) ) > 0.00000001  :
    print( ' Fehler bei supportCartesianCentered function test ')
    print(beta)


Trapez = [  [ 1 , 1] , [ 0 , 1] , [-3 , - 1] , [ 0 , -1 ] ]

a = Trapez[3][0] - Trapez[2][0]
b = Trapez[0][0] - Trapez[1][0]
h = Trapez[1][1] - Trapez[2][1]
e = Trapez[0][0] - Trapez[2][0]

center = rP.getCenter(Trapez)
# center form wikipedia formula for center of trapez
center2 = [ (a**2 - b**2 + e * ( a + 2.0 * b ) ) / ( 3.0 * ( a + b ) ) + Trapez[2][0] , h * ( a + 2.0 * b ) / ( 3.0 * ( a + b ) ) + Trapez[2][1]]
if M.dist( center, center2) > 0.00001:
    print('fehler bei getCenter for trapez')
    print( center )
    print( center2 )
    print(a)
    print(b)
    print(h)
    print(e)
    plotPoly( Trapez , 'r')

#rP.makeCentered( test )

#plotPoly( test , 'b')

#support_test = [ [ 0 , 1] , [ 0.5 * math.pi , 1 ] , [ 1.5 * math.pi , 1 ] ]

#u = [ 1 , 1 ]

#print( rP.supportFunction( support_test , u ) - 0.5 )

repeats = 1000

n = 0
good_triangles = 0
bad_triangles = 0
while( n < repeats ):
    P = rP.getRandomPolygon(3)

    if isConvex( P ) and isConvexRun( P ):
        good_triangles += 1
    else:
        bad_triangles += 1

    n += 1
#print( [ good_triangles , bad_triangles ])
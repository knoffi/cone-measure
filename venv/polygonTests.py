from cmath import sqrt
import math
import numpy
import random
import coneVol as cV
import matrix as M
import matplotlib.pyplot as mp
import randomPolygon3Centered as rP
import machEps as mE

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
        if( M.norm( M.subVek( point , orderedCartesianVertices ) ) <= test_eps ):
            return True
        else:
            return False

    if( len(orderedCartesianVertices) == 2 ):
        d = M.subVek(orderedCartesianVertices[0] , orderedCartesianVertices[1] )
        normal = [ d[1] , - d[0] ]
        if( abs( M.scal( point , normal ) - M.scal( orderedCartesianVertices[0] , normal ) <= test_eps) ):
            return True
        else:
            return False
    i=0
    n = len( orderedCartesianVertices )
    for v in orderedCartesianVertices:
        w = orderedCartesianVertices[ (i+1) % n ]
        d = M.subVek( w , v )
        normal = [ d[1] , -d[0] ]
        if ( M.scal( normal , point) + test_eps > M.scal( normal , v) ):
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
def isConvexRun( orderedCartesianVertices ):
    n = len(orderedCartesianVertices)
    for i in range( n ):
        v= orderedCartesianVertices[i % n]
        w= orderedCartesianVertices[ ( i + 1 ) % n]
        x = orderedCartesianVertices[ ( i + 2 ) % n]
        d= M.subVek( w , v )
        e = M.subVek( x , w )
        if( isNotStrictlyLeft( d , e) ):
            return False
    return True



def isNotStrictlyLeft( u ,v_notLeft ):
    alpha = rP.getAngle( u )

    beta = rP.getAngle( v_notLeft )

    if alpha < math.pi:
        if alpha < beta + convexRun_eps and beta < alpha + math.pi + convexRun_eps:
            return False
        else:
            return True
    else:
        if alpha < beta + convexRun_eps or beta < alpha - math.pi + convexRun_eps:
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


dist = 1.0 / 2.0
P_test = [ [ 2 , 5 ] , [ 0 , 5 ] , [ 0 , -3  ] , [ 2 , -3 ] ]

rP.makeCentered( P_test )






# eventuell wird ein tester besser, wenn ich die polarkoordinaten nicht komplett umwandle
# HIER SOLLTE ALS NÄCHSTES GEARBEITET WERDEN; UM TESTER ZU VERBESSERN ::: WELCHE SIND DIE SINNVOLLEN EPSILONS FÜR DIE TEST-SCHRANKENß SOLLTE MAN SCHNELL MAL MACHINE EPSILON ZUVOR BERECHNEN? UND KONDITIONIEREUNG BESTIMMEN?
# warum kriege ich so oft  benachbarte vertices mit fast gleichem Winkel (Zufall?) und (unnötig) fast gleichem Radius (Fehler bei minRadius, maxRadius) ?
test_eps = mE.getMachEps()
convexRun_eps = test_eps
test_repeats= 10000
test_verticeNumber = 7
#makeEdgeNumberTest( test_repeats , test_verticeNumber )
#print(makeTest(5))

test = [ [ 8 , 1 ] , [ 8 , 8 ] , [ 0 , 2 ] , [ -2 , -2 ] ]

#rP.makeCentered( test )

#plotPoly( test , 'b')

#support_test = [ [ 0 , 1] , [ 0.5 * math.pi , 1 ] , [ 1.5 * math.pi , 1 ] ]

#u = [ 1 , 1 ]

#print( rP.supportFunction( support_test , u ) - 0.5 )
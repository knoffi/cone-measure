import randomPolygon3Centered as rP
import matrix as M
import math
import fractions as f

convexRun_eps = 0.01
containsPoint_eps = 0.01

def getRandomRoundedPolygon( number, digits ):
    K = []
    while len(K) == 0 or not isConvexRun(K) or not isConvex(K):

        K = rP.getRandomPolygon( number )

        rP.roundVertices(K, digits)

    return K

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

def getRandomRationalPolygon( number , digits):
    K = getRandomRoundedPolygon( number, digits)
    result = getRationalWithDigits( K , digits )
    return result

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
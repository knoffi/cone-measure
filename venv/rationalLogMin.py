import randomPolygon3Centered as rP
import coneVol2 as cV
import random as r
import polygonTests as pT
import matrix as M
import math as math
import random as rd



def logMinTest( orderedCartPolygon_1 , orderedCartPolygon_2 , bringInPosition , eps_prod ):
    K = orderedCartPolygon_1
    L = orderedCartPolygon_2
    bringInPosition(K)
    bringInPosition(L)
    rP.makeUnitVolume(K)
    rP.makeUnitVolume(L)
    # not sure if makeUnitVolume destroys the centered property... it is clear that it needs it for the normalization of volume
    bringInPosition(K)
    bringInPosition(L)
    coneVol = cV.getConeVol( K )
    #print( coneVol )

    prod = 1

    for i in range( len( K ) ):
        # Reihenfolge der Normalenvektoren sollte jetzt egal sein
        u = [ coneVol[i-1][0] , coneVol[i-1][1] ]
        h_1 = rP.supportFunctionCartesianCentered( L , u )
        h_2 = rP.supportFunctionCartesianCentered( K , u )
        quotient = h_1 / h_2
        if quotient == 0:
            print( 'quotient von logMin war gleich 0...')
            print( u )
            print( L[ i-1 ])
            print( K[ i-1 ])
        #print( quotient)
        #print( K[ i - 1] )
        #print( u )
        if quotient >= 0:
            prod = prod * math.pow( quotient , coneVol[ i - 1 ][ 2 ] )
        else:
            print( ' quotient von logMin war negativ' )
            print( u )
            print( coneVol[ i - 1][2] )
            pT.plotPoly( K , 'r')
            pT.plotPoly( L , 'g')
    if prod + eps_prod >= 1:
        return True
    else:
        print( 1 - prod)
        return False


def logMinRationalAutoTest( repeats , roundMethod, digits , bringInPositionRational, getPositionRational , eps_positionRational , eps_volumeRational , eps_normalsRational , eps_logMinRational):
    n = 0
    logMinTrue = 0
    logMinFalse = 0
    logMinTriangleFalse = 0
    while n < repeats:
        K = []
        L = []
        # error possible or is len(K) == 0 checked first?
        while len(K) * len(L) == 0 or not pT.isConvex(K) or not pT.isConvex(L) or not pT.isConvexRun(K) or not pT.isConvexRun(L):
            P = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
            Q = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
            K = rP.getCartesian( P )
            L = rP.getCartesian( Q )
            roundMethod( K , digits )
            roundMethod( L , digits )

        M = rP.getRationalWithDigits(K , digits )
        N = rP.getRationalWithDigits(L , digits )

        if not logMinTestRational( M , N , bringInPositionRational , eps_logMinRational ):
            if M.norm( getPositionRational(M)) +  M.norm( getPositionRational(N) ) < eps_positionRational:
                if math.fabs( rP.getVolumeRational(M) - 1 ) + math.fabs( rP.getVolumeRational(N) - 1 ) < eps_volumeRational:
                    # which bound is sufficient for a stable logMin calculation? a bound of 0.00003 seems to be always fulfilled...
                    if pT.getWorsteNormalsDirectionRational(M) < eps_normals:
                        if len(K) == 3 or len(L) == 3:
                            #print( 'das darf nicht sein')
                            #print(K)
                            #print(L)
                            logMinTriangleFalse += 1
                        else:
                            print( 'logMin Gegenbeispiel' )
                            print( M )
                            print( N )
                            pT.plotPoly(K, 'r')
                            pT.plotPoly(L, 'g')
                            logMinFalse += 1
        else:
            logMinTrue += 1

        n = n+1
    print('done')
    print( [ logMinTrue , logMinFalse , logMinTriangleFalse ] )


def getAveragePolygonConditions( digits , eps  ):
    pT.containsPoint_eps = eps
    pT.convexRun_eps = eps
    repeats = 0
    result = 0
    while( repeats < 10000):

        n = 0

        K = []
        L = []
        # error possible or is len(K) == 0 checked first?
        while  len(K) * len(L) == 0 or not pT.isConvex(K) or not pT.isConvex(L) or not pT.isConvexRun(K) or not pT.isConvexRun(L):
            P = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
            Q = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
            K = rP.getCartesian( P )
            L = rP.getCartesian( Q )
            rP.roundVertices( K , digits )
            rP.roundVertices( L , digits )
            n += 1
        result += n
        repeats += 1
        if repeats % 1000 == 0:
            print(repeats)
            pT.plotPoly( K ,'r')
            pT.plotPoly( L , 'g')

    print( result / repeats)

#getAveragePolygonConditions( 1 , 0.05 )

result_for_2digits_0Point1AsEps = 2.68502
result_for_1digits_0Point1AsEps = 2.6558

result_for_1digits_0Point01AsEps = 1.29
result_for_1digits_0Point05AsEps = 1.7934



logMinRationalAutoTest( 1000 , rP.roundVertices , 1 ,  rP.makeCentered , rP.getCenter , 0.000000001 , 0.000000001 , 0.000000001, 0.5)
#logMinAutoTest( 1000 , rP.makeBaryCentered , rP.getBaryCenter , 0.000000001 , 0.000000001 , 0.5)


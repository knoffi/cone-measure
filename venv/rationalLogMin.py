import randomPolygon3Centered as rP
import coneVol2 as cV
import random as r
import polygonTests as pT
import matrix as M
import math as math
import random as rd
import fractions as f


def logMinTestRational( orderedCartPolygonRational_1 , orderedCartPolygonRational_2 , eps_prod ):
    A = orderedCartPolygonRational_1
    B = orderedCartPolygonRational_2
    # is volume formula false, if polytop is not centered?
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)
    rP.translateRational( A , info_A[0] )
    rP.translateRational( B , info_B[0] )
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)

    # does this work? A has rational points...
    coneVol = cV.getPseudoConeVolRational( A )


    # the -1 is necessary because python calculates ln( 1 + x )
    # dont worry, this meanas that exactly ln( (x - 1 ) + 1 ) = ln( x ) will be calculated, as long as x-1 > -1
    logMinBound = math.exp( info_A[1] * f.Fraction( 1 , 2 ) * math.log1p( info_B[1] / info_A[1]  - 1) )

    prod = 1

    for i in range( len( A ) ):
        # Reihenfolge der Normalenvektoren sollte jetzt egal sein
        u = [ coneVol[i-1][0] , coneVol[i-1][1] ]
        h_1 = rP.supportFunctionCartesianCenteredRational( B , u )
        h_2 = rP.supportFunctionCartesianCenteredRational( A , u )
        quotient = h_1 / h_2
        if quotient == 0:
            print( 'quotient von logMin war gleich 0...')
            print( u )
            print( A[ i-1 ])
            print( B[ i-1 ])
        #print( quotient)
        #print( K[ i - 1] )
        #print( u )
        if quotient >= 0:
            #print('maybe I can improve this, here starts the inexactness...')
            prod = prod * math.pow( quotient , coneVol[ i - 1 ][ 2 ] )
        else:
            print( ' quotient von logMin war negativ' )
            print( u )
            print( coneVol[ i - 1][2] )
            #pT.plotPoly( A , 'r')
            #pT.plotPoly( B , 'g')
    if prod + eps_prod >= logMinBound:
        return True
    else:
        print( 'log min is false')
        print( logMinBound - prod)
        return False

def logMinTestRational2(orderedCartPolygonRational_1, orderedCartPolygonRational_2, eps_prod):
    A = orderedCartPolygonRational_1
    B = orderedCartPolygonRational_2
    # is volume formula false, if polytop is not centered?
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)
    rP.translateRational(A, info_A[0])
    rP.translateRational(B, info_B[0])
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)

    # does this work? A has rational points...
    coneVol = cV.getPseudoConeVolRational(A)

    # the -1 is necessary because python calculates ln( 1 + x )
    # dont worry, this meanas that exactly ln( (x - 1 ) + 1 ) = ln( x ) will be calculated, as long as x-1 > -1
    logMinBound =  math.pow( info_B[1] / info_A[1] , info_A[1] * f.Fraction(1, 2) )

    prod = 1

    for i in range(len(A)):
        # Reihenfolge der Normalenvektoren sollte jetzt egal sein
        u = [coneVol[i - 1][0], coneVol[i - 1][1]]
        h_1 = rP.supportFunctionCartesianCenteredRational(B, u)
        h_2 = rP.supportFunctionCartesianCenteredRational(A, u)
        quotient = h_1 / h_2
        if quotient == 0:
            print('quotient von logMin war gleich 0...')
            print(u)
            print(A[i - 1])
            print(B[i - 1])
        # print( quotient)
        # print( K[ i - 1] )
        # print( u )
        if quotient >= 0:
            # print('maybe I can improve this, here starts the inexactness...')
            prod = prod * math.pow(quotient, coneVol[i - 1][2])
        else:
            print(' quotient von logMin war negativ')
            print(u)
            print(coneVol[i - 1][2])
            # pT.plotPoly( A , 'r')
            # pT.plotPoly( B , 'g')
    if prod + eps_prod >= logMinBound:
        return True
    else:
       # print('log min is false')
        #print(logMinBound - prod)
        return False


def logMinRationalAutoTest( repeats , roundMethod, digits ,  eps_positionRational , eps_volumeRational , eps_normalsRational , eps_logMinRational):
    n = 0
    logMinTrue = 0
    logMinFalse = 0
    logMinTriangleFalse = 0
    centeredFails = 0
    while n < repeats:
        K = []
        L = []
        # error possible or is len(K) == 0 checked first?
        while len(K) * len(L) == 0  or not pT.isConvexRun(K) or not pT.isConvexRun(L) or not pT.isConvex(K) or not pT.isConvex(L):
            P = rP.getRandomNoncenteredPolarPolygon( math.floor( 3 ) )
            Q = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
            K = rP.getCartesian( P )
            L = rP.getCartesian( Q )
            roundMethod( K , digits )
            roundMethod( L , digits )

        R = rP.getRationalWithDigits(K , digits )
        S = rP.getRationalWithDigits(L , digits )
        rP.makeCentered( K )
        rP.makeCentered( L )

        info_R = rP.getCenterAndVolumeRational(R)
        info_S = rP.getCenterAndVolumeRational(S)
        rP.translateRational(R, info_R[0])
        rP.translateRational(S, info_S[0])
        info_R = rP.getCenterAndVolumeRational(R)
        info_S = rP.getCenterAndVolumeRational(S)

        if not logMinTestRational2( R , S , eps_logMinRational) and not logMinTestRational( R , S , eps_logMinRational):
            centeredNorm = info_R[0][0] * info_R[0][0] + info_R[0][1] * info_R[0][1] + info_S[0][0] * info_S[0][0] + info_S[0][1] * info_S[0][1] * 1.0
            print( centeredNorm )
            if centeredNorm  < eps_positionRational:
                # I am not normalizing the volumes anymore.
                #if math.fabs( rP.getVolumeRational(M) - 1 ) + math.fabs( rP.getVolumeRational(N) - 1 ) < eps_volumeRational:
                if len(K) == 3 or len(L) == 3:
                    #print( 'das darf nicht sein')
                    #rP.printRationalPolygon(R)
                    #rP.printRationalPolygon(S)
                    #pT.plotPoly(K, 'r')
                    #pT.plotPoly(L, 'g')
                    logMinTriangleFalse += 1
                else:
                    print( 'logMin Gegenbeispiel' )
                    rP.printRationalPolygon( R )
                    rP.printRationalPolygon( S )
                    pT.plotPoly(K, 'r')
                    pT.plotPoly(L, 'g')
                    logMinFalse += 1
            else:
                centeredFails += 1
        else:
            logMinTrue += 1

        n = n+1
    print('done')
    print( [ logMinTrue , logMinFalse , centeredFails, logMinTriangleFalse  ] )


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


logMinRationalAutoTest( 1000 , rP.roundVertices , 1 , 0.000000001 , 0.000000001 , 0.000000001, 0.5)



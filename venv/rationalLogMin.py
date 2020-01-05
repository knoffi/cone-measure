import randomPolygon3Centered as rP
import coneVol2 as cV
import random as r
import polygonTests as pT
import matrix as M
import math as math
import random as rd
import fractions as f
import randomRationalPolygon as rRP


# this method should center both polygons and return array with elements [volume_A , volume_B , array_of_pairs[ quotient, volume_of_cone ]  ]
def logMinTestRationalPreparation(orderedCartPolygonRational_1, orderedCartPolygonRational_2):
    A = orderedCartPolygonRational_1
    B = orderedCartPolygonRational_2
    # is volume formula false, if polytop is not centered?
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)
    rP.translateRational(A, info_A[0])
    rP.translateRational(B, info_B[0])
    info_A = rP.getCenterAndVolumeRational(A)
    info_B = rP.getCenterAndVolumeRational(B)
    logMinValues = []

    cD_A = cV.getPseudoConeVolRational(A)

    for data in cD_A:
        u = [ data[0] , data[1] ]
        h_1 = rP.supportFunctionCartesianCenteredRational(B, u)
        h_2 = rP.supportFunctionCartesianCenteredRational(A, u)
        quotient = h_1 / h_2
        logMinValues.append( [ quotient , data[2] ])

    return [ info_A[1] , info_B[1] , logMinValues]



def logMinTestRational1( logMinValues , eps_prod ):
    volume_A = logMinValues[0]
    volume_B = logMinValues[1]

    # does this work? A has rational points...
    logMinBound = math.exp( volume_A * f.Fraction( 1 , 2 ) * math.log( volume_B / volume_A ) )

    prod = 1

    for i in range( len( logMinValues[2] ) ):
        # Reihenfolge der Normalenvektoren sollte jetzt egal sein
        quotient = logMinValues[2][i-1][0]
        if quotient == 0:
            print( 'quotient von logMin war gleich 0...')
        if quotient > 0:
            volume_of_cone = logMinValues[2][i-1][1]
            prod = prod * math.pow( quotient , volume_of_cone )
        else:
            print( ' quotient von logMin war negativ' )
    if prod + eps_prod >= logMinBound:
        return True
    else:
        #print( 'log min is false')
        #print( logMinBound - prod)
        return False




def logMinTestRational2(logMinValues , eps_prod):
    volume_A = logMinValues[0]
    volume_B = logMinValues[1]

    logMinBound =  math.pow( volume_B / volume_A , volume_A * f.Fraction(1, 2) )

    prod = 1

    for i in range( len(logMinValues[2]) ):
        quotient = logMinValues[2][i-1][0]
        if quotient == 0:
            print('quotient von logMin war gleich 0...')
        if quotient >= 0:
            volume_of_cone = logMinValues[2][i-1][1]
            # print('maybe I can improve this, here starts the inexactness...')
            prod = prod * math.pow(quotient, volume_of_cone)
        else:
            print(' quotient von logMin war negativ')
    if prod + eps_prod >= logMinBound:
        return True
    else:
        # print('log min is false')
        #print(logMinBound - prod)
        return False


def logMinTestRationalNormalizing( logMinValues , eps_prod):
    volume_A = logMinValues[0]
    volume_B = logMinValues[1]
    alpha = 1 / math.sqrt(logMinValues[0])
    beta = 1 / math.sqrt(logMinValues[1])

    logMinBound = 1

    prod = 1

    for i in range(len( logMinValues[2])):

        quotient = logMinValues[2][i - 1][0] * beta / alpha
        if quotient == 0:
            print('quotient von logMin war gleich 0...')
        if quotient >= 0:
            volume_of_cone = logMinValues[2][i - 1][1]
            # print('maybe I can improve this, here starts the inexactness...')
            prod = prod * math.pow(quotient, alpha * volume_of_cone)
        else:
            print(' quotient von logMin war negativ')
    if prod + eps_prod >= logMinBound:
        return True
    else:
        # print('log min is false')
        print(prod)
        #print( logMinValues )
        return False


def logMinRationalAutoTest( repeats , digits ,  eps_logMinRational):
    n = 0
    logMinTrue = 0
    logMinFalse = 0
    logMinTriangleFalse = 0
    centeredFails = 0
    while n < repeats:

        K = rRP.getRandomRoundedPolygon( 5 , digits)
        L = rRP.getRandomRoundedPolygon( 6 , digits)

        R = rRP.getRationalWithDigits( K , digits )
        S = rRP.getRationalWithDigits(L, digits)

        info_R = rP.getCenterAndVolumeRational(R)
        info_S = rP.getCenterAndVolumeRational(S)
        rP.translateRational(R, info_R[0])
        rP.translateRational(S, info_S[0])
        info_R = rP.getCenterAndVolumeRational(R)
        info_S = rP.getCenterAndVolumeRational(S)

        logMinValues = logMinTestRationalPreparation( R , S )
        while(True):
            try:
                if not logMinTestRational2( logMinValues , eps_logMinRational) and not logMinTestRational1( logMinValues , eps_logMinRational) and not logMinTestRationalNormalizing( logMinValues , eps_logMinRational ):
                    #centeredNorm = info_R[0][0] * info_R[0][0] + info_R[0][1] * info_R[0][1] + info_S[0][0] * info_S[0][0] + info_S[0][1] * info_S[0][1]
                    #print( centeredNorm )
                    if len(R) == 3 or len(S) == 3:
                        #print( 'das darf nicht sein')
                        rP.printRationalPolygon(R)
                        rP.printRationalPolygon(S)
                        pT.plotPoly(K, 'r')
                        pT.plotPoly(L, 'g')
                        logMinTriangleFalse += 1
                        break
                    else:
                        print( 'logMin Gegenbeispiel' )
                        rP.printRationalPolygon( R )
                        rP.printRationalPolygon( S )
                        pT.plotPoly(R, 'r')
                        pT.plotPoly(S, 'g')
                        logMinFalse += 1
                else:
                    logMinTrue += 1

                n = n+1
                break
            except OverflowError or ZeroDivisionError:
                break
    print('done')
    print( [ logMinTrue , logMinFalse , centeredFails, logMinTriangleFalse  ] )

#logMinRationalAutoTest( 1000 , 3 , 0)


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


#def testSupportFunctions( K , L):


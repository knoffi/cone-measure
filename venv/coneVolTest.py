import coneVol2 as cV
import matrix as M
import randomPolygon3Centered as rP
import polygonTests as pT
import math
import numpy
import random

def coneVolTest( polygon , coneVolume , eps):
    diff = 0
    coneData = cV.getConeVol( polygon )
    for i in range(len(polygon)):
        diff += M.dist( coneData[ i - 1 ] , coneVolume[ i - 1])**2
    if diff >= eps:
        print( 'coneVolTest failed with polygon')
        print( polygon )

def getSumOuterNormal( polygon ):
    coneVol = cV.getConeVol( polygon )
    normSum=0
    for v in coneVol:
        n = [ v[0] , v[1] ]
        normSum = normSum + M.norm( n )
    return normSum



#testPolar_cV= [[0.12052685145948405, 0.10016898087234995], [-1.601236023721101, 1.2234795179795563], [-0.007010810061268179, -0.6251879546153744], [0.8160675487328953, -0.4299405664378391], [0.6716524335899903, -0.26851997779869285]]
#test_cV = [[0.12052685145948405, 0.10016898087234995], [-1.601236023721101, 1.2234795179795563], [-0.007010810061268179, -0.6251879546153744], [0.8160675487328953, -0.4299405664378391], [0.6716524335899903, -0.26851997779869285]]

# rP.getCartesian( testPolar_cV )
#[ [ 1 ,0 ] , [ 0 , 1 ] , [ 0 , -1 ] ]
# print( coneVolTest( test_cV , [ 0.5 , 0 , 0.5 ] ) )

def normalsTest( repeats , verticesNumber, test_eps):
    for i in range(repeats):
        P = rP.getRandomPolygon(verticesNumber)
        coneVol = cV.getConeVol( P )
        diff = math.fabs( getSumOuterNormal(P) - verticesNumber )
        if diff > test_eps :
            #print( P )
            #print( diff )
            normalNormsPrint(coneVol)
            print(' normalsTest failed at')
            print(len(coneVol))
            print(len(P))
            pT.plotPoly( P , 'r')
            break

def normalNormsPrint( coneVolume ):
    for data in coneVolume:
        print( M.norm( [ data[0] , data[1] ] ) )


def polyFunctionalTest( repeats , edgeNumber, eps ):
    turns = 0
    while turns < repeats:
        P = rP.getRandomPolygon( edgeNumber )
        cD = cV.getConeVol( P )
        vertex = P[0]
        value_vertex = cV.polyFunctional( cD , vertex )
        grad_vertex = cV.gradPolyFunctional( cD , vertex )

        # getConeVol, polyFunctional oder gradPolyFunctional könnte falsch sein...

        if( value_vertex >= eps):
            print( ' Fehler bei polyFunctionalTest , value zu hoch')
            print( P )
            print( value_vertex )

        if( M.norm(grad_vertex) >= eps):
            print(' Fehler bei polyFunctionalTest , grad zu hoch')
            print( P )
            print( grad_vertex )
        turns  += 1

def hMethodPolyFunctionalCheck( cD , point , eps , delta ):

    while( True ):
        try:
            value = cV.polyFunctional( cD , point )
            break

        except ZeroDivisionError:
            return True
            break

    e_1 = [ 1 , 0 ]
    e_2 = [ 0 , 1 ]
    grad_point = cV.gradPolyFunctional( cD , point )

    # diff is a row, grad is a column...
    diff = M.Matrix( [ grad_point ])
    diff_1 = diff.image( e_1  )[0]
    diff_2 = diff.image( e_2  )[0]

    pointTrans_1 = M.addVek( point , M.scaleVek( delta , e_1 ) )
    pointTrans_2 = M.addVek( point , M.scaleVek( delta , e_2 ) )

    diffQuotient_1 = cV.polyFunctional( cD , point ) - cV.polyFunctional( cD , pointTrans_1)
    diffQuotient_1 = delta * diffQuotient_1

    diffQuotient_2 = cV.polyFunctional( cD , point ) - cV.polyFunctional( cD , pointTrans_2)
    diffQuotient_2 = delta * diffQuotient_2

    dist_1 = math.fabs(diff_1 - diffQuotient_1 )
    dist_2 = math.fabs(diff_2 - diffQuotient_2 )

    if dist_1 > eps or dist_2 > eps:
        print(' difference equals ')
        print( dist_1 )
        print( dist_2 )
        return False

    return True

def gradPolyFunctionalTest( repeats ,  edgeNumber ,  eps , delta ):

    n = 0

    while n <= repeats:

        P = rP.getRandomPolygon( edgeNumber )

        cD_test = cV.getConeVol( P )
        point_test = [ random.random() * 10 , random.random() * 10 ]

        # could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...

        check = hMethodPolyFunctionalCheck( cD_test , point_test , eps , delta )

        if not check:
            print( ' gradPolyFunctionalTest failed with polytop and point_test : ')
            print( P )
            print( point_test )

        n = n + 1

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...
gradPolyFunctionalTest( 3 , 5 , 0.1 , 0.000001 )
#polyFunctionalTest( 20 , 5 , 0.000001 )

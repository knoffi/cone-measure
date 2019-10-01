import coneVol2 as cV
import matrix as M
import randomPolygon3Centered as rP
import polygonTests as pT
import math
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

polyFunctionalTest( 20 , 5 , 0.000001 )

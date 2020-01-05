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

coneVolTest( [[1,1] , [ -1 , 1] , [ -1 , -1 ] , [1 , -1]] , [ [ 1 , 0 , 1 ] , [ 0 , 1 , 1] , [-1 , 0 , 1 ] , [ 0 , -1 , 1 ]] , 0.1)

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

def hMethodPolyFunctionalCheck( cD , point , eps , h ):

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

    pointTrans_1 = M.addVek( point , M.scaleVek( h , e_1 ) )
    pointTrans_2 = M.addVek( point , M.scaleVek( h , e_2 ) )

    diffQuotient_1 = cV.polyFunctional( cD , pointTrans_1) / h
    substraction_1 = cV.polyFunctional( cD , point ) / h
    diffQuotient_1 = diffQuotient_1 - substraction_1

    diffQuotient_2 = cV.polyFunctional(cD, pointTrans_2) / h
    substraction_2 = cV.polyFunctional(cD, point) / h
    diffQuotient_2 = diffQuotient_2 - substraction_2

    dist_1 = math.fabs(diff_1 - diffQuotient_1 )
    dist_2 = math.fabs(diff_2 - diffQuotient_2 )

    if dist_1 > eps or dist_2 > eps:
        print(' difference equals ')
        print( dist_1 )
        print( dist_2 )
        print(' calculated diff and approximated diff:')
        diff.list()
        print( [diffQuotient_1 , diffQuotient_2])
        return False

    return True

def gradPolyFunctionalTest_random( repeats ,  edgeNumber ,  eps , h , delta , nearby_option):

    n = 1
    if not nearby_option:
        while n <= repeats:

            P = rP.getRandomPolygon( edgeNumber )

            cD_test = cV.getConeVol( P )
            point_test = [ random.random() * 10 , random.random() * 10 ]

            # could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...

            check = hMethodPolyFunctionalCheck( cD_test , point_test , eps , h )

            if not check:
                print(' gradPolyFunctionalTest failed with polytop , point_test , value : ')
                print(P)
                print(point_test)
                print(cV.polyFunctional(cD_test, point_test))

            n = n + 1
    else:
        while n <= repeats:

            P = rP.getRandomPolygon(edgeNumber)

            cD_test = cV.getConeVol(P)
            point_test = [random.random() , random.random() ]
            norm = M.norm(point_test)
            point_test = M.scaleVek( delta / norm , point_test )
            print(' here is the translation')
            print(point_test)
            point_test = M.addVek( point_test , P[0])
            # could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...

            check = hMethodPolyFunctionalCheck(cD_test, point_test, eps, h)

            if not check:
                print(' gradPolyFunctionalTest failed with polytop , point_test , value : ')
                print(P)
                print(point_test)
                print(cV.polyFunctional(cD_test, point_test))

            n = n + 1

def gradPolyFunctionalTest_array( polygon_list , point_list , eps , h ):
    index = 0
    for P in polygon_list:

        cD_test = cV.getConeVol(P)
        point_test = point_list[index]

        # could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...

        check = hMethodPolyFunctionalCheck(cD_test, point_test, eps, h)

        if not check:
            print(' gradPolyFunctionalTest failed with polytop , point_test , value : ')
            print(P)
            print(point_test)
            print( cV.polyFunctional( cD_test , point_test ) )

        index += 1

def hMethodNextVertexCheck( u , V, point , eps , h):
    nextVektor = cV.getNextVertex(u, V, point)
    while(True):
        try:
            nextVektor = cV.getNextVertex( u , V , point)
            break
        except ZeroDivisionError:
            return True
            break

    diff = cV.diffGetNextVertex( u , V , point)
    e_1 = [ 1 , 0 ]
    e_2 = [ 0 , 1 ]

    diff_1 = diff.image( e_1 )
    diff_2 = diff.image( e_2 )

    point_diff1 = M.addVek( point , M.scaleVek( h , e_1 ))
    quotient_diff1 = M.subVek( M.scaleVek( 1 / h , cV.getNextVertex( u , V , point_diff1 ) ) , M.scaleVek( 1 / h , cV.getNextVertex( u , V, point ) ) )

    point_diff2 = M.addVek(point, M.scaleVek(h, e_2))
    quotient_diff2 = M.subVek(M.scaleVek(1 / h, cV.getNextVertex(u, V, point_diff2)), M.scaleVek(1 / h, cV.getNextVertex(u, V, point)))

    difference_1 = M.dist( quotient_diff1 , diff_1)
    difference_2 = M.dist( quotient_diff2 , diff_2)

    print(' the differences equals')
    print(difference_1)
    print(difference_2)

    if difference_1 >= eps or difference_2 >= eps :
        print( ' the differences euqal' )
        print( difference_1 )
        print( difference_2 )
        return False
    return True

def hMethodConeVolIteratorCheck( cD ,  point , eps , h):

    while(True):
        try:
            vektor = cV.coneVolIterator( cD , point)
            break
        except ZeroDivisionError:
            return True
            break

    diff = cV.diffConeVolumeIterator( cD , point)

    e_1 = [ 1 , 0 ]
    e_2 = [ 0 , 1 ]

    diff_1 = diff.image( e_1 )
    diff_2 = diff.image( e_2 )

    point_diff11 = M.addVek( point , M.scaleVek( h , e_1 ))
    point_diff12 = M.subVek(point, M.scaleVek(h, e_1))
    quotient_diff1 = M.subVek( M.scaleVek( 1 / (2*h) , cV.coneVolIterator( cD , point_diff11 ) ) , M.scaleVek( 1 / (2*h) , cV.coneVolIterator( cD , point_diff12 ) ) )

    point_diff21 = M.addVek(point, M.scaleVek(h, e_2))
    point_diff22 = M.subVek(point, M.scaleVek(h, e_2))
    quotient_diff2 = M.subVek(M.scaleVek(1 / (2*h), cV.coneVolIterator( cD , point_diff21 ) ), M.scaleVek(1 / (2*h), cV.coneVolIterator( cD ,  point_diff22)))

    difference_1 = M.dist( quotient_diff1 , diff_1)
    difference_2 = M.dist( quotient_diff2 , diff_2)

    print(' the differences equals')
    print(diff_1)
    print(quotient_diff1)
    print(diff_2)
    print(quotient_diff2)

    if difference_1 >= eps or difference_2 >= eps :
        print( ' the differences equals' )
        print( difference_1 )
        print( difference_2 )
        return False
    return True

def coneVoliteratorTest ( cD , firstPointOfPolygon , eps):
    vector = cV.getConeVolIterator( cD , firstPointOfPolygon )
    if M.dist( vector , firstPointOfPolygon ) >= eps :
        print( ' error at coneVolIterator Test with first polygon vertex and cD :')
        print(firstPointOfPolygon)
        print(cD)
        return False
    return True

# has to be taken with diff from coneVolumeIterator because gamma changes determinant
def getDistOfDet( params , cD ):
    try:
        diff = cV.diffConeVolumeIterator( cD , params)
        det = diff.rows[0][0] * diff.rows[1][1] - diff.rows[1][0] * diff.rows[0][1]
        return math.fabs( det - 1 )
    except ZeroDivisionError or OverflowError:
        return 0

def checkDetOfDiff( params , cD , eps ):
    if getDistOfDet( params , cD) < eps:
        return True
    else:
        return False

def getDistOfDetApprox( params , cD ):
    diff = cV.diffConeVolIteratorApprox( params , cD , 0.001)
    det = diff.rows[0][0] * diff.rows[1][1] - diff.rows[1][0] * diff.rows[0][1]
    return math.fabs( det - 1 )

def checkDetOfDiffApprox( params , cD , eps ):
    if getDistOfDetApprox( params , cD) < eps:
        return True
    else:
        return False

def testerDetOfDiff( repeats, range_test, eps ):
    complete_fail = 0
    approx_fail = 0
    analytic_fail = 0
    for i in range( repeats ):
        P = rP.getRandomPolygon( 5)#random.random() * 5 + 3 )
        cD = cV.getConeVol( P )
        point = [ random.random() * 2 * range_test - range_test , random.random() * 2 * range_test - range_test ]
        if not checkDetOfDiff( point , cD , eps):
            if not checkDetOfDiffApprox( point , cD , eps ):
                complete_fail += 1
            else:
                analytic_fail += 1
        else:
            if not checkDetOfDiffApprox( point , cD , eps):
                approx_fail += 1

    print( [ repeats , complete_fail , analytic_fail , approx_fail ] )

testerDetOfDiff( 10000 , 1 , 0.00000001)



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# could be improved by making delta smaller if check is false ... orby setting point_test in the near of an vertex...
#gradPolyFunctionalTest_random( 1 , 5 , 0.1 , 0.000000000001  , 0.0001 , True )
#polyFunctionalTest( 20 , 5 , 0.000001 )

#polygons_test = [
#    [[9.203330824769969, 5.091434487244416], [-5.6003542709433995, 2.0960543313403526], [-2.469763362739936, -1.9552999854419217], [-1.689083779038116, -2.2524032885794893], [0.5558705879514816, -2.979785544563356]]
#]
#cD_test = cV.getConeVol(polygons_test[0])
#points_test = [ [9.202238149608856, 5.089797764821817] ]

#polygon_easy = [ [ 1, 0 ] , [ 0 , 1 ] , [ -1 , 0 ] , [ 0 , -1 ] ]
#cD_testeasy = cV.getConeVol(polygon_easy)
#point_testeasy = [ 1 , 1.001 ]

#pT.plotPoly( polygons_test[0] , 'r')
#gradPolyFunctionalTest_array( polygons_test , points_test , 0.1 , 0.0001 )
#cD_test = cV.getConeVol( polygons_test[0])
# print( cV.getConeVolIterator( cD_testeasy , point_testeasy ))
#coneVoliteratorTest( cD_testeasy , [ 1 , 1.001 ] , 0.1)
#hMethodConeVolIteratorCheck( cD_testeasy , [1 , 0.000001] , 0.1 , 0.01)




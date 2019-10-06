import matrix as M
import randomPolygon3Centered as rP
import coneVol as cV
import polygonTests as pT
import math
import matplotlib as mp
import machEps as mE

def psi(cD , params):
    return M.norm( cV.f(cD,params))


#print(coneVolumeNewton(cD,startPoint()[0]))

# interesting result                 cD                                          result
#            [[1, 1, 1], [0, -1, 1], [-1, 0, 1],[-1,-1,1],[-1,-2,1]]    [255253327.86799613, 357518606.484159]
#

def getClockRotation(vector):
    return[ vector[ 1 ] , -vector[ 0 ] ]


# exspect vertices to be ordered counter clockwise
def getConeVol(vertices):
    result=[]
    n=len(vertices)
    for i in range( n ):
        v = M.subVek( vertices[ i % n] , vertices[ i - 1 ])
        u_tilde = getClockRotation( v )
        divisor = M.norm(v)
        if u_tilde[0] == 0 and u_tilde[1] == 0:
            divisor = mE.getMachEps()
        u = M.scaleVek( 1.0 / divisor , u_tilde )

        height= M.scal(u,vertices[i-1])
        if height < 0:
            print('height is negativ at ')
            print(i)
        facetVolume= M.norm(v)
        coneVolume= 0.5 * height * facetVolume

        result.append([u[0],u[1],coneVolume])
    return result


cD_test = [[0, 1, 1], [1, 0, 1], [ 0 , -1, 1],[-1,0,1] ]
polygon_hardcoreTest = rP.getRandomPolygon( 5 )
cD_hardcoreTest = getConeVol( polygon_hardcoreTest )
p=[1,0.5]

def getGrid(n):
    result=[]
    l=0
    k=0
    while(l<=n**2):
        while(k<=n**2):
            if(k+l>0):
                result.append([k*1/n,l*1/n])
            k=k+1
        k=0
        l=l+1
    return result

# result will be like [[argmin_x, argmin_y], min]!!!
def min_psi(cD, grid):
    return min( grid , cV.phi )

def min(grid , cD , f ):
    result= [[0,0],float("inf")]
    for pair in grid:
        while(True):
            try:
              if(result[1]>f( cD , pair)):
                result=[pair, f( cD ,pair)]
              break
            except ZeroDivisionError:
                #print('zero division at',pair)
                break
    return result

def scalGrad( cD , params ):
    return M.scal( cV.gradPhi( params , cD ) , params )

def scalAbsGrad( cD , params ):
    return math.fabs(M.scal( cV.gradPhi( params , cD ) , params ))

def getOrthGrad( cD ):

    result = [ [ 0 , 0 ] , math.inf ]
    n = 2

    while( result[1] > eps_orthGrad):
        temp = min( getGrid( n ) , cD , scalAbsGrad )
        if temp[1] < result[1]:
            result = temp
        n = n + 1

    return result





cD_random = cV.getConeVol(rP.getRandomPolygon( 5 ))

eps_orthGrad = 0.01
orthGrad_test =  getOrthGrad( cD_test )

# print(cV.coneVolumeNewton( cD_test , [ orthGrad_test[0][0] , orthGrad_test[0][1] ] ) )

# newton läuft gegen [ 0 , inf ] mit orthGrad_test-startpunkt bei [5.25, 7.0]. Der presolver ist nicht geeignet. Baue lieber Quadrat im 1. Quadranten, der eine Aufsteigsrichtung enthält. Verfeinere dann dieses Quadrat, ohne es zu vergrößern

def startPoint():
    n=4
    result=[[0,0],float("inf")]
    while(result[1]>min_grid):
        grid=getGrid(n)
        result=min_psi(grid)
        if(result[1]==0):
            break
        n=n+1
    return result
# returns empty array, if there is not point with positiv grad. returns point with smallest, positiv grad
def getPositiveGrad( grid , cD ):
    result = []
    gradNorm_result = math.inf
    for point in grid:
        while (True):
            try:
                scal_point = M.scal(cV.gradPhi( point , cD ), point)
                gradNorm_point = M.norm(cV.gradPhi(point, cD))
                if scal_point >= 0 and gradNorm_point < gradNorm_result:
                    result = point
                    gradNorm_result = gradNorm_point
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    return result

def getZeroGradOnLine( cD , point_neg , point_pos  , eps):
    f = scalGrad
    midPoint = M.getMidPoint( point_neg , point_pos )
    while math.fabs( f( cD , midPoint )) > eps:
        if f( cD , midPoint ) > 0:
            point_pos = midPoint
            midPoint = M.getMidPoint( point_neg , point_pos)
        else:
            point_neg = midPoint
            midPoint = M.getMidPoint(point_neg, point_pos)
        print( midPoint )
        print( f( cD ,midPoint) )
    return midPoint

#point_zero = [ 0.1 , 0.001 ]
#print( scalGrad( cD_test , point_zero ))
#print( cV.gradPhi( point_zero , cD_test))
#print( getZeroGradOnLine( cD_test  , [ 1.25 , - 1.25 ] , [0 , 0 ] , 0.001 ))


# is a traingular grid!
# ordering could be improved by starting in the middle and then go to the diagonalStart/diagonalEnd
# but when you start in the middle, you get an edge case every uneventh diagonal (has an even number of points, no middle)
def getDistOrderedGrid( steps , stepSize ):
    grid = []
    n = math.floor( steps / stepSize )
    #range( n ) = [ 0 , 1 , ... , n - 1 ], so ...
    for k in range( n + 1):
        diagonalStart = [ k * stepSize , 0 ]
        grid.append(diagonalStart)
        for l in range(k):
            diagonalStart = M.addVek( diagonalStart , [ - stepSize , stepSize ] )
            grid.append( diagonalStart )

    return grid

# print( getDistOrderedGrid( 3 , 1.0 / 2) )


# endPoint lies in the upper right corner
def transGridEnd( grid , newEndPoint):
    for i in range( len( grid ) ):
        grid[i] = M.subVek( newEndPoint , grid[i] )

# cuts away parts of the grid that are not in the first (positive) quadrant
# hier einen tester schreiben, falls gridOrdering verändert wird
# 'mid' means the mid of the boundary edge (hypothenuse)
def transGridMid( orderedGrid , newMidPoint):

    length = len(orderedGrid)
    start = orderedGrid[0]
    end_leftup = orderedGrid[ length -1 ]
    end_rightdown = [ end_leftup[1] , start[1] ]
    midOfGrid = M.getMidPoint( end_leftup , end_rightdown )
    translation = M.subVek( newMidPoint , midOfGrid )
    i = 0
    # delete every point in orderedGrid, that is not in first (postive) quadrant. Because of 'index out of bounce' error, i and length may have to be reduced
    while i < length:
        orderedGrid[i] = M.addVek( orderedGrid[i] , translation )
        if orderedGrid[i][0] < 0 or orderedGrid[i][1] < 0:
            orderedGrid.pop(i)
            i = i - 1
            length = length - 1
        i = i + 1




grid_test = getDistOrderedGrid( 2 , 1 )
midPoint_test = [ 1 , 1 ]
transGridMid( grid_test , midPoint_test )

if not grid_test == [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [2.0, 0.0], [1.0, 1.0], [0.0, 2.0]]:
    print( 'error at transGridMid' )

def getNearestNegGrad( orderedGrid , cD , params_posGrad):
    transGridMid( orderedGrid , params_posGrad )
    result = [ ]
    for point in orderedGrid:
        while (True):
            try:
                if scalGrad( cD , point ) <= 0:
                    if len(result) == 0 or M.dist( point , params_posGrad ) < M.dist( result , params_posGrad ) :
                        result = point
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    return result

grid_test2 = getDistOrderedGrid( 1 , 0.2)
point_test2 = [ 1.1 , 1.001 ]
result_test2 = getNearestNegGrad( grid_test2 , cD_test , [ 1 , 1 ])
if scalGrad( cD_test , result_test2 ) > 0:
    print( ' test von getNearestNegGrad failed. Result was empty or has wrong scalar product with grad')


def getPosGradStart( cD ):
    posGrad = []
    n = 4
    # dieser part wächst mit n^4 ... kann man zu n^2 verbessern, indem man das grid in vorher berechnete Schranken baut
    while len(posGrad) == 0:
        grid_pos = getGrid(n)
        posGrad = getPositiveGrad(grid_pos, cD)
        n = n + 1
    return posGrad

def getStart_PolytopRetrieval( cD , eps ):

    posGrad = getPosGradStart( cD )
    print( 'habe pos grad gefunden')
    print(posGrad)

    steps = 3
    stepSize = 0.5
    negGrad = []


    while len(negGrad) ==   0:
        grid_neg = getDistOrderedGrid( steps , stepSize )
        transGridMid( grid_neg , posGrad )
        negGrad = getNearestNegGrad( grid_neg , cD , posGrad )
        stepSize = stepSize / 2.0
        if stepSize / 2.0 <= mE.getMachEps():
            steps = steps + 1
            stepSize = 0.5

    print( 'habe neg grad gefunden ')
    print( negGrad )

    return getZeroGradOnLine( cD , negGrad , posGrad , eps )

#pT.plotPoly( polygon_hardcoreTest , 'b' )
#result_test3 = getStart_PolytopRetrieval( cD_hardcoreTest , 0.0001 )
#print( polygon_hardcoreTest )
#print( cV.gamma( cD_hardcoreTest , result_test3 ) )
#print( 'jetzt kommt das phi vom Ergebnis ' )
#print( cV.phi( result_test3 , cD_hardcoreTest ))

#posGrad_test = getPositiveGrad( grid_test , cD_test)
#orderedGrid_test = getDistOrderedGrid( 4 , 0.1 )
#print(posGrad_test)
#negGrad_test = getNearestNegGrad( orderedGrid_test , cD_test , posGrad_test )
#print( negGrad_test )
#if len( negGrad_test ) > 0:
#    print( scalGrad( cD_test , negGrad_test ) )

#print(scalGrad( cD_test , [ 3.0 , 0.65 ] ) )
#zeroGrad_test =  getZeroGradOnLine( cD_test , negGrad_test , posGrad_test , 0.001)

#print( M.norm( cV.gradPhi(  result_test , cD_test)) )
#print( cV.gamma( cD_test , getPositiveGrad( grid_test , cD_test )))
#print( cV.gradPhi( [ 1 , 1 ] , cD_test))
# wie können meine verfahren ( newton , descend ) verwenden, dass das polytop centered ist? dann können sie die antwort vielleicht schneller berechnen...



min_grid=0.2
cD_test = [[0, 1, 1], [1, 0, 1], [ 0 , -1, 1],[-1,0,1],[-1, 2,1]]
p=[1,0.5]

def preSolverPhi( params , cD ):
    v = cV.gamma( cD , params )
    t = M.norm(v)
    n = len( cD )
    if t >= 1 :
        return ( t ** ( 2 * n ) * cV.phi( params , cD ) )
    else:
        return cV.phi( params , cD )

def preSolverPlot( cD , params ):
    x = cV.domain( 0.22 , 0.223 , 0.0001 )
    v = cV.gamma( cD , params )
    for x_value in x:
        stretchedParams = M.scaleVek( x_value , v )
        y_value = preSolverPhi( stretchedParams  , cD )
        mp.plot( x_value , y_value , 'ro' )
    mp.show()

# preSolverPlot( cD_test , p )

# preSolver( cD )


def lineSolution( cD , params ):
    v = cV.gamma( cD , params )
    stepsize = 0.5

def boxGradPreSolver( cD , upperBound)
    #-5.551115123125783e-17
    #habe
    #pos
    #grad
    #gefunden
    #[0.0, 2.25]
    #habe
    #neg
    #grad
    #gefunden
    #[0.0, 1.875]
    #[[46.180908729875426, 82.72003028915721], [43.912175071500265, 82.67680172607484],
    # [-197.32750935609224, -324.58496864389923], [53.78747564167159, 78.59642944297092],
    # [53.44694991304499, 80.59170718569625]]
    #[-1.9793352079767308, 0.5797744254998762]
    #jetzt
    #kommt
    #das
    #phi
    #vom
    #Ergebnis
    #3.9678031744037927e+37


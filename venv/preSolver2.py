import matrix as M
import randomPolygon3Centered as rP
import coneVol2 as cV
import coneVolTest as cVT
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
        v = M.subVek( vertices[ i ] , vertices[ ( i - 1 + n ) % n ])
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
              if(result[1]>f( pair , cD)):
                result=[pair, f( pair , cD )]
              break
            except ZeroDivisionError:
                #print('zero division at',pair)
                break
    return result

def scalGrad( params , cD  ):
    return M.scal( cV.gradPhi( params , cD ) , params )

def scalGradNormed( params , cD ):
    v = params
    w = cV.gradPhi( params , cD )
    result = M.scal( w , v )
    norm_1 = M.norm( v )
    norm_2 = M.norm( w )

    #if norm_2 * norm_1 == 0:
    #    return 0
    return result / ( norm_1 * norm_2 )

def scalGradNormedApprox( params , cD ):
    v = params
    w = cV.gradPhiApprox( params , cD , 0.001)
    result = M.scal( v , w )

    result = result / ( M.norm( v ) * M.norm( w ) )

    return result

def scalAbsGrad( params , cD ):
    return math.fabs(M.scal( cV.gradPhi( params , cD ) , params ))

# exspect epsOrthGrad_array to contain upper bound for scalar product of pont and gradient
def getOrthGrad( cD , epsOrthGrad_array ):

    epsOrthGrad = epsOrthGrad_array[0]
    result = [ 0 , 0 ]
    value_result = math.inf
    n = 2

    while( result[1] > epsOrthGrad):
        temp = min( getGrid( n ) , cD , scalAbsGrad )
        if temp[1] < result[1]:
            result = temp
        n = n + 1

    return result

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

def getLowPointOnGrid( grid , cD , upperBound ):
    result = []
    value_result = upperBound
    for point in grid:
        while (True):
            try:
                value = cV.phi( point , cD )
                if  value_result >= value:
                    result = point
                    value_result = value
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    return result


# exspect array to contain upperBound
def getGoodLowNearPosGrad( cD , upperBound_array ):
    upperBound = upperBound_array[0]
    result = []
    steps = 3
    stepSize = 0.5
    posGrad = getPosGradStart( cD )
    print(posGrad)

    while len(result) == 0:
        grid_neg = getDistOrderedGrid(steps, stepSize)
        transGridMid(grid_neg, posGrad)
        result = getLowPointOnGrid(grid_neg, cD, upperBound)
        stepSize = stepSize / 2.0
        if stepSize + 2 == 2:
            print( ' wow das ist klein ')
            steps = steps + 1
            stepSize = 0.5
    return result




def getZeroGradOnLine( cD , point_neg , point_pos  , eps):
    f = scalGrad
    midPoint = M.getMidPoint( point_neg , point_pos )
    while math.fabs( f( midPoint, cD )) > eps:
        if f( midPoint , cD ) > 0:
            point_pos = midPoint
            midPoint = M.getMidPoint( point_neg , point_pos)
        else:
            point_neg = midPoint
            midPoint = M.getMidPoint(point_neg, point_pos)
        print( midPoint )
        print( f( midPoint , cD ) )
    return midPoint

#point_zero = [ 0.1 , 0.001 ]
#print( scalGrad( cD_test , point_zero ))
#print( cV.gradPhi( point_zero , cD_test))
#print( getZeroGradOnLine( cD_test  , [ 1.25 , - 1.25 ] , [0 , 0 ] , 0.001 ))


# is a traingular grid!
# ordering could be improved by starting in the middle and then go to the diagonalStart/diagonalEnd
# but when you start in the middle, you get an edge case every odd diagonal (has an even number of points, no middle)
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
                if scalGrad( point, cD ) <= 0:
                    if len(result) == 0 or M.dist( point , params_posGrad ) < M.dist( result , params_posGrad ) :
                        result = point
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    return result



def getPosGradStart( cD ):
    posGrad = []
    n = 4
    while len(posGrad) == 0:
        grid_pos = getGrid(n)
        posGrad = getPositiveGrad(grid_pos, cD)
        n = n + 1
    return posGrad


# exspects array to contain [ eps ] whereas eps is upper bound for scalar product of point and gradient
def getOrthogonalGradStart( cD , eps_array ):

    eps = eps_array[0]
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
cD_test = [[1, 0, 1], [ 0 , 1 , 1] , [ -1 , 0, 1] , [0, -1 ,1] ]
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

#print( 'the low point is')
#result_test = getGoodLowPoint( cD_test , 1 )
#print(result_test)
#print( cV.phi( result_test , cD_test ))
#print( cV.phi( [ 1 , 1 ] , cD_test ))
# gthe problem is: low points will come far away which will not lead to a minimum... so I shoudl better search behind the posGrad points...


def getMinOnGrid( grid , f , cD ):
    result = [ -1 ,  -1 ]
    value_result = math.inf

    for point in grid:
        while(True):
            try:
                value_point = f( point, cD )
                if value_point < value_result:
                    value_result = value_point
                    result[0] = point[0]
                    result[1] = point[1]
                break
            except ZeroDivisionError:
                break
    return [ result , value_result ]

def getQuadraticGrid( size, stepSize , midPoint):
    grid = []
    n = math.floor(size * 1.0 / stepSize)
    for i in range(n):
        for j in range(n):
            point = [i * stepSize - size / 2.0 + midPoint[0], j * stepSize - size / 2.0 + midPoint[1]]
            if point[0] >= 0 and point[1] >= 0:
                grid.append( point )

    return grid


# exspects array to contain [ f , upperBound ] whereas f is a non negativ function f(params,coneData )
def quadraticMinSearcher( cD , F_Bound_array ):
    f = F_Bound_array[0]
    upperBound = F_Bound_array[1]
    size = 1000
    stepSize = size / 100.0
    midPoint = [ size / 2.0 , size / 2.0 ]
    grid = getQuadraticGrid( size , stepSize , midPoint )
    minData = getMinOnGrid( grid , f , cD )
    result = minData[0]
    value_result = minData[1]

    while( value_result > upperBound ):
        size = size / 2
        if( size + 128 == 128 ):
            break
        stepSize = size / 100.0
        midPoint = result
        grid = getQuadraticGrid(size, stepSize, midPoint)
        minData = getMinOnGrid(grid, f, cD)
        result = minData[0]
        value_result = minData[1]
    return [ result , value_result ]


#P_hardTest = rP.getRandomPolygon( 5 )
#cD_hardTest = cV.getConeVol( P_hardTest )
#print( 'here comes the polygon:')
#print( P_hardTest )
#print( ' here comes minData of phi with gamma point : ')
#minData_phi = quadraticMinSearcher( cD_hardTest , [ cV.phi , 0.001 ] )
#print( minData_phi )
#print( cV.gamma( cD_hardTest , minData_phi[0] ) )
#print( ' here comes minData of sigma with gamma point : ')
#minData_sigma = quadraticMinSearcher( cD_hardTest , [ cV.sigma , 0.001 ] )
#print( minData_sigma )
#print( cV.gamma( cD_hardTest , minData_sigma[0] ) )

def quadraticMinSearcherTest( repeats , f ):
    eps = 0.00001
    fails = 0
    while repeats > 0:
        print( repeats )
        repeats -= 1
        P = rP.getRandomPolygon( 5 )
        cD = cV.getConeVol( P )

        result = quadraticMinSearcher( cD , [ f , eps ] )[0]
        if M.dist( cV.gamma( cD , result) , P[0]) > 0.1:
            if cV.sigma( result , cD ) < eps:
                print( 'not unique solution' )
                print( P )
                print(result)
                print( cV.gamma( cD , result ) )
                fails += 1
                break
            else:
                print('quadratic searcher failed by by')
                print( cV.sigma( result , cD ) )
                print(P)
                print(result)
                print(cV.gamma(cD, result))
                fails += 1
    return fails

def quadraticMinSearcherTestWithoutGamma( repeats , f ):
    eps = 0.01
    notUnique = 0
    singDiffMatrix = 0
    fails = 0
    while repeats > 0:
        repeats -= 1
        P = rP.getRandomPolygon( 5 )
        cD = cV.getConeVol( P )

        result = quadraticMinSearcher( cD , [ f , eps ] )[0]
        if M.dist( result , P[0]) > 0.1:
            value = M.dist(cV.getConeVolIterator( cD ,  result ) , result)
            if  value < eps:
                print( 'not unique solution' )
                print( P )
                print(result)
                notUnique += 1
                break
            else:
                if f( result , cD ) < 0.01:
                    print(" singular diff matrix ")
                    print(P)
                    print(result)
                    print( f( result , cD ) )
                    singDiffMatrix += 1
                else:
                    print('quadratic searcher failed by')
                    print(P)
                    print(result)
                    fails += 1
    return [notUnique , singDiffMatrix , fails]

def isNilpotent( point , cD ):
    A = cV.diffConeVolumeIterator( cD , point )
    A.scale(-1)
    A.add(M.idMatrix(2))
    A.makeMatrixNormalized()
    #result = diff.rows[0][0] * diff.rows[1][1] - diff.rows[0][1] * diff.rows[1][0]
    result = A.mult(A)
    return result.rows[0][0]**2 + result.rows[1][1]**2 + result.rows[1][0]**2 + result.rows[0][1]**2


#print( quadraticMinSearcherTestWithoutGamma( 100, isNilpotent) )



def nonCenteredCleverLowValue( cD ):
    return True






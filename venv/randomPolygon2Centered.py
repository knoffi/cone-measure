from cmath import sqrt
import math
import numpy
import random
import coneVol as cV
import matrix as M
import matplotlib.pyplot as mp

def getRandomPolygon(n):
    result = [ ]
    angles = randomAngles(n)
    steps = 0
    while( steps < n - 1 ):
        result.append( getNextPoint2( result , angles[steps] , n - 1 - steps ) )
        steps = steps + 1
    # adding a point to make the set result centered may destroy polytop property...
    makeCentered( result )

    return result

def makeCentered( vertices ):
    lastPoint = [ 0 , 0 ]
    for v in vertices:
        lastPoint = M.subVek( lastPoint , v )
    vertices.append( lastPoint )

# better distribution: do not set angle of first point as 0
def getNextPoint2( vertices , angle , stepsAfter):
    if len(vertices) == 0:
        return polar( 0 , randomRadius( 0 , math.inf ) )

    max = maxRadius(vertices, angle , stepsAfter )
    min = minRadius(vertices , angle)

    if min > max:
        print('oh shit')
    r = randomRadius( min , max )
    # print( beta / math.pi * 180 )
    # print(max)
    # print(min)
    result=polar( angle , r )
    #while(noCondition(vertices,result)):
    #    r=randomRadius()
    #    beta=randomAngle(alpha,gamma)
    #    result=polar(beta,r)
    return result

def getAngle(point):
    x=point[0]
    y=point[1]
    if x  == 0 and y > 0:
        return math.pi
    if x == 0 and y < 0:
        return 1.5*math.pi
    if x > 0 and y > 0:
        return math.atan(y/x)
    if x > 0 and y < 0:
        return math.atan(y/x) + 2* math.pi
    if x < 0:

        return math.atan(y/x) + math.pi
    return 0

#better distribution: select random n random angles between 0 and pi, then order them and return an array with every angle
# berücksichtigt das hier die periodizität, also pi = 3*pi

def isContained( array , value ):
    for entry in array:
        if entry == value:
            return True
    return False

def randomAngles( n ):
    result = [ 0 ]
    steps = 1

    while steps < n :
        value = random.random() * 2 * math.pi
        while( isContained( result , value ) ):
            value = random.random() * 2 * math.pi
        result.append( value )
        steps = steps + 1

    result.sort()
    return result

print( randomAngles( 5 ) )


# has to be improved, so that it will give a good distributed number between 0 and +infty
def randomRadius(min, max):
    #randomRadius will always be below 10 + min . Is that good?
    if( max == math.inf ):
        return random.random() * ( 50 - min ) + min

    return min + ( max - min ) * random.random()

def polar(angle, radius):
    x= radius * math.cos(angle)
    y= radius * math.sin(angle)
    return[x,y]

def maxRadius(vertices,angle , steptsAfter):
    result= math.inf
    # improvable by noticing: it is not necessary to check all previous vertices
    others = vertices.copy()
    for v in vertices :
        # remove(x) removes an/any entry which equals x
        others.remove(v)
        for w in others:
            if(getMax(v,w,angle)<result):
                result=getMax(v,w,angle)

    if steptsAfter == 0:
        if( getMax( vertices[0] , vertices[1] , angle ) < result ):
            result = getMax( vertices[0] , vertices[1] , angle)
        print('hey')

    return result

def minRadius(vertices,angle):
    result=0
    # improvable by noticing: it is not necessary to check all previous vertices
    others = vertices.copy()
    for v in vertices :
        # remove(x) removes an/any entry which equals x
        others.remove(v)
        for w in others:
            if( getMin(v,w,angle) > result ):
                result=getMin(v,w,angle)
    return result

def getMax(v,w,angle):
    u_1 = polar(angle,1)
    u_2 = M.subVek(v,w)
    U= M.Matrix([u_1,u_2]).copyTrans()
    result=U.lgsSolve(v)
    # my lgsSolver returns a 'solution' even if v is not in im(M) . Therefore the if-check needs to get adjusted
    if(len(result)==0):
        return math.inf
    if( result[0] <= 0 ):
        return math.inf
    return result[0]

def getMin(v,w,angle):
    u_1 = polar(angle,1)
    u_2 = M.subVek(v,w)
    U = M.Matrix([u_1,u_2]).copyTrans()
    result=U.lgsSolve(v)
    if( len(result) == 0 ):
        return 0
    if ( result[0] <= 0):
        return 0
    return result[0]



v=[0,1]
w=[1,0]
angle= math.pi/4
vertices=[ [ 1, 0 ],  [ 0 , -1 ] , [ -1 , 0 ] ]
# ein result mit zwei gleichen Punkten: [[1.3506563545534096, 0.0], [1.8164857151592047, 0.0], [1.3506563545534096, 0.0], [1.4946652011162949, 0.0]]
P = getRandomPolygon(5)
# - minRadius testen und einbauen in getNextPoint von getRandomPolygon()
# - minRadius ist notwendig, aber vllt muss man bei maxRadius und minRadius nur die Vorgänger überprüfen ... -> mache 2
# - habe einmal [[44.44310714415645, 0.0], [-2.868745831991136, 18.486864345855214], [-14.100910361720999, -23.4225736149311], [-13.38842882794231, -23.35396107847429], [-14.085022122502007, 28.289670347550175]]
#   rausbekommen. Das ist kein 5-Eck, zwei Punkte sind in der konvexen Hülle der anderen enthalten...

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

def plotPoly(vertices):
    x=[]
    y=[]

    for v in vertices:
        x.append(v[0])
        y.append(v[1])
    # make a line between last point and first point
    x.append(vertices[0][0])
    y.append(vertices[0][1])
    mp.plot( x , y , 'r--' )
    mp.plot( x , y , 'ro' )
    mp.axis(getPolyDomain(vertices))
    print(vertices)
    mp.show()

plotPoly(P)
import coneVol2 as cV
import machEps as mE
import matrix as M
import polygonTests as pT
import preSolver2 as pS
import randomPolygon3Centered as rP
import matplotlib.pyplot as mp
import math

def vizPosGrad( cD , n ):
    grid = pS.getGrid( n )
    x = []
    y= []
    b = 0
    d = 0

    for point in grid:
        if pS.scalGrad( cD , point ) >= 0:
            x.append( point[0] )
            y.append( point[1] )
        if point[0] > b :
            b = point[0]
        if point[1] > d:
            d = point[1]

    mp.plot(x, y, 'ro')
    mp.axis( 0 , b , 0 , d )
    # print(vertices)
    mp.show()

def vizLowValue( cD , n , first_bound, second_bound , third_bound):
    grid = pS.getGrid(n)
    x_1 = []
    y_1 = []
    x_2 = []
    y_2 = []
    x_3 = []
    y_3 = []
    b = 0
    d = 0

    for point in grid:
        while (True):
            try:
                value = cV.phi(point, cD)
                if value <= first_bound:
                    x_1.append(point[0])
                    y_1.append(point[1])
                else:
                    if value <= second_bound:
                        x_2.append(point[0])
                        y_2.append(point[1])
                    else:
                        if value <= third_bound:
                            x_3.append(point[0])
                            y_3.append(point[1])

                if point[0] > b:
                    b = point[0]
                if point[1] > d:
                    d = point[1]
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    mp.plot(x_1 , y_1 , 'go')
    mp.plot(x_2 , y_2 , 'yo')
    mp.plot(x_3, y_3, 'ro')
    mp.axis( [ 0, b, 0, d ] )
    # print(vertices)
    mp.show()
    if len(x_1) > 0 :
        params = [ x_1[0] , y_1[0] ]
        vertex = cV.gamma(cD, params)
        print( ' Point near to result , gamma params , phi-value , gradient ')
        print(vertex)
        print( params)
        print( cV.phi( params , cD ) )
        print( cV.gradPhi( params , cD) )
    if len(x_1) > 1 :
        lastIndex = len( x_1 ) - 1
        params = [x_1[lastIndex],y_1[lastIndex]]
        print(' ANOTHER one is  ')
        vertex = cV.gamma(cD, params)
        print(vertex)
        print(params)
        print(cV.phi(params, cD))
        print(cV.gradPhi(params, cD))

def visiualizeLowValueOnGrid( size , stepSize , cD , param , value_1 , value_2 , value_3 ):
    grid = []
    n = math.floor(size * 1.0/stepSize)
    for i in range(n):
        for j in range(n):
            grid.append( [ i * stepSize - size / 2.0 + param[0], j *stepSize - size / 2.0 + param[1]])
    print('grid ist fertig')
    min = math.inf
    point_min = param
    x_1 = []
    y_1 = []
    x_2 = []
    y_2 = []
    x_3 = []
    y_3 = []
    x_4 = []
    y_4 = []
    for point in grid:
        while(True):
            try:
                value= cV.phi(point , cD )

                if value < min:
                    point_min = point
                    min = value
                if value < value_1:
                    x_1.append(point[0])
                    y_1.append([point[1]])
                else:
                    if value < value_2:
                        x_2.append(point[0])
                        y_2.append([point[1]])
                    else:
                        if value < value_3:
                            x_3.append(point[0])
                            y_3.append([point[1]])
                        else:
                            x_4.append(point[0])
                            y_4.append(point[1])

                break
            except ZeroDivisionError:
                break

    a_x = param[0] - size / 2.0
    b_x = param[0] + size / 2.0 + 5

    a_y = param[1] - size / 2.0
    b_y = param[1] + size / 2.0 + 5

    mp.plot( x_1 , y_1 , 'go')
    mp.plot( x_2 , y_2 , 'yo')
    mp.plot( x_3 , y_3 , 'ro')
    mp.plot( x_4 , y_4 , 'bo')
    mp.axes( [a_x , b_x , a_y , b_y])

    print('here comes minimal point, gamma , value, gradient:')
    print( point_min )
    print( cV.gamma( cD , point_min ))
    print( cV.phi( point_min , cD ))
    print( cV.gradPhi( point_min , cD))

    mp.show()

def vizLowGrad( cD , n , first_bound, second_bound , third_bound):
    grid = pS.getGrid(n)
    x_1 = []
    y_1 = []
    x_2 = []
    y_2 = []
    x_3 = []
    y_3 = []
    b = 0
    d = 0

    for point in grid:
        while (True):
            try:
                value = M.norm(cV.gradPhi(point, cD))
                if value <= first_bound:
                    x_1.append(point[0])
                    y_1.append(point[1])
                else:
                    if value <= second_bound:
                        x_2.append(point[0])
                        y_2.append(point[1])
                    else:
                        if value <= third_bound:
                            x_3.append(point[0])
                            y_3.append(point[1])

                if point[0] > b:
                    b = point[0]
                if point[1] > d:
                    d = point[1]
                break
            except ZeroDivisionError:
                # print('zero division at',pair)
                break
    mp.plot(x_1 , y_1 , 'go')
    mp.plot(x_2 , y_2 , 'yo')
    mp.plot(x_3, y_3, 'ro')
    mp.axis( [ 0, b, 0, d ] )
    # print(vertices)
    mp.show()

# besser kontrollierbaren Grid-Maker schreiben


#polygon_ohShit = [[0.9440772944478957, -0.5997771631098889], [-16.804067260622325, -40.38884325863358], [-17.999708289662934, -10.829240656802343], [-32.34730498625317, 88.7483757416121]]
#pT.plotPoly( polygon_ohShit , 'r')
polygon_test = [[2.673368179682499, 3.09152986544487], [1.2086453601351808, 4.28111986768648], [-1.1761317014903958, -0.022433820601322707], [-3.4952312190856785, -4.881491593765966], [0.789349380758395, -2.4687243187640626]]
#pT.plotPoly( polygon_test , 'r')

#print( polygon_test )
cD_test = cV.getConeVol( polygon_test )#[ [ 1 , 0 , 1  ] , [ 0 , 1 , 1 ] , [ -1 , 0 , 1 ] , [ 0 , -1 , 1 ] ]

#vizLowValue( cD_test , 10 , 0.1 , 1 , 500)
#vizLowGrad( cD_test , 10 , 2 , 1 , 1)

visiualizeLowValueOnGrid( 0.01 , 0.0001 , cD_test , [4.3575, 1.6424999999999996] , 0.00001 , 0.0001 , 0.001)
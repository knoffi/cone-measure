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
    mp.plot([param[0]], [param[1]], 'wo')
    mp.axes( [a_x , b_x , a_y , b_y])

    #print('here comes minimal point, gamma , value, gradient:')
    #print( point_min )
    #print( cV.gamma( cD , point_min ))
    #print( cV.phi( point_min , cD ))
    #print( cV.gradPhi( point_min , cD))

    mp.show()

def vizNiveauGradOnGrid( size , stepSize , cD , param , grad , eps):
    grid = []
    value_param = cV.phi( param , cD )
    n = math.floor(size * 1.0/stepSize)
    for i in range(n):
        for j in range(n):
            grid.append( [ i * stepSize - size / 2.0 + param[0], j *stepSize - size / 2.0 + param[1]])
    print('grid ist fertig')

    x_1 = []
    y_1 = []
    for point in grid:
        while(True):
            try:
                value= cV.phi(point , cD )

                if math.fabs(value - value_param) < eps:
                    x_1.append(point[0])
                    y_1.append([point[1]])
                break
            except ZeroDivisionError:
                break

    d= M.scaleVek( 1.0 / M.norm(grad) , grad )
    x_d = []
    y_d = []

    stepBound = 0.71 * size
    step = - stepBound
    while( step <= stepBound ):
        point = M.addVek( param , M.scaleVek( step , d ) )
        x_d.append( point[0] )
        y_d.append( point[1] )
        step = step + stepSize



    a_x = param[0] - size / 2.0
    b_x = param[0] + size / 2.0 + 5

    a_y = param[1] - size / 2.0
    b_y = param[1] + size / 2.0 + 5

    mp.plot( x_1 , y_1 , 'go')
    mp.plot( x_d , y_d , 'yo')

    mp.plot([param[0]], [param[1]], 'wo')
    mp.axes( [a_x , b_x , a_y , b_y])

    #print('here comes minimal point, gamma , value, gradient:')
    #print( point_min )
    #print( cV.gamma( cD , point_min ))
    #print( cV.phi( point_min , cD ))
    #print( cV.gradPhi( point_min , cD))

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

def vizConeDataParamFunction(size, stepSize, function , cD, param, value_1, value_2, value_3):

    grid = []
    n = math.floor(size * 1.0 / stepSize)
    for i in range(n):
        for j in range(n):
            grid.append([i * stepSize - size / 2.0 + param[0], j * stepSize - size / 2.0 + param[1]])
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
        while (True):
            try:
                value = function(point, cD)

                if value < min:
                    point_min = point
                    min = value
                if value < value_1:
                    x_1.append(point[0])
                    y_1.append(point[1])
                else:
                    if value < value_2:
                        x_2.append(point[0])
                        y_2.append(point[1])
                    else:
                        if value < value_3:
                            x_3.append(point[0])
                            y_3.append(point[1])
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

    mp.plot(x_1, y_1, 'go')
    mp.plot(x_2, y_2, 'yo')
    mp.plot(x_3, y_3, 'ro')
    mp.plot(x_4, y_4, 'bo')
    mp.plot([point_min[0]], [point_min[1]], 'wo')
    mp.axes([a_x, b_x, a_y, b_y])

    print('here comes minimal point, gamma , value, gradient:')
    print( point_min )
    #print( cV.gamma( cD , point_min ))
    #print( cV.phi( point_min , cD ))
    #print( cV.gradPhi( point_min , cD))

    mp.show()

def getDistOfDet2( params , cD ):
    while(True):
        try:
            diff = cV.diffConeVolumeIterator( cD , params)
            det = diff.rows[0][0] * diff.rows[1][1] - diff.rows[1][0] * diff.rows[0][1]
            return math.fabs( det - 1 )
        except ZeroDivisionError or OverflowError:
            return 0

def vizDiffDet(size, stepSize, cD, param, value_1, value_2):

    grid = []
    n = math.floor(size * 1.0 / stepSize)
    for i in range(n):
        for j in range(n):
            grid.append([i * stepSize - size / 2.0 + param[0], j * stepSize - size / 2.0 + param[1]])
    x_1 = []
    y_1 = []

    for point in grid:
        try:
            value = getDistOfDet2(param, cD)
        except OverflowError or ZeroDivisionError:
            value = 1

        if value_1 <= value <= value_2:
            x_1.append(point[0])
            y_1.append(point[1])



    a_x = param[0] - size / 2.0
    b_x = param[0] + size / 2.0 + 5

    a_y = param[1] - size / 2.0
    b_y = param[1] + size / 2.0 + 5

    mp.plot(x_1, y_1, 'go')
    mp.axes([a_x, b_x, a_y, b_y])

    mp.show()


# besser kontrollierbaren Grid-Maker schreiben


#polygon_ohShit = [[0.9440772944478957, -0.5997771631098889], [-16.804067260622325, -40.38884325863358], [-17.999708289662934, -10.829240656802343], [-32.34730498625317, 88.7483757416121]]
#pT.plotPoly( polygon_ohShit , 'r')

vertexNumber = 5
polygon_test = [[2.5867398322535426, 0.35370196134820064], [2.816765578870172, 0.6295822977163863], [2.876297688138016, 0.7085658788323118], [3.0156770817979806, 0.920133822103927], [3.2394931255109474, 1.5199632701112769], [3.264789304173254, 1.595191542976061], [3.724878042101699, 3.0170215939392877], [3.814877212031399, 3.562393608394033], [3.805267457792396, 3.625719686521581], [3.7323733429202934, 3.7628157979801506], [3.4945376684513922, 4.032118077767379], [3.4021986658415564, 3.962113846872688], [2.4505649979872577, 3.2406425887674666], [2.385287887899589, 3.1911533801671967], [1.9963758828748048, 2.896301904722165], [1.7526950672949497, 2.7115560857408734], [1.4785918040258232, 2.5037455673666718], [1.1702316222363573, 2.269963267234045], [0.5302116894870117, 1.7847341606554012], [0.44507910318623856, 1.720191157098469], [0.43071118097515304, 1.709298161893054], [0.42830602333192025, 1.7074746991618641], [0.4036475938508899, 1.6887799881124186], [0.3229526946722635, 1.6276014031503094], [0.24413354710762858, 1.5678449133941952], [0.04640950546993311, 1.4179410507701142], [-0.17957128234478759, 1.2466144245712263], [-0.32300368032403437, 1.1378716014708878], [-0.5880100510817052, 0.9369578516354184], [-0.6355138006568612, 0.9009430326182588], [-0.6665890121590617, 0.8773834563844612], [-0.7208415097330012, 0.8362520894116748], [-0.7300663954692099, 0.8292582699055915], [-0.7824545990073234, 0.7895403090423702], [-0.7900257292263145, 0.7838002780334713], [-0.8279311661190941, 0.7550623780209548], [-0.8431717474537053, 0.7435077738836858], [-0.8774692505965414, 0.7175052174523378], [-0.8946005706825306, 0.7045171547510899], [-0.9027183852317148, 0.6983626561604043], [-0.980975673677574, 0.6390321055035905], [-1.1437591255835025, 0.5156182712988058], [-1.1571655604007427, 0.5054542080680189], [-1.2121263435464593, 0.4637855506049331], [-1.2391048516264667, 0.44333170809264905], [-1.2947322885531918, 0.40115750588408106], [-1.4065446449744032, 0.31638605977647105], [-1.8161441053769063, 0.0058434199872439], [-1.823014096743991, 0.0006348551607779501], [-1.8288016930344846, -0.0037530791232921606], [-1.9145610626624578, -0.06877255557914909], [-1.9936022040773425, -0.12869933204554046], [-2.2038634659290217, -0.28811493361000384], [-2.814872197671972, -0.7514658864113118], [-3.5227364675839805, -1.2884736665985472], [-4.6384411646339885, -2.134891875319733], [-5.179392634218635, -2.5453042166976307], [-4.575135177671799, -2.5491832973702873], [-4.406993908926179, -2.5484823712208424], [-3.89424447437459, -2.542094311062124], [-3.6642999583764757, -2.534673161773632], [-3.0202764273616363, -2.5086776327746], [-2.0470511252220422, -2.4543626947998436], [-1.7612691106880378, -2.429745430121149], [-1.5404746093384005, -2.4103618645634293], [-1.5287149810424265, -2.4091034912847946], [-1.3894247703113551, -2.393782991805184], [-1.2994005389089274, -2.38363505034361], [-1.206691313666303, -2.3721780133196804], [-0.8711180767800243, -2.3042102161462186], [-0.7550131627905892, -2.267194799492451], [-0.6843971511040434, -2.221638072237676], [-0.5845851531938262, -2.1569017289347756], [-0.2994920181681502, -1.9672474240024926], [-0.2192620692964681, -1.9120937245860663], [-0.16420466742060374, -1.8742446690311845], [-0.08707743585254532, -1.821223712225415], [-0.08086383376133327, -1.8169519967919534], [-0.027151202685557208, -1.7797807755444641], [-0.007573798549373401, -1.766205938865122], [0.1395859706232025, -1.6640090617139216], [0.2652031839025275, -1.576470287551836], [0.3462086463897328, -1.5193205447610185], [0.3754295420740775, -1.4985757427793513], [0.4865241418674051, -1.4128001450344456], [0.8805607448622572, -1.0820321270396074], [1.1339464919009248, -0.8690378214383169], [1.2525537079152953, -0.7692638914303036], [1.3318867589631302, -0.7025274548171956], [1.3753089931526346, -0.6659995318706317], [1.4370282196278574, -0.6140795033168567], [1.549044208701108, -0.5198480589690921], [1.7401003076406498, -0.3591120554768841], [1.8107863165430262, -0.2996412267221694], [1.8593681193512923, -0.25875992631402345], [1.9162321726508709, -0.21090661037899877], [1.9204859833256185, -0.20732680009425675], [2.1090823744184988, -0.048606704424120004], [2.1211746417837225, -0.03842949148610582], [2.1369155998954, -0.02517399121459922]]
polygon_test = [ [ 1 , 1 ] , [ -1 , 1 ] , [-1 , -1 ] , [ 1 , -1 ] ]#rP.getRandomPolygon(vertexNumber)
pT.plotPoly( polygon_test , 'r')
print( 'here comes the polygon')
print( polygon_test )
cD_test = cV.getConeVol( polygon_test )#[ [ 1 , 0 , 1  ] , [ 0 , 1 , 1 ] , [ -1 , 0 , 1 ] , [ 0 , -1 , 1 ] ]

#vizLowValue( cD_test , 10 , 1 , 10 , 100)
#vizLowGrad( cD_test , 10 , 2 , 1 , 1)

a_plot = 1 * 5.0
b_plot = 10 * 5.0
c_plot = 100 * 5.0

size_plot = 100
midPoint_plot= [ size_plot / 2.0 , size_plot / 2.0 ]
stepSize_plot = size_plot / 10.0
#vizConeDataParamFunction( size_plot , stepSize_plot , cV.phi , cD_test , midPoint_plot , a_plot , b_plot , c_plot )
#vizConeDataParamFunction( size_plot , stepSize_plot , cV.sigma , cD_test , midPoint_plot , a_plot , b_plot , c_plot )
#vizDiffDet( size_plot , stepSize_plot , cD_test , midPoint_plot, 0 , 0.0001 )
#vizConeDataParamFunction( size_plot , stepSize_plot , pS.scalGradNormed , cD_test , midPoint_plot , 0 , 0.5 , 0.8)
#vizConeDataParamFunction( 10 , 0.1 , pS.scalGradNormedApprox , cD_test , [ 5, 5 ] , 0 , 0.5 , 0.8)
#vizConeDataParamFunction( 10 , 0.1 , pS.scalGradNormedApprox , cD_test , [ 5, 5 ] , 0 , 0.5 , 0.8)

n = 1

def scalarDifferenceGradient( point , cD):
    v = cV.getConeVolIterator( cD , point )
    d = M.subVek( v , point )

    diff = cV.diffConeVolumeIterator( cD , point )
    diff.scale(-1)
    diff.add(M.idMatrix(2))

    u = M.Matrix([d]).mult(diff)

    # equal to scalar product of u and d :
    return u.image(d)[0]

vizConeDataParamFunction( 10 , 0.1 , scalarDifferenceGradient , cD_test , [ 5 , 5  ] , -0.5 , 0 , 0.5)

def phiNormScaled( params , cD):
    return cV.phi( params , cD ) * ( ( M.norm( params ) ) ** n + 1 )
#vizConeDataParamFunction( size_plot , stepSize_plot , phiNormScaled , cD_test , midPoint_plot , a_plot , b_plot , c_plot)
#n = 2
#vizConeDataParamFunction( size_plot , stepSize_plot , phiNormScaled , cD_test , midPoint_plot , a_plot , b_plot , c_plot)
#n = 3
#vizConeDataParamFunction( size_plot , stepSize_plot , phiNormScaled , cD_test , midPoint_plot , a_plot , b_plot , c_plot)


quadratic_result = [[10.078125, 7.96875], 0.002654158755822429]
print( cV.gamma( cD_test , quadratic_result[0]))

def bla():
    vizDiffDet(size_plot, stepSize_plot, cD_test, midPoint_plot, 0 , 0.00001 )

#bla()

# Interessant: bei polygon_test = [[2.673368179682499, 3.09152986544487], [1.2086453601351808, 4.28111986768648], [-1.1761317014903958, -0.022433820601322707], [-3.4952312190856785, -4.881491593765966], [0.789349380758395, -2.4687243187640626]]
# gibt es den Punkt = [0.0, 63.0] ? bei dem gradPhi eine Singularität hat, obwohl phi keine Singularität hat ! Wegen Numerik?
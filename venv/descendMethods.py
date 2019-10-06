import coneVol2 as cV
import machEps as mE
import matrix as M
import polygonTests as pT
import preSolver2 as pS
import randomPolygon3Centered as rP
import matplotlib.pyplot as mp
import BFGS

eps_descend = mE.getMachEps()


def stepSize_posGrad(cD, point, d , n):
    alpha = 0.5
    nextPoint = M.addVek(point, M.scaleVek(alpha, d))
    while ( cV.phi(point, cD ) < cV.phi( nextPoint , cD ) ): #or pS.scalGrad( cD , nextPoint ) < 0 ):
        alpha = alpha * 0.5
        nextPoint = M.addVek(point, M.scaleVek(alpha, d))
    if n % 300 == 0:
        print( alpha )
    return alpha

def notEfficient( cD , point , d , stepSize , crease ):
    v = M.addVek( point , M.scaleVek( stepSize , d ) )
    upper_bound = cV.phi( point , cD ) + crease * stepSize * M.scal( cV.gradPhi( point , cD ) , d )
    return cV.phi( v , cD) > upper_bound

def stepSizeArmijo( cD, point , d , crease , efficiency , beta_1 , beta_2):
    result = - efficiency * ( M.scal( cV.gradPhi( point , cD) , d ) / ( M.norm(d) ) **2 )
    while notEfficient( cD , point , d , result , crease ):
        result = beta_2 * result
    return result

# works only if global minimum is zero and if its the only extremum
def coneVolumeDescend( cD , start ):

    result = start
    n=0

    print( 'gehe in while schleife ')
    while ( cV.phi( result , cD ) > eps_descend ):
        d = M.scaleVek( -1 , cV.gradPhi( result , cD ) )
        d = M.scaleVek( 1 / M.norm(d) , d )
        alpha = stepSize_posGrad( cD, result , d , n )
        result = M.addVek( result , M.scaleVek( alpha, d ) )
        if n % 300 == 0:
            print('here comes result')
            print(result)
            print( cV.gradPhi( result , cD ) )
            print( cV.phi( result , cD ))
        n = n + 1
    return result


def coneVolumeDescendArmijo( cD , start , crease , efficiency , beta_1 , beta_2 ):

    result = start
    n=0

    previous = [ -1 , -1 ]

    print( 'gehe in while schleife ')
    while ( cV.phi( result , cD ) > eps_descend ):
        if M.dist( previous , result) > 0.1:
            previous[0] = result[0]
            previous[1] = result[1]
            d = M.scaleVek( -1 , cV.gradPhi( result , cD ) )
            #d = M.scaleVek( 1 / M.norm(d) , d )
            alpha = stepSizeArmijo( cD, result , d , crease , efficiency , beta_1 , beta_2 )
            result = M.addVek( result , M.scaleVek( alpha, d ) )
        else:
            print('here com previous and result')
            print(result)
            print(previous)
            result = coneVolumeDescendBFGSApprox( cD , result , previous , crease , efficiency , beta_1 , beta_2 )

        n = n + 1
    return result

def coneVolumeDescendBFGS( cD , params_new, params_prev , crease , efficiency , beta_1 , beta_2 ):
    A_k = M.idMatrix(2)
    n = 0
    while( cV.phi( params_new , cD) >  eps_descend):
        A_k = BFGS.getBFGS( cD , params_new , params_prev , A_k)
        d = A_k.image( M.scaleVek( -1 , cV.gradPhi( params_new , cD ) ) )
        d = M.scaleVek( 1.0 / M.norm(d) , d)

        alpha = stepSizeArmijo( cD , params_new , d ,  crease , efficiency , beta_1 , beta_2 )

        params_prev = [params_new[0] , params_new[1]]
        d=  M.scaleVek( alpha , d )
        params_new = M.addVek( params_new , d )

        if n % 1 == 0:
            print(n)
            print('here comes params_new and params_prev BFGS result')
            print(params_new)
            print(params_prev)
            print(cV.gamma( cD , params_new))
            print( cV.gradPhi( params_new , cD ) )
            print( cV.phi( params_new , cD ))

        n = n + 1
    return params_new

def coneVolumeDescendBFGSApprox( cD , params_new, params_prev , crease , efficiency , beta_1 , beta_2 ):
    A_k = M.idMatrix(2)
    n = 0
    while( cV.phi( params_new , cD) >  eps_descend):
        A_k = BFGS.getBFGS( cD , params_new , params_prev , A_k)
        d = A_k.image( M.scaleVek( -1 , cV.gradPhiApprox( params_new , cD , 0.0001 ) ) )
        d = M.scaleVek( 1.0 / M.norm(d) , d)

        alpha = stepSize_posGrad( cD , params_new , d ,  n ) #crease , efficiency , beta_1 , beta_2 )

        params_prev = [params_new[0] , params_new[1]]
        d=  M.scaleVek( alpha , d )
        params_new = M.addVek( params_new , d )

        if n % 1 == 0:
            print(n)
            print('here comes params_new and params_prev BFGS result')
            print(params_new)
            print(params_prev)
            print(cV.gamma( cD , params_new))
            print( cV.gradPhi( params_new , cD ) )
            print( cV.phi( params_new , cD ))

        n = n + 1
    return params_new




def getPolygonReversedOrder( polygon ):
    result = []
    n = len(polygon)

    for i in range(n):
        result.append( polygon[n-1 - i] )
    return result
# es könnte sein, dass eines der Probleme die Verschiebung der Quadranten ist... damals in ConeVol hatte gamma den ersten und letzten Teil des coneDatas genommen... muss jetzt nochmal alles verändert werden?
cD_test = [ [ 1 , 0 , 1  ] , [ 0 , 1 , 1 ] , [ -1 , 0 , 1 ] , [ 0 , -1 , 1 ] ]
#cD_test = getPolygonReversedOrder( cD_test )

#print( pS.scalGrad( cD_test , start_test ) )
#print( start_test )
#result_test = coneVolumeDescendArmijo( cD_test , start_test , crease_test , efficiency_test , beta_1test , beta_2test )
#result_coordTest = cV.gamma( cD_test , result_test )
#print( result_coordTest )
#print( result_coordTest[0] * result_coordTest[1])

beta_1test = 0.5
beta_2test = 0.8
crease_test = 0.2
efficiency_test = 100
start_test = [ 0.9 , 1.1 ] #pS.getPosGradStart( cD_test )

polygon_test2 = [[2.673368179682499, 3.09152986544487], [1.2086453601351808, 4.28111986768648], [-1.1761317014903958, -0.022433820601322707], [-3.4952312190856785, -4.881491593765966], [0.789349380758395, -2.4687243187640626]]
 #rP.getRandomPolygon( 5 )

print( ' hier kommt das polytop mit coneData')
print(polygon_test2)

cD_test2 = cV.getConeVol( polygon_test2 )

print( cD_test2 )

#pT.plotPoly( polygon_test2 , 'r')
start_test2 = [4.7, 1.6]

print( 'hier kommt der start' )

print( start_test2 )
result_test2 = coneVolumeDescendArmijo( cD_test2 , start_test2 , crease_test , efficiency_test , beta_1test , beta_2test )

print( 'hier kommt result und sein Abstand zum ersten Eckpunkt')
print( cV.gamma( cD_test2 , result_test2 ) )
print( M.dist( cV.gamma(result_test2) , polygon_test2[0] ) )

# stepSize scheint lange zu brauchen bei P = [[88.6803578610117, 8.551930964805006], [88.5180013132855, 8.760717405689423], [85.52161100268587, 10.739874317313332], [-351.45828922627163, -33.498541339396326], [88.73831904928856, 5.4460186515885685]]
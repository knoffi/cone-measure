import coneVol2 as cV
import machEps as mE
import matrix as M
import polygonTests as pT
import preSolver2 as pS
import randomPolygon3Centered as rP
import matplotlib.pyplot as mp
import stepSize as sS
import BFGS

eps_descend = 0.001 #mE.getMachEps()


# TODO dafür sorgen, dass ich für jede verwendete Schrittweiten-methode eine andere Methode schreibe und alles auch mit Approx!

# works only if global minimum is zero and if its the only extremum
def coneVolumeDescend( cD , start ):

    result = start
    previous = [ -1 , -1 ]
    n=0

    print( 'gehe in while schleife ')
    while ( cV.phi( result , cD ) > eps_descend ):
        if M.dist(previous, result) > 0.1:
            previous[0] = result[0]
            previous[1] = result[1]
            d = M.scaleVek(-1, cV.gradPhiApprox(result, cD , 0.0001))
            # d = M.scaleVek( 1 / M.norm(d) , d )
            alpha = sS.stepSize_posGrad(cD, result , d , n)
            result = M.addVek(result, M.scaleVek(alpha, d))
            n = n+1

        else:
            print('versuche es mit BFGS und BFGSApprox, Start bei:')
            print( cV.gamma( cD , result ) )

            result_2 = coneVolumeDescendBFGS(cD, result, previous, crease_test, efficiency_test, beta_1test, beta_2test)
            print('BFGS ist fertig mit:')
            print( cV.gamma( cD , result_2 ))
            result_1 = coneVolumeDescendBFGSApprox(cD, result, previous, crease_test, efficiency_test, beta_1test, beta_2test)
            print( 'BFGS Approx ist fertig!')

            if cV.phi(result_1, cD) < cV.phi(result_2, cD):
                return result_1
            else:
                return result_2

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
            alpha = sS.stepSize_posGrad( cD, result , d , crease , efficiency , beta_1 , beta_2 )
            result = M.addVek( result , M.scaleVek( alpha, d ) )
        else:
            print('versuche es mit BFGS und BFGSApprox')
            result_1 = coneVolumeDescendBFGSApprox( cD , result , previous , crease , efficiency , beta_1 , beta_2 )
            result_2 = coneVolumeDescendBFGS(cD, result, previous, crease, efficiency, beta_1, beta_2)

            if cV.phi( result_1 , cD ) < cV.phi( result_2 , cD ):
                return result_1
            else:
                return result_2

        n = n + 1
    return result

def coneVolumeDescendBFGS( cD , params_new, params_prev , crease , efficiency , beta_1 , beta_2 ):
    A_k = M.idMatrix(2)
    n = 0

    while( cV.phi( params_new , cD) >  eps_descend):

        # ICH VERÄNDERE HIER MEIN CREASE...
        crease = 0.000001 #cV.phi( params_new , cD) * 0.00000001

        while (True):
            try:
                A_k = BFGS.getBFGS(cD, params_new, params_prev, A_k)
                break
            except ZeroDivisionError:
                print('es geht nicht besser')
                return params_new
                break

        antiGrad = M.scaleVek(-1, cV.gradPhi( params_new , cD ))
        d = A_k.lgsSolve(antiGrad)
        d = M.scaleVek( 1.0 / M.norm(d) , d)

        alpha = sS.stepSize_posGrad( cD , params_new , d , n ) # crease , efficiency , beta_1 , beta_2 )
        d = M.scaleVek(alpha, d)

        params_prev = [params_new[0] , params_new[1]]
        params_new = M.addVek( params_new , d )

        if( n % 300 == 0):
            print(cV.gamma( cD , params_new ))

        n = n + 1

    return params_new

def coneVolumeDescendBFGSApprox( cD , params_new, params_prev , crease , efficiency , beta_1 , beta_2 ):
    A_k = M.idMatrix(2)
    n = 0
    while( cV.phi( params_new , cD) >  eps_descend):
        while(True):
            try:
                A_k = BFGS.getBFGSApprox( cD , params_new , params_prev , A_k)
                break
            except ZeroDivisionError:
                print('es geht nicht besser')
                return params_new
                break

        antiGrad = M.scaleVek( -1 , cV.gradPhiApprox( params_new , cD , 0.0001 ) )
        d =  A_k.lgsSolve(antiGrad)
        d = M.scaleVek( 1.0 / M.norm(d) , d)

        alpha = sS.stepSizeArmijo( cD , params_new , d , crease , efficiency , beta_1 , beta_2 )
        d = M.scaleVek(alpha, d)

        params_prev = [params_new[0], params_new[1]]
        params_new = M.addVek(params_new, d)
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
beta_2test = 0.5
crease_test = 0.5
efficiency_test = 50
start_test = [ 0.9 , 1.1 ]
#start_test = pS.getPosGradStart( cD_test )

polygon_test2 = [[2.673368179682499, 3.09152986544487], [1.2086453601351808, 4.28111986768648], [-1.1761317014903958, -0.022433820601322707], [-3.4952312190856785, -4.881491593765966], [0.789349380758395, -2.4687243187640626]]
 #rP.getRandomPolygon( 5 )

print( ' hier kommt das polytop mit coneData')
print(polygon_test2)

cD_test2 = cV.getConeVol( polygon_test2 )

print( cD_test2 )

#pT.plotPoly( polygon_test2 , 'r')
start_test2 = [4.7, 1.6]
start_test2 = pS.getLowPoint( cD_test2 , 0.001)
print( 'hier kommt der start' )

print( start_test2 )
result_test2 = coneVolumeDescend( cD_test2 , start_test2 )

print( 'hier kommt result und sein Abstand zum ersten Eckpunkt')
print( cV.gamma( cD_test2 , result_test2 ) )
print( M.dist( cV.gamma( cD_test2 , result_test2 ) , polygon_test2[0] ) )

# stepSize scheint lange zu brauchen bei P = [[88.6803578610117, 8.551930964805006], [88.5180013132855, 8.760717405689423], [85.52161100268587, 10.739874317313332], [-351.45828922627163, -33.498541339396326], [88.73831904928856, 5.4460186515885685]]
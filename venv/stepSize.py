import matrix as M
import math
import coneVol2 as cV
import machEps as mE

eps = 100 * mE.getMachEps()

def stepSize_posGrad(cD, point, d , n):
    alpha = 0.5
    nextPoint = M.addVek(point, M.scaleVek(alpha, d))
    while ( cV.phi(point, cD ) + eps < cV.phi( nextPoint , cD ) ): #or pS.scalGrad( cD , nextPoint ) < 0 ):
        alpha = alpha * 0.5
        nextPoint = M.addVek(point, M.scaleVek(alpha, d))
    if n % 300 == 0:
        print( alpha )
    if alpha < eps == 1:
        return 0
    return alpha

def notEfficient( cD , point , d , stepSize , crease ):
    v = M.addVek( point , M.scaleVek( stepSize , d ) )
    upper_bound = cV.phi( point , cD ) + crease * stepSize * M.scal( cV.gradPhi( point , cD ) , d ) / M.norm(d)
    return cV.phi( v , cD) > upper_bound

def stepSizeArmijo( cD, point , d , crease , efficiency , beta_1 , beta_2):
    result = - efficiency * ( M.scal( cV.gradPhi( point , cD) , d ) / ( M.norm(d)**2 ) )
    while notEfficient( cD , point , d , result , crease ):
        result = beta_2 * result
    return result

def powellFunction_1( cD , params , d , step ):
    dividend = cV.phi( M.addScaleVek( params , step , d )   , cD ) - cV.phi( params , cD )
    divisor = step * M.scal( cV.gradPhi( params , cD ) , d)

    if( step <= eps):
        return 1

    return dividend / divisor


def powellFunction_2( cD, params , d , step):
    v = cV.gradPhi( M.addScaleVek( params , step , d ))
    dividend = M.scal( v , d )
    divisor = M.scal( cV.gradPhi( params , cD ) , d )

    return dividend / divisor


# step_start > 0 , gamma > 1 , 0 < alpha < ro < 1
def startStepSizePowell( cD , point , d , step_start , gamma , alpha , ro):
    return 0


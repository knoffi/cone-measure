import coneVol2 as cV
import machEps as mE
import matrix as M
import polygonTests as pT
import preSolver2 as pS
import randomPolygon3Centered as rP
import matplotlib.pyplot as mp


# I use copyTrans() two times on s_Matrix and y_Matrix, thus speed of this can be improved
def getBFGS( cD , params_new , params_k , A_k ):
    s_vector = M.subVek(params_new, params_k)
    y_vector = M.subVek( cV.gradPhi( params_new , cD ) , cV.gradPhi( params_k , cD ) )

    # s_Matrix and y_Matix is a column, not a row, therefor copyTrans
    s_Matrix = M.Matrix( [ s_vector ] ).copyTrans()
    y_Matrix = M.Matrix( [ y_vector ] ).copyTrans()

    dividend = A_k.mult(s_Matrix)
    dividend = dividend.mult(dividend.copyTrans())
    divisor = M.scal( s_vector , A_k.image( s_vector))
    quotient = dividend
    quotient.scale( 1.0/divisor )

    rankOneMod = M.subMatrix( A_k , quotient)

    dividend2 = y_Matrix.mult( y_Matrix.copyTrans() )

    # here could be a division through zero. If this occurence, then I should just translate params_new a liiiitle bit...
    divisor2 = M.scal(y_vector , s_vector )

    quotient = dividend2
    quotient.scale( 1.0 / divisor2 )

    rankTwoMod = M.addMatrix( rankOneMod , quotient )

    return rankTwoMod





    return rankOneMod

cD_test = [ [ 1 , 0 , 1 ] ,  [ 0 , 1 , 1 ] , [ -1 , 0 , 1 ] , [ 0 , -1 , 1 ] ]
params_newTest = [ 1.0 , 1.5 ]
params_kTest = [ 1.2 , 1.0 ]
s_vectorTest = M.subVek(params_newTest, params_kTest)
A_kTest= M.Matrix([ [ 9 , 2 ] , [ 2 , 7 ] ])
result_test = getBFGS( cD_test , params_newTest , params_kTest , A_kTest  )
result_test.list()

print(result_test.image(s_vectorTest))
print('here comes the anti gradient of the new params and the bfgs direction:')
grad_test = cV.gradPhi( params_kTest , cD_test )
direction_test = A_kTest.image( grad_test )
print(M.scaleVek( 1.0 / M.norm(grad_test) , grad_test ) )
print( M.scaleVek( 1.0 / M.norm(direction_test) , direction_test ))

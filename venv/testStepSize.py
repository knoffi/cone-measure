import coneVol2 as cV
import matrix as M
import machEps as mE
import vizualization as v


polygon_test2 = [[2.673368179682499, 3.09152986544487], [1.2086453601351808, 4.28111986768648], [-1.1761317014903958, -0.022433820601322707], [-3.4952312190856785, -4.881491593765966], [0.789349380758395, -2.4687243187640626]]
cD = cV.getConeVol( polygon_test2 )

params_now = [4.7, 1.6152821997297826]
print(cV.gamma( cD , params_now))
print(cV.phi( params_now , cD ))
d = [-0.2083940408151545, -0.9780449497608644]
grad = M.scaleVek(1.0 / M.norm( cV.gradPhi( params_now , cD)) , cV.gradPhi( params_now , cD) )
grad_approx = cV.gradPhiApprox( params_now , cD , 0.0000001)
stepSize = 1
params_next = M.addVek( params_now , M.scaleVek( stepSize , d))
v.visiualizeLowValueOnGrid(0.001 , 0.00001 , cD , params_now , 0.02405 , 0.024059 , 0.02406)
v.vizNiveauGradOnGrid(0.001 , 0.00001 , cD , params_now , d , 0.000001)
v.vizNiveauGradOnGrid(0.001 , 0.00001 , cD , params_now , grad , 0.000001)
v.vizNiveauGradOnGrid(0.001 , 0.00001 , cD , params_now , grad_approx , 0.000001)
n = 0
diff = ( 1 * cV.phi(params_now , cD )) - ( 1 * cV.phi( params_next , cD ))
d = [ -1 , 2]
while( diff < 0  ):
    stepSize = stepSize * 0.5
    params_next = M.addVek(params_now, M.scaleVek(stepSize, d))
    diff = cV.phi(params_now , cD ) - cV.phi( params_next , cD )
    if ( n +1 ) % 10 == 0:
        print(diff)


print(stepSize)

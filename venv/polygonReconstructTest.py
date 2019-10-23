import polygonReconstruct as pR
import randomPolygon3Centered as rP
import coneVol2 as cV
import math
import random
import matrix as M

repeats = 0
noFail = 0
fail_soft = 0
fail_medium = 0
fail_hard = 0

while repeats < 50:
    eps_test =  0.00001
    P_test = rP.getRandomPolygon( math.floor( 10 * random.random() + 3 ) )
    cD_test = cV.getConeVol( P_test )

    result = pR.polygonReconstructCentered( cD_test , 0.0001 , 0.5 , 0.9 , 0.5 , 0.5 , eps_test)

    if M.dist( cV.gamma( cD_test , result ) , P_test[0]) > 0.01:
        sigma_result = cV.sigma(result , cD_test )

        if( sigma_result < eps_test):
            noFail += 1
        else:
            if( sigma_result < eps_test * 100 ):
                fail_soft += 1
            else:
                if (sigma_result < eps_test * 100 * 100):
                    fail_medium += 1
                else:
                    fail_hard += 1


        if sigma_result < eps_test:
            print(' not unique!!! polygon, result, gamma result , phi value, sigma value:')
            print( P_test )
            print( result )
            print( cV.gamma( cD_test , result ))
            print( cV.phi( result , cD_test ) )
            print( cV.sigma( result , cD_test ) )
        else:
            print(' FAIL!!  polygon, result, gamma result , phi value, sigma value:')
            print(P_test)
            print(result)
            print(cV.gamma(cD_test, result))
            print(cV.phi(result, cD_test))
            print(cV.sigma(result, cD_test))

    repeats += 1
    print(repeats)

print( [ noFail , fail_soft , fail_medium , fail_hard ] )

# results where everything worked
P_1 = [[3.50838286589081, -1.5107444243803627], [2.325398721911787, -0.21029997610538254], [-1.485558813118278, 3.5717000639909195], [-2.8198717746014843, 2.0978600608178417], [-1.5283510000828333, -3.9485157243230167]]
result_1 = [1.6455078125, 3.017578125]
P_2 = [[2.7625155828289434, 3.5462931123540065], [-0.7791440004189139, 0.22586884369998636], [-1.0551126144569645, -0.4313604423010734], [-1.1614789148144171, -1.416368117005431], [0.23321994686135245, -1.9244333967474883]]
result_2 = [1.8603515625, 2.71484375]
results_good = [ [ P_1 , result_1] , [ P_2 , result_2 ] ]
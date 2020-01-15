import coneVol2 as cV
import polygonTests as pT
import matrix as M
#import vizualization

def polygonReconstruct( cD , vertex_first ):
    result = [ vertex_first ]
    otherVertices = cV.getConeVolIteratedVertices( cD , vertex_first )
    n = len(otherVertices)
    index = 0

    while index < n - 1 :
        result.append(otherVertices[ index ])
        index = index + 1



    return result


def compareResults( P , alternativeParams):
    cD = cV.getConeVol( P )
    alternatePoint = cV.gamma( cD , alternativeParams )

    P_alternative = polygonReconstruct( cD , alternatePoint )
    cD_alternative = cV.getConeVol( P_alternative )

    pT.plotPoly( P , 'b')
    pT.plotPoly( P_alternative , 'r')
    print( len(P_alternative ))
    print( 'difference in coneVolume matrices')
    print( M.distMatrix( M.Matrix(cD) , M.Matrix(cD_alternative) ))
    print( ' sigma value of alternative params')
    print( cV.sigma(alternativeParams , cD ) )

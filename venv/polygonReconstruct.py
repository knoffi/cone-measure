import preSolver2 as pS
import descendMethods as dM
import randomPolygon3Centered as rP
import coneVol2 as cV
def polygonReconstructCentered( cD_centered , start_bound , crease ,  efficiency , beta_1 , beta_2 , eps):
    cD = cD_centered
    dM.eps_descend = eps
    # result will be like [ start, value_start ]
    start = pS.quadraticMinSearcher( cD , [ cV.sigma , start_bound ] )[0]

    result = dM.coneVolumeDescendArmijo( cD , start , crease , efficiency , beta_1 , beta_2 )
    return result
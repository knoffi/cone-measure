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

# examples where qudratic searcher failed
P_1 = [[57.20537767280718, 42.28766842082867], [57.91318199804691, 44.27186911171974], [58.42079511943523, 45.74267073041455], [-229.8723195447645, -173.55223874552073], [56.33296475447521, 41.25003048255775]]
alternateParams_1 = [10.594346662890045, 14.447841263504436]

P_2 = [[40.82526023155419, -45.43410551553221], [-167.10154616440508, 196.01157780595076], [40.587628233827246, -49.85362199195141], [42.67453006415001, -52.17960094319114], [43.01412763487367, -48.54424935527602]]
alternateParams_2 = [92.62605192020703, 146.51875351224746]

P_3 = [[9.050900470991463, 2.3487168724296517], [8.826425787439392, 2.9108511362725995], [-35.04288737494956, -2.7940254463476935], [7.102866044533535, -2.0051110311852893], [10.062695071985175, -0.4604315311692686]]
alternateParams_3 = [86.40624999999979, 87.00433280135637]


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



# hier haben beide quadratischen Suchverfahren verkackt...
#[[128.2477454625004, 209.14570676021197], [-514.9895573465076, -829.7657853974156], [128.98912168765278, 206.23349880950573], [128.93562085930907, 206.76735407370927], [128.8170693370453, 207.61922575398827]]
 #here comes minData of phi with gamma point :
#[[1479.9999999999634, 918.9214957890601], 1.517361875548798]
#[-33.45394684653445, 2167.9854804271117]
 #here comes minData of sigma with gamma point :
#[[62.202550635158715, 0.0], 5023.109934924761]
#[-21.736755404708155, 58.28096404483522]

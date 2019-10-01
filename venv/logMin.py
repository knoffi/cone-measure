import randomPolygon3Centered as rP
import coneVol as cV
import random as r
import polygonTests as pT
import matrix as M
import math as math
import random as rd
#int log h_l/h_l>= L

#P = rP.getRandomPolygon(5)
#Q = rP.getRandomPolygon(5)
# !!!! error in randomPolytop3???!!!! this is no polygon with 6 vertices: [[1.1401725238813245, 0.0], [-1.5758520160340928, -1.2311345299857555], [-118.96995447650002, -10.308197516035044], [-0.4625773828624123, 2.7606389546275225]]

#pT.plotPoly( P , 'r')
#rP.makeUnitVolume( P )
#pT.plotPoly( P , 'r' )

#P = [[132.36947145112424, 220.43230303575757], [-31.257699600609765, -44.27262366537498], [-33.561775200998504, -56.46240797834192], [-34.14982146100093, -59.83445736682929], [-33.40017518851507, -59.86281402521137]]
#Q = [[3.805955300147715, 3.554616667948144], [-0.6020247375518881, 0.8375692779673509], [-2.198015207864401, -0.1537467139313896], [-3.4184482679185253, -0.9317998766221556], [2.4125329131870994, -3.3066393553619498]]
#False
# ist das ein Gegenbeispiel? Es war keinmal ein domain-Fehler bei math.pow( quotient , ... ) dabei.

# integrate over dV_K ( cone volume measure of K ), which means: log h_L / log h_K ...



def logMinTest( orderedPolarPolygon_1 , orderedPolarPolygon_2 ):
    K = rP.getCartesian(orderedPolarPolygon_1)
    L = rP.getCartesian(orderedPolarPolygon_2)
    rP.makeCentered(K)
    rP.makeCentered(L)
    rP.makeUnitVolume(K)
    rP.makeUnitVolume(L)
    coneVol = cV.getConeVol( K )
    #print( coneVol )

    prod = 1

    for i in range( len( K ) ):
        u = [ coneVol[i-1][0] , coneVol[i-1][1] ]
        quotient = rP.supportFunction( orderedPolarPolygon_2 , u ) / M.scal( K[ i - 1 ] , u )
        if quotient == 0:
            print( u )
            print( L[ i-1 ])
            print( K[ i-1 ])
        #print( quotient)
        #print( K[ i - 1] )
        #print( u )
        if quotient >= 0:
            prod = prod * math.pow( quotient , coneVol[ i - 1 ][ 2 ])
        else:
            print( ' bang' )
            print( u )
            print( coneVol[ i - 1][2] )
            pT.plotPoly( K , 'r')
            pT.plotPoly( L , 'g')
    if prod >= 1:
        return True
    else:
        print(K)
        print(L)
        print( 1 - prod)
        return False

#quotients = [ 0.05260231151731214,
#0.8504298083447889,
#0.2878130680022156,
#0.09475793078867864,
#0.1372928559122144 ]
#check = 1
#for number in quotients:
#    check = check * number
#print( 1 - check )

# K = [ [ 1 , 0 ] , [ 0 , 1 ] , [ 0 , -1 ] ]
# L = [ [ 1 , 0 ] , [ 0 , 1 ] , [ 0 , -1 ] ]

def logMinAutoTest( repeats ):
    n = 0
    while n < repeats:
        P = rP.getRandomNoncenteredPolarPolygon( math.floor( rd.random() * 10 ) + 3 )
        Q = rP.getRandomNoncenteredPolarPolygon( math.floor( rd.random() * 10 ) + 3 )

        if not logMinTest( P , Q ):
            print( P )
            print( Q )
            print( 'oh shit' )
            print( 'sollte tester für polytop-eigenschaft benutzen' )
        n = n+1
        print('done')

logMinAutoTest( 20 )

# mögliche Fehler: h wird negativ (support Function), coneVolume wird negativ, u hat norm null... solche Dinge.

# Für P = [[0.43327452367750274, 1.1615678212687457], [2.6048180628406046, 1.2742308939562415], [4.893940117032688, 1.6736464013784849], [5.973152059553502, 362.6363399104863], [5.973152059553502, 362.6363399104863], [5.994497208615119, 107.67142130068721]]
# und Q = [[1.5227088527724186, 1.1028604092237395], [1.5544981710465458, 1.1020959884645107], [2.5678271294458743, 1.0380503276779611], [5.415598275274818, 1.538541840649036], [5.948540515870024, 180.7707730015919]]
# ist aufeinmal ein Normalenvektor 'u' mit Norm 0 ... wie geht das? u wird durch das coneVolume bestimmt



# Für P = [[0.12052685145948405, 0.10016898087234995], [-1.601236023721101, 1.2234795179795563], [-0.007010810061268179, -0.6251879546153744], [0.8160675487328953, -0.4299405664378391], [0.6716524335899903, -0.26851997779869285]]
# und Q = [[1.40777956825808, -3.635020850376281], [1.3914049764315826, -3.5826546263421952], [-5.539665271747164, 14.623890722053702], [1.269343333690422, -3.546059420119178], [1.4711373933670784, -3.860155825216049]]
# in Polarkoordinaten (winkel ist negativ...) scheint logMin falsch zu sein. Aber es gibt manchmal immernoch diesen norm( u ) == 0-Fehler... vielleicht läuft etwas nicht richtig am Programm
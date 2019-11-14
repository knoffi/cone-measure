import randomPolygon3Centered as rP
import coneVol2 as cV
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



def logMinTest( orderedCartPolygon_1 , orderedCartPolygon_2 , bringInPosition , eps_prod ):
    K = orderedCartPolygon_1
    L = orderedCartPolygon_2
    bringInPosition(K)
    bringInPosition(L)
    rP.makeUnitVolume(K)
    rP.makeUnitVolume(L)
    # not sure if makeUnitVolume destroys the centered property... it is clear that it needs it for the normalization of volume
    bringInPosition(K)
    bringInPosition(L)
    coneVol = cV.getConeVol( K )
    #print( coneVol )

    prod = 1

    for i in range( len( K ) ):
        # Reihenfolge der Normalenvektoren sollte jetzt egal sein
        u = [ coneVol[i-1][0] , coneVol[i-1][1] ]
        h_1 = rP.supportFunctionCartesianCentered( L , u )
        h_2 = rP.supportFunctionCartesianCentered( K , u )
        quotient = h_1 / h_2
        if quotient == 0:
            print( 'quotient von logMin war gleich 0...')
            print( u )
            print( L[ i-1 ])
            print( K[ i-1 ])
        #print( quotient)
        #print( K[ i - 1] )
        #print( u )
        if quotient >= 0:
            prod = prod * math.pow( quotient , coneVol[ i - 1 ][ 2 ] )
        else:
            print( ' quotient von logMin war negativ' )
            print( u )
            print( coneVol[ i - 1][2] )
            pT.plotPoly( K , 'r')
            pT.plotPoly( L , 'g')
    if prod + eps_prod >= 1:
        return True
    else:
        #print( 1 - prod)
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

def logMinAutoTest( repeats , bringInPosition, getPosition , eps_position , eps_volume , eps_normals , eps_logMin):
    n = 0
    logMinTrue = 0
    logMinFalse = 0
    logMinTriangleFalse = 0
    while n < repeats:
        P = rP.getRandomNoncenteredPolarPolygon( math.floor( 3 ) )
        Q = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
        K = rP.getCartesian( P )
        L = rP.getCartesian( Q )
        if not logMinTest( K , L , bringInPosition , eps_logMin ):

            if pT.isConvex(K) and  pT.isConvex(L) and pT.isConvexRun(K) and  pT.isConvexRun(L):
                if M.norm( getPosition(K)) +  M.norm( getPosition(L) ) < eps_position:
                    if math.fabs( rP.getVolume(K) - 1 ) + math.fabs( rP.getVolume(L) - 1 ) < eps_volume:
                        # which bound is sufficient for a stable logMin calculation? a bound of 0.00003 seems to be always fulfilled...
                        if pT.getWorsteNormalsDirection(K) < eps_normals:
                            if len(K) == 3 or len(L) == 3:
                                print( 'das darf nicht sein')
                                print(K)
                                print(L)
                                pT.plotPoly(K, 'r')
                                pT.plotPoly(L, 'g')
                                logMinTriangleFalse += 1
                            else:
                                print( 'logMin Gegenbeispiel' )
                                print( K )
                                print( L )
                                pT.plotPoly(K, 'r')
                                pT.plotPoly(L, 'g')
                                logMinFalse += 1
        else:
            logMinTrue += 1

        n = n+1
    print('done')
    print( [ logMinTrue , logMinFalse , logMinTriangleFalse ] )

def logMinRoundAutoTest( repeats , roundMethod, digits , bringInPosition, getPosition , eps_position , eps_volume , eps_normals , eps_logMin):
    n = 0
    logMinTrue = 0
    logMinFalse = 0
    logMinTriangleFalse = 0
    while n < repeats:
        P = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
        Q = rP.getRandomNoncenteredPolarPolygon( math.floor( 4 ) )
        K = rP.getCartesian( P )
        L = rP.getCartesian( Q )
        roundMethod( K , digits )
        roundMethod( L , digits )

        #because of roundig, convexity can be lost and support function can get negativ
        # because of rounding, two edges can be the same.
        # and then norm vector will be [] or something ? Because a type error ouccurs
        while(True):
            try:
                if not logMinTest( K , L , bringInPosition , eps_logMin ):

                    if pT.isConvex(K) and  pT.isConvex(L) and pT.isConvexRun(K) and  pT.isConvexRun(L):
                        if M.norm( getPosition(K)) +  M.norm( getPosition(L) ) < eps_position:
                            if math.fabs( rP.getVolume(K) - 1 ) + math.fabs( rP.getVolume(L) - 1 ) < eps_volume:
                                # which bound is sufficient for a stable logMin calculation? a bound of 0.00003 seems to be always fulfilled...
                                if pT.getWorsteNormalsDirection(K) < eps_normals:
                                    if len(K) == 3 or len(L) == 3:
                                        print( 'das darf nicht sein')
                                        print(K)
                                        print(L)
                                        pT.plotPoly(K, 'r')
                                        pT.plotPoly(L, 'g')
                                        logMinTriangleFalse += 1
                                    else:
                                        print( 'logMin Gegenbeispiel' )
                                        print( K )
                                        print( L )
                                        pT.plotPoly(K, 'r')
                                        pT.plotPoly(L, 'g')
                                        logMinFalse += 1
                else:
                    logMinTrue += 1
                break
            except TypeError:
                break
        n = n+1
    print('done')
    print( [ logMinTrue , logMinFalse , logMinTriangleFalse ] )



logMinAutoTest( 10000 , rP.makeCentered , rP.getCenter , 0.000000001 , 0.000000001 , 0.000000001, 0.5)
#logMinAutoTest( 1000 , rP.makeBaryCentered , rP.getBaryCenter , 0.000000001 , 0.000000001 , 0.5)

# mögliche Fehler: h wird negativ (support Function), coneVolume wird negativ, u hat norm null... solche Dinge.

# Für P = [[0.43327452367750274, 1.1615678212687457], [2.6048180628406046, 1.2742308939562415], [4.893940117032688, 1.6736464013784849], [5.973152059553502, 362.6363399104863], [5.973152059553502, 362.6363399104863], [5.994497208615119, 107.67142130068721]]
# und Q = [[1.5227088527724186, 1.1028604092237395], [1.5544981710465458, 1.1020959884645107], [2.5678271294458743, 1.0380503276779611], [5.415598275274818, 1.538541840649036], [5.948540515870024, 180.7707730015919]]
# ist aufeinmal ein Normalenvektor 'u' mit Norm 0 ... wie geht das? u wird durch das coneVolume bestimmt
P = [[0.6382449192404283, 1.6421863649812602], [0.6861978846456941, 1.6369216410302712], [1.5834756164177117, 2.324465507422115], [1.753926046773701, 2.6213448009357405], [2.6374198407024996, 2.4416720850944116], [4.539079562681944, 2.145404233227454], [4.9663691753868955, 2.045820722537232], [5.202494901347077, 1.9840737483293533], [5.790025528402832, 2.121973857368612], [5.832478556954785, 2.1416589259974934], [6.173265688074134, 2.119135219790364], [6.200373785109342, 2.1034826562411713]]
Q  =[[2.7395254834422404, 1.9239457880688033], [3.3822025830862747, 3.5844757752769514], [5.970781090091909, 1.062412424236289]]

#print( logMinTest( P , Q ) )

K = [[0.3807152635484133, 0.23855106591659164], [0.3653212298547383, 0.25578485447326504], [-0.014746505847977457, 0.6332808078166909], [-0.1461081359300426, 0.7075506993916125], [-0.6331093337976781, 0.2975421240087332], [-0.11460707240921618, -0.6681895893215861], [0.1446550024065254, -0.6291607922348909], [0.2678094494117747, -0.5617998939135442], [0.542085002599915, -0.3430237125832624], [0.559291909951667, -0.3220078981049258], [0.6116594392314314, -0.1165771003047147], [0.6087054994910835, -0.09942801543260288]]
L = [[-0.20084726556302285, 0.5038871568097025], [-1.163452816696621, -0.40039248009544], [1.3643000822596438, -0.10349467671426259]]

#pT.plotPoly( K , 'r')
#pT.plotPoly( L , 'b' )

# Für P = [[0.12052685145948405, 0.10016898087234995], [-1.601236023721101, 1.2234795179795563], [-0.007010810061268179, -0.6251879546153744], [0.8160675487328953, -0.4299405664378391], [0.6716524335899903, -0.26851997779869285]]
# und Q = [[1.40777956825808, -3.635020850376281], [1.3914049764315826, -3.5826546263421952], [-5.539665271747164, 14.623890722053702], [1.269343333690422, -3.546059420119178], [1.4711373933670784, -3.860155825216049]]
# in Polarkoordinaten (winkel ist negativ...) scheint logMin falsch zu sein. Aber es gibt manchmal immernoch diesen norm( u ) == 0-Fehler... vielleicht läuft etwas nicht richtig am Programm
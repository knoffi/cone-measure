import math
import matrix as M
def getMachEps():
    result = 1
    while( 1 + result > 1 ):
        result = result / 2
    return result

def splitAdd( bigNumber , smallNumber):
    big_exp = math.floor( math.log( bigNumber , 2 ) )
    small_exp = math.floor(math.log(smallNumber, 2))
    k = small_exp - big_exp

    result = smallNumber + ( 2**k * bigNumber )

    while( k < 0 ):
        result += ( 2**k * bigNumber )
        k +=  1

    return result

def splitVekAdd( bigVector , smallVector ):
    result = []
    for i in range(len(bigVector)):
        result.append( splitAdd(bigVector[i] , smallVector[i]))
    return result


test_number = 11.6
print(getMachEps())
print( splitAdd( test_number , getMachEps() ) - ( test_number + getMachEps()))
print( splitVekAdd( [ test_number , test_number ] , [ getMachEps() , getMachEps()]))
print( M.addVek( [ test_number , test_number ] , [ getMachEps() , getMachEps() ] ) )
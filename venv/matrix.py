import math

class Matrix:
    rows = [[0]]

    def __init__(self, a):
        self.rows = a

    def size(self):
        m = len(self.rows)
        n = len(self.rows[0])
        return [m, n]

    def list(self):
        i = 0
        if self.rows == []:
            print("empty matrix")
        else:
            while i < len(self.rows):
                print(self.rows[i])
                i += 1

    def scale(self, t):
        i = 0
        j = 0
        m = self.size()[0]
        n = self.size()[1]
        while i < m:
            while j < n:
                self.rows[i][j] *= t
                j += 1
            j = 0
            i += 1

    def copyRow(self, i):
        a = []
        j = 0
        n = self.size()[1]
        while j < n:
            a.append(self.rows[i][j])
            j += 1
        return a

    def copyCol(self, j):
        a = []
        i = 0
        m = self.size()[0]
        while i < m:
            a.append(self.rows[i][j])
            i += 1
        return a

    def copy(self):
        a = []
        if self.rows == []:
            return Matrix([])
        i = 0
        m = self.size()[0]
        while i < m:
            a.append(self.copyRow(i))
            i += 1
        C = Matrix(a)
        return C

    def add(self, B):
        i = 0
        j = 0
        m = self.size()[0]
        n = self.size()[1]
        while i < m:
            while j < n:
                self.rows[i][j] += B.rows[i][j]
                j += 1
            j = 0
            i += 1

    def plus(self, B):
        C = self.copy()
        C.add(B)
        return C

    def mult(self, B):
        c = []
        i = 0
        j = 0
        k = 0
        m = self.size()[0]
        n = B.size()[1]
        while i < m:
            d = []
            while j < n:
                entry = 0
                while k < self.size()[1]:
                    entry += self.rows[i][k] * B.rows[k][j]
                    k += 1
                k = 0
                d.append(entry)
                j += 1
            j = 0
            c.append(d)
            i += 1
        C = Matrix(c)
        return C

    def copyTrans(self):
        a = []
        j = 0
        n = self.size()[1]
        while j < n:
            a.append(self.copyCol(j))
            j += 1
        C = Matrix(a)
        return C

    def switchRow(self, i, j):
        a = self.copyRow(i)
        self.rows[i] = self.rows[j]
        self.rows[j] = a

    def addScaleRow(self, i, j, alpha):
        k = 0
        while k < self.size()[1]:
            self.rows[i][k] = self.rows[i][k] + alpha * self.rows[j][k]
            k += 1

    def scaleRow(self, i, alpha):
        k = 0
        while k < self.size()[1]:
            self.rows[i][k] *= alpha
            k += 1

    def appendRow(self, row):
        self.rows.append(row)

    def deleteRow(self, index):
        temp = []
        i = 0
        while i < self.size()[0]:
            if i != index:
                temp.append(self.rows[i])
            i += 1
        self.rows = temp

    def zsf(self):
        i = 0
        j = 0
        k = 1
        count = 0
        while i < self.size()[0] and j < self.size()[1]:
            count = 0
            k = 1
            if self.rows[i][j] != 0:
                while (k < self.size()[0] - i):
                    self.addScaleRow(i + k, i, -self.rows[i + k][j] / self.rows[i][j])
                    k += 1
                j += 1
                i += 1
            else:
                while (k < self.size()[0] - i and count == 0):
                    if self.rows[i + k][j] != 0:
                        count += 1
                        self.switchRow(i, i + k)
                    k += 1
                if count == 0:
                    j += 1

    def nzsf(self):
        i = 0
        j = 0
        k = 1
        l = 1
        count = 0
        while i < self.size()[0] and j < self.size()[1]:
            count = 0
            k = 1
            if self.rows[i][j] != 0:
                self.scaleRow(i, 1 / self.rows[i][j])
                while (k < self.size()[0] - i):
                    self.addScaleRow(i + k, i, -self.rows[i + k][j])
                    k += 1
                while (l <= i):
                    self.addScaleRow(i - l, i, -self.rows[i - l][j])
                    l += 1
                l = 1
                j += 1
                i += 1
            else:
                while (k < self.size()[0] - i and count == 0):
                    if self.rows[i + k][j] != 0:
                        count += 1
                        self.switchRow(i, i + k)
                    k += 1
                k = 1
                if count == 0:
                    j += 1

    def nzsfHeads(self):
        i = 0
        j = 0
        k = 1
        l = 1
        count = 0
        a = []
        while i < self.size()[0] and j < self.size()[1]:
            count = 0
            k = 1
            if self.rows[i][j] != 0:
                a.append([i, j])
                self.scaleRow(i, 1 / self.rows[i][j])
                while (k < self.size()[0] - i):
                    self.addScaleRow(i + k, i, -self.rows[i + k][j])
                    k += 1
                while (l <= i):
                    self.addScaleRow(i - l, i, -self.rows[i - l][j])
                    l += 1
                l = 1
                j += 1
                i += 1
            else:
                while (k < self.size()[0] - i and count == 0):
                    if self.rows[i + k][j] != 0:
                        count += 1
                        self.switchRow(i, i + k)
                    k += 1
                k = 1
                if count == 0:
                    j += 1
        return a

    def kernel(self):
        Copy = self.copy()
        heads = Copy.nzsfHeads()
        j = 0
        k = 0
        l = 0
        count = 0
        vector = []
        kernel = []
        temp = 0
        while j < Copy.size()[1]:
            vector = []
            count = 0
            # check whether jth column is headless
            while count == 0 and k <= j and k < len(heads):
                if j == heads[k][1]:
                    count += 1
                k += 1
            k = 0
            # if there is no head append vector as a kernel column, otherwise increment j
            if count == 0:
                # set jth component of vector as 1
                while k < Copy.size()[1]:
                    if k != j:
                        vector.append(0)
                    else:
                        vector.append(1)
                    k += 1
                k = 0
                # if a_kj not equal 0, find head in kth row
                while k < Copy.size()[0]:
                    if Copy.rows[k][j] != 0:
                        l = 0
                        temp = 0
                        # If head is in a_kl, override lth component of vector
                        while l < len(heads):  # and l<j
                            if heads[l][0] == k:
                                temp = heads[l][1]
                            l += 1
                        vector[temp] = -Copy.rows[k][j]
                    k += 1
                k = 0
                kernel.append(vector)
            j += 1
        if (kernel == []):
            return Matrix([])
        else:
            kernelMatrix = Matrix(kernel).copyTrans()
            return kernelMatrix

    def lgs(self, vector):
        i = 0
        copy = self.copy()
        while i < copy.size()[0]:
            copy.rows[i].append(vector[i])
            i += 1
        return copy

    # If vector is in Im(self)     then lgsSolve will  return a correct solution

    # if vector is NOT in Im(self) then lgsSolve will compute a solution for every  component (of vector) <= range(self)
    def lgsSolve(self, vector):
        copy = self.lgs(vector)
        i = 0
        j = 0
        k = 1
        l = 1
        count = 0
        headColumns = []
        while i < self.size()[0] and j < self.size()[1]:
            count = 0
            k = 1
            if copy.rows[i][j] != 0:
                headColumns.append(j)
                copy.scaleRow(i, 1 / copy.rows[i][j])
                while (k < copy.size()[0] - i):
                    copy.addScaleRow(i + k, i, -copy.rows[i + k][j])
                    k += 1
                while (l <= i):
                    copy.addScaleRow(i - l, i, -copy.rows[i - l][j])
                    l += 1
                l = 1
                j += 1
                i += 1
            else:
                while (k < copy.size()[0] - i and count == 0):
                    if copy.rows[i + k][j] != 0:
                        count += 1
                        copy.switchRow(i, i + k)
                    k += 1
                k = 1
                if count == 0:
                    j += 1
        coeff = copy.copyCol(copy.size()[1] - 1)
        k = 0
        solver = []
        while k < self.size()[1]:
            solver.append(0)
            k += 1
        k = 0
        while k < len(headColumns):
            solver[headColumns[k]] = coeff[k]
            k += 1
        return solver

    def image(self, b):
        im = []
        i = 0
        j = 0
        temp = 0
        while i < self.size()[0]:
            j = 0
            temp = 0
            while j < self.size()[1]:
                temp += self.rows[i][j] * b[j]
                j += 1
            im.append(temp)
            i += 1
        return im

    def appended(self, B, b, v, indices):
        i = 0
        Matrix = self.copy()
        vektor = []
        while i < len(v):
            vektor.append(v[i])
            i += 1
        i = 0
        while i < len(indices):
            Matrix.rows.append(B.rows[indices[i]])
            vektor.append(b[indices[i]])
            i += 1
        return [Matrix, vektor]


a = [[0, 1], [1, 0], [2, 4]]
b = [[2, 1, 2, 2], [0, 1, 2, 2], [0, 0, 5, 1]]
E = Matrix([[2, 1], [1, 2], [2, 1]])
D = Matrix([[0, 0, 0, 0, 0, 0], [0, 1, 1, 1, 1, 1], [0, 2, 2, 2, 2, 3], [0, 0, 0, 1, 0, 1]])


def scal(v, w):
    i = 0
    scalar = 0
    while i < len(v):
        scalar += v[i] * w[i]
        i += 1
    return scalar


def checkZero(d):
    i = 0
    check = True
    while i < len(d):
        if d[i] > 0.001 or d[i] < -0.001:
            check = False
        i += 1
    return check


def addVek(v, w):
    i = 0
    d = []
    while i < len(v):
        d.append(v[i] + w[i])
        i += 1
    return d


def subVek(v, w):
    i = 0
    d = []
    while i < len(v):
        d.append(v[i] - w[i])
        i += 1
    return d


# only works, if both vectors have the same length?
def polarAngleSub( v_polar , w_polar ):
    alpha = v_polar[0]
    beta = w_polar[0]
    r_1 = v_polar[1]
    r_2 = w_polar[1]
    if alpha - beta >= 0:
        return alpha - beta
    else:
        return alpha - beta + 2 * math.pi

def scaleVek(alpha, v):
    i = 0
    w = []
    while i < len(v):
        w.append(alpha * v[i])
        i += 1
    return w

def getMidPoint( v , w ):
    result = []
    i = 0
    while i < len( v ):
        result.append( ( v[i] + w[i] ) / 2.0 )
        i = i + 1
    return result

def addScaleVek( v , alpha , w ):
    return addVek( v , scaleVek( alpha , w ) )


# recently, this method calculated the norm squared (without root), and coneVolumeNewton worked
def norm(vector):
    return math.sqrt(vector[0] ** 2 + vector[1] ** 2)

def dist( v , w ):
    diff = subVek( v , w )
    return norm( diff )

def subMatrix( A , B):
    n = A.size()[0]
    result_rows = []
    for i in range(n):
        v = A.rows[i]
        w = B.rows[i]
        result_rows.append( subVek( v , w ) )
    return Matrix(result_rows)

def addMatrix( A , B):
    n = A.size()[0]
    result_rows = []
    for i in range(n):
        v = A.rows[i]
        w = B.rows[i]
        result_rows.append( addVek( v , w ) )
    return Matrix(result_rows)


def distMatrix( A , B ):
    result = 0
    n = A.size()[0]
    m = B.size()[1]
    for i in range(n):
        for j in range(m):
            entry_A = A.rows[i][j]
            entry_B = B.rows[i][j]
            result += ( entry_A - entry_B ) ** 2
    return math.sqrt( result )

A_test = Matrix( [ [1 , 0] , [ 0 , 1] , [ 0 , 0 ] ] )
B_test = Matrix( [ [ 1 , 1 ] , [ -1 , 1 ] , [ 0 , 0 ] ] )

if ( distMatrix( A_test , B_test ) - math.sqrt(2) ) > 0.00001:
    print( ' Fehler bei distMatrix ')

def minLamb(lamb, j):
    temp = lamb[j]
    i = j
    index = j
    while i < len(lamb):
        if lamb[i] < temp:
            temp = lamb[i]
            index = i
        i += 1
    return index


A = Matrix(a)
B = Matrix(b)


# A.list()
# B=A.copyTrans()
# A=A.mult(A) #warum gibt es heir keine exception oder so? A*A ist doch gar nicht defineirt????
# A.list()
# B.list()
# A=A.mult(A)
# x=E.lgsSolve([3,3,3])
# print(x)
# print(E.image([1,0]))
class LEQ:
    eMatrix = Matrix([[1]])
    eVektor = [1]
    qMatrix = Matrix([1])
    qVektor = [1]

    # LEQ is meant for the min problem with function '' 1/2 x^T qMatrix x + q^T x''

    def __init__(self, A, b, Q, q):
        self.eMatrix = A
        self.eVektor = b
        self.qMatrix = Q
        self.qVektor = q

    # min method may not work if eMatrix.kernel() is the empty-matrix.
    # min method gives the correct result if qMatrix is positiv definite. If Q has a kernel, Q.lgsSolve(q) may not work.
    def min(self):

        if self.eMatrix.rows == []:
            minimizer = self.qMatrix.lgsSolve(scaleVek(-1, self.qVektor))
            return [minimizer, []]

        M = self.eMatrix
        y = M.lgsSolve(self.eVektor)
        Z = self.eMatrix.kernel()

        # if Z is empty because eMatrix is injective
        if Z.rows == []:
            minimizer = y
        else:
            tZ = Z.copyTrans()
            q = []
            i = 0
            # construction of -q^~
            while i < Z.size()[0]:
                q.append(-y[i] - self.qVektor[i])
                i += 1
            i = 0
            q = tZ.image(q)
            Q = tZ.mult(self.qMatrix)
            Q = Q.mult(Z)
            z = Q.lgsSolve(q)
            z = Z.image(z)
            minimizer = []
            while i < Z.size()[0]:
                minimizer.append(z[i] + y[i])
                i += 1
            i = 0
        At = self.eMatrix.copyTrans()
        lamb = At.lgsSolve(addVek(self.qMatrix.image(minimizer), self.qVektor))  # warum +??? das muesste fehler ergeben

        return [minimizer, lamb]


# test=LEQ( Matrix([[1,1,0],[0,0,1]]), [2,1] , Matrix([[1,0,0],[0,1,0],[0,0,1]]),[0,0,0])
# print(test.min())


class LIQ:
    eMatrix = Matrix([[1]])
    iMatrix = Matrix([1])
    eVektor = [1]
    iVektor = [1]
    qMatrix = Matrix([1])
    qVektor = [1]

    # aktive=[1]

    def __init__(self, A, G, b, h, Q, q):
        self.eMatrix = A
        self.iMatrix = G
        self.eVektor = b
        self.iVektor = h
        self.qMatrix = Q
        self.qVektor = q

    # self.aktive=[]
    # active indices lieber als SET ausgeben lassen??
    def activI(self, x):
        activx = []
        i = 0
        j = 0
        temp = 0
        while i < self.iMatrix.size()[0]:
            j = 0
            temp = 0
            while j < self.iMatrix.size()[1]:
                temp += self.iMatrix.rows[i][j] * x[j]
                j += 1
            if temp == self.iVektor[i]:
                activx.append(i)
            i += 1
        return activx

    def lambdCheck(self, lambd):
        i = len(self.eVektor)
        check = True
        while i < len(lambd):
            if lambd[i] < 0:
                check = False
            i += 1
        return check

    def stepSize(self, x, d):
        size = 1
        i = 0
        while i < self.iMatrix.size()[0]:
            v = self.iMatrix.rows[i]
            # k=0
            # scal=0
            # while k<len(d):
            #    scal+=d[k]*self.iMatrix.rows[i][k]
            #    k+=1
            temp1 = scal(d, v)
            if temp1 < 0:
                temp2 = (self.iVektor[i] - scal(x, v)) / temp1
                if size > temp2:
                    size = temp2

            i += 1

        return size

    def activSet(self, start):
        if self.eVektor == [] and self.iVektor == []:
            return self.qMatrix.lgsSolve(scaleVek(-1, self.qVektor))
        x = start
        active = self.activI(x)  # active falsch fuer start=(0.5,0.5)
        C = self.eMatrix.appended(self.iMatrix, self.iVektor, self.eVektor, active)  # laeuft
        subProb = LEQ(C[0], C[1], self.qMatrix, self.qVektor)
        subResult = subProb.min()
        temp = subResult[0]
        d = subVek(temp, x)  # d=(0,0) wenn start=(0.5 , 0.5)
        lamb = subResult[1]
        while (not self.lambdCheck(lamb)) or (not checkZero(d)):
            print(x)
            if not checkZero(d):
                step = self.stepSize(x, d)
                x = addVek(x, scaleVek(step, d))
                active = self.activI(x)

                C = self.eMatrix.appended(self.iMatrix, self.iVektor, self.eVektor, active)
                subProb = LEQ(C[0], C[1], self.qMatrix, self.qVektor)
                subResult = subProb.min()
                temp = subResult[0]
                d = subVek(temp, x)
                lamb = subResult[1]
            else:
                index = minLamb(lamb, len(self.eVektor))
                subProb.eMatrix.deleteRow(index)
                subResult = subProb.min()
                temp = subResult[0]
                d = subVek(temp, x)
                lamb = subResult[1]
        print(x)
        return x


# THERE WERE SOME COMBINATION WHERE I GOT AN ENDLESS LOPE, BUT IAM NOT SURE IF I ENTERED WRONG VALUES OR IF IT IS A MISTAKE IN THE CODE/ALGORITHM!!!!!


# test2=LIQ( Matrix([[1,1,0],[1,0,0]]) , Matrix([[1,1,0],[0,0,1],[0,2,0]]), [1,1],  [2,1,2] , Matrix([[1,0,0],[0,1,0],[0,0,1]]),[0,0,0] )
# print(test2.activI([1,1,0]))
A = Matrix([[0, 0]])
G = Matrix([[1, 1], [1, -1], [-1, -1], [-1, 1]])
b = [0]
h = [3, -3, -9, -3]
# A=Matrix([])
# G=Matrix([[-1,-1],[-1,1],[1,-1],[-1,-1]])
# b=[]
# h=[-1,-1,-1,-1]
# Q=Matrix([[1,0],[0,1]])
# q=[-2,3]
Q = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
q = [-2, -3, -4]


# test2=LIQ(A,G,b,h,Q,q)

# minimizer=test2.activSet([1,0,0])#scheint zu laufen, ausser wenn A eine empty matrix ist und Q=Id und start bei (0.5,0.5)
# print(minimizer)
# print(G.lgsSolve(h))


def idMatrix(n):
    i = 0
    j = 0
    Id = []
    temp = []
    while i < n:
        j = 0
        while j < n:
            if i == j:
                temp.append(1)
            else:
                temp.append(0)
            j += 1
        Id.append(temp)
        temp = []
        i += 1
    Id = Matrix(Id)
    return Id


# get nearestpPoint to p on the polytop defined by the equalities (A,b) and inequailities (G,h), start is a point of the polytop
def nearestPoint(A, b, G, h, p, start):
    if b == [] and h == []:
        return p

    l = len(p)
    Id = idMatrix(l)
    condMatrix = []
    condVektor = []
    q = []
    bTranslated = []
    hTranslated = []
    i = 0
    while i < l:
        q.append(0)
        i += 1
    i = 0
    if b != []:
        while i < len(b):
            bTranslated.append(b[i] - scal(A.rows[i], p))
            condVektor.append(bTranslated[i])
            condMatrix.append(A.rows[i])
            i += 1
    i = 0
    if h != []:
        while i < len(h):
            hTranslated.append(h[i] - scal(G.rows[i], p))
            condVektor.append(hTranslated[i])
            condMatrix.append(G.rows[i])
            i += 1
    i = 0
    # start may work even when the lgs is not solvable. lgsSolve gives solution for all linearly independent components
    start = addVek(start, scaleVek(-1, p))  # Matrix(condMatrix).lgsSolve(condVektor)
    problem = LIQ(A, G, bTranslated, hTranslated, Id, q)
    nearPoint = problem.activSet(start)
    nearPoint = addVek(nearPoint, p)
    return nearPoint


def nearestPoint2(A, b, G, h, p, start):
    # start punkt einfach finden, indem man den normalenvektor einer facette an den Randpunkt einer Facette bildet und so lange verkuerzt, bis alle Bedingung erfuellt sind? Aber dann muesste man ja schon einen Randpunkt kennen!
    Q = idMatrix(len(p))
    q = scaleVek(-1, p)
    problem = LIQ(A, G, b, h, Q, q)
    return problem.activSet(start)


# p = [-2, 5]
# print('The closest point of the given intersection of half spaces to', p, 'ist', nearestPoint2(A, b, G, h, p, [6, 3]), '.')

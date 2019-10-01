# sphinx_gallery_thumbnail_number = 3
from cmath import sqrt
import math
import numpy
import random


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


def scaleVek(alpha, v):
    i = 0
    w = []
    while i < len(v):
        w.append(alpha * v[i])
        i += 1
    return w

# recently, this method calculated the norm squared (without root), and coneVolumeNewton worked
def norm(vector):
    return math.sqrt(vector[0] ** 2 + vector[1] ** 2)


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

# only works in R^2
def nextVertex(u, volume, point):
    orth = []
    # im uhrzeigersinn um 90 grad drehen
    orth.append(u[1])
    orth.append(-u[0])

    # assumes that 0 lies in K (i.e. point is not orthogonal)
    d = scaleVek(2 * volume / scal(point, u), orth)
    temp = addVek(point, d)
    point.clear()
    point.append(temp[0])
    point.append(temp[1])


def getNextVertex(u, volume, point):
    result = []
    nextVertex(u, volume, point)
    for x in point:
        result.append(x)
    return result

def domain(start, end, step):
    result = []
    i = 0
    while (start + i * step <= end):
        result.append(start + i * step)
        i += 1
    return result

# IMPORTANT: does nextVertex work with local variable rules of python?
def coneVolIterator(coneData, point):
    result = []
    for facet in coneData:
        nextVertex([facet[0], facet[1]], facet[2], point)
    for x in point:
        result.append(x)
    return result


def getConeVolIterator(coneData, point):
    result = []
    for x in point:
        result.append(x)
    coneVolIterator(coneData, result)
    return result


# berechnet abstand zum quadrat
def polyFunctional(coneData, point):
    return norm(subVek(point, getConeVolIterator(coneData, point))) ** 2


def arrayPolyFunctional1(coneData, point, param):
    result = []
    for i in param:
        result.append(polyFunctional(coneData, scaleVek(i, point)))
        if polyFunctional(coneData, scaleVek(i, point)) < eps_graph:
            print(i)
    return result


def arrayPolyFunctional2(coneData, param):
    result = []
    for i in param:
        result.append(polyFunctional(coneData, [1, i]))
    return result

def gradGetNextVertex(u, volume, point):
    B = Matrix([[u[1]], [-u[0]]])
    A = Matrix([[-u[0], -u[1]]])
    M = B.mult(A)
    c = scal(point, u)
    c =2*volume / (c**2)
    M.scale(c)
    # Matrix.plus oder Matrix.add verwenden?
    M.add(idMatrix(2))
    return M.copyTrans()


def diffGetNextVertex(u, volume, point):
    B = Matrix([[u[1]], [-u[0]]])
    A = Matrix([[-u[0], -u[1]]])
    M = B.mult(A)
    c = scal(point, u)
    c =2*volume / (c**2)
    M.scale(c)
    # Matrix.plus oder Matrix.add verwenden?
    M.add(idMatrix(2))
    return M
#gradGetNextVertex([1,1],1,[1,1]).list()


def gradConeVolumeIterator(coneData, point):
    result = idMatrix(2)
    for facet in coneData:
        result = gradGetNextVertex([facet[0], facet[1]], facet[2], point).mult(result)
        nextVertex([facet[0], facet[1]], facet[2], point)
    return result.copyTrans()

# Does the for-loop start with the correct facet? It depends on coneData...
def diffConeVolumeIterator(coneData, point):
    result = idMatrix(2)
    for facet in coneData:
        result = diffGetNextVertex([facet[0], facet[1]], facet[2], point).mult(result)
        nextVertex([facet[0], facet[1]], facet[2], point)
    return result


def gradPolyFunctional(coneData, point):
    M = idMatrix(2)
    A = gradConeVolumeIterator(coneData, point)
    A.scale(-1)
    M.add(A)

    b = subVek(point, getConeVolIterator(coneData, point))
    b = scaleVek(2, b)

    result = M.image(b)
    return result


def gamma(coneData, params):
    d_0 = [coneData[len(coneData) - 1][1], -coneData[len(coneData) - 1][0]]
    d_1 = [-coneData[0][1], coneData[0][0]]
    return addVek(scaleVek(params[0], d_0), scaleVek(params[1], d_1))


def gradGamma(coneData, params):
    d_0 = [coneData[len(coneData) - 1][1], -coneData[len(coneData) - 1][0]]
    d_1 = [-coneData[0][1], coneData[0][0]]
    return Matrix([d_0, d_1]).copyTrans()


def diffGamma(coneData, params):
    d_0 = [coneData[len(coneData) - 1][1], -coneData[len(coneData) - 1][0]]
    d_1 = [-coneData[0][1], coneData[0][0]]
    return Matrix([d_0, d_1])


def phi(params,cD):
    return polyFunctional(cD, gamma(cD,params))

# is this the correct gradient?
def gradPhi(params,cD):
    M = gradGamma(cD, params)
    v = gradPolyFunctional(cD, gamma(cD, params))
    return M.image(v)


def stepSize(coneData, point, d):
    alpha = const_step
    while (phi(point, coneData) <= phi(addVek(point, scaleVek(alpha, d)),coneData)):
        alpha = alpha * const_step
    return alpha


# works only if global minimum is zero and if its the only extremum
def coneVolumeDescend(coneData, start):
    result = start
    while (phi(result,coneData) > eps_descend):
        d = scaleVek(-1, gradPhi(result,coneData))
        d= scaleVek(1/norm(d),d)
        alpha = stepSize(coneData, result, d)
        result = addVek(result, scaleVek(alpha, d))
        print(d)
    return result


def g(a, x):
    return a * x[0] ** 2


def gradG(x):
    return [2 * x[0]]

def f(cD,params):
    return subVek(getConeVolIterator(cD,gamma(cD,params)),gamma(cD,params))

def diff(cD, params):
    A=diffGamma(cD,params)
    B=diffConeVolumeIterator(cD,gamma(cD,params))
    C=idMatrix(2)
    C.scale(-1)
    B.add(C)
    return B.mult(A)

def regulator(params):
    return norm(f(cD,params))*norm(params)

def coneVolumeNewton(cD,params):
    value=float("inf")
    while(regulator(params)>eps_newton):
        print(regulator(params))
        v=f(cD,params)
        M=diff(cD,params)
        params=subVek(params,M.lgsSolve(v))
        #print(value)

    return params

# es ist wichtig, dass u_1 in der richtigen reihenfolge ist und als erster rechts von p steht
# p=[1,10] hat eine sehr interessanten Minimizer auf spann([1,10])

#x = domain(0.7, 19, 0.03)
#y = arrayPolyFunctional2(cD, x)
# IMPORTANT: algorithm throws exception if one of the p_i is orthogonal to a u_i.
#print(coneVolumeNewton(cD,p))

def getGrid(n):
    result=[]
    l=0
    k=0
    while(l<=n**2):
        while(k<=n**2):
            if(k+l>0):
                result.append([k*1/n,l*1/n])
            k=k+1
        k=0
        l=l+1
    return result

# result will be like [[argmin_x, argmin_y], min]!!!
def min(grid):
    result= [[0,0],float("inf")]
    for pair in grid:
        while(True):
            try:
              if(result[1]>psi(pair)):
                result=[pair, psi(pair)]
              break
            except ZeroDivisionError:
                #print('zero division at',pair)
                break
    return result


def startPoint():
    n=4
    result=[[0,0],float("inf")]
    while(result[1]>min_grid):
        grid=getGrid(n)
        result=min(grid)
        if(result[1]==0):
            break
        n=n+1
    print(result)
    return result


def psi(params):
    return norm(f(cD,params))


#print(coneVolumeNewton(cD,startPoint()[0]))

# interesting result                 cD                                          result
#            [[1, 1, 1], [0, -1, 1], [-1, 0, 1],[-1,-1,1],[-1,-2,1]]    [255253327.86799613, 357518606.484159]
#

def getClockRotation(vector):
    return[vector[1],-vector[0]]


# exspect vertices to be ordered clockwise
def getConeVol(vertices):
    result=[]
    n=len(vertices)
    a=0
    b=1
    v= subVek(vertices[a],vertices[b])
    for i in range(n):
        u_tilde=getClockRotation(v)
        u=scaleVek(1/norm(u_tilde),u_tilde)

        height= scal(u,vertices[a])
        facetVolume= norm(v)
        coneVolume= 0.5*height * facetVolume

        result.append([u[0],u[1],coneVolume])

        a = (a+1) % n
        b =( b + 1 ) % n
        v = subVek(vertices[a], vertices[b])
    print('check_getConeVolume')
    return result
# ziemlich unsymmetrisches fünfeck, 5. und 3. facette sind parallel
min_grid=0.2
eps_newton=0.0000000000001
eps_graph= 0.0000001
const_step= 0.1
eps_descend = 0.1
cD = [[1, 1, 1], [0, -1, 1], [-1, 0, 1],[-1,-1,1],[-1,-2,1]]
p=[1,0.7]
# cD2 wird gar nicht bei start verwendet, sondern cD....
# cD2=getConeVol([[1,1],[1,-1],[-1,-1],[-1,1]])
# print(coneVolumeDescend(cD,startPoint()[0]))

def getRandomPolygon(n):
    result=[]
    steps = 0
    while( steps < n ):
        result.append(getNextPoint(result))
        steps = steps + 1
    return result

# better distribution: choose all angles randomly in [0,pi] then sort them. After that choose radii accordingly
# better distribution: do not set angle of first point as 0
def getNextPoint(vertices):
    if(len(vertices) == 0):
        beta = 0
        radius = randomRadius( 0 , math.inf )
        return polar( beta , radius )

    alpha= getAngle(vertices[len(vertices)-1])
    gamma= getAngle(vertices[0])
    beta= randomAngle(alpha,gamma)

    if( len(vertices) == 1 ):
        beta = 2 * math.pi * ( 1 - random.random() )

    max = maxRadius(vertices, beta)
    min = minRadius(vertices , beta)
    r = randomRadius( min , max )
    # print( beta / math.pi * 180 )
    # print(max)
    # print(min)
    result=polar(beta,r)
    #while(noCondition(vertices,result)):
    #    r=randomRadius()
    #    beta=randomAngle(alpha,gamma)
    #    result=polar(beta,r)
    return result

def getAngle(point):
    x=point[0]
    y=point[1]
    if x  == 0 and y > 0:
        return math.pi
    if x == 0 and y < 0:
        return 1.5*math.pi
    if x > 0 and y > 0:
        return math.atan(y/x)
    if x > 0 and y < 0:
        return math.atan(y/x) + 2* math.pi
    if x < 0:
        return math.atan(y/x) + math.pi
    return 0

#better distribution: select random n random angles between 0 and pi, then order them and return an array with every angle
# berücksichtigt das hier die periodizität, also pi = 3*pi
def randomAngle( a , b):
    return a + random.random()*(b-a)

# has to be improved, so that it will give a good distributed number between 0 and +infty
def randomRadius(min, max):
    #randomRadius will always be above 1. Is that good?
    if( max == math.inf ):
        return 1/random.random() + min

    return min + ( max - min ) * random.random()

def polar(angle, radius):
    x= radius * math.cos(angle)
    y= radius * math.sin(angle)
    return[x,y]

def maxRadius(vertices,angle):
    result= math.inf
    # improvable by noticing: it is not necessary to check all previous vertices
    others = vertices.copy()
    for v in vertices :
        # remove(x) removes an/any entry which equals x
        others.remove(v)
        for w in others:
            if(getMax(v,w,angle)<result):
                result=getMax(v,w,angle)
    return result

def minRadius(vertices,angle):
    result=0
    # improvable by noticing: it is not necessary to check all previous vertices
    others = vertices.copy()
    for v in vertices :
        # remove(x) removes an/any entry which equals x
        others.remove(v)
        for w in others:
            if( getMin(v,w,angle) > result ):
                result=getMin(v,w,angle)
    return result

def getMax(v,w,angle):
    u_1 = polar(angle,1)
    u_2 = subVek(v,w)
    U= Matrix([u_1,u_2]).copyTrans()
    result=U.lgsSolve(v)
    # my lgsSolver returns a 'solution' even if v is not in im(M) . Therefore the if-check needs to get adjusted
    if(len(result)==0):
        return math.inf
    if( result[0] <= 0 ):
        return math.inf
    return result[0]

def getMin(v,w,angle):
    u_1= polar(angle,1)
    u_2= subVek(v,w)
    U=Matrix([u_1,u_2]).copyTrans()
    result=U.lgsSolve(v)
    if( len(result)==0 ):
        return 0
    if ( result[0] <= 0):
        return 0
    return result[0]



v=[0,1]
w=[1,0]
angle= math.pi/4
vertices=[ [ 1, 0 ],  [ 0 , -1 ] , [ -1 , 0 ] ]
# ein result mit zwei gleichen Punkten: [[1.3506563545534096, 0.0], [1.8164857151592047, 0.0], [1.3506563545534096, 0.0], [1.4946652011162949, 0.0]]
print(getRandomPolygon(4))
# minRadius testen und einbauen in getNextPoint von getRandomPolygon()
# minRadius ist notwendig, aber vllt muss man bei maxRadius und minRadius nur die Vorgänger überprüfen ... -> mache 2 versionen




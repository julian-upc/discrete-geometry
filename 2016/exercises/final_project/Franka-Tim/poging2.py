import itertools
import numpy
import copy
import random


def checkSigns(signList): #This function scans the input line from left to right and checks how many times the sign switches
    switch = 0
    for i in [1..5]: #Standard length 6
        if signList[i] != signList[i-1]: #Whenever the sign is not the same as the previous one, add 1 to the variable 'switch'
            switch += 1
    return switch



def getComplement(i,j,k,l):
    compList = [0,1,2,3,4,5,6,7]
    compList.remove(i)
    compList.remove(j)
    compList.remove(k)
    compList.remove(l)
    return compList

##### Checks if two line segments intersect #####
def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

##### Checks if point is in triangle #####
def sign(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def PointInAABB(pt, c1, c2):
    return c2[0] <= pt[0] <= c1[0] and c2[1] <= pt[1] <= c1[1]

def PointInTriangle(pt, v1, v2, v3):
    b1 = sign(pt, v1, v2) <= 0
    b2 = sign(pt, v2, v3) <= 0
    b3 = sign(pt, v3, v1) <= 0

    return ((b1 == b2) and (b2 == b3)) and PointInAABB(pt, map(max, v1, v2, v3), map(min, v1, v2, v3))

##### Checks if Gale Diagram corresponds to neighborly polytope #####
#input: galeDiagram = [[x0,y0,s0],[x1,y1,s1],...,[x7,y7,s7]] where xi, yi, si are resp. x-coordinate, y-coordinate, sign of the i-th point

def checkNB(galeDiagram):
    #Remove every possible combination of 4 points and check the amount of +'s
    NB = True
    counter = 0
    for i, j, k, l in itertools.product([0..4], [1..5], [2..6], [3..7]):
        if i != j and i!= k and i != l and j != k and j != l and k != l:
            pointList1 = copy.copy(galeDiagram)  #Copy the list, otherwise the input list changes as well
            del pointList1[l]                    #pointList1 is the pointlist with 4 points removed from the input
            del pointList1[k]
            del pointList1[j]
            del pointList1[i]
            #res = checkSigns(pointList1)
            res = countPositives(pointList1)[0]
            if res == 0 and res == 4:
                NB = false
                #compList = getComplement(i,j,k,l)
                #print(compList)
            elif res == 1:
                locPlus = countPositives(pointList1)[1]
                locMin = countPositives(pointList1)[2]
                NB = PointInTriangle(pointList1[locPlus[0]],pointList1[locMin[0]],pointList1[locMin[1]],pointList1[locMin[2]])
                #print(pointList1[locPlus[0]],pointList1[locMin[0]],pointList1[locMin[1]],pointList1[locMin[2]])
            elif res == 2:
                locPlus = countPositives(pointList1)[1]
                locMin = countPositives(pointList1)[2]
                NB = intersect(pointList1[locPlus[0]],pointList1[locPlus[1]],pointList1[locMin[0]],pointList1[locMin[1]])
                #print(pointList1[locPlus[0]],pointList1[locPlus[1]],pointList1[locMin[0]],pointList1[locMin[1]])
    #print(counter)
    return NB



def countPositives(SubgaleDiagram): #counts the number of points with postive sign & allocates the positive and negative signes
    counter = 0
    plusList = list()
    minList = list()
    if len(SubgaleDiagram) != 4:
        print("Input length does not equal 4.")
    else:
        for i in [0..len(SubgaleDiagram)-1]:
            dummy = SubgaleDiagram[i]
            if dummy[2] == '+':
                plusList.append(i)
                counter += 1
            else: minList.append(i)
    return counter, plusList, minList

signs = ['+','-']
allCombinations = [p for p in itertools.product(signs, repeat=8)] #Generate all combinations of + and - of length 8

allCombinations

#len(allCombinations) #Check that there are indeed 2^8 = 256 combinations

#def checkTekens(): #For all possible combinations, check if they are neighbourly
#    for i in [0..255]:
#        print(i)
#        checkNB(list(allCombinations[i]))

#checkTekens()

##### Checks if Gale Diagram is Simplicial #####

def checkSimplicial(galeDiagram):
    isSimplicial = True
    for i in [0..len(galeDiagram)-3]:
        for j in [i+1..len(galeDiagram)-2]:
            for k in [j+1..len(galeDiagram)-1]:
                #print(galeDiagram[i],galeDiagram[j],galeDiagram[k]) #Check that the functions indeed checks all combinations
                print(galeDiagram[i],galeDiagram[j],galeDiagram[k])
                if checkCollinear(galeDiagram[i],galeDiagram[j],galeDiagram[k]) == True:
                    isSimplicial = False
                    #We want to break the loop here but Python doesn't allow us to do so :'(
                    #This means that we have to check all 8 choose 3 = 56 combinations
    return isSimplicial

def checkCollinear(i,j,k):  #Check for every combination of 3 points whether or not they are collinear
    a = i[0]
    b = i[1]
    c = j[0]
    d = j[1]
    e = k[0]
    f = k[1]
    M = matrix(3, 3, [1, a, b , 1,c,d, 1,e,f])
    #print(M)
    if numpy.linalg.det(M) == 0:
        return True
    else: return False


##### testlijsten #####
lijst1 = [[1,1,'-'],[2,2,'+'],[3,4,'+'],[4,5,'+'],[5,6,'-'],[4,0,'+'],[6,1,'-'],[3,5,'+']]

lijst2 = [[1,1,'-'],[2,2,'+'],[3,4,'+'],[4,5,'+']]

##### Generate random Gale Diagrams #####
#Generate all combinations of + and - of length 8
signs = ['+','-']
allCombinations = [p for p in itertools.product(signs, repeat=8)]

# Count the number of '+' signs in a given combination
# A bit unnecessary because we already had a function for this made for input length 4..

def combCounter(combo):
    counter = 0
    if len(combo) != 8:
        print("Input length does not equal 8.")
    else:
        for i in [0..len(combo)-1]:
            dummy = combo[i]
            if dummy == '+':
                counter += 1
    return counter

#Filter out all indices of combinations that have 3 or 4 '+' signs, because we know a Gale Diagram corresponding to a neighborly polytope has 3 or 4 '+' (or '-') signes

def checkCombinations():
    combList = list()
    for i in [0..255]:
        cur = allCombinations[i]
        if (combCounter(cur) == 3) or (combCounter(cur) == 4):
            combList.append(i)
    return combList


goodCombs = checkCombinations()
len(goodCombs)

#Generate the coordinates at random and append a correct signlist (index as input) to generate a Gale Diagram

def generateRandomGale(integer):
    randomGD = list()
    goodCombs = checkCombinations()
    randomGD.append([random.randrange(-5, 5) for _ in range(0, 8)])
    randomGD.append([random.randrange(-5, 5) for _ in range(0, 8)])
    k = goodCombs[integer]
    randomGD.append(list(allCombinations[k]))
    return randomGD

#Transpose matrix
uitkomst = generateRandomGale(4) #Example to check if it works
randomGDt = map(list,map(None,*uitkomst))

#generate a random coordinate and a random correct signlist to generate a random gale diagram

def mainFunction():
    neighbourlyGD = list()
    for i in [0..5]: #Generate 500 random points
        rand = random.randint(0,125) #Pick a random sign configuration
        uitkomst = generateRandomGale(rand) #Generate corresponding matrix
        randomGDt = map(list,map(None,*uitkomst)) #Transpose matrix
        print(randomGDt)
        if checkSimplicial(randomGDt) == True: #Somehow we managed to mess up the input form for checkCollinear, which is used in checkSimplicial. It's probably something with () instead of [] (or the other way around) but we haven't been able to fix this error in time.
            if checkNB(randomGDt) == True:
                neighbourlyGD.append(randomGDt)
            else:
                break
        else:
            break
    return neighbourlyGD



##### Vertex Facet Structure from Gale Diagram & Check equivalence #####
# Here we would construct the vertex facet structure from the set of Gale Diagrams corresponding to neighborly and simplicial polytopes. 
# With this Vertex Facet Structure we could check which Gale Diagrams correspond to polytope which are combinatorically equivalent.
# Now we would be able to find 3 classes of polytopes: N_(8), N^(*)_(8), C_(4)(8). (See http://ac.els-cdn.com/S0021980067800553/1-s2.0-S0021980067800553-main.pdf?_tid=d64f88c4-defc-11e6-96d8-00000aab0f02&acdnat=1484908999_94ff4a6a8058689c6eaa70d77ef58ee4)

##### From Gale Diagram to Polytope #####
# Here we would construct the point configuration from the 3 gale diagrams we got in the previous step.
# We started building this inverse Gale transformation.

#def insert_row(M,k,row): #Insert a new row in a matrix
#    return matrix(M.rows()[:k]+[row]+M.rows()[k:])

#def makeGaleMatrix(pointList): #Given the list of points, calculate the Gale diagram
#    L = list()
#    for i in [0..len(pointList)-1]:
#        L.append(pointList[i].c) #Get the coordinates and add them
#    tussen0 = insert_row(matrix(L).transpose(),4,[1,1,1,1,1,1,1,1]) #Add the all one row to get homogeneous coordinates
#    tussen = tussen0.rref() #Apply Gaussian elimination such that we get I|B
#    res = tussen[:,[5,6,7]] #Select B
#    res1 = insert_row(res,5,[-1,0,0]) #Add rows such that we get [B|-I]^T
#    res2 = insert_row(res1,6,[0,-1,0])
#    res3 = insert_row(res2,7,[0,0,-1])
#    return res3











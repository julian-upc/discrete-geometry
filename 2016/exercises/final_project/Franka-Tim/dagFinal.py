# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Dear Julian,
#
# We had a hard time with this assignment. Even though Simon was kind enough to help us a little bit out with Python,
# we unfortunately did not succeed to complete the assignment, although we put much effort and time in trying to do so. 
# We feel sorry for handing in such an incomplete assignment. Hopefully we were at least working in the right direction.
# Please let us know if there is anything else we could do... 
#
# Mit freundlichem gruß,
#
# Franka & Tim


import itertools

class Point: #We use this in the functions to make the Gale diagram
    def __init__(self,coordinates,sign):
        self.c = coordinates
        self.s = sign

def insert_row(M,k,row): #Insert a new row in a matrix
    return matrix(M.rows()[:k]+[row]+M.rows()[k:])

def makeGaleMatrix(pointList): #Given the list of points, calculate the Gale diagram
    L = list()
    for i in [0..len(pointList)-1]:
        L.append(pointList[i].c) #Get the coordinates and add them
    tussen0 = insert_row(matrix(L).transpose(),4,[1,1,1,1,1,1,1,1]) #Add the all one row to get homogeneous coordinates
    tussen = tussen0.rref() #Apply Gaussian elimination such that we get I|B
    res = tussen[:,[5,6,7]] #Select B
    res1 = insert_row(res,5,[-1,0,0]) #Add rows such that we get [B|-I]^T
    res2 = insert_row(res1,6,[0,-1,0])
    res3 = insert_row(res2,7,[0,0,-1])
    return res3

pointList = [Point([-1,0,0,2],0),Point([1,0,0,3],0),Point([0,-1,0,4],0),Point([0,1,0,5],0),Point([0,0,1,8],0),Point([0,0,-1,5],0),Point([1,0,1,1],0),Point([1,0,1,4],0)] #Example of a point list
test = makeGaleMatrix(pointList) #Check if the function works

def getSign(pointList,row): #Given a certain row of a point matrix, output the sign vector that corresponds to it
    S = list()
    for i in [0..7]:
        if pointList[i].c[row] > 0:
            S.append('+')
        elif pointList[i].c[row] == 0:
            S.append('0')
        else:
            S.append('-')
    return S

test2 = getSign(pointList,0)

def checkSigns(signList): #This function scans the input line from left to right and checks how many times the sign switches
    switch = 0
    for i in [1..5]: #Standard length 6
        if signList[i] != signList[i-1]: #Whenever the sign is not the same as the previous one, add 1 to the variable 'switch'
            switch += 1
    return switch

import copy

import itertools

def checkNB(signList):
    #Remove every possible pair of points and check if the number of sign switches per subsequence is 2
    NB = true
    count = 0
    while (NB == true and count < 49): #It has to stop after at most 7x7 iterations
        for i, j in itertools.product([0..6], [1..7]):
            if j != i:
                pointList1 = copy.copy(signList) #Copy the list, otherwise the input list changes as well
                del pointList1[j]
                del pointList1[i]
                res = checkSigns(pointList1)
                if res != 2:
                    NB = false
                    break
        count += 1
    print(NB)
    return NB

signs = ['+','-']
allCombinations = [p for p in itertools.product(signs, repeat=8)] #Generate all combinations of + and - of length 8

allCombinations

len(allCombinations) #Check that there are indeed 2^8 = 256 combinations

def checkTekens(): #For all possible combinations, check if they are neighbourly
    for i in [0..255]:
        print(i)
        checkNB(list(allCombinations[i]))

checkTekens()

# # # # # # # 
# What we wanted to do after this:
# Delete two more points from every list with +'s and -'s. Check if the intersection of the positive and negative convex hull 
# empty is or not. Now you know which lists are neighborly. Get the corresponding pointconfiguration. 
# This is possible by calculating the kernels and consider the different signs. You can check which points form a face. 
# Use the Altshuler determinant to determine which polytopes are equivalent up to isomorphy.
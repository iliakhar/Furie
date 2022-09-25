import math
import matplotlib.pyplot as plt

def dft(coord, isForward):
    origCoord = [[]]
    origCoord.append([])
    divider = len(coord[0])
    numSign = -1
    if isForward == False:
        divider = 1
        numSign = 1
    for k in range(len(coord[0])):

        re = 0.0
        im = 0.0
        for n in range(len(coord[0])):
            x = coord[0][n]
            y = coord[1][n]
            arg = 2*math.pi*k*n/len(coord[0])
            reCos = math.cos(arg)
            imSin = math.sin(arg)
            re += (x*reCos-y*imSin)
            im += (x*imSin + y*reCos)*numSign
        origCoord[0].append(re/divider)
        origCoord[1].append(im/divider)
    return origCoord

###############################################
def dftA1(data, p1, p2, isForward):
    furieA1=[[]]
    furieA1.append([])
    divider = p1
    numSign = -1
    if isForward == False:
        divider = 1
        numSign = 1
    for j2 in range(p2):
        for k1 in range(p1):
            re = 0.0
            im = 0.0
            for j1 in range(p1):
                x = data[0][j2 + p2*j1]
                y = data[1][j2 + p2*j1]
                arg = 2 * math.pi *j1*k1/p1
                reCos = math.cos(arg)
                imSin = math.sin(arg)
                re += (x*reCos-y*imSin)
                im += (x*imSin + y*reCos)*numSign
            furieA1[0].append(re/divider)
            furieA1[1].append(im/divider)
    return furieA1


def dftA2(data, p1, p2, isForward):
    furieA2 = [[]]
    furieA2.append([])
    furieA1 = dftA1(data, p1, p2, isForward)
    divider = p2
    numSign = -1
    if isForward == False:
        divider = 1
        numSign = 1
    for k2 in range(p2):
        for k1 in range(p1):
            re = 0.0
            im = 0.0
            for j2 in range(p2):
                x = furieA1[0][k1 + j2*p1]
                y = furieA1[1][k1 + j2*p1]
                arg = 2*math.pi*((j2/(p1*p2))*(k1+p1*k2))
                reCos = math.cos(arg)
                imSin = math.sin(arg)
                re += (x*reCos-y*imSin)
                im += (x*imSin + y*reCos)*numSign
            furieA2[0].append(re/divider)
            furieA2[1].append(im/divider)
    return furieA2



###############################################

data = [[3,8,4,7,15,11], [0,0,0,0,0,0]]
data.append([])

furieList = dft(data, True)

origCoord = dft(furieList, False)

print(data[0])

print("\nRe and Im:")
for i in range(len(furieList[0])):
    print(round(furieList[0][i], 5), "\t", round(furieList[1][i], 5), "\n" )

print("\n Answer:")
for i in range(len(origCoord[0])):
    print(round(origCoord[0][i], 5), "\t", round(origCoord[1][i], 5), "\n" )

print("\n\n\t\t\t\t\tSecond Furie\nRe and Im:")
furieA2 = dftA2(data, 2, 3, True)
for i in range(len(origCoord[0])):
    print(round(furieA2[0][i], 5), "\t", round(furieA2[1][i], 5), "\n" )

origCoord.clear()
origCoord = dftA2(furieA2, 2, 3, False)
print("\n Answer:")
for i in range(len(origCoord[0])):
    print(round(origCoord[0][i], 5), "\t", round(origCoord[1][i], 5), "\n" )
import math
import matplotlib.pyplot as plt

def CalculateFunc(x):
    #return math.sin(x)
    return x*x*x-math.sin(x)

def CreateCoordList(start, stop, step):
    allCoord = [[]]
    allCoord.append([])
    #pointCoord = []
    curPos = start
    while curPos < stop:
        allCoord[0].append(curPos)
        allCoord[1].append(CalculateFunc(curPos))
        curPos+=step

    return allCoord

def dft(data):
    furieList = [[]]
    furieList.append([])
    for k in range(len(data[0])):
        re = 0.0
        im = 0.0
        for n in range(len(data[0])):
            y = data[1][n]
            arg = 2 * math.pi *k*n/len(data[0])
            reCos = math.cos(arg)
            imSin = math.sin(arg)
            re += y * reCos
            im -= y * imSin
            
        furieList[0].append(re/len(data[0]))
        furieList[1].append(im/len(data[0]))


    return furieList


def idft(furieCoord):
    origCoord = [[]]
    origCoord.append([])
    for k in range(len(furieCoord[0])):

        re = 0.0
        im = 0.0
        for n in range(len(furieCoord[0])):
            x = furieCoord[0][n]
            y = furieCoord[1][n]
            arg = 2*math.pi*k*n/len(furieCoord[0])
            reCos = math.cos(arg)
            imSin = math.sin(arg)
            re += (x*reCos-y*imSin)
            im += (x*imSin + y*reCos)
        origCoord[0].append(re)
        origCoord[1].append(im)
    return origCoord

data = [[0,0,0,0], [1,0,1,0]]
data.append([])


#data = CreateCoordList(-3, 7, 0.1)


furieList = dft(data)

origCoord = idft(furieList)

print(data[1])

print("\nRe and Im:")
for i in range(len(furieList[0])):
    print(round(furieList[0][i], 5), "  ", round(furieList[1][i], 5), "\n" )

print("\n Answer:")
for i in range(len(origCoord[0])):
    print(round(origCoord[0][i], 5), "  ", round(origCoord[1][i], 5), "\n" )

#plt.plot(data[0], data[1], 'b-')
#plt.plot(data[0], furieList[0], 'k-')
#plt.plot(data[0], furieList[1], 'g-')

#plt.grid()
#plt.show()

#plt.plot(data[0], origCoord[0], 'r-')
#plt.plot(data[0], furieList[0], 'k-')
#plt.plot(data[0], furieList[1], 'g-')

#plt.grid()
#plt.show()

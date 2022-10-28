import math
import matplotlib.pyplot as plt
N1 =0
def dft(coord, isForward):
    global N1
    N1 = 0
    origCoord = [[]]
    origCoord.append([])
    divider = len(coord[0])
    numSign = -1
    insCount = 0
    if isForward == False:
        divider = 1
        numSign = 1
    for k in range(len(coord[0])):

        re = 0.0
        im = 0.0
        for n in range(len(coord[0])):
            N1+=5
            x = coord[0][n]
            y = coord[1][n]
            arg = 2*math.pi*k*n/len(coord[0])*numSign
            reCos = math.cos(arg)
            imSin = math.sin(arg)
            re += (x*reCos-y*imSin)
            im += (x*imSin + y*reCos)
        origCoord[0].append(re/divider)
        origCoord[1].append(im/divider)
    return origCoord


###############################################
def dftA1(data, p1, p2, isForward):
    global N1
    furieA1 = [[None]*p1*p2]
    furieA1.append([None]*p1*p2)
    divider = p1
    numSign = -1
    if isForward == False:
        divider = 1
        numSign = 1
    for k1 in range(p1):
        for j2 in range(p2):
            re = 0.0
            im = 0.0
            for j1 in range(p1):
                N1+=5
                x = data[0][j2 + p2*j1]
                y = data[1][j2 + p2*j1]
                arg = 2 * math.pi *j1*k1/p1*numSign
                reCos = math.cos(arg)
                imSin = math.sin(arg)
                re += (x*reCos-y*imSin)
                im += (x*imSin + y*reCos)
            furieA1[0][k1*p2 + j2] = re/divider
            furieA1[1][k1*p2 + j2] = im/divider
            
    return furieA1



def dftA2(data, p1, p2, isForward):
    
    furieA2 = [[None]*p1*p2]
    furieA2.append([None]*p1*p2)
    global N1
    N1 = 0
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
                N1+=5
                x = furieA1[0][k1*p2 + j2]
                y = furieA1[1][k1*p2 + j2]
                arg = 2*math.pi*((j2*(k1+p1*k2)/(p1*p2)))*numSign
                reCos = math.cos(arg)
                imSin = math.sin(arg)
                re += (x*reCos-y*imSin)
                im += (x*imSin + y*reCos)
            furieA2[0][k2*p1 + k1] = re/divider
            furieA2[1][k2*p1 + k1] = im/divider
    return furieA2

def DecToRevBinary(binar, num):
    powerOf2 = 0
    while(num != 0):
        binar.append(int(num%2))
        num=int(num/2)
        powerOf2+=1
    return powerOf2

def DecToBinary(binar, num):
    powerOf2 = 0
    while(num != 0):
        binar.insert(0, int(num%2))
        num=int(num/2)
        powerOf2+=1
    return powerOf2

def BinaryToDec(binar):
    num = 0
    powerOf2 = 1
    for i in range(len(binar)-1, -1,-1):
        num += binar[i]*powerOf2
        powerOf2*=2
    return num

def RevBinaryToDec(binar):
    num = 0
    powerOf2 = 1
    for i in range(len(binar)):
        num += binar[i]*powerOf2
        powerOf2*=2
    return num

def fftAn(r, dt, binarMaxSize, isForward):
    global N1
    answ = [[]]
    answ.append([])
    divider = 2
    numSign = -1
    if isForward == False:
        divider = 1
        numSign = 1
    for i in range(len(dt[0])):
        N1+=5
        binar = []
        DecToBinary(binar, i)
        for j in range(binarMaxSize - len(binar)):
            binar.insert(0, 0)
        
        ksum = 0
        powerOf2 = 1
        for j in range(r):
            ksum+=powerOf2*binar[j]
            powerOf2*=2

        binar[r-1] = 0
        ind = BinaryToDec(binar)
        re = dt[0][ind]
        im = dt[1][ind]

        arg = 2*math.pi*ksum/powerOf2*numSign
        binar[r-1] = 1
        ind = BinaryToDec(binar)
        x = dt[0][ind]
        y = dt[1][ind]
        reCos = math.cos(arg)
        imSin = math.sin(arg)
        re += (x*reCos-y*imSin)
        im += (x*imSin + y*reCos)

        

        answ[0].append(re/divider)
        answ[1].append(im/divider)
    return answ


def fft(data, isForward):
    global N1
    N1 = 0
    k = []
    dt = data[:]
    answ = [[]]
    answ.append([])
    binarMaxSize = DecToRevBinary(k, len(data[0]) - 1)
    for i in range(binarMaxSize):
        dt = fftAn(i+1, dt, binarMaxSize, isForward)

    for i in range(len(dt[0])):
        k = []
        DecToRevBinary(k, i)
        for j in range(binarMaxSize - len(k)):
            k.append(0)
        ind = BinaryToDec(k)
        answ[0].append(dt[0][ind])
        answ[1].append(dt[1][ind])
    dt = answ
    return dt
    

###############################################

#data = [[1,2,3,4], [0,0,0,0]]
data = [[1,2,3,4,5,6,7,8], [0,0,3,0,0,0,0,0]]
#data = [[1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

size = len(data[0])
for i in range(1, math.floor(math.sqrt(len(data[0])))+1):
    if(size%i==0):
        mult1 = i
        mult2 = int(size/i)

print("\n\nMult: ",mult1," ",mult2,"\n\n")
data.append([])

furieList= dft(data, True)
kCount = N1
origCoord = dft(furieList, False)
origCount = N1

print(data[0])

print("\nRe and Im:")
for i in range(len(furieList[0])):
    print(round(furieList[0][i], 5), "\t", round(furieList[1][i], 5), "\n" )



print("\nDir and Rev = ", kCount, "  ", origCount)
print("\n\n Answer:")
for i in range(len(origCoord[0])):
    print(round(origCoord[0][i], 5), "\t", round(origCoord[1][i], 5), "\n" )

print("\n\n\t\t\t\t\tSecond Furie\nRe and Im:")
furieA2 = dftA2(data, mult1, mult2, True)
kCount = N1
for i in range(len(furieA2[0])):
    print(round(furieA2[0][i], 5), "\t", round(furieA2[1][i], 5), "\n" )

origCoord.clear()
origCoord = dftA2(furieA2, mult1, mult2, False)
origCount = N1
print("\n Answer:")
for i in range(len(origCoord[0])):
    print(round(origCoord[0][i], 5), "\t", round(origCoord[1][i], 5), "\n" )

print("\nDir and Rev = ", kCount, "  ", origCount)



print("\n\n\t\t\t\t\tThird Furie\nRe and Im:")
furieList= fft(data, True)

for i in range(len(furieList[0])):
    print(round(furieList[0][i], 5), "\t", round(furieList[1][i], 5), "\n" )
kCount = N1
furieA3 = fft(furieList, False)
print("\n\nAnswer:")
origCount = N1
for i in range(len(furieA3[0])):
    print(round(furieA3[0][i], 5), "\t", round(furieA3[1][i], 5), "\n" )
print("\nDir and Rev = ", kCount, "  ", origCount)


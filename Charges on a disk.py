from numpy import zeros, random, pi, sqrt, roll, exp, degrees, arctan, linspace
import pylab as pl
from  copy import deepcopy as cpy
import time
import itertools as it

N = 5 # number of particals
D=100
R=D/2  # Radus of disk
delta = 0.5  # stepsize
M = 1000*N  # number of steps per iteration for T
T = (10**-1, 10**-7)  # temperature
test_temp = 1


def dist(P_1, P_2):  # get distance between [x1,y1], [x2,y2]
    x1 = P_1[0]
    y1 = P_1[1]
    x2 = P_2[0]
    y2 = P_2[1]
    return sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def E(P_1, P_2):  # get potenital between [x1,y1], [x2,y2], for adding diff charges later
    return 1 / dist(P_1, P_2)

def all_energy(co_ords):
    lst = []
    for i in range(N):
        for j in range(1, N - i):
            lst.append([co_ords[i], co_ords[i + j]])
        #		print(lst)
            
    en = it.starmap(E, lst)
    return en

def zeroth_energy(co_ords, full_energy): # find energy of 0th co-ord to all others, 
    en_new = cpy(full_energy)
    lst = []
    for j in range(1, N):
        lst.append([co_ords[0], co_ords[j]])
    for x, val in enumerate(it.starmap(E, lst)):
        en_new[x] = val
    return en_new

def rand_fill(): #create N 2dimentional arrays
    p = zeros((N, 2))
    for n in range(N):
        p[n] = [random.randint(-R,R), random.randint(-R,R)]
        while sum(p[n] ** 2) > R ** 2:
            p[n] = [random.randint(-R,R), random.randint(-R,R)]
    return p


def rand_select(co_ords): #shifts the co-ords for a random 0th co_ord
    temp = cpy(co_ords)
    shift = random.randint(N)
    temp = cpy(roll(temp, shift, axis=0))
    return temp


def mover(co_ords, tempature):
    m = 0  # counter
    co_work = cpy(co_ords)
    while m <= M:
#        print(m)
        m += 1
        co_olds = rand_select(co_work)
        co_new = cpy(co_olds)

        axis = random.choice((0, 1))
        sign = random.choice((+1, -1))
        co_new[0, axis] = co_new[0, axis] + delta * sign
        if sum(co_new[0] ** 2) > R ** 2:
            m -= 1
            continue

        e_old = list(all_energy(co_olds))
        e_new = zeroth_energy(co_new, e_old)
        change_e = abs(sum(e_old)) - abs(sum(e_new))
        #		print (change_e)
#        if change_e <= 0:             
#          print('change_e=%f\t change/t=%f exp=%f'%(change_e, change_e / tempature, exp(change_e / tempature) ))
        if change_e >= 0:
            co_work = cpy(co_new)
            continue
        
        elif random.rand() < exp(change_e / tempature):

            co_work = cpy(co_new)
            continue
        else:
            continue
    return co_work


def worker():
    start_time = time.time()
    a = rand_fill()
#    for i in a:
#        print(sum(i ** 2))
#    print(sum(all_energy(a)))
    plott(a)
    for t in linspace(T[0], T[1], 20):
        a = mover(a, t)
        plott(a)
    print("--- %s seconds ---" % (time.time() - start_time))
#    for i in a:
#        print(sum(i ** 2))
#    print(sum(all_energy(a)))
    plott_sav(a)
    return True

def plott(co_ords):
    pl.clf()
    for x, y in iter(co_ords):
#        print('%f \t %f'%(x,y))
        pl.scatter(x, y)
    circle1=pl.Circle((0,0),R,color='g',fill=False)
    fig = pl.gcf()
    fig.gca().add_artist(circle1)
    pl.xlim(-R-1,R+1)
    pl.ylim(-R-1,R+1)    
    pl.axes().set_aspect('equal')
    pl.show()
    return True


def plott_sav(co_ords):
    pl.clf()
    for x, y in iter(co_ords):
#        print('%f \t %f'%(x,y))
        pl.scatter(x, y)
    circle1=pl.Circle((0,0),R,color='g',fill=False)
    fig = pl.gcf()
    fig.gca().add_artist(circle1)
    pl.xlim(-R-1,R+1)
    pl.ylim(-R-1,R+1)   
    title='%d charges.'%N
    pl.title(title)
    pl.axes().set_aspect('equal')
#    pl.show()
    save='%d_charges_on_rad_%d_rand_%d'%(N,R,random.randint(1,100))
    pl.savefig(save)
    return True
    
def to_polar(co_ords):
    points=[]
    for x, y in co_ords:
        print(x,y)
        r=sqrt(x**2+y**2)
        theta=degrees(arctan(x/y))
        points.append[r,theta]
    return points

##~~~~~ Field bodge added random number to save for parelleisum 
for run in range(14,17):
    N=run
    M = 1000*N
    worker()

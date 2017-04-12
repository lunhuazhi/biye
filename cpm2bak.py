import numpy as np
import matplotlib.pyplot as pl
def gt(t,L, T):
    gt = []

    for index, value in enumerate(t):
        gt.append((1 / 2 / L / T) * (1 - np.cos(2 * np.pi * value / L / T)))
    return gt


def qt(t, gt):
    qt = []
    for index, value in enumerate(t):
        qt.append(sum(gt[0:index]) * 0.01)
    return qt


def soucrce(numble, m):
    temp = []
    yuan = 2 * np.random.randint(0, m, numble) - (m - 1)
    for item in yuan:
        temp.append(item)
    return temp

def modulation(qt,source,nsample):
    phase =[0,0]
    acc_ph =[]
    for index ,val in enumerate(source):
        if index >= L:
            phase.append(np.mod(phase[index-1]+np.pi*h*source[index-L],2*np.pi) )
        if index == 0:
            acc_ph.extend(2*np.pi*h*source[index]*qt[:nsample] + phase[index])
        else:
            acc_ph.extend(2*np.pi*h*source[index]*qt[:nsample] + 2*np.pi*h*source[index-1]*qt[nsample:] + phase[index])
    return acc_ph
def gen_tree(m, p):
    # gstate_from ={}
    to_state = {}  # 字典，key为当前态和新输入的，value为在当前态和新输入的数值到达的新状态
    phaslist = range(0, p)
    informationlist = [1, -1, 3, -3]

    for phas in phaslist:
        for information_old in informationlist:
            for information_new in informationlist:
                phas_new = np.mod(phas + information_old, p)
                to_state[(phas, information_old, information_new)] = (phas_new, information_new)

    return to_state



def bendixiangwei():
    pass
def corfilter():
    pass
def demod():
    pass




h = 0.25
m = 4
L = 2
T = 1
t = np.arange(0,2,200)
nsample = 100
nsymbol = 4
source =  soucrce(nsymbol,m)

gt = gt(t,L,T)
pl.plot(np.linspace(0,2,nsample*2),gt)
qt =np.array(qt(t,gt))
pl.plot(t,qt)
acc = modulation(qt,source,nsample)

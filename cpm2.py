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
            acc_ph.extend(np.mod(2*np.pi*h*source[index]*qt[:nsample] + phase[index],2*np.pi))
        else:
            acc_ph.extend(np.mod(2*np.pi*h*source[index]*qt[:nsample] + 2*np.pi*h*source[index-1]*qt[nsample:] + phase[index],2*np.pi))
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
#to_state 字典，当前为phas information_old 输入information_new  到达的新状态
    return to_state



def bendixiangwei():
    pass
def corfilter(qt,nsample):
    coor ={}
    fudu = [-1,+1,-3,+3]
    for info_old in fudu:
        for info_new in fudu:
           coor[(info_old,info_new)] =  np.exp(-1j*np.mod((2*np.pi*h*info_new*qt[:nsample]+2*np.pi*h*info_old*qt[nsample:]),2*np.pi))

    return coor
def demodmy(h,p,xinhao,fil,nsample,gentree):
    cnt = 0
    cnt2 = 0
    phaselist = range(0,p)
    val_old ={}
    val_new ={}
    survive={}
    survive2 = {}
    fudu =[-1,1,-3,3]
    number = int(len(xinhao)/nsample)
    print(number)
    demod_sig =[]

    for thenumber in range(0,number):#thenumber 表示第几个符号

        for phase in phaselist:
            for info_new in fudu:
                for state_from in gentree.keys():
                    if gentree[state_from] == (phase,info_new):
                        fenzhi = np.sum(xinhao[thenumber*nsample:(thenumber+1)*nsample]*fil[(state_from[1],state_from[2])]*1/nsample)*np.exp(-1j*state_from[0]*h*np.pi)
                        fenzhi = np.real(fenzhi)
                        survive2[(thenumber,state_from)] = fenzhi+val_old.get((state_from[0],state_from[1]),0)

                temp = sorted(survive2.values(),reverse=True)

                key_in_survive = list(survive2.keys())
                for key in key_in_survive:
                    if survive2[key]== temp[0]:
                        survive[key] = (phase, info_new)

                survive2 ={}
                #
                # for key in survive:
                #     if survive[key] != temp[0]:
                #         survive.pop(key)
                #     else:
                #         survive[key] = (phase,info_new)
                val_new[(phase,info_new)] = temp[0]
                #temp =[]


        cnt = cnt + 1
        if cnt == 10:
            cnt = 0
            cnt2 = cnt2 + 1
            demod_sig.extend(huisu(cnt2,survive,val_new))
            survive = {}

        val_old = val_new
        val_new = {}

    return demod_sig

def huisu(cnt2,survive,val):
    '''
    :param survive:幸存路径的集合 
    :param val: 是一个字典，key最新的状态，值是每个状态对应的似然值
    :return: 
    '''
    desig =[]
    maxval = max([v for v in val.values()])
    for state in val:
        if val[state] == maxval:
            start_point =  state
            break
    #for i in range(10*cnt2-1,-1,10*cnt2-11):
       # print(i)
    for symbol in range(10*cnt2-1,10*cnt2-11,-1):
        for branch in survive:
            if branch[0] == symbol and survive[branch] == start_point:
                desig.insert(0,branch[1][2])
                start_point = (branch[1][0],branch[1][1])
                break
    return desig





p =8
h = 0.25
m = 4
L = 2
T = 1
nsample = 100
t = np.arange(0,2,1/nsample)

nsymbol = 20
source =  soucrce(nsymbol,m)

gt = gt(t,L,T)
pl.plot(t,gt)
qt =np.array(qt(t,gt))
pl.plot(t,qt)
acc = np.array(modulation(qt,source,nsample))

xinhao = np.exp(1j*acc)
#######################demod####################

#生成相位网格
gentree = gen_tree(m,p)
#生成匹配滤波器
fil = corfilter(qt,nsample)
#开始解调，设置回溯深度为10，回溯深度不应该小于5*L
de = demodmy(h,p,xinhao,fil,nsample,gentree)
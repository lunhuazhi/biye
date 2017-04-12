import numpy as np
import matplotlib.pyplot as pl


# from scipy import *
def wgn(x, snr):
    snr = 10 ** (snr / 10.0)
    xpower = np.sum(x ** 2) / len(x)
    npower = xpower / snr
    return np.random.randn(len(x)) * np.sqrt(npower)


def gt(t, L, T):
    gt = []
    for index, value in enumerate(t):
        if value >= 0 and value <= L * T:
            gt.append((1 / 2 / L / T) * (1 - np.cos(2 * np.pi * value / L / T)))
        else:
            gt.append(0)
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


def moudlation(source, qt, nsample):
    phas = [0, 0]
    acc_ph = []
    for index, value in enumerate(source):
        if index >= 2:
            phas.append(np.mod(phas[index - 1] + np.pi * h * source[index - 2], 2 * np.pi))
        acc_ph.extend(np.mod(
            phas[index] + 2 * np.pi * h * qt[nsample:] * source[index - 1] + 2 * np.pi * h * qt[:nsample] * source[
                index], 2 * np.pi))

    return acc_ph


def demod(nsample, xinhao, m, p, corelate, state_tree):
    branch_metric = {}
    acc_old = {}
    acc_new = {}
    survive = {}
    phaselist = range(0, p)
    demod_sig =[]
    jishu = 0
    informationlist = [1, -1, 3, -3]
    for every_symbol in range(0, int(len(xinhao) / nsample)):
        #print(int(len(xinhao) / nsample))
        for phase in phaselist:
            for inforold in informationlist:
                for infornew in informationlist:
                    branch_metric[(phase, inforold, infornew)] = np.sum(
                        xinhao[every_symbol * nsample:(every_symbol + 1) * nsample] * np.conj(
                            np.exp(corelate[(inforold, infornew)] * 1j)))
                    branch_metric[(phase, inforold, infornew)] = np.real(
                        branch_metric[(phase, inforold, infornew)] * np.exp(-1j * np.pi * phase * h))
        for state_from, state_new in state_tree.items():
            if state_new in acc_new:
                temp = acc_old.get((state_from[0],state_from[1]),0) + (branch_metric[state_from])
                if acc_new[state_new] > temp:
                    pass
                else:
                    acc_new[state_new] = temp
                    survive[(every_symbol, state_new)] = state_from
            else:
                acc_new[state_new] = acc_old.get((state_from[0],state_from[1]), 0) + (branch_metric[state_from])
                survive[(every_symbol, state_new)] = state_from

        if jishu ==9:
            temp =[]
            jishu = 0
            maxval = max(acc_new.values())
            print(maxval)
            for state in acc_new:
                if acc_new[state] == maxval:
                    print("begin hui su")
                    index = state
                    break

            for fuhao in range(every_symbol,-1,-1):
                temp.insert(0, survive[(fuhao, index)][2])
                    #demod_sig.insert(0, survive[(fuhao, index)][2])
                index = (survive[(fuhao, index)][0], survive[(fuhao, index)][1])
            demod_sig.extend(temp)

        else:
            jishu = jishu +1


        acc_old = acc_new
        acc_new = {}
        branch_metric = {}

    return (survive), demod_sig


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


def corelate(nsample, m, h, qt):
    # phaselist = range(0,m)
    info_list = [-1, 1, -3, 3]
    bendixiangwei = {}
    for info_old in info_list:
        for info_new in info_list:
            bendixiangwei[(info_old, info_new)] = np.mod(
                2 * np.pi * h * qt[nsample:] * info_old + 2 * np.pi * h * qt[:nsample] * info_new, 2 * np.pi)
            # return np.mod(2*np.pi*h*qt[nsample:]*info_old+2*np.pi*h*qt[:nsample]*info_new,2*np.pi)

    return bendixiangwei


# def veitebi(survive, acc_old, nsymbol):
#     # print(acc_old)
#     # print('----------------------------')
#     demod_xinhao = []
#     maxva = max(acc_old.values())
#     # print(maxva)
#     for index in acc_old:
#         # 开始回溯
#
#         if acc_old[index] == maxva:
#             # print(index)
#             state_new = index
#             for symbol in range(nsymbol - 1, -1, -1):
#                 demod_xinhao.insert(0, survive[(symbol, state_new)][2])
#                 state_new = (survive[(symbol, state_new)][0], survive[(symbol, state_new)][1])
#
#     # print(demod_xinhao)
#     return demod_xinhao


def error(source, demod_xinhao):
    cnt = 0
    for index in range(0, len(source)):
        if source[index] != demod_xinhao[index]:
            cnt = cnt + 1
    return cnt / len(source)


m = 4
nsample = 100 # 1s种用100个点代替
L = 2  # 关联长度为2
nsymbol = 20
h = 0.25  # 调制指数为0.2
p = 8
T = 1  # 一个符号为1s
t = np.linspace(0, 2, 2 * nsample)  # 时间
gt = np.array(gt(t, L, T))  # 成型脉冲
qt = np.array(qt(t, gt))  # qt
# pl.plot(t,qt)
EbN0 = np.linspace(0, 15, 15)
snr = [i - 10 * np.log10(nsample) + 10 * np.log10(np.log2(m)) for i in EbN0]
state_tree = gen_tree(m, p)

corelate = corelate(nsample, m, h, qt)
source = soucrce(nsymbol, m)
# t_time = np.linspace(0,len(source)*T,len(source)*T*nsample)

acc_ph = np.array(moudlation(source, qt, nsample))
xinhao = np.exp(acc_ph * 1j)

error_rate = []

survive, de = demod(nsample, xinhao, m, p, corelate, state_tree)
# for sn in snr:
#     xinhao_noise = xinhao + wgn(xinhao, sn)
#     survive, acc_new = demod(nsample, xinhao_noise, m, p, corelate, state_tree)
#     demod_xinhao = veitebi(survive, acc_new, nsymbol)
#
#     error_rate.append(error(source, demod_xinhao))

# pl.semilogy(EbN0, error_rate)
# i1 = max(acc_new.values())
# for key in acc_new.keys():
#    if acc_new[key] == i1:
#        print(key)
#        break



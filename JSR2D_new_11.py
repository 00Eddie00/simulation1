import os
import numpy as np
from math import sqrt

global CCAJSR, NEWCCA, NEWF
global CTRL_AREA, TOTAL_AREA, All_of_det, All_of_area
global K1, K2, F
global abcMatrix, nlMatrix
global Gn
global B_INNER, B_WALL, B_INFLOW, B_OUTFLOW, B_SYMMETRY, B_SOURCE
global ICARYR, ICAFSR, AVG_CA_JSR, AVG_Gn_JSR, CURSTEP, RELEASE_TIMES
global NE, AREA

'''
针对Fn的修改五
'''


def INITIAL_PARAMETER():
    """
    初始化参数
    """
    global CCAMYO, CCAFSR, BCSQ, KDCSQ, DCAJSR, DF, DCARYR, DCAFSR, DT, RELEASE_TIMES
    global B_INNER, B_WALL, B_INFLOW, B_OUTFLOW, B_SYMMETRY, B_SOURCE
    global CCAJSR, Gn, NEWCCA, NEWF
    global TOTAL_AREA, CTRL_AREA, All_of_area, All_of_det
    global K1, K2, F
    global nmax

    # Jdiffusion扩散项
    DCAJSR = 3.5 * 10 ** 8  # 终池jSR 自由Ca离子扩散系数
    DF = 2 * 10 ** 7  # 荧光指示剂染料 扩散系数

    # Jbuff
    BCSQ = 14.0  # 方程中的B，所有肌钙集蛋白（肌浆网SR上调控钙储存的主要结合蛋白，调节SR钙释放）的量
    KDCSQ = 0.63  # 方程中的K，肌钙集蛋白的Ca离子解离常数

    # Jrefill,Jryr两项
    # 扩散系数
    DCARYR = 6.5 * 10 ** 7  # RyR通道 Ca离子扩散系数
    DCAFSR = 0.7854 * 10 ** 6  # 连接的fSR通道 Ca离子扩散系数
    # 浓度
    CCAMYO = 0.0001  # [Ca2+]cyto  胞浆的Ca浓度
    CCAFSR = 1.0  # [Ca2+]fSR  fSR通道中Ca离子浓度

    # Jdye染料
    F = 0.1  # [F]T   fluo-3总浓度
    K1 = 48800  # KF+  Ca离子和fluo-3的结合速率
    K2 = 19520  # KF-  Ca离子和fluo-3的解离速率
    # K1 = 0
    # K2 = 0

    # 关闭RyR通道用到
    DT = 2 * 10 ** -6  # dt
    RELEASE_TIMES = 2 * 10 ** -2  # 0.02s,这个就是ryr通道开放的时间，后80毫秒是恢复的

    # 点的属性
    B_INNER = 0  # 内部
    B_INFLOW = 2  # 入流
    B_WALL = 1  # 壁面
    B_OUTFLOW = 4  # 出流
    B_SYMMETRY = 8
    B_SOURCE = 16

    # 一些初始化
    nmax = 0
    CCAJSR = np.empty(PMAX, float)  # fn,肌质网每一点钙的浓度
    Gn = np.ones(PMAX, float) * (1 / 14)  # gn初始值 1/14  固定
    NEWCCA = np.empty(PMAX, float)  # 存每一次迭代的中间结果fn
    NEWF = np.empty(PMAX, float)  # 存每一次迭代的中间结果gn
    CTRL_AREA = np.zeros(PMAX, float)  # 每一个点对应的三角形单元面积
    TOTAL_AREA = 0.0  # TOTAL_AREA区域的面积

    for I in range(0, NP):
        CCAJSR[I] = 1.0  # 为肌质网每一点的初始钙浓度赋值为1

    All_of_det = np.zeros(NE, float)
    All_of_area = np.zeros(NE, float)
    MATRIX1 = np.zeros((3, 3))
    count = np.zeros(NP, float)

    TOTAL_CTRL_AREA = 0
    DET = 0.0
    AREA = 0.0
    nmax = 0
    for E in range(0, NE):  # NE为三角形单元的数量,计算D和V
        N1 = NOD[0, E]
        N2 = NOD[1, E]
        N3 = NOD[2, E]
        count[N1] = count[N1] + 1
        count[N2] = count[N2] + 1
        count[N3] = count[N3] + 1  # 计算每个点相关的三角形的个数
        for I in range(0, 3):
            MATRIX1[I, 0] = 1.0
            MATRIX1[I, 1] = PX[NOD[I, E]]
            MATRIX1[I, 2] = PY[NOD[I, E]]
        DET = np.linalg.det(MATRIX1)  # D 行列式
        AREA = DET / 2.0  # A 三角形面积
        TOTAL_AREA = TOTAL_AREA + AREA

        All_of_det[E] = DET  # 存储每个三角形的D
        All_of_area[E] = AREA  # 存储每个三角形的A

        for i in range(0, 3):
            CTRL_AREA[NOD[i, E]] = CTRL_AREA[NOD[i, E]] + AREA  # 每一点的控制面积

    for i in range(0, NP):
        TOTAL_CTRL_AREA = TOTAL_CTRL_AREA + CTRL_AREA[i]
    # print('TOTAL_CTRL_AREA', TOTAL_CTRL_AREA)
    # print('TOTAL_TRI_AREA', TOTAL_AREA)

    # save = open("DATA3\\TRI_AREA.dat", "w")
    # for E in range(0, NE):
    #     save.write(str(All_of_area[E]) + "\n")
    # print("写入每个三角形的面积")
    #
    # save = open("DATA3\\CTRL_AREA.dat", "w")
    # for I in range(0, NP):
    #     save.write(str(CTRL_AREA[I]) + "\n")
    # print("写入每个点的控制面积")

    # save_area = open("area.dat", "w")
    # save_ctrl_area = open("ctrl_area.dat","w")
    # for i in range(0,NE):
    #     save_area.write(str(All_of_area[i])+ "\n")
    # for i in range(0,NP):
    #     save_ctrl_area.write(str(CTRL_AREA[i])+ "\n")
    # save_ctrl_area.close()
    # save_area.close()

    nmax = int(max(count))
    # print("NMAX ", nmax)
    # print("NP ", NP)
    # print("NE ", NE)

    print("    执行cal_ABCNL()函数")
    cal_ABCNL()
    print("    执行search_Triangle()函数")
    search_Triangle()
    # print_triangle()
    print("INITIAL_PARAMETER()函数执行完成")


def DOLOOP_DIFFUSION():
    """
    循环函数
    """
    global DCARYR, DSAVE, NSTEP, CURSTEP, RELEASE_TIMES, STEP1, STEP2
    global AVG_CA_JSR, ICARYR, ICAFSR, AVG_Gn_JSR

    CURSTEP = 1  # 从多少步开始，默认从第1步开始
    DSAVE = 1  # 每做多少步保存一次
    NSTEP = 50000  # 程序到多少步结束
    ITERATION = 10  # 迭代次数
    PATH = "DATA\\DATA_blink"  # 保存路径

    STEP1 = -1  # 不读
    STEP2 = -1  # 不读
    pathsave = PATH + '\\SAVE'
    pathgn = PATH + '\\Gn'
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    if not os.path.exists(pathgn):
        os.makedirs(pathgn)

    if CURSTEP == 1:  # 从第一步开始，保存SAVE00000000.dat
        # 存储SAVE00000000文件
        file_f0 = pathsave + "\\SAVE00000000.dat"
        save = open(file_f0, "w")
        for I in range(0, NP):
            save.write(str(CCAJSR[I]) + "\n")  # 第0步的值在CCAJSR中，赋值为1
        AVERAGE_CCA()  # 计算平均浓度 没问题
        CA_CURRENT()  # 计算0时刻电流 没问题
        save.write(str(AVG_CA_JSR) + "\n")  # NP + 1
        save.write(str(ICARYR) + " ")  # NP + 2
        save.write(str(ICAFSR) + "\n")  # NP + 2
        save.write(str(CURSTEP) + "\n")  # NP + 3
        print("数据写入  " + file_f0)
        save.close()

        file_g0 = pathgn + "\\Gn00000000.dat"
        save2 = open(file_g0, "w")
        for I in range(0, NP):
            save2.write(str(Gn[I]) + "\n")  # 第0步的值在Gn中，赋值为(1/14)
        AVERAGE_CCA_Gn()
        save2.write(str(AVG_Gn_JSR) + "\n")  # NP + 1
        save2.write(str(CURSTEP) + "\n")  # NP + 2
        print("数据写入  " + file_g0)
        save2.close()
    else:  # 不是第一步开始，先读取先前保存的Save、Gn
        STRN1 = str(CURSTEP - 1).zfill(8)
        read_fn = pathsave + "\\SAVE" + STRN1 + '.dat'
        FIL4 = open(read_fn, "r")
        read_gn = pathgn + "\\Gn" + STRN1 + '.dat'
        FIL5 = open(read_gn, "r")

        I1 = 0
        I2 = 0
        for line in FIL4.readlines():
            if I1 < NP:
                CCAJSR[I1] = float(line)  # 存在CCAJSR中
            else:
                if I1 == NP:
                    AVG_CA_JSR = float(line)
                elif I1 == NP + 1:
                    tmp = line.split()
                    ICARYR = float(tmp[0])
                    ICAFSR = float(tmp[1])
                else:
                    STEP1 = int(line)
                    print("读取原有的Save", STEP1)
            I1 = I1 + 1
        FIL4.close()

        for line in FIL5.readlines():
            if I2 < NP:
                Gn[I2] = float(line)  # 存在Gn中
                NEWF[I2] = float(line)
            else:
                if I2 == NP:
                    AVG_Gn_JSR = float(line)
                else:
                    STEP2 = int(line)
                    print("读取原有的Gn", STEP2)
            I2 = I2 + 1
        FIL5.close()

    RELEASE_STEP = int((RELEASE_TIMES + (10 ** -6)) / DT)  # 释放时间的迭代次数
    if CURSTEP >= RELEASE_STEP:
        DCARYR = 0
    # 释放的点就停止释放

    for J in range(CURSTEP, NSTEP + 1):

        # 每一步需要迭代 5-10 次
        CORRECTION_EQUATION(ITERATION)
        GN_EQUATION(ITERATION)
        for I in range(0, NP):
            CCAJSR[I] = NEWCCA[I]

        if (J % DSAVE == 0) or (J == NSTEP):
            # 需要输出的时候再计算这三个值
            AVERAGE_CCA()
            AVERAGE_CCA_Gn()
            CA_CURRENT()

            STRN = str(J).zfill(8)
            FILENAME1 = pathsave + "\\SAVE" + STRN + '.dat'
            FIL1 = open(FILENAME1, "w")
            for I in range(0, NP):
                FIL1.write(str(CCAJSR[I]) + "\n")
            print("数据写入 ", FILENAME1)
            FIL1.write(str(AVG_CA_JSR) + "\n")
            FIL1.write(str(ICARYR) + " ")
            FIL1.write(str(ICAFSR) + "\n")
            FIL1.write(str(J) + "\n")

            FILENAME2 = pathgn + "\\Gn" + STRN + '.dat'
            FIL2 = open(FILENAME2, "w")
            for I in range(0, NP):
                FIL2.write(str(Gn[I]) + "\n")
            print("数据写入 ", FILENAME2)
            FIL2.write(str(AVG_Gn_JSR) + "\n")  # 计算平均值
            FIL2.write(str(J) + "\n")

            FIL1.close()
            FIL2.close()

        if J == RELEASE_STEP:
            DCARYR = 0
        #
        # if J >= 50:
        #     DCARYR = 0

        # print("DCARYR ", DCARYR)
        for i in range(0, NP):
            NEWCCA[i] = 0.
            NEWF[i] = 0.
    print("执行成功")


# ****************************************************************************
def cal_ABCNL():
    """
    计算a,b,c,N,L
    """
    global All_of_det, nlMatrix, abcMatrix
    nlMatrix = np.zeros((NE, 3, 3))
    abcMatrix = np.zeros((NE, 3, 3))
    for i in range(0, NE):
        for j in range(0, 3):
            p2 = NOD[(j + 1) % 3, i]
            p3 = NOD[(j + 2) % 3, i]
            D = All_of_det[i]

            # 第i个三角形，每个点为中心的a,b,c
            abcMatrix[i, j, 0] = (PX[p2] * PY[p3] - PX[p3] * PY[p2]) / D
            abcMatrix[i, j, 1] = (PY[p2] - PY[p3]) / D
            abcMatrix[i, j, 2] = (PX[p3] - PX[p2]) / D

            # 第i个三角形，每个点为中心的L，Nix,Niy
            nlMatrix[i, j, 0] = sqrt((PX[p2] - PX[p3]) ** 2 + (PY[p2] - PY[p3]) ** 2)
            nlMatrix[i, j, 1] = (PY[p3] - PY[p2]) / nlMatrix[i, j, 0]
            nlMatrix[i, j, 2] = (PX[p2] - PX[p3]) / nlMatrix[i, j, 0]


def search_Triangle():
    """
    存储某一点相邻的所有三角形单元
    """
    global triangleNumber, nmax
    triangleNumber = np.full([NP, nmax, 2], -1)  # 存放围绕此点的三角形单元编号

    if not os.path.exists("TRIANGLE_NUMBER"):
        os.mkdir("TRIANGLE_NUMBER")
    for N in range(0, NP):  # 遍历所有的点
        index = 0
        # 先检查有没有文件
        # 1.没有文件->写
        if not os.path.exists("TRIANGLE_NUMBER\\POINT" + str(N + 1) + ".dat"):
            with open("TRIANGLE_NUMBER\\POINT" + str(N + 1) + ".dat", "w") as file_object:
                for E in range(0, NE):  # 遍历每个三角形
                    for i in range(0, 3):  # 遍历三角形的三个顶点
                        if (N == NOD[i][E]):
                            triangleNumber[N][index][0] = E
                            triangleNumber[N][index][1] = i
                            file_object.write(str(E) + " ")
                            file_object.write(str(i) + "\n")
                            index = index + 1
                            break
            file_object.close()
        # 2.有文件->直接读取
        else:
            with open("TRIANGLE_NUMBER\\POINT" + str(N + 1) + ".dat", "r") as file_object:
                for line in file_object.readlines():
                    current_line = list(filter(not_empty, line.strip("\n").split(" ")))
                    triangleNumber[N][index][0] = current_line[0]
                    triangleNumber[N][index][1] = current_line[1]
                    index = index + 1
            file_object.close()


# def print_triangle():
#     path = "PARAMETER"
#     if not os.path.exists(path):
#         os.mkdir(path)
#
#     with open(path+"\\triangle_area.dat", "w") as dat_file:
#         for i in range(0,NE):
#             dat_file.write(str(All_of_area[i]) + "\n")
#         dat_file.write(str(TOTAL_AREA))
#
#     with open(path+"\\ctrl_area.dat", "w") as dat_file:
#         for i in range(0, NP):
#             dat_file.write(str(CTRL_AREA[i]) + "\n")
#
#     with open(path+"\\b.dat", "w") as dat_file:
#         for i in range(0 ,NE):
#             for j in range(0,3):
#                 dat_file.write(str(abcMatrix[i][j][1]) + " ")
#             dat_file.write("\n")
#
#     with open(path+"\\c.dat", "w") as dat_file:
#         for i in range(0 ,NE):
#             for j in range(0,3):
#                 dat_file.write(str(abcMatrix[i][j][2]) + " ")
#             dat_file.write("\n")
#
#     with open(path+"\\l.dat", "w") as dat_file:
#         for i in range(0 ,NE):
#             for j in range(0,3):
#                 dat_file.write(str(nlMatrix[i][j][0]) + " ")
#             dat_file.write("\n")
#
#     with open(path+"\\nix.dat", "w") as dat_file:
#         for i in range(0 ,NE):
#             for j in range(0,3):
#                 dat_file.write(str(nlMatrix[i][j][1]) + " ")
#             dat_file.write("\n")
#
#     with open(path+"\\niy.dat", "w") as dat_file:
#         for i in range(0 ,NE):
#             for j in range(0,3):
#                 dat_file.write(str(nlMatrix[i][j][2]) + " ")
#             dat_file.write("\n")


# def calmax_point_type():
#     outboundary = np.zeros(NE,int)
#     inboundary = np.zeros(NE,int)
#     save8 = open("in.dat","w")
#     save9 = open("out.dat", "w")
#     for i in range(0,NE):
#         for j in range(0,3):
#             N = NOD[j, i]
#             if NPOCH[N] == B_INFLOW:
#                 inboundary[i] = inboundary[i] + 1
#             elif NPOCH[N] == B_OUTFLOW:
#                 outboundary[i] = outboundary[i] + 1
#
#         save8.write(str(inboundary[i]) + "\n")
#         save9.write(str(outboundary[i]) + "\n")


def point_type(num_of_tri):
    """
    判断三角形不同点的类型并返回数组
    :param num_of_tri: 三角形序号
    :return:
    """
    global NP, NE, NPOCH, NOD, NOE, PX, PY
    global B_INNER, B_WALL, B_INFLOW, B_OUTFLOW
    inner_point = []
    out_boundary = []
    in_boundary = []
    out_L = 0.
    in_L = 0.
    for i in range(0, 3):
        N = NOD[i, num_of_tri]
        if NPOCH[N] == B_INNER:
            inner_point.append(N)  # 存的点位置
        elif NPOCH[N] == B_INFLOW:
            in_boundary.append(N)
        elif NPOCH[N] == B_OUTFLOW:
            out_boundary.append(N)

    if (len(in_boundary) == 2):
        for i in range(0, 3):
            N = NOD[i, num_of_tri]
            if (N != in_boundary[0]) & (N != in_boundary[1]):
                in_L = nlMatrix[num_of_tri][i][0]

    if (len(out_boundary) == 2):
        for i in range(0, 3):
            N = NOD[i, num_of_tri]
            if (N != out_boundary[0]) & (N != out_boundary[1]):
                out_L = nlMatrix[num_of_tri][i][0]

    return inner_point, out_boundary, in_boundary, in_L, out_L


# ****************************************************************************
def sort(N1, N2, N3, i):
    """
    :param N1:
    :param N2:
    :param N3:
    :param i: 中心点
    :return:
    """
    global Nod_2, Nod_3, Column_2, Column_3
    if (i == N1):
        Nod_2, Nod_3 = N2, N3
        Column_2, Column_3 = 1, 2
    elif (i == N2):
        Nod_2, Nod_3 = N3, N1
        Column_2, Column_3 = 2, 0
    elif (i == N3):
        Nod_2, Nod_3 = N1, N2
        Column_2, Column_3 = 0, 1
    return Nod_2, Nod_3, Column_2, Column_3


def CORRECTION_EQUATION(iteration):  # 求fn
    TEMPCCA = np.zeros(PMAX, float)  # 保存每一次迭代后的结果
    for it in range(0, iteration):
        for p in range(0, NP):
            m1, m5, m6 = 0.0, 0.0, 0.0  # 分别对应(1)(5)(6)式
            in_boundary_list = []  # 保存入流边界三角形编号
            out_boundary_list = []  # 保存出流边界三角形编号
            for o in range(0, nmax):  # 计算相邻三角形
                n = triangleNumber[p][o][0]  # n为三角形编号
                if n == -1:
                    break
                i = triangleNumber[p][o][1]  # 该点为三角形的第几个点编号
                inner_point, out_boundary, in_boundary, in_L, out_L = point_type(n)
                nod_2, nod_3, column_2, column_3 = sort(NOD[0, n], NOD[1, n], NOD[2, n], p)  # 其余两点的编号及在三角形中的序号
                avg = (1 / 3) * (CCAJSR[p] + CCAJSR[nod_2] + CCAJSR[nod_3])  # fni为第 i 个三角形单元的三个顶点 n 时刻的平均值
                m5 = m5 + (1.0 + ((BCSQ * KDCSQ) / (KDCSQ + avg) ** 2)) * All_of_area[n]  # 计算（5）式
                if (NPOCH[nod_2] != B_INNER) and (NPOCH[nod_3] != B_INNER):  # f2i和f3i对应的点都是边界点，则跳过该三角形
                    if len(in_boundary) == 2:  # 若为入流边界
                        in_boundary_list.append(o)  # 记录该点的triangleNumber的下标
                    elif len(out_boundary) == 2:  # 若为出流边界
                        out_boundary_list.append(o)  # 记录该点的triangleNumber的下标
                    continue
                f2i, f3i = 0.0, 0.0
                if it == 0:  # 判断是否为第一次迭代
                    f2i, f3i = CCAJSR[nod_2], CCAJSR[nod_3]  # f2i, f3i取n时刻的值
                elif it > 0:
                    f2i, f3i = (CCAJSR[nod_2] + NEWCCA[nod_2]) / 2, (CCAJSR[nod_3] + NEWCCA[nod_3]) / 2
                fb2 = f2i * abcMatrix[n, column_2, 1]  # f2i*b2i
                fb3 = f3i * abcMatrix[n, column_3, 1]  # f3i*b3i
                fc2 = f2i * abcMatrix[n, column_2, 2]  # f2i*c2i
                fc3 = f3i * abcMatrix[n, column_3, 2]  # f3i*c3i
                Li = nlMatrix[n, i, 0]  # 控制边界长度
                m1 = m1 + (nlMatrix[n, i, 1] * (fb2 + fb3) + nlMatrix[n, i, 2] * (fc2 + fc3)) * Li
                m6 = m6 + (nlMatrix[n, i, 1] * abcMatrix[n, i, 1] + nlMatrix[n, i, 2] * abcMatrix[n, i, 2]) * Li
            m2 = DT * (K2 * Gn[p] - K1 * CCAJSR[p] * (F - Gn[p])) * CTRL_AREA[p]  # 计算（2）式
            n7, n8, n9, n10 = 0.0, 0.0, 0.0, 0.0  # 分别对应(7)(8)(9)(10)式
            if len(in_boundary_list) > 0:  # 入流边界数目大于0
                for q in range(0, len(in_boundary_list)):
                    o1 = in_boundary_list[q]
                    in_n = triangleNumber[p][o1][0]  # 取出三角形编号
                    inner_point1, out_boundary1, in_boundary1, in_L1, out_L1 = point_type(in_n)
                    in_Nod_2, in_Nod_3, in_Column_2, in_Column_3 = sort(NOD[0, in_n], NOD[1, in_n], NOD[2, in_n], p)
                    Lj = in_L1  # 入流边界长度
                    f2j, f3j = CCAJSR[in_Nod_2], CCAJSR[in_Nod_3]
                    if it > 0:  # 判断是否为第一次迭代
                        f2j, f3j = NEWCCA[in_Nod_2], NEWCCA[in_Nod_3]
                    n7 = n7 + DCAFSR * (CCAFSR - (f2j + f3j) * (1 / 3)) * Lj
                    n8 = n8 + DCAFSR * Lj
                TEMPCCA[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6 + n7 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 + n8 * DT * (1 / 3))  # 存入临时数组保存
                continue
            elif len(out_boundary_list) > 0:  # 出流边界数目大于0
                for r in range(0, len(out_boundary_list)):
                    o2 = out_boundary_list[r]
                    out_n = triangleNumber[p][o2][0]  # 三角形编号
                    inner_point2, out_boundary2, in_boundary2, in_L2, out_L2 = point_type(out_n)
                    out_Nod_2, out_Nod_3, out_Column_2, out_Column_3 = sort(NOD[0, out_n], NOD[1, out_n], NOD[2, out_n],p)
                    Lk = out_L2
                    f2k, f3k = CCAJSR[out_Nod_2], CCAJSR[out_Nod_3]
                    if it > 0:
                        f2k, f3k = NEWCCA[out_Nod_2], NEWCCA[out_Nod_3]
                    n9 = n9 + DCARYR * (CCAMYO - (f2k + f3k) * (1 / 3)) * Lk
                    n10 = n10 + DCARYR * Lk
                TEMPCCA[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6 + n9 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 - n10 * DT * (1 / 3))  # 存入临时数组保存
                continue
            TEMPCCA[p] = (DT * DCAJSR * m1 + m2 + CCAJSR[p] * m5 + 0.5 * CCAJSR[p] * DT * DCAJSR * m6) / (
                    m5 - 0.5 * DT * DCAJSR * m6)  # 若无出流入流边则存入临时数组保存
        for t in range(0, NP):
            NEWCCA[t] = TEMPCCA[t]  # 保存上一次迭代的结果
            TEMPCCA[t] = 0.0


def GN_EQUATION(iteration):  # 求gn
    TEMPF = np.empty(PMAX, float)  # 保存每一次迭代后的结果
    for it in range(0, iteration):
        for p in range(0, NP):
            m1, m2, m3, m4, m5, m6 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  # 分别对应(1)(2)(3)(4)(5)(6)式
            for o in range(0, nmax):
                n = triangleNumber[p][o][0]  # 三角形编号
                if n == -1:
                    break
                i = triangleNumber[p][o][1]  # 该点为三角形的第几个点编号
                inner_point, out_boundary, in_boundary, in_L, out_L = point_type(n)
                Nod_2, Nod_3, Column_2, Column_3 = sort(NOD[0, n], NOD[1, n], NOD[2, n], p)
                if (NPOCH[Nod_2] != B_INNER) and (NPOCH[Nod_3] != B_INNER):  # f2i和f3i对应的点都是边界点，则跳过该三角形
                    continue
                Li = nlMatrix[n, i, 0]
                g2i, g3i = Gn[Nod_2], Gn[Nod_3]
                if it > 0:
                    g2i, g3i = (Gn[Nod_2] + NEWF[Nod_2]) * 0.5, (Gn[Nod_3] + NEWF[Nod_3]) * 0.5
                gb2 = g2i * abcMatrix[n, Column_2, 1]  # g2i*b2i
                gb3 = g3i * abcMatrix[n, Column_3, 1]  # g2i*b2i
                gc2 = g2i * abcMatrix[n, Column_2, 2]  # g2i*b2i
                gc3 = g3i * abcMatrix[n, Column_3, 2]  # g2i*b2i
                m1 = m1 + (nlMatrix[n, i, 1] * (gb2 + gb3) + nlMatrix[n, i, 2] * (gc2 + gc3)) * Li
                m6 = m6 + (nlMatrix[n, i, 1] * abcMatrix[n, i, 1] + nlMatrix[n, i, 2] * abcMatrix[n, i, 2]) * Li
            m2 = DT * (K1 * CCAJSR[p] * (F - Gn[p]) - K2 * Gn[p]) * CTRL_AREA[p]
            m3 = Gn[p] * CTRL_AREA[p]
            m5 = CTRL_AREA[p]
            m1 = m1 * DT * DF
            m4 = m6 * 0.5 * DT * Gn[p] * DF
            m6 = m6 * 0.5 * DT * DF
            TEMPF[p] = (m1 + m2 + m3 + m4) / (m5 - m6)
        for t in range(0, NP):
            NEWF[t] = TEMPF[t]
            TEMPF[t] = 0.0


# ****************************************************************************
def AVERAGE_CCA():
    """
    求平均钙离子浓度
    """
    global TOTAL_CA, AVG_CA_JSR

    TOTAL_CA = 0.0
    for I in range(0, NP):
        TOTAL_CA = TOTAL_CA + CCAJSR[I] * CTRL_AREA[I]  # 每一点Ca的浓度X每一个点对应的三角形单元面积=每一点钙的总面积
    TOTAL_CA = TOTAL_CA / 3
    if TOTAL_AREA > 0.00000000000001:
        AVG_CA_JSR = TOTAL_CA / TOTAL_AREA  # JSR平均浓度


def AVERAGE_CCA_Gn():
    """
    求平均荧光钙浓度
    """
    global TOTAL_Gn, AVG_Gn_JSR

    TOTAL_Gn = 0.0
    for I in range(0, NP):
        TOTAL_Gn = TOTAL_Gn + Gn[I] * CTRL_AREA[I]  # 每一点Ca的浓度X每一个点对应的三角形单元面积=每一点钙的总面积
    TOTAL_Gn = TOTAL_Gn / 3
    if TOTAL_AREA > 0.00000000000001:
        AVG_Gn_JSR = TOTAL_Gn / TOTAL_AREA  # JSR平均浓度


def CA_CURRENT():
    """
    求钙离子电流
    """
    global ICARYR, ICAFSR, AVG_CA_JSR, AVG_Gn_JSR, LENGTH
    global B_INFLOW, B_OUTFLOW
    H_JSR = 30.0  # JSR的高度
    UNITEC = 1.610217733  # 一个电荷的常数
    MOLNUM = 6.0221367  # 摩尔常数，单位
    ICARYR = 0.0
    ICAFSR = 0.0
    CA_IN = 0.0
    CA_OUT = 0.0
    LENGTH = 0.0
    for I in range(0, NE):
        N1 = NOD[0, I]
        N2 = NOD[1, I]
        N3 = NOD[2, I]
        #  OUT CURRENT
        if (NPOCH[N1] == B_OUTFLOW) and (NPOCH[N2] == B_OUTFLOW):
            LENGTH = sqrt((PX[N1] - PX[N2]) ** 2 + (PY[N1] - PY[N2]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (CCAJSR[N1] + CCAJSR[N2]) / 2
        if (NPOCH[N3] == B_OUTFLOW) and (NPOCH[N2] == B_OUTFLOW):
            LENGTH = sqrt((PX[N3] - PX[N2]) ** 2 + (PY[N3] - PY[N2]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (CCAJSR[N3] + CCAJSR[N2]) / 2
        if (NPOCH[N1] == B_OUTFLOW) and (NPOCH[N3] == B_OUTFLOW):
            LENGTH = sqrt((PX[N1] - PX[N3]) ** 2 + (PY[N1] - PY[N3]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (CCAJSR[N1] + CCAJSR[N3]) / 2
        #  IN CURRENT
        if NPOCH[N1] == B_INFLOW and NPOCH[N2] == B_INFLOW:
            LENGTH = sqrt((PX[N1] - PX[N2]) ** 2 + (PY[N1] - PY[N2]) ** 2)
            CA_IN = CA_IN + LENGTH * (CCAFSR - (CCAJSR[N1] + CCAJSR[N2]) / 2)
        if NPOCH[N3] == B_INFLOW and NPOCH[N2] == B_INFLOW:
            LENGTH = sqrt((PX[N3] - PX[N2]) ** 2 + (PY[N3] - PY[N2]) ** 2)
            CA_IN = CA_IN + LENGTH * (CCAFSR - (CCAJSR[N3] + CCAJSR[N2]) / 2)
        if NPOCH[N1] == B_INFLOW and NPOCH[N3] == B_INFLOW:
            LENGTH = sqrt((PX[N1] - PX[N3]) ** 2 + (PY[N1] - PY[N3]) ** 2)
            CA_IN = CA_IN + LENGTH * (CCAFSR - (CCAJSR[N1] + CCAJSR[N3]) / 2)

    CA_OUT = CA_OUT * DCARYR * H_JSR
    CA_IN = CA_IN * DCAFSR * H_JSR
    ICARYR = CA_OUT * UNITEC * 2 * MOLNUM * (10.0 ** -11)
    ICAFSR = CA_IN * UNITEC * 2 * MOLNUM * (10.0 ** -11)


# **********************************************************************************


def not_empty(s):
    return s and s.strip()


def LOAD_GRIDINFO():
    """
    加载网格信息
    """
    global NP, NE, PMAX, EMAX, NVEX, count
    global PX, PY, NOD, NOE, NPOCH

    NP = 0  # 节点个数
    NE = 0  # 三角形单元的数量
    PMAX = 50000  # 点数
    EMAX = 100000  # 三角形单元数
    NVEX = 3  # 三角形顶点数

    NPOCH = np.empty(PMAX, int)  # 每个节点的属性
    NOD = np.empty([NVEX, EMAX], int)  # 这些节点组成的三角形单元的信息
    NOE = np.empty([NVEX, EMAX], int)  # 这个是每一个三角形单元的三个邻居单元的编号
    PX = np.empty(PMAX, float)  # 点的横坐标
    PY = np.empty(PMAX, float)  # 点的纵坐标

    count = 0
    gridt_file = open('GRID/gridt.dat', 'r')  # gridt.dat，这个是空间离散化的每一个离散的节点的坐标，第一行是节点的数量
    for line in gridt_file.readlines():  # 把空间里面离散的每一点的坐标的横坐标和纵坐标分别加载到PX数组和PY数组，把结点的数量存入到NP中
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        PX[count] = current_line[0]
        PY[count] = current_line[1]
        # print(PX[count], PY[count])
        count = count + 1
    NP = count  # NP为点的数量
    # print(NP)
    gridt_file.close()

    count = 0
    gridt_file = open('GRID/npoch.dat', 'r')  # npoch.dat为每个节点的属性：0代表内部节点，2代表入流边界或远场的节点，1代表壁面的节点，4代表出流边界节点
    for line in gridt_file.readlines():  # 把每个结点的属性存入到NPOCH数组中
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NPOCH[count] = current_line[0]
        # print(NPOCH[count])
        count = count + 1
    gridt_file.close()

    count = 0
    gridt_file = open('GRID/nod.dat', 'r')  # nod.dat，这个是gridt.dat中这些节点组成的三角形单元的信息，第一行是三角形单元的数量，后面是每一个三角元对应的三个节点的编号；
    for line in gridt_file.readlines():  # 把每一个三角形单元所对应的三角形单元的编号分别存入到数组NOD中，即NOD是每个三角形单元所对应的节点编号的数组
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NOD[0, count] = int(current_line[0]) - 1
        NOD[1, count] = int(current_line[1]) - 1
        NOD[2, count] = int(current_line[2]) - 1
        # print(NOD[0, count], NOD[1, count], NOD[2, count])
        count = count + 1
    NE = count  # NE为三角形单元的数量
    gridt_file.close()

    count = 0
    gridt_file = open('GRID/noe.dat',
                      'r')  # noe.dat，这个是每一个三角形单元的三个邻居单元的编号（所谓邻居指的是和这个单元有一条边重合的三角形单元），编号0表示这条边是区域边界，所以没有邻居
    for line in gridt_file.readlines():  # 把每一个三角形单元所对应的邻居单元的编号存入到数组NOE中
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NOE[0, count] = current_line[0]
        NOE[1, count] = current_line[1]
        NOE[2, count] = current_line[2]
        # print(NOE[0, count], NOE[1, count], NOE[2, count])
        count = count + 1
    gridt_file.close()


def main():
    print("执行LOAD_GRIDINFO()函数")
    LOAD_GRIDINFO()  # 网格的信息
    print("执行INITIAL_PARAMETER()函数")
    INITIAL_PARAMETER()  # 初始化参数
    print("执行DOLOOP_DIFFUSION()函数")
    DOLOOP_DIFFUSION()  # 循环，一步一步往前推


main()

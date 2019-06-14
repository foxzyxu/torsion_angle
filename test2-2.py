import os
import pandas as pd
import numpy as np
import math
import seaborn as sns

keep_point = ['C_0j', 'O_0i', 'C_09', 'C_07',
              'N_0k', 'C_0j', 'O_0i', 'C_09',
              'C_0s', 'N_0r', 'C_0l', 'N_0k',
              'C_0t', 'C_0s', 'N_0r', 'C_0l',
              'S_13', 'C_12', 'C_0v', 'C_0u',
              'O_14', 'S_13', 'C_12', 'C_0v',
              'C_1a', 'N_18', 'S_13', 'C_12',
              'C_1b', 'C_1a', 'N_18', 'S_13',
              'N_1f', 'C_1b', 'C_1a', 'N_18',
              'C_1g', 'N_1f', 'C_1b', 'C_1a' ]
reserve = list(set(keep_point))

#，读取数据，存在两个分子的结构，拆成两个处理
info = os.listdir()
info.remove('test2-2.py')
compound = pd.DataFrame()
for cif in info:
    temp = pd.read_table(cif, delim_whitespace=True, 
                         header = None, skiprows = 26)
    temp.index = temp[0].tolist()
    temp.drop([1,5], axis=1,inplace=True)
    temp = temp.loc[reserve]
    temp = temp.round(decimals = 5)
    n = temp.shape[0]
    if n == 38:
        temp = temp.reset_index(drop=True)
        temp1 = temp.loc[0::2]
        temp2 = temp.loc[1::2]
        temp1.index = temp1[0].tolist()
        temp2.index = temp2[0].tolist()
        temp1.drop([0], axis=1,inplace=True)
        temp2.drop([0], axis=1,inplace=True)
        compound = compound.append(temp1)
        compound = compound.append(temp2)
    if n == 19:
        temp.drop([0], axis=1,inplace=True)
        compound = compound.append(temp)
    #一共459*10个角需要计算


def torsion_angle(a,b,c,d):#计算二面角
    #平面向量
    Rba = [round(a[0]-b[0],5), round(a[1]-b[1],5), round(a[2]-b[2],5)]
    Rbc = [round(c[0]-b[0],5), round(c[1]-b[1],5), round(c[2]-b[2],5)]
    Rcb = [round(b[0]-c[0],5), round(b[1]-c[1],5), round(b[2]-c[2],5)]
    Rcd = [round(d[0]-c[0],5), round(d[1]-c[1],5), round(d[2]-c[2],5)]
    #法向量
    nabc = [round(Rba[1]*Rbc[2]-Rba[2]*Rbc[1],5),
            round(Rba[2]*Rbc[0]-Rba[0]*Rbc[2],5),
            round(Rba[0]*Rbc[1]-Rba[1]*Rbc[0],5)]
    nbcd = [round(Rcb[1]*Rcd[2]-Rcb[2]*Rcd[1],5),
            round(Rcb[2]*Rcd[0]-Rcb[0]*Rcd[2],5),
            round(Rcb[0]*Rcd[1]-Rcb[1]*Rcd[0],5)]
    #点积
    dot_abcd = round(round(nabc[0]*nbcd[0],5) + round(nabc[1]*nbcd[1],5) 
    + round(nabc[2]*nbcd[2],5),5)
    dot_abc  = round(round(nabc[0]*nabc[0],5) + round(nabc[1]*nabc[1],5) 
    + round(nabc[2]*nabc[2],5),5)
    dot_bcd  = round(round(nbcd[0]*nbcd[0],5) + round(nbcd[1]*nbcd[1],5) 
    + round(nbcd[2]*nbcd[2],5),5)
    #反三角函数
    val = round(math.sqrt(dot_abc*dot_bcd),5)
    cos = round(dot_abcd/val,5)
    #防止数据溢出
    if cos > 1:
        cos = 1
    if cos < -1:
        cos = -1
    angle = math.acos(cos)
    return round(angle,5)

#柔性角循环计算
res = []
for i in range(0,459):
    j = i*19
    single = compound.iloc[j:j+19]
    for k in range(0,10):
        l = k*4
        a = list(single.loc[keep_point[l]  ])
        b = list(single.loc[keep_point[l+1]])
        c = list(single.loc[keep_point[l+2]])
        d = list(single.loc[keep_point[l+3]])
        res.append(torsion_angle(a,b,c,d))
angle = np.reshape(res,(459,10))
#np.savetxt('angle_list2-2.csv', angle, delimiter = ',')
sns.distplot(angle[0], kde=False,axlabel='angle0')
sns.distplot(angle[1], kde=False,axlabel='angle1')
sns.distplot(angle[2], kde=False,axlabel='angle2')
sns.distplot(angle[3], kde=False,axlabel='angle3')
sns.distplot(angle[4], kde=False,axlabel='angle4')
sns.distplot(angle[5], kde=False,axlabel='angle5')
sns.distplot(angle[6], kde=False,axlabel='angle6')
sns.distplot(angle[7], kde=False,axlabel='angle7')
sns.distplot(angle[8], kde=False,axlabel='angle8')
sns.distplot(angle[9], kde=False,axlabel='angle9')
        
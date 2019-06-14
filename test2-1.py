import os
import pandas as pd
import numpy as np
import math
import seaborn as sns

#确定有用坐标
compound0 = pd.read_table('0.xyz', delim_whitespace=True, header = None, 
                     skiprows = 2)#读取文件
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
drop_point = list(set(compound0[0]) - set(keep_point))

#批量读取
info = os.listdir()
info.remove('test2-1.py')
compound = pd.DataFrame()
for xyz in info:
    temp = pd.read_table(xyz, delim_whitespace=True, 
                         header = None, skiprows = 2)
    temp.index = temp[0].tolist()
    temp.drop([0], axis=1,inplace=True)
    temp.drop(index=drop_point, inplace=True)
    temp = temp.round(decimals = 2)
    compound = compound.append(temp)
#print(compound1.info())


def torsion_angle(a,b,c,d):#计算二面角
    #平面向量
    Rba = [round(a[0]-b[0],2), round(a[1]-b[1],2), round(a[2]-b[2],2)]
    Rbc = [round(c[0]-b[0],2), round(c[1]-b[1],2), round(c[2]-b[2],2)]
    Rcb = [round(b[0]-c[0],2), round(b[1]-c[1],2), round(b[2]-c[2],2)]
    Rcd = [round(d[0]-c[0],2), round(d[1]-c[1],2), round(d[2]-c[2],2)]
    #法向量
    nabc = [round(Rba[1]*Rbc[2]-Rba[2]*Rbc[1],2),
            round(Rba[2]*Rbc[0]-Rba[0]*Rbc[2],2),
            round(Rba[0]*Rbc[1]-Rba[1]*Rbc[0],2)]
    nbcd = [round(Rcb[1]*Rcd[2]-Rcb[2]*Rcd[1],2),
            round(Rcb[2]*Rcd[0]-Rcb[0]*Rcd[2],2),
            round(Rcb[0]*Rcd[1]-Rcb[1]*Rcd[0],2)]
    #点积
    dot_abcd = round(round(nabc[0]*nbcd[0],2) + round(nabc[1]*nbcd[1],2) 
    + round(nabc[2]*nbcd[2],2),2)
    dot_abc  = round(round(nabc[0]*nabc[0],2) + round(nabc[1]*nabc[1],2) 
    + round(nabc[2]*nabc[2],2),2)
    dot_bcd  = round(round(nbcd[0]*nbcd[0],2) + round(nbcd[1]*nbcd[1],2) 
    + round(nbcd[2]*nbcd[2],2),2)
    #反三角函数
    val = round(math.sqrt(dot_abc*dot_bcd),2)
    cos = round(dot_abcd/val,2)
    #防止数据溢出
    if cos > 1:
        cos = 1
    if cos < -1:
        cos = -1
    angle = math.acos(cos)
    return round(angle,2)

#柔性角循环计算
res = []
for i in range(0,700):
    j = i*19
    single = compound.iloc[j:j+19]
    for k in range(0,10):
        l = k*4
        a = list(single.loc[keep_point[l]  ])
        b = list(single.loc[keep_point[l+1]])
        c = list(single.loc[keep_point[l+2]])
        d = list(single.loc[keep_point[l+3]])
        res.append(torsion_angle(a,b,c,d))
angle = np.reshape(res,(700,10))
#依次画图，当然也可以循环加正则表达作图
#np.savetxt('angle_list.csv', angle_list, delimiter = ',')
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
    
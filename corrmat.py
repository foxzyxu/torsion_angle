import pandas as pd

angle1 = pd.read_csv('angle_list2-1.csv', header = None)
angle2 = pd.read_csv('angle_list2-2.csv', header = None)
#依次画图，当然也可以循环作图
angle = pd.concat([angle1, angle2], axis=1)
corrmat = angle.corr()

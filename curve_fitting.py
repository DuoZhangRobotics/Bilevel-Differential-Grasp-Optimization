import numpy as np
import os
# BarrettHand data
data = [1.03552, 0.332847,
        1.05205, 0.33816,
        0.991532, 0.318707,
        1.52492, 0.490154,
        0.33841, 0.108775,
        0,0,
        1.16725, 0.375188,
        0.0776315, 0.024953,
        ]
data = []
file_list = os.listdir('./shadowhand_parameters')
cnt = 0
for file in file_list:
	with open('./shadowhand_parameters/'+file,'r') as f:
		initial_paramters =[*map(float, f.readline().split())]
		data.append(initial_paramters[8])
		data.append(initial_paramters[9])
		data.append(initial_paramters[12])
		data.append(initial_paramters[13])
		data.append(initial_paramters[16])
		data.append(initial_paramters[17])
		data.append(initial_paramters[21])
		data.append(initial_paramters[22])
	cnt += 1
	if cnt == 3:
		break

data = np.array(data).reshape((-1, 2))
print(data)
x = data[:, 0]
y = data[:, 1]
# print(x, y)
A = np.vstack([x, np.ones((x.shape[0]))]).T
m, c = np.linalg.lstsq(A, y, rcond=None)[0]
print("m = ", m)
print(f'c = {c}')
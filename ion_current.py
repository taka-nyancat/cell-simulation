import math
import numpy as np
from matplotlib import pyplot as plt
gcabar = 0.005
eca = 120
t    = np.linspace(0, 4, 40000) 
freq = 0.5

# v_in = 1 * (np.sin(2 * np.pi * freq * (t)))
# a = []
# for i in v_in:
#     if i>0:
#         a.append(i)
#     else:
#         a.append(0)
# v_in = a


v_in = -50 + 20* (np.sin(2 * np.pi * freq * (t - 0.001)))
# print(type(v_in))
# v_in = []
# s = 1
# e = 1.5
# m = 1.25
# for i in t:
#     an = 40 * np.exp(-((i-m)**2)/(2*0.1**2)) -60
#     v_in.append(an)
# v_in = np.array(v_in)

s_inf = 1/(1+np.exp(-(v_in+10)/2))
ica = gcabar * s_inf * (v_in-eca)

# plt.plot(t, v_in, label="$V_{in}$")
plt.plot(t, ica, label="$I_{ca}$")
plt.xlabel('Time[s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.show()
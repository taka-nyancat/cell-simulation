import neuron
from neuron import h
import numpy as np
from pprint import pprint
import tqdm
import matplotlib.pyplot as plt

tstop = 3500.0
v_init = -65.0

# retinafile = "swc"

# myRetina = np.loadtxt(retinafile)
# myRetina[:,2:6]=0.001*8*myRetina[:,2:6]
# interim=myRetina.shape
# ret_nofcomps=interim[0]
# compdiam=myRetina[:,5]*2.0

# Retina = {}
# print("start loading medulafile")
# for i in tqdm.tqdm(range(1, ret_nofcomps, 1)):
#     aind=int(myRetina[i,0]-1)
#     bind=int(myRetina[i,6]-1)
#     axyz=myRetina[aind,2:5]
#     bxyz=myRetina[bind,2:5]
#     # print("aind:{}, bind:{}".format(aind, bind))

#     complength=np.sqrt(np.sum((axyz-bxyz)**2))
#     meandiam=(compdiam[aind]+compdiam[bind])*0.5
#     # print("aind:{}, bind:{}, axyz:{}, bxyz:{}, complength:{}, meandiam:{}".format(aind, bind, axyz, bxyz, complength, meandiam))
#     Retina[str(aind)] = h.Section()
#     Retina[str(aind)].diam = meandiam
#     Retina[str(aind)].L = complength
#     Retina[str(aind)].insert ("mole")
#     Retina[str(aind)].cm = 0.6
#     Retina[str(aind)].Ra = 0.1

#     if bind == 0:
#         continue
    
#     Retina[str(aind)].connect(Retina[str(bind)], 1)
# print("end loading medulafile")


laminafile = "swc/L1 home_3529395.swc"

myLamina = np.loadtxt(laminafile)
myLamina[:,2:6]=0.001*8*myLamina[:,2:6]
interim=myLamina.shape
la_nofcomps=interim[0]
compdiam=myLamina[:,5]*2.0

Lamina = {}
print("start loading medulafile")
for i in tqdm.tqdm(range(1, la_nofcomps, 1)):
    aind=int(myLamina[i,0]-1)
    bind=int(myLamina[i,6]-1)
    axyz=myLamina[aind,2:5]
    bxyz=myLamina[bind,2:5]
    # print("aind:{}, bind:{}".format(aind, bind))

    complength=np.sqrt(np.sum((axyz-bxyz)**2))
    meandiam=(compdiam[aind]+compdiam[bind])*0.5
    # print("aind:{}, bind:{}, axyz:{}, bxyz:{}, complength:{}, meandiam:{}".format(aind, bind, axyz, bxyz, complength, meandiam))
    Lamina[str(aind)] = h.Section()
    Lamina[str(aind)].diam = meandiam
    Lamina[str(aind)].L = complength
    Lamina[str(aind)].insert ("mole")
    Lamina[str(aind)].cm = 0.6
    Lamina[str(aind)].Ra = 0.1

    if bind == 0:
        continue
    
    Lamina[str(aind)].connect(Lamina[str(bind)], 1)
lamina_post_synapse = h.gsyn2(0.5, sec = Lamina[str(10)])
lamina_post_synapse.k = 0.05
lamina_post_synapse.numsyn = 200
lamina_post_synapse.vre = -80
lamina_post_synapse.vth = -48
print("end loading medulafile")

medulafile = "swc/Tm1 home_226.swc"

myMedulla = np.loadtxt(medulafile)
myMedulla[:,2:6]=0.001*8*myMedulla[:,2:6]
interim=myMedulla.shape
med_nofcomps=interim[0]
compdiam=myMedulla[:,5]*2.0

Medulla = {}
print("start loading medulafile")
for i in tqdm.tqdm(range(1, med_nofcomps, 1)):
    aind=int(myMedulla[i,0]-1)
    bind=int(myMedulla[i,6]-1)
    axyz=myMedulla[aind,2:5]
    bxyz=myMedulla[bind,2:5]
    # print("aind:{}, bind:{}".format(aind, bind))

    complength=np.sqrt(np.sum((axyz-bxyz)**2))
    meandiam=(compdiam[aind]+compdiam[bind])*0.5
    # print("aind:{}, bind:{}, axyz:{}, bxyz:{}, complength:{}, meandiam:{}".format(aind, bind, axyz, bxyz, complength, meandiam))
    Medulla[str(aind)] = h.Section()
    Medulla[str(aind)].diam = meandiam
    Medulla[str(aind)].L = complength
    Medulla[str(aind)].insert ("mole")
    Medulla[str(aind)].cm = 0.6
    Medulla[str(aind)].Ra = 0.1

    if bind == 0:
        continue
    
    Medulla[str(aind)].connect(Medulla[str(bind)], 1)
medula_post_synapse = h.gsyn2(0.5, sec = Medulla[str(10)])
# medula_post_synapse.tau1 = 0.5
# medula_post_synapse.tau2 = 1
# medula_post_synapse.e = 0
# medula_post_synapse.k = 0.05
# medula_post_synapse.numsyn = 200
# medula_post_synapse.vre = -80
# medula_post_synapse.vth = -48
print("end loading medulafile")

nclist = []
nclist.append(h.NetCon(Lamina[str(la_nofcomps-1)](0.5)._ref_v, None, Lamina[str(la_nofcomps-1)](0.5)))
synlist = []


ic = h.SinCurrent(0.5, sec=Lamina[str(1)])
ic.st = 0
ic.en = 3500 # ms
ic.amp = 50 # nA
ic.offset = 50.0
ic.delay = 0.0
ic.freqency = 0.001

cvode = h.CVode ()
cvode.active (1)
cvode.atol (1.0e-5)

vv = h.Vector() # membrane potential vector
vv1 = h.Vector() # membrane potential vector
vv2 = h.Vector() # membrane potential vector
cv = h.Vector ()
tv = h.Vector() # time stamp vector

vv.record(Lamina[str(la_nofcomps-1)](0.5)._ref_v)
vv1.record(Medulla[str(10)](0.5)._ref_v)
vv2.record(Medulla[str(med_nofcomps-1)](0.5)._ref_v)
cv.record(ic._ref_i)
tv.record(h._ref_t)

h.finitialize(v_init)
h.fcurrent()
neuron.run(tstop)
fig = plt.figure(figsize=(18, 18), dpi=60)
ax = fig.add_subplot(2, 2, 1)
ax.set_ylim([-80,40])
ax.plot(tv.as_numpy(), vv.as_numpy())
ax.set_xlabel("Time(ms)")
ax.set_ylabel("Potential(mV)")
ax.set_title("end")

ax1 = fig.add_subplot(2, 2, 2)
ax1.set_ylim([-80,40])
ax1.plot(tv.as_numpy(), vv1.as_numpy())
ax1.set_xlabel("Time(ms)")
ax1.set_ylabel("Potential(mV)")
ax1.set_title("start")

ax2 = fig.add_subplot(2, 2, 3)
ax2.set_ylim([-80,40])
ax2.plot(tv.as_numpy(), vv2.as_numpy())
ax2.set_xlabel("Time(ms)")
ax2.set_ylabel("Potential(mV)")
ax2.set_title("middle")

sx = fig.add_subplot(2, 2, 4)
sx.set_ylim([-80,100])
sx.plot(tv.as_numpy(), cv.as_numpy())
sx.set_xlabel("Time(ms)")
sx.set_ylabel("Potential(mV)")
sx.set_title("input")

plt.savefig("figs/n02.png")
plt.show()
self.pc.target_var(syn_obj,getattr(syn_obj,"_ref_" + con["target_synapse"]["value"]),index)

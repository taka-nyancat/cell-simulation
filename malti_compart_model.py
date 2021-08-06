import neuron
from neuron import h
import numpy as np
from pprint import pprint
from numpy.core.fromnumeric import argmax
import tqdm
import matplotlib.pyplot as plt


def insertMod(sec):
    r = 4
    sec.insert("hh")
    sec.gnabar_hh = 0.190 * r
    sec.gkbar_hh = 0.060 * r 
    sec.gl_hh = 0.0001 * 0.0001 * 0.3 * r
    sec.el_hh = -67
    sec.ena = 50
    sec.ek = -100

    sec.insert("IAHP2")
    sec.gahpbar_IAHP2 = 0.004 * 8
    sec.eahp_IAHP2 = -140
    sec.a_IAHP2 = 1.0 
    sec.b_IAHP2 = 0.0006 
    sec.co_IAHP2 = 2.0 	

    sec.insert("cad")
    sec.ctau_cad = 0.0124 * 7
    sec.cainf_cad = 2.4e-6
    sec.a_cad = 0.002 * 0.166
    sec.Bt_cad = 6.0e-3 * 5
    sec.kBin_cad = 5000
    sec.kBout_cad = 0.04 * 0.05
    sec.bbrini_cad = 0.01

    sec.insert("ICA")
    sec.eca = 120
    sec.gcabar_ICA = 0.005 * 1.0

    sec.Ra = 0.001

    print(sec.cai)

tstop = 3500.0
v_init = -65.0

swcfile = "swc/Tm1 home_226.swc"
cellname = swcfile[4:-4]
para_list = {"model": "mole", "cm": 0.6, "Ra": 0.1, "gcabar": 1.1}
filename = cellname + ",model:{},cm:{},Ra:{}, gca:{}".format(para_list["model"], para_list["cm"], para_list["Ra"], para_list["gcabar"])

mycell = np.loadtxt(swcfile)
mycell[:,2:6]=0.001*8*mycell[:,2:6]
interim=mycell.shape
nofcomps=interim[0]
compdiam=mycell[:,5]*2.0
print("nofcomps:{}".format(nofcomps))

cell = {}
print("start loading swcfile")
for i in tqdm.tqdm(range(1, nofcomps, 1)):
    aind=int(mycell[i,0]-1)
    bind=int(mycell[i,6]-1)
    axyz=mycell[aind,2:5]
    bxyz=mycell[bind,2:5]
    print("aind:{}, bind:{}".format(aind, bind))

    complength=np.sqrt(np.sum((axyz-bxyz)**2))
    meandiam=(compdiam[aind]+compdiam[bind])*0.5
    if complength == 0:
        print("ng")
    print(complength)
    # print("aind:{}, bind:{}, axyz:{}, bxyz:{}, complength:{}, meandiam:{}".format(aind, bind, axyz, bxyz, complength, meandiam))
    cell[str(aind)] = h.Section()
    cell[str(aind)].diam = meandiam
    cell[str(aind)].L = complength
    # cell[str(aind)].insert (para_list["model"])
    # cell[str(aind)].cm = para_list["cm"]
    # cell[str(aind)].Ra = para_list["Ra"]
    insertMod(cell[str(aind)])

    if bind == 0:
        continue
    
    cell[str(aind)].connect(cell[str(bind)], 1)

print("end loading swcfile")

ic = h.SinCurrent(0.5, sec=cell[str(1)])
ic.st = 0
ic.en = 3500 # ms
ic.amp = 0.5 # nA
ic.offset = 0.5
ic.delay = 0.0
ic.freqency = 0.001

cvode = h.CVode ()
cvode.active (1)
cvode.atol (1.0e-5)

vend= h.Vector() # membrane potential vector
vst = h.Vector() # membrane potential vector
vmid = h.Vector() # membrane potential vector
inj = h.Vector ()
t = h.Vector() # time stamp vector
cain = h.Vector()

vend.record(cell[str(nofcomps-1)](0.5)._ref_v)
vst.record(cell[str(10)](0.5)._ref_v)
vmid.record(cell[str(nofcomps//2)](0.5)._ref_v)
inj.record(ic._ref_i)
t.record(h._ref_t)
cain.record(cell[str(nofcomps-1)](0.5)._ref_cai)



h.finitialize(v_init)
h.fcurrent()
neuron.run(tstop)

vend_l = vend.as_numpy()
vst_l = vst.as_numpy()
vmid_l = vmid.as_numpy()
inj_l = inj.as_numpy()
t_l = t.as_numpy()
cain_l = cain.as_numpy()

recored_list = {"inj point":{"v_list":inj_l, "maxtime": t_l[inj_l.argmax()], "maxv": max(inj_l), "minv":min(inj_l)},
    "start point":{"v_list":vst_l, "maxtime": t_l[vst_l.argmax()], "maxv": max(vst_l), "minv": min(vst_l)},
    "middle point":{"v_list":vmid_l, "maxtime": t_l[vmid_l.argmax()], "maxv": max(vmid_l), "minv": min(vmid_l)},
    "end point":{"v_list":vend_l, "maxtime": t_l[vend_l.argmax()], "maxv": max(vend_l), "minv": min(vend_l)},
    "cal":{"v_list":cain_l}
}

fig = plt.figure(figsize=(18, 18), dpi=60)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t_l, cain_l)

ax = []
plot_pos = 0
time_text_list = []
text = ""

fig = plt.figure(figsize=(18, 18), dpi=60)
for recored_point, data in recored_list.items():
    ax.append(fig.add_subplot(3, 2, plot_pos+1))
    # ax[plot_pos].set_ylim([data["minv"] -5, data["maxv"] + 5])
    ax[plot_pos].plot(t, data["v_list"])
    ax[plot_pos].set_xlabel("Time(ms)")
    ax[plot_pos].set_ylabel("Potential(mV)")
    ax[plot_pos].set_title(recored_point)
    # time_text_list.append("{} max time : {} msec\n".format(recored_point, data["maxtime"]))
    plot_pos += 1

# for com in time_text_list:
#     text += com
# deley = recored_list["end point"]["maxtime"] - recored_list["start point"]["maxtime"]
# add_text = "deley time : {}".format(deley)
# text += add_text
# fig.text(0.4, 0.03, text)

plt.savefig("figs/" + filename + ".png")
plt.show()

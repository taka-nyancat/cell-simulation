import neuron
from neuron import h
import numpy as np
import matplotlib.pyplot as plt

def insertMod(sec):
    # r = 1
    # sec.insert("hh")
    # print("default gna:{}".format(sec.gnabar_hh))
    # sec.gnabar_hh = 0.190 * r
    # print("fixed gna:{}".format(sec.gnabar_hh))
    # print("------------------------------------")
    # print("default gk:{}".format(sec.gkbar_hh))
    # sec.gkbar_hh = 0.060 * r 
    # print("fixed gk:{}".format(sec.gkbar_hh))
    # print("------------------------------------")
    # print("default gl:{}".format(sec.gl_hh))
    # sec.gl_hh = 0.0001 * 0.0001 * 0.3 * r
    # print("fixed gl:{}".format(sec.gl_hh))
    # print("------------------------------------")
    # print("default el:{}".format(sec.el_hh))
    # sec.el_hh = -67
    # print("fixed el:{}".format(sec.el_hh))
    # print("------------------------------------")
    # print("default ena:{}".format(sec.ena))
    # sec.ena = 50
    # print("fixed ena:{}".format(sec.ena))
    # print("------------------------------------")
    # print("default ek:{}".format(sec.ek))
    # sec.ek = -100
    # print("fixed ek:{}".format(sec.ek))
    # print("------------------------------------")

    sec.insert("mole")
    sec.cm = 0.8
    sec.Ra = 0.01

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

    # sec.insert("ICA")
    # sec.eca = 120
    # sec.gcabar_ICA = 0.005 * 1.0

# list_ = [i for i in range(10)]
# print("list_:{}".format(list_))
# print("list_ type:{}".format(type(list_)))
# list_.as_numpy()
# print("list_:{}".format(list_))
# print("list_ type:{}".format(type(list_)))


tstop = 1200.0
v_init = -65.0

soma = h.Section()
soma.diam = 0.4
soma.L = 140.0
soma.nseg = 29
insertMod(soma)

ic = h.IClamp(0.5, sec=soma)
ic.delay = 1000.0 # ms
ic.dur = 1 # ms
ic.amp = 50 # nA

# ic = h.SinCurrent(0.1, sec=soma)
# ic.st = 0
# ic.en = 3500 # ms
# inp = 50
# ic.amp = inp # nA
# ic.offset = inp-20
# ic.delay = 0.0
# ic.freqency = 0.001

cvode = h.CVode ()
cvode.active (1)
cvode.atol (1.0e-5)

vv = h.Vector() # membrane potential vector
tv = h.Vector() # time stamp vector
cai = h.Vector()
inj = h.Vector()
vv.record(soma(0.5)._ref_v)
vv_list = vv.as_numpy()
tv.record(h._ref_t)
cai.record(soma(0.5)._ref_cai)
inj.record(ic._ref_i)
print("original type : {}".format(type(vv)))


h.finitialize(v_init)
h.fcurrent()
neuron.run(tstop)

st = 500

fig = plt.figure(figsize=(18, 18), dpi=60)
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(tv[st:], inj[st:])
ax1.set_xlabel("Time(ms)")
ax1.set_ylabel("Potential(nA)")
ax1.set_title("input current")
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(tv[st:], vv[st:])
ax2.set_xlabel("Time(ms)")
ax2.set_ylabel("Potential(mV)")
ax2.set_title("point0.5 potential")

ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(tv[st:], cai[st:])
ax3.set_xlabel("Time(ms)")
ax3.set_ylabel("molar(mM)")
ax3.set_title("calcium molar")

plt.savefig("figs/calmod.png")
plt.show()
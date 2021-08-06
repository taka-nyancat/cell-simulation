import neuron
from neuron import h
import numpy as np
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

tstop = 1000.0
v_init = -65.0

soma = h.Section()
soma.diam = 0.4
soma.L = 140.0
soma.nseg = 29
insertMod(soma)

# ic = h.IClamp(0.5, sec=soma)
# ic.delay = 100.0 # ms
# ic.dur = 100 # ms
# ic.amp = 0 # nA

cvode = h.CVode ()
cvode.active (1)
cvode.atol (1.0e-5)

vv = h.Vector() # membrane potential vector
tv = h.Vector() # time stamp vector
cai = h.Vector()
# inj = h.Vector()
vv.record(soma(0.5)._ref_v)
tv.record(h._ref_t)
cai.record(soma(0.5)._ref_cai)
# inj.record(ic._ref_i)

h.finitialize(v_init)
h.fcurrent()
neuron.run(tstop)

fig = plt.figure(figsize=(18, 18), dpi=60)
# ax1 = fig.add_subplot(2, 2, 1)
# ax1.plot(tv.as_numpy(), inj.as_numpy())
# ax1.set_xlabel("Time(ms)")
# ax1.set_ylabel("Potential(mV)")
# ax1.set_title("ic")

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(tv.as_numpy(), vv.as_numpy())
ax2.set_xlabel("Time(ms)")
ax2.set_ylabel("Potential(mV)")
ax2.set_title("vv")

ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(tv.as_numpy(), cai.as_numpy())
ax3.set_xlabel("Time(ms)")
ax3.set_ylabel("Potential(mV)")
ax3.set_title("cai")

plt.savefig("figs/n02.pdf")
plt.show()
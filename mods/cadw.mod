TITLE Fast mechanism for submembranal Ca++ concentration (cai)
:
: Takes into account:
:
:	- increase of cai due to calcium currents
:	- extrusion of calcium with a simple first order equation
:
: This mechanism is compatible with the calcium pump "cad" and has the 
: same name and parameters; however the parameters specific to the pump
: are dummy here.
:
: Parameters:
:
:	- depth: depth of the shell just beneath the membran (in um)
:	- cainf: equilibrium concentration of calcium (2e-4 mM)
:	- taur: time constant of calcium extrusion (must be fast)
:	- kt,kd: dummy parameters
:
: Written by Alain Destexhe, Salk Institute, 1995
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica WRITE cai
      :	RANGE depth,kt,kd,cainf,taur
	RANGE depth,kt,kd,taur
	RANGE ctau, cainf, a, Bt, kBin, kBout, bbrini
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT {
	FARADAY = 96489		(coul)		: moles do not appear in units
:	FARADAY = 96.489	(k-coul)	: moles do not appear in units
}

PARAMETER {
	depth	= .1	(um)		: depth of shell
	:taur	= 5	(ms)		: rate of calcium removal
	ctau    = 0.0000125 (1/ms)
	cainf	= 2.4e-4	(mM)
	kt	= 0	(mM/ms)		: dummy
	kd	= 0	(mM)		: dummy
	a	= 0.002	
    Bt      = 1.2e-2  (mM)           :total ca buffer =50(binding ratio)x cainf
    kBout  = 1.02e    (mM/ms)           :determin ca increase from B
    kBin   =0.05      (mM/ms/mM)    : time const =20ms order
    bbrini = 0.99

}
 
STATE {
	cai		(mM)
    bbr             :bind ratio in ca-buffer 

        
       
}

INITIAL {
	cai = cainf
    :bbr = 0.98
	: 49/50
    bbr=bbrini

}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
    pump_ch (mM/ms)
}
	
BREAKPOINT {
    SOLVE state METHOD derivimplicit    
}

DERIVATIVE state { 

	:drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	drive_channel = - a * ica
	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward  > not pump
    pump_ch= (cai-cainf)*ctau*(-1.0)
    if (pump_ch >=0.) {pump_ch = 0 } : cannnot pump inward  
	:cai' = drive_channel + (cainf-cai)*ctau
    cai' = drive_channel + pump_ch + 2*kBout*Bt*bbr - cai*cai*kBin*Bt*(1-bbr) :new ca dynamics
    bbr'= -kBout*bbr+cai*cai*kBin*(1-bbr)   :new buffer dynamics
}
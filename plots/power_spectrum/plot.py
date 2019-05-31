import analytics.shot_noise_driven.lif_neuron as ana
import analytics.gaussian_white_noise_driven.lif_neuron as gwnana

import numpy as np
import param_scan.io
import param_scan.simulation_run
from param_scan.io import InDirectory
import os
from latex_param_values import LatexParamValues
from grace_plot import GracePlot

def with_tau_m(tau_m, prms):
    """ convert back from dimnesionless units """
    p = dict(prms)
    p["tau_m"] = tau_m
    p["rin_e"] = prms["rin_e"] / tau_m
    p["tr"] = prms["tr"] * tau_m
    p["df"] = prms["df"] / tau_m
    p["dt"] = prms["dt"] * tau_m
    p["f_c"] = prms["f_c"] / tau_m
    p["f_max"] = prms["f_max"] / tau_m
    p["f_sig"] = prms["f_sig"] / tau_m
    p["r_sample"] = prms["r_sample"] / tau_m
    return p

tau_m = 0.02 # s

with InDirectory(os.path.abspath(__file__)):

    rfn, r = param_scan.io.load_newest_run()
    prms = r["parameters"]
    gr = GracePlot("plot")
    try:
        s = np.loadtxt(param_scan.simulation_run.get_output_path(r, prms, "sxx")).T
        gr.plot(s[0][1:]/tau_m, s[1][1:]/tau_m)
    except IOError:
        gr.plot([0],[0])
    fs = np.linspace(1e-4, 10., 1000)
    gr.plot(fs/tau_m, ana.powspec(prms, fs=fs)/tau_m)
    gr.plot(fs/tau_m, gwnana.powspec(prms, mu=prms["mu"]+prms["a_e"]*prms["rin_e"], D=prms["rin_e"]*prms["a_e"]**2, fs=fs)/tau_m)
    #gr.plot(fs, ana.r0(prms) * (1+2*np.real([ana.st_rate_d(prms, f=f, tau=1, a_i=0, rin_i=0) for f in fs])))
    gr.save() 
    LatexParamValues().write("paramvalues.tex", with_tau_m(tau_m * 1000, prms))

    

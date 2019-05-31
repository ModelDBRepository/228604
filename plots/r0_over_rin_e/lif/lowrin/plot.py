import analytics.shot_noise_driven.if_neuron as ana
import analytics.gaussian_white_noise_driven.if_neuron as gwnana
import numpy as np
import param_scan.io
import param_scan.parameter_sets
import param_scan.simulation_run
from param_scan.io import InDirectory
import os
from latex_param_values import LatexParamValues
from grace_plot import GracePlot
from math import log10

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
    for pp in param_scan.parameter_sets.unroll(prms, axes="expweights"):
        for p in param_scan.parameter_sets.unroll(pp, axes="a_e_cutoff"):
            re, r0 = param_scan.simulation_run.read_values(r, param_scan.parameter_sets.unroll(p), ["rin_e", "r0"], "stdout", ignore_errors=True)
            re, r0 = map(np.array, [re, r0])
            gr.plot(re/tau_m, r0/tau_m)
    theo_res = np.logspace(log10(min(re)), log10(max(re)), 50)
    models = {"pif": ana.PIF(), "lif": ana.LIF(), "qif": ana.QIF()}
    gwnmodels = {"pif": gwnana.PIF(), "lif": gwnana.LIF(), "qif": gwnana.QIF()}
    gr.plot(theo_res/tau_m, [ana.r0(prms, model=models[prms["model"]], rin_e=re)/tau_m for re in theo_res])
    gr.plot(theo_res/tau_m, [gwnana.r0(prms, model=gwnmodels[prms["model"]], mu=prms["mu"]+prms["a_e"]*re, D=re*prms["a_e"]**2)/tau_m for re in theo_res])

    gr.plot(theo_res/tau_m, np.minimum(theo_res, 1./prms["tr"])/tau_m) 
    gr.plot(theo_res/tau_m, theo_res * np.exp(-(prms["vt"]-prms["mu"])/prms["a_e"])/tau_m)

    gr.save()
    LatexParamValues().write("paramvalues.tex", with_tau_m(tau_m*1000, p))

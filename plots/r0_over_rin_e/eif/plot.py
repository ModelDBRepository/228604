import analytics.shot_noise_driven.if_neuron as ana
import analytics.gaussian_white_noise_driven.if_neuron as gwnana
import pylab as pl
import numpy as np
import param_scan.io
import param_scan.parameter_sets
import param_scan.simulation_run
from param_scan.io import InDirectory
import os
from latex_param_values import LatexParamValues
from grace_plot import GracePlot
import cProfile

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
    gr = GracePlot("plot", rows=3)
    gr.focus(row=1)
    gr.plot([-5,5],[0,0])
    #prms = r["parameters"]
    for prms in param_scan.parameter_sets.unroll(r["parameters"], axes="d"):
        gr.focus(row=0)
        re, r0 = param_scan.simulation_run.read_values(r, param_scan.parameter_sets.unroll(prms), ["rin_e", "r0"], "stdout")
        re, r0 = map(np.array, [re, r0])
        gr.plot(re/tau_m, r0/tau_m, legend=r"simulation (\xD\f{} = %g)" % prms["d"])
        theo_res = np.linspace(min(re), max(re), 50)
        models = {"pif": ana.PIF(), "lif": ana.LIF(), "qif": ana.QIF(), "eif": ana.EIF(d=prms["d"], vtb=prms["vtb"])}
        gwnmodels = {"pif": gwnana.PIF(), "lif": gwnana.LIF(), "qif": gwnana.QIF(), "eif": gwnana.EIF(d=prms["d"], vtb=prms["vtb"])}
        anar0s = []
        

        gr.plot(theo_res/tau_m, [ana.r0(prms, model=models[prms["model"]], rin_e=re, vt=prms["vtb"]+100*prms["d"])/tau_m for re in theo_res], legend=r"theory (\xD\f{} = %g)" % prms["d"])
        gr.plot(theo_res/tau_m, [gwnana.r0(prms, model=gwnmodels[prms["model"]], mu=prms["mu"]+prms["a_e"]*re, D=re*prms["a_e"]**2, vt=prms["vtb"]+100*prms["d"])/tau_m  for re in theo_res], legend=r"diffusion approximation (\xD\f{} = %g)" % prms["d"])
        
        gr.focus(row=1)
        vs = np.linspace(-5,5,5000)
        gr.plot(vs, models[prms["model"]].f(vs, prms["mu"]))
        gr.plot([models[prms["model"]].intervals(prms["mu"], prms["vr"], prms["vt"])[-1][0]], [0])
        gr.focus(row=2)
        gr.plot(vs, gwnmodels[prms["model"]].U(vs, prms["mu"]))
        gr.text(r"\xD\f{} = %g" % prms["d"])
    gr.save()
    LatexParamValues().write("paramvalues.tex", with_tau_m(tau_m*1000, prms))


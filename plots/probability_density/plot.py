import analytics.shot_noise_driven.if_neuron as ana
import analytics.gaussian_white_noise_driven.if_neuron as gwnana
import pylab as pl
import numpy as np
import param_scan
import param_scan.io
import param_scan.simulation_run
from param_scan.io import InDirectory
import os
from grace_plot import GracePlot
from latex_param_values import LatexParamValues

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


with InDirectory(os.path.abspath(__file__)):

    rfn, r = param_scan.io.load_newest_run()
    gr = GracePlot("plot")

    #prms = r["parameters"]
    for prms in param_scan.parameter_sets.unroll(r["parameters"]):
        print r
        h = np.loadtxt(param_scan.simulation_run.get_output_path(r, prms, "vhist")).T
        gr.plot(h[0]-(h[0][1]-h[0][0])/2, h[1])
        vs = np.linspace(prms["vhist_l"], prms["vhist_r"], 500)
        models = {"pif": ana.PIF(), "lif": ana.LIF(), "qif": ana.QIF(), "eif": ana.EIF(d=prms["d"], vtb=prms["vtb"])}
        gr.plot(vs, [ana.P0(prms, model=models[prms["model"]], v=v, vt=prms["vtb"]+100*prms["d"]) for v in vs])
        gwnmodels = {"pif": gwnana.PIF(), "lif": gwnana.LIF(), "qif": gwnana.QIF(), "eif": gwnana.EIF(d=prms["d"], vtb=prms["vtb"])} 
        gr.plot(vs, [gwnana.P0(prms, mu=prms["mu"]+prms["a_e"]*prms["rin_e"], D=prms["a_e"]**2*prms["rin_e"], model=gwnmodels[prms["model"]], v=v, vt=prms["vtb"]+50*prms["d"]) for v in vs])
        print "r0 =", ana.r0(prms, model=models[prms["model"]], vt=prms["vtb"]+50*prms["d"])
        #print "da_r0 =", gwnana.r0(prms, mu=prms["mu"]+prms["a_e"]*prms["rin_e"], D=prms["a_e"]**2*prms["rin_e"], model=gwnmodels[prms["model"]])
    gr.save()
    LatexParamValues().write("paramvalues.tex", with_tau_m(20, prms))



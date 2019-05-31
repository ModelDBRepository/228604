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

with InDirectory(os.path.abspath(__file__)):

    rfn, r = param_scan.io.load_newest_run()
    prms = r["parameters"]
    gr = GracePlot("plot")    

    re, r0 = param_scan.simulation_run.read_values(r, param_scan.parameter_sets.unroll(prms), ["rin_e", "r0"], "stdout")
    gr.plot(re, r0)
    theo_res = np.linspace(min(re), max(re), 50)
    models = {"pif": ana.PIF(), "lif": ana.LIF(), "qif": ana.QIF()}
    gwnmodels = {"pif": gwnana.PIF(), "lif": gwnana.LIF(), "qif": gwnana.QIF()}
    print [ana.r0(prms, model=models[prms["model"]], rin_e=re, vt=prms["vtb"]+100*prms["d"]) for re in theo_res]
    #gr.plot(theo_res, [ana.r0(prms, model=models[prms["model"]], rin_e=re) for re in theo_res])
    #gr.plot(theo_res, [gwnana.r0(prms, model=gwnmodels[prms["model"]], mu=prms["mu"]+prms["a_e"]*re, D=re*prms["a_e"]**2) for re in theo_res])

    gr.save()
    LatexParamValues().write("paramvalues.tex", prms)


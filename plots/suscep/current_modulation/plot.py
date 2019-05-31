import analytics.shot_noise_driven.lif_neuron as ana
import analytics.gaussian_white_noise_driven.lif_neuron as gwnana
import pylab as pl
import numpy as np
import param_scan.io
import param_scan.simulation_run
from param_scan.io import InDirectory
import os
from grace_plot import GracePlot
from latex_param_values import LatexParamValues

with InDirectory(os.path.abspath(__file__)):

    rfn, r = param_scan.io.load_newest_run()
    prms = r["parameters"]

    sx = np.loadtxt(param_scan.simulation_run.get_output_path(r, r["parameters"], "sx.spec")).T
    ss = np.loadtxt(param_scan.simulation_run.get_output_path(r, r["parameters"], "ss.spec")).T
    
    fs = np.linspace(1e-4, 10., 100)

    gr = GracePlot("plot", rows=2)
    gr.focus(row=0)
    simsus = (sx[1]+1j*sx[2])/ss[1];
    gr.plot(sx[0], np.abs(simsus))
    gr.plot(fs, prms["eps_v"] * np.abs(ana.suscep(prms, fs=fs)))
    gr.plot(fs, prms["eps_v"] * np.abs(ana.suscep_highf(prms, fs=fs)))
    gr.plot(fs, prms["eps_v"] * np.abs(gwnana.suscep(prms, mu=prms["mu"]*prms["a_e"]*prms["rin_e"], D=prms["rin_e"]*prms["a_e"]**2, fs=fs)))
    gr.focus(row=1)
    gr.plot(sx[0], np.angle(simsus))
    gr.plot(fs, np.angle(ana.suscep(prms, fs=fs)))
    gr.plot(fs, np.angle(ana.suscep_highf(prms, fs=fs)))
    gr.plot(fs, np.angle(gwnana.suscep(prms, mu=prms["mu"]*prms["a_e"]*prms["rin_e"], D=prms["rin_e"]*prms["a_e"]**2, fs=fs)))

    LatexParamValues().write("paramvalues.tex", prms)
    gr.save()

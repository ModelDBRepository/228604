import analytics.shot_noise_driven.lif_neuron as ana
import analytics.gaussian_white_noise_driven.lif_neuron as gwnana
import pylab as pl
import numpy as np
import param_scan.io
import param_scan.simulation_run
import param_scan.parameter_sets as ps
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
    #p["f_sig"] = prms["f_sig"] / tau_m
    p["r_sample"] = prms["r_sample"] / tau_m
    return p

tau_m = 0.02 # s
with InDirectory(os.path.abspath(__file__)):

    rfn, r = param_scan.io.load_newest_run()
    prms = r["parameters"]

    fs = np.logspace(-1, 2, 1000)
    
    fs_sim, resus, imsus = param_scan.simulation_run.read_values(r, ps.unroll(r["parameters"]), ["f_sig_c", "resus", "imsus"], "stdout", ignore_errors=True)
    fs_sim, resus, imsus = map(np.array, [fs_sim, resus, imsus])
    simsus = (resus+1j*imsus);

    gr = GracePlot("plot", rows=2)
    gr.focus(row=0)
    gr.plot(fs_sim/tau_m, np.abs(simsus/tau_m))
    gr.plot(fs/tau_m, np.abs(ana.suscep(prms, fs=fs)/tau_m))
    gr.plot(fs/tau_m, np.abs(ana.suscep_highf(prms, fs=fs)/tau_m))
    gr.plot(fs/tau_m, np.abs(gwnana.suscep(prms, mu=prms["mu"]*prms["a_e"]*prms["rin_e"], D=prms["rin_e"]*prms["a_e"]**2, fs=fs)/tau_m))
    gr.focus(row=1)
    gr.plot(fs_sim/tau_m, np.angle(simsus))
    gr.plot(fs/tau_m, np.angle(ana.suscep(prms, fs=fs)/tau_m))
    gr.plot(fs/tau_m, np.angle(ana.suscep_highf(prms, fs=fs)/tau_m))
    gr.plot(fs/tau_m, np.angle(gwnana.suscep(prms, mu=prms["mu"]*prms["a_e"]*prms["rin_e"], D=prms["rin_e"]*prms["a_e"]**2, fs=fs)/tau_m))

    LatexParamValues().write("paramvalues.tex", with_tau_m(tau_m * 1000, prms))
    gr.save()

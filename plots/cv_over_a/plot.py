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

    rfn, r = param_scan.io.load_newest_run(d="./sim")
    prms = r["parameters"]
    rfn_theo, r_theo = param_scan.io.load_newest_run_matching({}, d="./theo")
    rfn_da, r_da = param_scan.io.load_newest_run_matching({}, d="./diffapp")
    gr = GracePlot("plot")    

    a, rin_e, cv = param_scan.simulation_run.read_values(r, param_scan.parameter_sets.unroll(prms), ["a_e", "rin_e", "cv"], "stdout", ignore_errors=True)
    gr.plot(a, cv)
    
    a_theo, rin_e_theo, cv_theo = param_scan.simulation_run.read_values(r_theo, param_scan.parameter_sets.unroll(r_theo["parameters"]), ["a_e", "rin_e", "cv"], "stdout", ignore_errors=True)
    gr.plot(a_theo, cv_theo)
    
    a_da, rin_e_da, cv_da = param_scan.simulation_run.read_values(r_da, param_scan.parameter_sets.unroll(r_da["parameters"]), ["a_e", "rin_e", "cv"], "stdout", ignore_errors=True)
    gr.plot(a_da, cv_da)
    
    gr.plot(a_theo, [1 for a in a_theo])
    gr.plot(a_theo, [1./3**0.5 for a in a_theo])

    gr.save()
    LatexParamValues().write("paramvalues.tex", with_tau_m(tau_m * 1000, prms))


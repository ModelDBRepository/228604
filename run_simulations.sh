set -e

echo "Running simulations for probability density (EIF)"
(cd plots/probability_density/eif &&
../../../run_simulation sn_sig.sim)

echo "Running simulations for probability density (LIF)"
(cd plots/probability_density/lif &&
../../../run_simulation sn_spont.sim)

(echo "Generating theory for CV (QIF)"
cd plots/cv_over_a/qif/theo &&
    ../../../../run_simulation sn_fpt_theo.sim)

(echo "Generating diffusion approximation theory for CV (QIF)"
cd plots/cv_over_a/qif/diffapp &&
 ../../../../run_simulation qif_diff_app.sim)

(echo "Running simulations for CV (QIF)"
cd plots/cv_over_a/qif/sim &&
 ../../../../run_simulation sn_spont.sim)

echo "Running simulations for firing rate (EIF)"
(cd plots/r0_over_rin_e/eif &&
../../../run_simulation sn_sig.sim)

echo "Running simulations for firing rate (LIF)"
(cd plots/r0_over_rin_e/lif/lowrin &&
../../../../run_simulation sn_spont.sim)

echo "Running simulations for power spectrum (LIF)"
(cd plots/power_spectrum &&
../../run_simulation sn_spont.sim)

echo "Running simulations for susceptibility (LIF)"
(cd plots/suscep/current_modulation/cosstim &&
../../../../run_simulation sn_sig.sim)

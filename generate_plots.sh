set -e

python -m plots.probability_density.eif.plot
python -m plots.probability_density.lif.plot
python -m plots.r0_over_rin_e.lif.lowrin.plot
python -m plots.cv_over_a.qif.plot
python -m plots.power_spectrum.plot
echo "plotting r0 over rin for EIF... this may take some time"
python -m plots.r0_over_rin_e.eif.plot
python -m plots.suscep.current_modulation/cosstim/plot

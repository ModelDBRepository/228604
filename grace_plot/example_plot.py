import numpy as np
from grace_plot import GracePlot

gr=GracePlot("exampleplot", rows=1, cols=2)

x = np.linspace(0,10,100)

gr.focus(col=0) # erste spalte
gr.plot(x, x**3)

gr.focus(col=1) # zweite spalte
gr.plot(x, np.cos(x), legend="ein cosinus")

gr.save()

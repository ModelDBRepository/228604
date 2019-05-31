#! /usr/bin/env python
from numpy import *
from pygrace import grace

p=grace()
#p=GracePlot()
x=linspace(1,10,20)
y=sin(x)

p.plot(x,y)
p.multi(2,1)

p.focus(1,0)
p.plot(2*x,y)

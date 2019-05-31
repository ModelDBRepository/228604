
from pygrace.gracePlot import gracePlot
import os
import time
from subprocess import check_call

class GracePlot:
    def __init__(self, name, rows=1, cols=1, nogui=True):
        self.name = name
        self.numlines = 0
        if os.path.exists(name + ".par"):
            self.param_file = name + ".par"
        elif os.path.exists("grace.par"):
            self.param_file = "grace.par"
        else:
            self.param_file = None
        
        self.gr = gracePlot(param_file=None, nogui=nogui)#self.param_file)
        self.strings = []
        self.legends = {}
        for i in range(rows*cols):
            self.legends[i] = []
        if rows*cols > 1:
            self.gr.multi(rows, cols)
            #if self.param_file:
            #    self.gr.grace('getp "' + self.param_file + '"')

    def focus(self, row=0, col=0):
        self.gr.focus(row, col)

    def plot(self, x, y=None, dy=None, legend=None):
        self.gr.plot(x, Y=y, dy=dy)
        self.legends[self.gr.curr_graph.gID].append(legend)

    def line(self, x1, y1, x2, y2):
        self.numlines += 1
        self.command("with line %i" % (self.numlines - 1))
        self.command("    line on")
        self.command("    line loctype world")
        self.command("    line g0")
        self.command("    line %f, %f, %f, %f" % (x1, y1, x2, y2))
        self.command("    line linewidth 1.8")
        self.command("    line linestyle 1")
        self.command("    line color 1")
        self.command("    line arrow 3")
        self.command("    line arrow type 2")
        self.command("    line arrow length 1.000000")
        self.command("    line arrow layout 1.000000, 1.000000")
        self.command("line def")


    def text(self, text):
        self.strings.append(text)

    def command(self, cmd):
        self.gr.grace(cmd)

    def save(self):
        if os.path.exists(self.name + ".eps"):
            os.remove(self.name + ".eps")
        if self.param_file:
            self.gr.grace('getp "' + self.param_file + '"')

        # now that the parameters are loaded, print the strings
        i = 0
        for text in self.strings:
            self.command("with string %i" % i)
            i += 1
            self.command("string on")
            self.command('string def "' + str(text) + '"')

        #if len(self.legends) > 0:
        for g in self.gr.g:
            for i in range(len(self.legends[g.gID])):
                if self.legends[g.gID][i] is not None:
                    self.command( ('g%s.s%s legend "' % (g.gID, i)) + self.legends[g.gID][i] + '"' )
                    self.command('with g%s; legend on' % g.gID)
                    

        self.gr.grace('redraw')
        self.gr.grace('hardcopy device "EPS"')
        self.gr.grace('print to "' + self.name + '.eps"')
        self.gr.grace('print')
        self.gr.grace('saveall "' + self.name + '.agr"')
        self.gr.grace('exit')

        # check if eps is complete, otherwise wait 1 sec
        while not os.path.exists(self.name + ".eps"):
            time.sleep(.5)
        fh = open(self.name + ".eps", "r")
        epsstr = fh.read()
        fh.close()
        if len(epsstr) > 5 and not (epsstr[-6:] == "%%EOF\n"):
            time.sleep(1)
        check_call(['epstopdf', self.name + '.eps'])
            


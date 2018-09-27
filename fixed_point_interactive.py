import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import matplotlib.colors as colors
import matplotlib.cm as cmx
import warnings

def poly(coeffs):
  def calc(x):
    result = 0
    for i,c in enumerate(reversed(coeffs)):
      result += c*(x**i)
    return result
  return calc


def get_fp_iteration_display_points(Gfun, x0, maxits=6, atol=10**(-5), rtol=10**(-5)):
    # first iteration: xaxis (x0, 0) -> y=g (x0, Gfun(x0)) -> y=x (Gfun(x0), Gfun(x0)) -> 
    xs = [x0]
    ys = [0]

    guesses = [x0]

    for i in range(maxits):
        #y=g (yk, Gfun(xk))
        xk0 = xs[-1]
        yk0 = Gfun(xk0)

        #y=x (Gfun(xk), Gfun(xk))
        xk1 = Gfun(xk0)
        yk1 = Gfun(xk0)

        #y=0, x=xk+1 (Gfun(xk), 0)
        xk2 = Gfun(xk0)
        yk2 = 0

        xs += [xk0, xk1] #, xk2]
        ys += [yk0, yk1] #, yk2]

        guesses += [xk2]

        s = xs[-1] - Gfun(xs[-1])
        if (abs(s) < atol + rtol*s):
            return xs, ys, guesses, i+1, True # return number of iterations and wether or not there was convergence

    return xs, ys, guesses, i+1, False #  return number of iterations and wether or not there was convergence


def poly_to_Gfun(p):
    # Function to convert to root finding problem given f(x). 'f(x*) = 0' --> 'g(x*) = x*' 
    #Gfun = lambda Ffun: lambda x: p(x) + x
    Gfun = poly(p.c)
    return Gfun

Path = mpath.Path

fig, ax = plt.subplots()

pathdata = [
    (Path.MOVETO, (-2.37, 0.7)),
    (Path.MOVETO, (0.0, 0.35)),
    (Path.MOVETO, (2.4, -2.2))
    ]

domain = ((-3,3), (-3,3))

codes, verts = zip(*pathdata)
path = mpath.Path(verts, codes)
patch = mpatches.PathPatch(path, facecolor=None, edgecolor=None, alpha=0.5, visible=False)
ax.add_patch(patch)


class InteractivePolynomial(object):

    epsilon = 25  # max pixel distance to count as a vertex hit

    def plot_iterations(self, axsi, xs, ys, guesses, its, convergence):

        axsi.plot(xs, ys, color='m', marker=None, linestyle='dashed', linewidth=0.5)


        ylim = axsi.get_ylim()
        #axsi.text(0, ylim[0],text, color=text_color, fontsize=8, ha='center', va='bottom', alpha=0.8)

        cm = plt.get_cmap("hot") 
        cNorm = colors.Normalize(vmin=0, vmax=len(guesses))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        for idx,g in enumerate(guesses):
            colorVal = scalarMap.to_rgba(idx)
            alpha=1 #1/(idx+1)
            markersize= 20/(idx+1)
            axsi.plot(g,0, marker='x', color=colorVal, alpha=alpha,markersize=markersize)

        cNorm = colors.Normalize(vmin=0, vmax=len(xs))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        for idx in range(len(xs)):
            colorVal = scalarMap.to_rgba(idx)
            alpha=0.75 #1/(idx+1)
            markersize= 10/(idx+1)
            axsi.plot(xs[idx],ys[idx], marker='.', color=colorVal, alpha=alpha, markersize=markersize)


    def __init__(self, pathpatch, domain):

        self.ax = pathpatch.axes
        canvas = self.ax.figure.canvas
        self.pathpatch = pathpatch
        self.pathpatch.set_animated(True)

        x,y = self.get_verts()

        self.domain = domain
        self.num_points = 300 # number of points to plot on graph each function


        # make pretty
        #self.ax.set_title(title, fontsize=8)
        self.ax.grid(True, which='both')
        self.ax.axhline(y=0, color='k')
        self.ax.axvline(x=0, color='k')
        self.ax.set_aspect('equal') # equal
        self.ax.set_ylim([self.domain[0][0],self.domain[1][1]]) # set y lim

        # draw polynomial
        xs, ys, self.f = self.poly_points()
        #self.function_latex = r"$cos(x)$"
        self.polynomial, = self.ax.plot(xs, ys, animated=True) #, label=self.function_latex)



        # plot y=x line
        target_line_xs = np.linspace(self.domain[0][0],self.domain[1][1],self.num_points)
        target_line_ys = target_line_xs # compute data points for plot y = x
        # create plot for this function
        self.traget_line = self.ax.plot( target_line_xs, target_line_ys, 'r', animated=False) #, label=latex)     # y = g(x) in red


        # polynomial plot
        self.line, = self.ax.plot(x, y, marker='.', markerfacecolor='r', animated=True)
        self.line.set_linestyle('None')
        self.line.set_markerfacecolor('None')
        self.line.set_markeredgecolor('b')
        self.line.set_markersize(20)
        self._ind = None  # the active vert


        # draw x0
        # initial guess
        self.x0 = 2
        #self.initial_guess, = ax.axvline(x=self.x0,animated=True) #ax.plot(self.x0, 0, marker='X', animated=True)
        self.initial_guess, = self.ax.plot(self.x0, 0, marker='X', animated=True)
         #, marker='X', animated=True)
        self.initial_guess.set_markersize(10)
        self.initial_updated = True
        self.released = True


        # fixed point plot
        Gfun = poly_to_Gfun(self.f)
        xs, ys, guesses, its, convergence = get_fp_iteration_display_points(Gfun, self.x0)
        self.fixed_point_plot, = self.ax.plot(xs, ys, color='black', marker=None, linestyle='dashed', linewidth=1)


        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas

    def poly_points(self):

        x,y = self.get_verts()

        # calculate polynomial
        with np.testing.suppress_warnings() as sup:
            sup.filter(np.RankWarning, "Polyfit may be poorly conditioned")
            z = np.polyfit(x, y, len(x))
            f = np.poly1d(z)

        # calculate new x's and y's
        #m = 300 # number of points to plot on graph
        x_plt_points = np.linspace(self.domain[0][0], self.domain[0][1], self.num_points)
        y_plt_points = f(x_plt_points)

        return x_plt_points, y_plt_points, f

    def update_polynomial(self):
        x_plt_points, y_plt_points, self.f = self.poly_points()
        #self.function_latex = r"$cos(x)$"
        self.polynomial.set_data(x_plt_points, y_plt_points)

    def get_verts(self):
        x, y = zip(*self.pathpatch.get_path().vertices)
        return list(x), list(y)

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.polynomial)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.initial_guess)
        self.ax.draw_artist(self.fixed_point_plot)
        self.canvas.blit(self.ax.bbox)


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # moving polynomial strechers
        # display coords
        xy = np.asarray(self.pathpatch.get_path().vertices)
        xyt = self.pathpatch.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
        ind = d.argmin()

        if d[ind] >= self.epsilon:
            ind = None

        return ind

    def check_for_x0_position(self, event):
        'get the index of the vertex under point if within epsilon tolerance of x0'

        # moving x0
        # display coords
        xyt = self.pathpatch.get_transform().transform(np.asarray([self.x0, 0]))
        xt, yt = xyt[0], xyt[1]

        d = np.sqrt((xyt[0] - event.x)**2 + (xyt[1] - event.y)**2)

  
        if d >= self.epsilon+10:
            self.initial_updated = False
        else:
            self.initial_updated = True


    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)
        
        self.check_for_x0_position(event)
        self.released = False
 

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if event.button != 1:
            return
        self._ind = None
        self.released = True


    def motion_notify_callback(self, event):
        'on mouse movement'
        #if self._ind is None:
        #    return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        # self_ind refers to polynomial strech markers
        # update polynomial stretch markers
        if self._ind is not None:
            vertices = self.pathpatch.get_path().vertices
            vertices[self._ind] = x, y

            self.canvas.restore_region(self.background)

            self.line.set_data(zip(*vertices))
            self.ax.draw_artist(self.line)

            self.ax.draw_artist(self.initial_guess)

            self.update_polynomial()
            self.ax.draw_artist(self.polynomial)

            self.canvas.blit(self.ax.bbox)

        # update x0
        if self.initial_updated: # and self.released:
            self.canvas.restore_region(self.background)

            self.ax.draw_artist(self.line)
            self.ax.draw_artist(self.polynomial)

            x, y = event.inaxes.transData.inverted().transform((event.x, event.y))
            self.x0 = x
   
            self.initial_guess.set_data(self.x0, 0)
            self.ax.draw_artist(self.initial_guess)

            if self.released:
                self.initial_updated = False

            self.canvas.blit(self.ax.bbox)

        # call fixed point iteration
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.polynomial)
        self.ax.draw_artist(self.initial_guess)
        Gfun = poly_to_Gfun(self.f)
        xs, ys, guesses, its, convergence = get_fp_iteration_display_points(Gfun, self.x0)
        self.fixed_point_plot.set_data(xs, ys)
        self.ax.draw_artist(self.fixed_point_plot)



interactor = InteractivePolynomial(patch, domain)
ax.set_title('drag points to update fixed point iteration')
ax.set_xlim(domain[0][0], domain[0][1])
ax.set_ylim(domain[1][0], domain[1][1])

plt.show()

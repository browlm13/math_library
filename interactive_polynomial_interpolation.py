import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


Path = mpath.Path

fig, ax = plt.subplots()

pathdata = [
    (Path.MOVETO, (-2.37, 0.7)),
    (Path.MOVETO, (0.0, 0.35)),
    (Path.MOVETO, (2.4, -2.2))
    ]

domain = ((-3,3), (-3,3))

codes, verts = zip(*pathdata)
path = mpath.Path(verts, codes) #, _interpolation_steps=2, closed=False)
patch = mpatches.PathPatch(path, facecolor=None, edgecolor=None, alpha=0.5, visible=False)
ax.add_patch(patch)

class InteractivePolynomial(object):

    epsilon = 15  # max pixel distance to count as a vertex hit

    def __init__(self, pathpatch, domain):

        self.ax = pathpatch.axes
        canvas = self.ax.figure.canvas
        self.pathpatch = pathpatch
        self.pathpatch.set_animated(True)

        x,y = self.get_verts()

        self.domain = domain

        # make pretty
        #self.ax.set_title(title, fontsize=8)
        self.ax.grid(True, which='both')
        self.ax.axhline(y=0, color='k')
        self.ax.axvline(x=0, color='k')
        self.ax.set_aspect('equal') # equal
        self.ax.set_ylim([self.domain[0][0],self.domain[1][1]]) # set y lim
        self.ax.legend(loc=9, bbox_to_anchor=(0.5, -0.2))

        xs, ys, self.f = self.poly_points()
        self.polynomial, = ax.plot(xs, ys, animated=True)

        self.line, = ax.plot(x, y, marker='.', markerfacecolor='r', animated=True)
        self.line.set_linestyle('None')
        self.line.set_markerfacecolor('None')
        self.line.set_markeredgecolor('b')
        self.line.set_markersize(20)

        self._ind = None  # the active vert

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
        m = 300 # number of points to plot on graph
        x_plt_points = np.linspace(self.domain[0][0], self.domain[0][1], m)
        y_plt_points = f(x_plt_points)

        return x_plt_points, y_plt_points, f

    def update_polynomial(self):
        x_plt_points, y_plt_points, self.f = self.poly_points()
        self.polynomial.set_data(x_plt_points, y_plt_points)

    def get_verts(self):
        x, y = zip(*self.pathpatch.get_path().vertices)
        return list(x), list(y)

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.polynomial)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # display coords
        xy = np.asarray(self.pathpatch.get_path().vertices)
        xyt = self.pathpatch.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
        ind = d.argmin()

        if d[ind] >= self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if event.button != 1:
            return
        self._ind = None

    def motion_notify_callback(self, event):
        'on mouse movement'
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        vertices = self.pathpatch.get_path().vertices
        vertices[self._ind] = x, y

        self.canvas.restore_region(self.background)

        self.line.set_data(zip(*vertices))
        self.ax.draw_artist(self.line)

        self.update_polynomial()
        self.ax.draw_artist(self.polynomial)

        self.canvas.blit(self.ax.bbox)


interactor = InteractivePolynomial(patch, domain)
ax.set_title('drag points to update polynomial')
ax.set_xlim(domain[0][0], domain[0][1])
ax.set_ylim(domain[1][0], domain[1][1])

plt.show()

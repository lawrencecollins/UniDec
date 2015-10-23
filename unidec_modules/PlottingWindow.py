# base class for all plotting windows used in UniDec
# contains basic setup functionality

import wx
from matplotlib import interactive
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import numpy as np

from unidec_modules.isolated_packages.ZoomSpan import ZoomSpan
from unidec_modules.isolated_packages.ZoomBox import ZoomBox
from unidec_modules.isolated_packages import FileDialogs

interactive(True)


class PlottingWindow(wx.Window):
    """
    Class for wx window with embedded matplotlib plots
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize plot window parameters.

        Optional keywords:
        figsize: size of figure in inches
        integrate: 0 or 1 of whether the plot should send the integrate pubsub when a right click is activated.
        smash: 0 or 1 of whether the plot should send the integrate pubsub when a right click is activated.

        :param args: Arguments
        :param kwargs: Keywords
        :return:
        """
        if "figsize" in kwargs:
            self.figure = Figure(figsize=kwargs["figsize"])  # , dpi=200)
            del kwargs["figsize"]
            self._axes = [0.2, 0.2, 0.7, 0.7]
        else:
            self.figure = Figure()
            self._axes = [0.1, 0.1, 0.8, 0.8]

        if "integrate" in kwargs:
            self.int = kwargs["integrate"]
            del kwargs["integrate"]
        else:
            self.int = 0
        if "smash" in kwargs:
            self.smash = kwargs["smash"]
            del kwargs["smash"]
        else:
            self.smash = 0
        wx.Window.__init__(self, *args, **kwargs)

        self.subplot1 = None
        self.zoom = None
        self.subplot1 = None
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.Bind(wx.EVT_SIZE, self.size_handler)
        self.resize = 1
        self.flag = False
        self.kda = False
        self.kdnorm = 1.
        self.normalticks = True
        self.nativez = []
        self.text = []
        self.lines = []
        self.cbar = None
        self.datalims = None
        self.cmap = None
        self.set_color()
        self.xlabel = ""
        self.ylabel = ""
        self.zoomtype = "box"
        self.tickcolor = "black"
        self.canvas.mpl_connect('button_release_event', self.on_release)

    def on_release(self, event):
        """
        Function triggered on button release event from plot.
        Currently wired to trigger on_save_figure_dialog on middle button.
        :param event: wx.Event
        :return: None
        """
        if event.button == 2:
            self.on_save_fig_dialog(event)

    def on_save_fig_dialog(self, evt):
        """
        Open a save figure dialog for specified plot.
        :param evt: wx.Event (unused)
        :return: None
        """
        path = FileDialogs.save_file_dialog()
        if path is not None:
            self.save_figure(path)

    def on_save_fig(self, evt, path, **kwargs):
        """
        Save figure to path.
        :param evt: wx.Event (unused)
        :param path: Path to save figure to
        :param kwargs: keywords passed to save_figure
        :return: None
        """
        if path is not None:
            self.save_figure(path, **kwargs)

    def save_figure(self, path, **kwargs):
        """
        Saves Figure to path.
        :param path: Path to save figure at.
        :param kwargs: Keywords passed to matplotlib.figure.savefig
        :return: None
        """
        if "transparent" in kwargs:
            t = kwargs["transparent"]
        else:
            t = True
        self.figure.savefig(path, transparent=t)
        print "Saved Figure: ", path

    def kda_test(self, xvals):
        """
        Test whether the axis should be normalized to convert mass units from Da to kDa.
        Will use kDa if: xvals[int(len(xvals) / 2)] > 100000 or xvals[len(xvals) - 1] > 1000000

        If kDa is used, self.kda=True and self.kdnorm=1000. Otherwise, self.kda=False and self.kdnorm=1.
        :param xvals: mass axis
        :return: None
        """
        try:
            if xvals[int(len(xvals) / 2)] > 100000 or xvals[len(xvals) - 1] > 1000000:
                self.kdnorm = 1000.
                self.xlabel = "Mass (kDa)"
                self.kda = True
            else:
                self.xlabel = "Mass (Da)"
                self.kda = False
                self.kdnorm = 1.
        except TypeError or ValueError:
            self.xlabel = "Mass (Da)"
            self.kdnorm = 1.
            self.kda = False

    def plotadddot(self, x, y, colval, markval):
        """
        Adds a scatter plot to the figure. May be one or more.
        :param x: x values
        :param y: y values
        :param colval: Color
        :param markval: Marker
        :return: None
        """
        self.subplot1.plot(np.array(x) / self.kdnorm, y, color=colval, marker=markval, linestyle='None', clip_on=False)

    def repaint(self):
        """
        Redraw and refresh the plot.
        :return: None
        """
        self.canvas.draw()

    def clear_plot(self, *args):
        """
        Clear the plot and rest some of the parameters.
        :param args: Arguments
        :return:
        """
        self.figure.clear()
        self.flag = False
        self.nativez = []
        self.text = []
        self.lines = []
        self.kda = False
        self.kdnorm = 1.
        if "nopaint" not in args:
            self.repaint()

    def set_nticks(self, bins):
        """
        Set the number of ticks in the x-axis.
        :param bins: Number of ticks in the x-axis
        :return: None
        """
        if self.normalticks:
            self.subplot1.tick_params(axis="x", labelsize=12)
            self.subplot1.tick_params(axis="y", labelsize=12)
            self.subplot1.xaxis.set_major_locator(MaxNLocator(nbins=bins))
        self.repaint()

    def add_legend(self, location=1):
        """
        Adds a legend to the plot.
        :param location: Integer code for location
        :return: None
        """
        handles, labels = self.subplot1.get_legend_handles_labels()
        if location == 1:
            self.subplot1.legend(handles, labels, loc=location, bbox_to_anchor=(1, 1))
        else:
            self.subplot1.legend(handles, labels, loc=location)
        self.repaint()

    def set_color(self, rgbtuple=None):
        """
        Sets background color
        :param rgbtuple: background color
        :return:
        """
        # Set figure and canvas colours to be the same
        if not rgbtuple:
            rgbtuple = [255., 255., 255.]
        col = [c / 255.0 for c in rgbtuple]
        self.figure.set_facecolor(col)
        self.figure.set_edgecolor(col)
        self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))

    def set_tickcolor(self):
        """
        Sets tick colors based on the colormap set at self.cmap
        :return: None
        """
        # TODO: Better determination of cmaps and Viridis

        output = cm.ScalarMappable(norm=None, cmap=self.cmap).to_rgba(0)
        if sum(output[:2]) > 0.9:
            self.tickcolor = "black"
        else:
            self.tickcolor = "white"
        '''
        if self.cmap[-1] == "r":
            self.tickcolor = "black"
        else:
            self.tickcolor = "white"
        '''

    def size_handler(self, *args, **kwargs):
        """
        Resizes the plots
        :param args:
        :param kwargs:
        :return: None
        """
        if self.resize == 1:
            self.canvas.SetSize(self.GetSize())

    def setup_zoom(self, plots, zoom, data_lims=None):
        """
        Set up zoom on axes.
        :param plots: Axes objects to setup
        :param zoom: Type of zoom ('span' or 'box')
        :param data_lims: Optional manual description of the data limits (where to go when fully zoomed out)
        :return: None
        """
        # setup for zoom box
        if zoom == 'span':
            self.zoom = ZoomSpan(
                plots,
                None,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))
        if zoom == 'box':
            self.zoom = ZoomBox(
                plots,
                None,
                drawtype='box',
                useblit=True,
                button=1,  # so zoombox is left button
                onmove_callback=None,
                spancoords='data',
                rectprops=dict(alpha=0.2, facecolor='yellow'),
                data_lims=data_lims,
                integrate=self.int, smash=self.smash)

"""

power_fit_gui.py

Purpose: Copied from fit_example_power_norm.py and edited to use in
         power_normalization_v2.py. 

Revisions:
        fit_example_power_norm.py
    2018 May 09 - gsteranka - Original version
        power_fit_gui.py
    2018 May 10 - gsteranka - Copied to use in GitHub code
"""

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
    NavigationToolbar2TkAgg)
from matplotlib.figure import Figure
import numpy as np

import tkinter
from tkinter import Tk, IntVar, StringVar
from tkinter.ttk import Frame, Combobox, Label, Button, Entry, LabelFrame

class PowerFitGui(Frame):
    """
    Class for a GUI to make a spline fit to power. Includes buttons to adjust
    parameters of the fit, including fit order, freespace locations, and knot
    locations. This is not intended to be called by a user directly, but rather
    to be used inside of the Normalization class when the USE_GUI keyword is
    set to True (which is the default).

    Args:
        parent (tkinter.Tk): Instance of tkinter.Tk class. This is the parent
            that the GUI is based on
        norm_inst: Instance of the Normalization class
        spm_vals_down (np.ndarray): Array of SPM values for the power values to
            fit. Downsampling is done using scipy.signal.resample_poly in
            the Normalization class. Downsampled to the spacing "dt_down" in
            the Normalization class (1/2 second by default).
        p_obs_down (np.ndarray): Array of downsampled power values to fit.
            Gotten from downsampling frequency-corrected complex samples
        spm_fit (np.ndarray): SPM values to evaluate the fit for

    Attributes:
        ax: Axes instance corresponding to the matplotlib.pyplot figure on
            the GUI
        canvas: Canvas on which the matplotlib.pyplot and the widgets are put
        cvar: Text variable in the combo box for selecting fit order
        f: Matplotlib.pyplot figure on the GUI
        fit_deg (int): Degree of the spline fit. Starts off at 1
        ind (np.ndarray): Index numbers to the "x" attribute of the regions
            being fit (specified by the "xlim" attribute)
        knots_entry: Entry box for specifying the knots used in the fit.
            Changing the entry to this box changes the "knots_spm" attribute
            (see below)
        knots_spm (list): Knots of the fit. Adjustable in the GUI in the box
            labeled "Knots". Each knot must be separated by a comma. For
            example, if you wanted to use the SPM knots 31000, 36500, and
            39000, then you would just enter "31000, 36500, 39000" (without
            the quotation marks). I usually use one knot per region being fit,
            which in this case is the regions in the "xlim" attribute (below)
        lvar: Variable for the label next to the combo box for selecting fit
            order
        norm_inst: The input instance of the Normalization class
        not_ind (np.ndarray): Index numbers of the "x" attribute of the regions
            not being fit (specified by the "xlim" attribute)
        parent (tkinter.Tk): Input parent
        toolbar: Matplotlib.pyplot toolbar that appears below the
            matplotlib.pyplot plot
        x (np.ndarray): The input spm_vals_down values
        xfit (np.ndarray): The input spm_fit values
        xlim (list): Set of regions to make the fit over. Starts off as the
            full range of SPM. Adjustable in the GUI, in the box labeled
            "Fit ranges". To specify, separate minimum and maximum values by a
            comma, and separate ranges by a semicolon. For example, if you want
            to make a fit using the SPM ranges [30500, 33000], [36000, 37000],
            and [38000, 40000], then you would enter into the box:
            "30500, 33000 ; 36000, 37000 ; 38000, 40000" (without the quotation
            marks)
        xlim_entry: Entry box for specifying regions to make the fit over.
            Entering into this box changes the "xlim" attribute (see above)
        y (np.ndarray): The input p_obs_down values
        yfit (np.ndarray): Resulting fit from the GUI. This is the attribute
            used to access the fit after the program is closed
        
    """

    def __init__(self, parent, norm_inst, spm_vals_down, p_obs_down, spm_fit):

        Frame.__init__(self, parent)

        self.norm_inst = norm_inst
        self.x = spm_vals_down
        self.y = p_obs_down
        self.xfit = spm_fit

        self.fit_deg = norm_inst._k
        self.xlim = norm_inst._freespace_spm
        self.knots_spm = norm_inst._knots_spm
        self.yfit = self._get_fit()

        self.ind = np.arange(len(self.x))
        self.not_ind = None

        self.parent = parent
        self.initUI()


    def initUI(self):
        """
        Initialize the user interface. Called by __init__()
        """

        self.parent.title('fit_example')
        self.pack(fill=tkinter.BOTH, expand=1)

        self.f = Figure(figsize=(5, 4))
        self.ax = self.f.add_subplot(111)

        self.plot_data()

        self.cvar = IntVar()

        fit_order_frame = LabelFrame(self, text='Fit Order')

        # Box to select fit order
        combo = Combobox(fit_order_frame, textvariable=self.cvar)
        combo.bind('<<ComboboxSelected>>', self.adjust_deg)
        combo['values'] = [1, 2, 3, 4, 5]
        combo.current(0)
        combo.grid(row=0, column=0)

        # Variable outside of box that says current fit order
        self.lvar = IntVar()
        self.lvar.set('1')
        lbl = Label(fit_order_frame, textvariable=self.lvar)
        lbl.grid(row=0, column=1)

        # Frame for fit order entry
        fit_order_frame.pack(side=tkinter.LEFT, padx=15)

        # Frame for freespace and knot entry
        fit_info_frame = LabelFrame(self, text='Fit ranges and knots (SPM)')

        # Label for free space region entry
        xlim_entry_lbl = Label(fit_info_frame, text='Fit ranges:')
        xlim_entry_lbl.grid(row=0, column=0)

        # Entry box for regions to fit
        evar1 = IntVar()
        self.xlim_entry = Entry(fit_info_frame, textvariable=evar1, width=50)
        self.xlim_entry.grid(row=0, column=1)

        # Label for knot entry
        knots_entry_lbl = Label(fit_info_frame, text='Knots:')
        knots_entry_lbl.grid(row=1, column=0)

        # Entry box for knots
        evar2 = IntVar()
        self.knots_entry = Entry(fit_info_frame, textvariable=evar2, width=50)
        self.knots_entry.grid(row=1, column=1)

        fit_info_frame.pack(pady=15)

        # Button to use the x range entered in the box
        set_xrange_btn = Button(self, text='Set fit range and knots',
            command=self.adjust_range_and_knots)
        set_xrange_btn.pack()

        # Button to revert to original x range
        revert_xrange_btn = Button(self, text='Reset fit range and knots',
            command=self.revert_range)
        revert_xrange_btn.pack()

        # Button to press when you're content with the fit and want to continue
        ok_btn = Button(self, text='OK', command=self.quit)
        ok_btn.pack(side=tkinter.LEFT, padx=120)


    def _get_fit(self):
        """
        Get the spline fit to power. Called by __init__(), adjust_deg(),
        adjust_range_and_knots(), and revert_range()
        """

        spm_fit, spline_fit = self.norm_inst.get_spline_fit(spm_fit=self.xfit,
            k=self.fit_deg, knots_spm=self.knots_spm, freespace_spm=self.xlim,
            USE_GUI=False)
        return spline_fit

    def plot_data(self):
        """
        Make the original matplotlib.pyplot plot on the canvas. Called by
        initUI()
        """

        self.ax.plot(self.x, self.y, color='b')
        self.ax.plot(self.xfit, self.yfit, color='r')
        self.ax.set_title('Power')
        self.ax.set_xlabel('SPM')
        self.canvas = FigureCanvasTkAgg(self.f, master=self.parent)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH,
            expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH,
            expand=1)


    def update_plot(self):
        """
        Update the plot if an event happened (namely, you pressed a button or
        entered something new into a box). Called by adjust_deg(),
        adjust_range_and_knots(), and revert_range()
        """

        ind = self.ind
        not_ind = self.not_ind

        self.ax.cla()
        self.ax.plot(self.x[ind], self.y[ind], color='b')
        self.ax.plot(self.xfit, self.yfit, color='r')
        if not_ind is not None:
            self.ax.plot(self.x[not_ind], self.y[not_ind], color='b', alpha=0.1)
        self.ax.set_title('Power')
        self.ax.set_xlabel('SPM')
        self.canvas.show()


    def adjust_deg(self, e):
        """
        Adjust the degree of the fit if you adjusted the "Fit order" combobox.
        Bound to the "combo" variable in initUI()
        """

        new_deg = e.widget
        self.lvar.set(new_deg.get())
        self.fit_deg = int(new_deg.get())

        self.yfit = self._get_fit()
        self.update_plot()


    def adjust_range_and_knots(self):
        """
        Adjust ranges to fit over and the knots used if you changed what's
        entered in the text boxes and press the "Set fit range and knots"
        button. Bound to the "set_xrange_btn" variable in initUI()
        """

        xlim_entry = self.xlim_entry.get()
        xlim = []
        split1 = xlim_entry.split(';')
        for i in range(len(split1)):
            split2 = (split1[i]).split(',')
            _xlim = [float(split2[0]), float(split2[1])]
            xlim.append(_xlim)
        self.xlim = xlim

        ind = []
        for i in range(len(xlim)):
            ind.append(np.argwhere((self.x >= xlim[i][0]) &
                (self.x <= xlim[i][1])))
        ind = np.reshape(np.concatenate(ind), -1)
        self.ind = ind

        not_ind = []
        not_ind.append(np.argwhere(self.x < xlim[0][0]))
        for i in range(len(xlim)-1):
            not_ind.append(np.argwhere((self.x > xlim[i][1]) &
                (self.x < xlim[i+1][0])))
        not_ind = np.reshape(np.concatenate(not_ind), -1)
        self.not_ind = not_ind

        knots_entry = self.knots_entry.get()
        knots_entry_split = knots_entry.split(',')
        self.knots_spm = [float(knot) for knot in knots_entry_split]

        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


    def revert_range(self):
        """
        Revert to the original range and knots to fit over if you press the
        "Reset fit range and knots" button. Bound to the
        "revert_xrange_btn" in initUI()
        """

        self.xlim = [[min(self.x), max(self.x)]]
        self.ind = np.arange(len(self.x))
        self.not_ind = None
        self.knots_spm = [min(self.x) + 2.0, max(self.x) - 2.0]
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


if __name__ == '__main__':
    pass

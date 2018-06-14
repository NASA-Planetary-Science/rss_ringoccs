"""

f_resid_fit_gui.py

Purpose: Copied from fit_example_resid_fit.py and edited to use residual
         frequency. Edited to use inside of freq_offset_fit.py

Revisions:
      fit_example.py
  2018 Mar 28 - gsteranka - Original version
      fit_example_v2.py
   2018 May 08 - gsteranka - Edited to be able to specify multiple x ranges
      f_resid_fit_gui.py
   2018 May 08 - gsteranka - copied from prev version
   2018 May 30 - gsteranka - Make decent initial guess at fit
"""

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
    NavigationToolbar2TkAgg)
from matplotlib.figure import Figure
import numpy as np
import pdb

import tkinter
from tkinter import Tk, IntVar, StringVar
from tkinter.ttk import Frame, Combobox, Label, Button, Entry, LabelFrame


class FResidFitGui(Frame):
    """
    Class for a GUI to make a polynomial fit to residual frequency. Includes
    buttons to adjust parameters of the fit, including fit order and freespace
    locations. This is not intended to be called by a user directly, but rather
    to be used inside of the FreqOffsetFit class when the USE_GUI keyword is
    set to True (which is the default).

    Args:
        parent (tkinter.Tk): Instance of tkinter.Tk class. This is the parent
            that the GUI is based on
        fit_inst: Instance of the FreqOffsetFit class
        x (np.ndarray): Array of SPM values for residual frequency. Inside of
            FreqOffsetFit, this is the "f_spm" variable inside the
            set_f_sky_resid_fit() method
        y (np.ndarray): Array of residual frequency values. Inside of
            FreqOffsetFit, this is the "f_sky_resid" variable in the
            set_f_sky_resid_fit() method

    Attributes:
        ax: Axes instance corresponding to the matplotlib.pyplot figure on
            the GUI
        canvas: Canvas on which the matplotlib.pyplot and the widgets are put
        cvar: Text variable in the combo box for selecting fit order
        f: Matplotlib.pyplot figure on the GUI
        fit_deg (int): Degree of the polynomial fit. Starts off at 3
        fit_inst: Instance of the FreqOffsetFit class
        ind (np.ndarray): Index numbers to the "x" attribute of the regions
            being fit (specified by the "xlim" attribute)
        knots_entry: Entry box for specifying the knots used in the fit
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
        not_ind (np.ndarray): Index numbers of the "x" attribute of the regions
            not being fit (specified by the "xlim" attribute)
        parent (tkinter.Tk): Input parent
        toolbar: Matplotlib.pyplot toolbar that appears below the
            matplotlib.pyplot plot
        x (np.ndarray): The input x values (SPM)
        xlim (list): Set of regions to make the fit over. Starts off as the
            full range of SPM. Adjustable in the GUI, in the box labeled
            "Fit range (SPM)". To specify, separate minimum and maximum values
            by a comma, and separate ranges by a semicolon. For example, if you
            want to make a fit using the SPM ranges [30500, 33000],
            [36000, 37000], and [38000, 40000], then you would enter into the
            box: "30500, 33000 ; 36000, 37000 ; 38000, 40000" (without the
            quotation marks)
        xlim_entry: Entry box for specifying regions to make the fit over.
            Entering into this box changes the "xlim" attribute (see above)
        y (np.ndarray): The input y values (residual frequency in Hz)
        yfit (np.ndarray): Resulting fit from the GUI. This is the attribute
            used to access the fit after the program is closed
    """

    def __init__(self, parent, fit_inst, x, y):

        Frame.__init__(self, parent)

        self.parent = parent
        self.fit_inst = fit_inst
        self.x = x
        self.y = y
        self.initUI()


    def initUI(self):
        """
        Initialize the user interface. Called by __init__()
        """

        self.parent.title('FResidFitGui')
        self.pack(fill=tkinter.BOTH, expand=1)

        self.fit_deg = 3
        self.xlim = self.fit_inst._spm_include
        dummy, self.yfit = self.fit_inst.get_f_sky_resid_fit()

        self.ind = np.arange(len(self.x))
        self.not_ind = None

        self.f = Figure(figsize=(5, 4))
        self.ax = self.f.add_subplot(111)

        self.plot_data()

        self.cvar = IntVar()

        # Frame for fit order entry
        fit_order_frame = LabelFrame(self, text='Fit order')

        # Box to select fit order
        combo = Combobox(fit_order_frame, textvariable=self.cvar)
        combo.bind('<<ComboboxSelected>>', self.adjust_deg)
        combo['values'] = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        combo.current(2)
        combo.grid(row=0, column=0)

        # Variable outside of box that says current fit order
        self.lvar = IntVar()
        self.lvar.set('3')
        lbl = Label(fit_order_frame, textvariable=self.lvar)
        lbl.grid(row=0, column=1)

        fit_order_frame.pack(side=tkinter.LEFT, padx=15)

        # Frame for fit range entry
        fit_range_frame = LabelFrame(self, text='Fit range (SPM)')

        # Box to enter x regions to fit
        evar1 = IntVar()
        self.xlim_entry = Entry(fit_range_frame, textvariable=evar1, width=50,
            justify='center')
        self.xlim_entry.grid(row=0, column=0)

        # Description of x range entry
        xrange_btn_lbl = Label(fit_range_frame, text='min1, max1 ; min2, max2...')
        xrange_btn_lbl.grid(row=1, column=0)

        fit_range_frame.pack()

        # Button to use the x range entered in the box
        set_xrange_btn = Button(self, text='Set fit range (SPM)',
            command=self.adjust_range)
        set_xrange_btn.pack()

        # Button to revert to original x range
        revert_xrange_btn = Button(self, text='Reset fit range',
            command=self.revert_range)
        revert_xrange_btn.pack()

        # Quit application if fit is okay
        quit_btn = Button(self, text='OK', command=self.quit)
        quit_btn.pack()


    def _get_fit(self):
        """
        Get the polynomial fit to residual frequency. Called by __init__(),
        adjust_deg(), adjust_range(), and revert_range()
        """

        self.fit_inst.set_f_sky_resid_fit(k=self.fit_deg,
            spm_include=self.xlim, USE_GUI=False)
        _x, f_sky_resid_fit = self.fit_inst.get_f_sky_resid_fit()
        return f_sky_resid_fit


    def plot_data(self):
        """
        Make the original matplotlib.pyplot plot on the canvas. Called by
        initUI()
        """

        self.ax.plot(self.x, self.y, color='g')
        self.ax.plot(self.x, self.yfit, color='r')
        self.ax.set_xlabel('SPM')
        self.ax.set_ylabel('Hz')
        self.ax.set_title('Residual Frequency')
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
        adjust_range(), and revert_range()
        """

        ind = self.ind
        not_ind = self.not_ind

        self.ax.cla()
        self.ax.plot(self.x[ind], self.y[ind], color='g')
        self.ax.plot(self.x, self.yfit, color='r')
        self.ax.set_ylim([min(self.y[ind]) - 1.0, max(self.y[ind]) + 1.0])

        y_arrow = np.mean(self.ax.get_ylim())
        for i in range(len(self.xlim)):
            self.ax.axvline(self.xlim[i][0], color='b')
            self.ax.axvline(self.xlim[i][1], color='b')
            _x_arrow = float(self.xlim[i][0])
            _dx_arrow = float(np.diff(self.xlim[i]))
            self.ax.arrow(_x_arrow, y_arrow, _dx_arrow, 0,
                color='b', head_width=0.1, head_length=50,
                length_includes_head=True)

        if not_ind is not None:
            self.ax.plot(self.x[not_ind], self.y[not_ind], color='g', alpha=0.1)

        self.ax.set_xlabel('SPM')
        self.ax.set_ylabel('Hz')
        self.ax.set_title('Residual Frequency')
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


    def adjust_range(self):
        """
        Adjust ranges to fit over if you changed what's entered in the text
        box and then press the "Set fit range (SPM)" button. Bound to the
        "set_xrange_btn" variable in initUI()
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

        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


    def revert_range(self):
        """
        Revert to the original range and knots to fit over if you press the
        "Reset fit range" button. Bound to the "revert_xrange_btn" in initUI()
        """

        self.xlim = [[min(self.x), max(self.x)]]
        self.ind = np.arange(len(self.x))
        self.not_ind = None
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


if __name__ == '__main__':
    pass

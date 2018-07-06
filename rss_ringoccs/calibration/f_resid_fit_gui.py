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
from scipy.interpolate import interp1d

import tkinter
from tkinter import Tk, IntVar, StringVar, messagebox
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
        is_chord (bool): Record if occultation is a chord occultation
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

    Notes:
        [1] Will spit out cryptic error messages at you after you press "OK"
            and then press either yes or no in the dialogue box. I don't know
            how to stop these, but they seem to be pretty harmless
    """

    def __init__(self, parent, fit_inst, x, x_rho, y):
        """
        Args:
            parent: tkinter.Tk instance
            fit_inst: Instance of the FreqOffsetFit class
            x: SPM values of the residual frequency (y)
            x_rho: Radius values corresponding to x
            y: Residual frequency values in Hz
        """

        Frame.__init__(self, parent)

        self.parent = parent
        self.fit_inst = fit_inst
        self.x = x
        self.x_rho = x_rho
        self.y = y
        self.initUI()


    def initUI(self):
        """
        Initialize the user interface. Called by __init__()
        """

        self.parent.title('FResidFitGui')
        self.pack(fill=tkinter.BOTH, expand=1)

        # Default fit
        self.fit_deg = 3
        self.xlim = self.fit_inst._spm_include
        dummy, self.yfit = self.fit_inst.get_f_sky_resid_fit()

        # Indices inside and outside of freespace. Starts off with all
        #     highlighted like they're inside of freespace regions, despite
        #     the original default
        self.ind = np.arange(len(self.x))
        self.not_ind = None

        # Figure to plot on
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

        # Put fit order box on the window
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

        # Put fit range box on the window
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
        quit_btn = Button(self, text='OK')
        quit_btn.bind('<Button-1>', self.quit_app)
        quit_btn.pack()


    def _get_fit(self):
        """
        Get the polynomial fit to residual frequency. Called by __init__(),
        adjust_deg(), adjust_range(), and revert_range()
        """

        self.fit_inst.set_f_sky_resid_fit(poly_order=self.fit_deg,
            spm_include=self.xlim, USE_GUI=False)
        _x, f_sky_resid_fit = self.fit_inst.get_f_sky_resid_fit()
        return f_sky_resid_fit


    def _get_rho_tick_labels(self, spm_tick_labels):
        """
        Given tick labels in SPM, get the new tick labels in radius

        Args:
            spm_tick_labels (list): SPM tick labels to get radius matches for
        """
        spm_to_rho_func = interp1d(self.x, self.x_rho)
        rho_tick_labels = spm_to_rho_func(spm_tick_labels)
        return ['%.1f' % rho_tick_label for rho_tick_label in rho_tick_labels]


    def plot_data(self):
        """
        Make the original matplotlib.pyplot plot on the canvas. Called by
        initUI()
        """

        # Plot residual frequency and fit
        self.ax.plot(self.x, self.y, color='g')
        self.ax.plot(self.x, self.yfit, color='r')
        self.ax.set_xlabel('SPM')
        self.ax.set_ylabel('Hz')

        # Make rho axis on upper x-axis
        self.ax_rho = self.ax.twiny()
        rho_diff = np.diff(self.x_rho)
        if (np.any(rho_diff > 0)) & (np.any(rho_diff < 0)):
            self.is_chord = True
            ax_rho_ticks = (self.ax.get_xticks())[
                (self.ax.get_xticks() >= min(self.x))
                & (self.ax.get_xticks() <= max(self.x))]
            self.ax_rho.set_xlim(self.ax.get_xlim())
            self.ax_rho.set_xticks(ax_rho_ticks)
            self.ax_rho.set_xticklabels(self._get_rho_tick_labels(
                ax_rho_ticks))
            self.ax_rho.set_xlabel('Rho (km)')
        else:
            self.is_chord = False
            self.ax_rho.plot(self.x_rho, self.y, alpha=0)
        self.ax_rho.set_xlabel('Rho (km)')

        # Show the canvas
        self.canvas = FigureCanvasTkAgg(self.f, master=self.parent)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH,
            expand=1)

        # Put the matplotlib toolbar on the screen
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

        # Indices to make highlighted or dimmed based on the freespace regions
        ind = self.ind
        not_ind = self.not_ind

        # Plot the residual frequency within the freespace regions and its fit
        self.ax.cla()
        self.ax.plot(self.x[ind], self.y[ind], color='g')
        self.ax.plot(self.x, self.yfit, color='r')
        self.ax.set_ylim([min(self.y[ind]) - 1.0, max(self.y[ind]) + 1.0])

        # Plot the radius scale on the upper x-axis
        if self.is_chord:
            ax_rho_ticks = (self.ax.get_xticks())[
                (self.ax.get_xticks() >= min(self.x))
                & (self.ax.get_xticks() <= max(self.x))]
            self.ax_rho.set_xlim(self.ax.get_xlim())
            self.ax_rho.set_xticks(ax_rho_ticks)
            self.ax_rho.set_xticklabels(self._get_rho_tick_labels(
                ax_rho_ticks))
        else:
            self.ax_rho.plot(self.x_rho, self.y, alpha=0)
        self.ax_rho.set_xlabel('Rho (km)')

        # Blue vertical lines to designate where the freespace regions are. Not
        #     doing anymore because I think it looks bad, and doesn't really
        #     help any since it's easy enough to look at the entered values
        #     for fit range and find where they are on the x-axis
#         y_arrow = np.mean(self.ax.get_ylim())
#         for i in range(len(self.xlim)):
#             self.ax.axvline(self.xlim[i][0], color='b')
#             self.ax.axvline(self.xlim[i][1], color='b')
#             _x_arrow = float(self.xlim[i][0])
#             _dx_arrow = float(np.diff(self.xlim[i]))
#             self.ax.arrow(_x_arrow, y_arrow, _dx_arrow, 0,
#                 color='b', head_width=0.1, head_length=50,
#                 length_includes_head=True)

        # Plot the residual frequency outside of the freespace regions
        if not_ind is not None:
            self.ax.plot(self.x[not_ind], self.y[not_ind], color='g', alpha=0.1)

        self.ax.set_xlabel('SPM')
        self.ax.set_ylabel('Hz')
        self.canvas.show()


    def adjust_deg(self, e):
        """
        Purpose:
        Adjust the degree of the fit if you adjusted the "Fit order" combobox.
        Bound to the "combo" variable in initUI()

        Args:
            e: Event received from selecting a fit order from the drop-down menu
        """

        # Set attribute for the new fit order
        new_deg = e.widget
        self.lvar.set(new_deg.get())
        self.fit_deg = int(new_deg.get())

        # Adjust the fit and re-plot
        self.yfit = self._get_fit()
        self.update_plot()


    def adjust_range(self):
        """
        Purpose:
        Adjust ranges to fit over if you changed what's entered in the text
        box and then press the "Set fit range (SPM)" button. Bound to the
        "set_xrange_btn" variable in initUI()
        """

        # Get the freespace region entry from the textbox, using semicolons to
        #     separate different regions, and commas to separate the minimum
        #     and maximum within a region
        try:
            xlim_entry = self.xlim_entry.get()
            xlim = []
            split1 = xlim_entry.split(';')
            for i in range(len(split1)):
                split2 = (split1[i]).split(',')
                _xlim = [float(split2[0]), float(split2[1])]
                xlim.append(_xlim)
            self.xlim = xlim
        except ValueError:
            print('WARNING (FResidFitGui.adjust_range): Illegal input in the '
                + '"Fit range (SPM)" text box. Returning to default fit ranges')
            self.revert_range()
            return

        # New indices of self.x falling inside of freespace regions
        ind = []
        for i in range(len(xlim)):
            ind.append(np.argwhere((self.x >= xlim[i][0]) &
                (self.x <= xlim[i][1])))
        ind = np.reshape(np.concatenate(ind), -1)
        self.ind = ind

        # New indices of self.x falling outside of freespace regions
        not_ind = []
        not_ind.append(np.argwhere(self.x < xlim[0][0]))
        for i in range(len(xlim)-1):
            not_ind.append(np.argwhere((self.x > xlim[i][1]) &
                (self.x < xlim[i+1][0])))
        not_ind = np.reshape(np.concatenate(not_ind), -1)
        self.not_ind = not_ind

        # Make a new fit and plot it
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


    def revert_range(self):
        """
        Purpose:
        Revert to the original range and knots to fit over if you press the
        "Reset fit range" button. Bound to the "revert_xrange_btn" in initUI()
        """

        # Get default radii exclusion zones from the instance of FreqOffsetFit
        #     and set the xlim, ind, and not_ind attributes accordingly
        rho_exclude = self.fit_inst._rho_exclude
        self.xlim = []
        for i in range(len(rho_exclude)-1):
            self.xlim.append(
                [self.x[np.argmin(abs(self.x_rho - rho_exclude[i][1]))],
                self.x[np.argmin(abs(self.x_rho - rho_exclude[i+1][0]))]])
        self.ind = np.arange(len(self.x))
        self.not_ind = None

        # Make the new fit and plot it
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()


    def quit_app(self, e):
        """
        Purpose:
        Called when you hit the "OK" button. Makes sure that you really want
        to exit the GUI
        """

        # Check if you really want to exit the GUI, and gives a firm reminder
        #     to set fit parameters before continuing
        ret = messagebox.askquestion('Question', 'Are you sure of this fit? '
            + 'If not, then click "No", adjust the fit parameters, and '
            + 'REMEMBER TO HIT THE "Set Fit Range (SPM)" BUTTON!!!')

        # Destroy GUI if yes, return to GUI if no
        if ret == 'yes':
            self.quit()
            self.parent.destroy()
        else:
            return


if __name__ == '__main__':
    pass

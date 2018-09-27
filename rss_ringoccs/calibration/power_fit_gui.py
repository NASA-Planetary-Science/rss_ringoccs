"""

power_fit_gui.py

Purpose: Copied from fit_example_power_norm.py and edited to use in
         power_normalization_v2.py.

Notes on text box entry:
    [1] In the "Fit ranges" text box, you need to separate freespace
        regions by a semicolon, and separate the minimum and maximum in each
        freespace region by a comma. For example, if you want to use the
        freespace regions as [30500, 33000], [36000, 37000], and
        [38000, 40000], then the entry into the textbox would be:
        "30500, 33000 ; 36000, 37000 ; 38000, 40000"
    [2] In the "Knots" text box, you need to separate knots by a comma. For
        example, if you want to use the SPM knots 31000, 36500, and 39000,
        then you would just enter "31000, 36500, 39000"
    [3] I usually have one knot per freespace region

Revisions:
        fit_example_power_norm.py
    2018 May 09 - gsteranka - Original version
        power_fit_gui.py
    2018 May 10 - gsteranka - Copied to use in GitHub code
    2018 Sep 19 - sflury - Updated GUI to have predicted free space regions
       highlighted using matplotlib.pyplot.fill_between. Changed plot colors
       to higher contrast and color-blind accessible.
    2018 Sep 26 - jfong - update to display fit read from file
    2018 Sep 26 - sflury - update to display correct fit order when GUI pops up
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


class PowerFitGui(Frame):
    """
    Class for a GUI to make a spline fit to power. Includes buttons to adjust
    parameters of the fit, including fit order, freespace locations, and knot
    locations. This is not intended to be called by a user directly, but rather
    to be used inside of the Normalization class when the USE_GUI keyword is
    set to True (which is the default).

    Args:
        parent (tkinter.Tk):
            Instance of tkinter.Tk class. This is the parent
            that the GUI is based on
        norm_inst:
            Instance of the Normalization class
        spm_vals_down (np.ndarray):
            Array of SPM values for the power values to
            fit. Downsampling is done using scipy.signal.resample_poly in
            the Normalization class. Downsampled to the spacing "dt_down" in
            the Normalization class (1/2 second by default).
        rho_km_vals_down (np.ndarray):
            Rho values corresponding to spm_vals_down
        p_obs_down (np.ndarray):
            Array of downsampled power values to fit.
            Gotten from downsampling frequency-corrected complex samples
        spm_fit (np.ndarray):
            SPM values to evaluate the fit for

    Attributes:
        ax:
            Axes instance corresponding to the matplotlib.pyplot figure on
            the GUI
        canvas:
            Canvas on which the matplotlib.pyplot and the widgets are put
        cvar:
            Text variable in the combo box for selecting fit order
        f:
            Matplotlib.pyplot figure on the GUI
        fit_deg (int):
            Degree of the spline fit. Starts off at 1
        ind (np.ndarray):
            Index numbers to the "x" attribute of the regions
            being fit (specified by the "xlim" attribute)
        is_chord (bool):
            True if chord occultation, False if not
        knots_entry:
            Entry box for specifying the knots used in the fit.
            Changing the entry to this box changes the "knots_spm" attribute
            (see below)
        knots_spm (list):
            Knots of the fit. Adjustable in the GUI in the box
            labeled "Knots". Each knot must be separated by a comma. For
            example, if you wanted to use the SPM knots 31000, 36500, and
            39000, then you would just enter "31000, 36500, 39000" (without
            the quotation marks). I usually use one knot per region being fit,
            which in this case is the regions in the "xlim" attribute (below)
        lvar:
            Variable for the label next to the combo box for selecting fit
            order
        norm_inst:
            The input instance of the Normalization class
        not_ind (np.ndarray):
            Index numbers of the "x" attribute of the regions
            not being fit (specified by the "xlim" attribute)
        parent (tkinter.Tk):
            Input parent
        toolbar:
            Matplotlib.pyplot toolbar that appears below the
            matplotlib.pyplot plot
        x (np.ndarray):
            The input spm_vals_down values
        x_rho (np.ndarray):
            The input rho_km_vals_down values
        xfit (np.ndarray):
            The input spm_fit values
        xlim (list):
            Set of regions to make the fit over. Starts off as the
            full range of SPM. Adjustable in the GUI, in the box labeled
            "Fit ranges". To specify, separate minimum and maximum values by a
            comma, and separate ranges by a semicolon. For example, if you want
            to make a fit using the SPM ranges [30500, 33000], [36000, 37000],
            and [38000, 40000], then you would enter into the box:
            "30500, 33000 ; 36000, 37000 ; 38000, 40000" (without the quotation
            marks)
        xlim_entry:
            Entry box for specifying regions to make the fit over.
            Entering into this box changes the "xlim" attribute (see above)
        y (np.ndarray):
            The input p_obs_down values
        yfit (np.ndarray):
            Resulting fit from the GUI. This is the attribute
            used to access the fit after the program is closed

    Notes:
        [1] Will spit out cryptic error messages at you after you press "OK"
            and then press either yes or no in the dialogue box. I don't know
            how to stop these, but they seem to be pretty harmless
    """

    def __init__(self, parent, norm_inst, spm_vals_down, rho_km_vals_down,
            p_obs_down, spm_fit):
        """
        Args:
            parent:
                tkinter.Tk instance
            norm_inst:
                Instance of the Normalization class
            spm_vals_down (np.ndarray):
                SPM values downsampled from raw resolution
            rho_km_vals_down (np.ndarray):
                Radius values corresponding to spm_vals_down
            p_obs_down (np.ndarray):
                Unnormalized power evaluated at
                spm_vals_down. Gotten from downsampling IQ_c and then
                calculating power as abs(IQ_c)**2
            spm_fit (np.ndarray):
                SPM values to evaluate the spline fit at
        """

        Frame.__init__(self, parent)

        self.norm_inst = norm_inst
        self.x = spm_vals_down
        self.x_rho = rho_km_vals_down
        self.y = p_obs_down
        self.xfit = spm_fit

        self.fit_deg = norm_inst.pnfp_splrep[2]
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

        self.parent.title('PowerFitGui')
        self.pack(fill=tkinter.BOTH, expand=1)

        self.f = Figure(figsize=(5, 4))
        self.ax = self.f.add_subplot(111)

        self.plot_data()

        self.cvar = IntVar()

        # Frame to put fit order box in
        fit_order_frame = LabelFrame(self, text='Fit Order')
        # Box to select fit order
        combo = Combobox(fit_order_frame, textvariable=self.cvar)
        combo.bind('<<ComboboxSelected>>', self.adjust_deg)
        combo['values'] = [i for i in range(1,6)]
        combo.current(int(self.fit_deg))
        combo.grid(row=0, column=0)

        # Variable outside of box that says current fit order
        self.lvar = IntVar()
        self.lvar.set(str(self.fit_deg))
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
        set_xrange_btn = Button(self, text='Set fit range and knots (SPM)',
            command=self.adjust_range_and_knots)
        set_xrange_btn.pack()

        # Button to revert to original x range
        revert_xrange_btn = Button(self, text='Reset fit range and knots',
            command=self.revert_range)
        revert_xrange_btn.pack()

        # Button to press when you're content with the fit and want to continue
        ok_btn = Button(self, text='OK')
        ok_btn.bind('<Button-1>', self.quit_app)
        ok_btn.pack(side=tkinter.LEFT, padx=120)

    def _get_fit(self):
        """
        Purpose:
        Get the spline fit to power. Called by __init__(), adjust_deg(),
        adjust_range_and_knots(), and revert_range()
        """

        spm_fit, spline_fit = self.norm_inst.get_spline_fit(
            spline_order=self.fit_deg, knots_spm=self.knots_spm,
            freespace_spm=self.xlim, USE_GUI=False, file_search=False)
        return spline_fit#, spline_rep

    def _get_rho_tick_labels(self, spm_tick_labels):
        """
        Purpose:
        Given tick labels in SPM, get the new tick labels in radius

        Args:
            spm_tick_labels (list):
                SPM tick labels to get corresponding radius values
        """
        spm_to_rho_func = interp1d(self.x, self.x_rho,
                                bounds_error=False,fill_value='extrapolate')
        rho_tick_labels = spm_to_rho_func(spm_tick_labels)
        return ['%.1f' % rho_tick_label for rho_tick_label in rho_tick_labels]

    def plot_data(self):
        """
        Purpose:
        Make the original matplotlib.pyplot plot on the canvas. Called by
        initUI()
        """

        # Shading of predicted free-space regions and normalization of power to
        #    max power
        self.rho_to_spm_func = interp1d(self.x_rho,self.x,
                                bounds_error=False,fill_value='extrapolate')
        self.rho_free = [[69100.0, 73500], [87350.0, 87450.0],
            [117690.0, 117780.0], [119850.0, 120020.0],
            [133500.0, 133650.0], [137000.0, 194400.0]]
        self.spm_free = [[self.rho_to_spm_func(xl[0]),self.rho_to_spm_func(xl[1])] for xl in self.rho_free]
        for xl in self.spm_free:
            self.ax.fill_between(xl,-1,2,color='0.75',zorder=0)

        self.ynorm = self.y/np.max(self.y)

        # Plot power and its default fit
        self.ax.plot(self.x, self.ynorm, color='0.0',zorder=1)
        self.ax.plot(self.xfit, self.yfit/np.max(self.y), color='r',lw=2,zorder=2)
        self.ax.set_xlabel('SPM')

        # Plot radius scale on upper x-axis
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
            self.ax_rho.plot(self.x_rho, self.ynorm, alpha=0)
            # Reverse x axis if ingress
            if np.all(rho_diff < 0):
                xlim_rho = self.ax_rho.get_xlim()
                self.ax_rho.set_xlim([xlim_rho[1], xlim_rho[0]])
        self.ax_rho.set_xlabel('Rho (km)')

        # set x and y limits
        self.ax.set_xlim(self.rho_to_spm_func(6e4),self.rho_to_spm_func(1.5e5))
        self.ax.set_ylim(-0.2,1.2)

        # Show the canvas
        self.canvas = FigureCanvasTkAgg(self.f, master=self.parent)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH,
            expand=1)

        # Matplotlib toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH,
            expand=1)

    def update_plot(self):
        """
        Purpose:
        Update the plot if an event happened (namely, you pressed a button or
        entered something new into a box). Called by adjust_deg(),
        adjust_range_and_knots(), and revert_range()
        """

        # Indices inside or outside of freespace range
        ind = self.ind
        not_ind = self.not_ind

        # plot power within freespace range and its fit
        self.ax.cla()

        # Shading of predicted free-space regions
        for xl in self.spm_free:
            self.ax.fill_between(xl,-1,2,color='0.75',zorder=0)
        # set x and y limits
        self.ax.set_xlim(self.rho_to_spm_func(6e4),self.rho_to_spm_func(1.5e5))
        self.ax.set_ylim(-0.2,1.2)
        # Plot power and its default fit
        self.ax.plot(self.x, self.ynorm, color='0.0',zorder=1)
        self.ax.plot(self.xfit, self.yfit/np.max(self.y), color='r',lw=2,zorder=2)

        if not_ind is not None:
            self.ax.plot(self.x[not_ind], self.y[not_ind]/np.max(self.y), color='k',
                alpha=0.1)

        # Plot radius scale on upper x-axis
        if self.is_chord:
            ax_rho_ticks = (self.ax.get_xticks())[
                (self.ax.get_xticks() >= min(self.x))
                & (self.ax.get_xticks() <= max(self.x))]
            self.ax_rho.set_xlim(self.ax.get_xlim())
            self.ax_rho.set_xticks(ax_rho_ticks)
            self.ax_rho.set_xticklabels(self._get_rho_tick_labels(
                ax_rho_ticks))
        else:
            self.ax_rho.plot(self.x_rho, self.ynorm, alpha=0)

        self.ax_rho.set_xlabel('Rho (km)')
        self.ax.set_xlabel('SPM')
        self.canvas.show()

    def adjust_deg(self, e):
        """
        Purpose:
        Adjust the degree of the fit if you adjusted the "Fit order" combobox.
        Bound to the "combo" variable in initUI()
        """

        # Get the new fit degree
        new_deg = e.widget
        self.lvar.set(new_deg.get())
        self.fit_deg = int(new_deg.get())

        # Get new fit and plot it
        self.yfit = self._get_fit()
        self.update_plot()

    def adjust_range_and_knots(self):
        """
        Purpose:
        Adjust ranges to fit over and the knots used if you changed what's
        entered in the text boxes and press the "Set fit range and knots"
        button. Bound to the "set_xrange_btn" variable in initUI()
        """

        # Get new freespace regions, using semicolons to separate regions
        #     and commas to separate minimum and maximum within a region.
        #     Revert to original fit if it fails
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
            print('WARNING (PowerFitGui.adjust_range_and_knots): Illegal '
                + 'input for "Fit ranges" text box. Reverting to original '
                + 'fit. Did you forget to enter a value before pressing '
                + 'the "Set fit ranges and knots (SPM)" button?')
            self.revert_range()
            return

        # Indices within freespace regions
        ind = []
        for i in range(len(xlim)):
            ind.append(np.argwhere((self.x >= xlim[i][0]) &
                (self.x <= xlim[i][1])))
        ind = np.reshape(np.concatenate(ind), -1)
        self.ind = ind

        # Indices outside of freespace regions
        not_ind = []
        not_ind.append(np.argwhere(self.x < xlim[0][0]))
        for i in range(len(xlim) - 1):
            not_ind.append(np.argwhere((self.x > xlim[i][1]) &
                (self.x < xlim[i + 1][0])))
        not_ind.append(np.argwhere(self.x > xlim[-1][1]))
        not_ind = np.reshape(np.concatenate(not_ind), -1)
        self.not_ind = not_ind

        # Get new knot values from user entry. Revert to original fit if it
        #     fails
        try:
            knots_entry = self.knots_entry.get()
            knots_entry_split = knots_entry.split(',')
            self.knots_spm = [float(knot) for knot in knots_entry_split]
        except ValueError:
            print('WARNING (PowerFitGui.adjust_range_and_knots): Illegal '
                + 'input for knots. Reverting to original fit. '
                + 'Did you forget to enter a value in the "Knots" text box '
                + 'before hitting the "Set fit range and knots (SPM)" button?')
            self.revert_range()
            return

        # Update fit and plot it
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()

    def revert_range(self):
        """
        Purpose:
        Revert to the original range and knots to fit over if you press the
        "Reset fit range and knots" button. Bound to the
        "revert_xrange_btn" in initUI()
        """

        # Get default freespace ranges from Normalization instance, and update
        #     xlim, ind, not_ind, and knots_spm attributes accordingly
        freespace_km = self.norm_inst._freespace_km
        self.xlim = []
        for _range in freespace_km:
            self.xlim.append([self.x[np.argmin(abs(self.x_rho - _range[0]))],
                self.x[np.argmin(abs(self.x_rho - _range[1]))]])
        self.knots_spm = []
        for _knot_km in self.norm_inst._knots_km:
            self.knots_spm.append(
                self.x[np.argmin(abs(self.x_rho - _knot_km))])
        self.ind = np.arange(len(self.x))
        self.not_ind = None

        # Update fit and plot it
        new_yfit = self._get_fit()
        self.yfit = new_yfit
        self.update_plot()

    def quit_app(self, e):
        """
        Called when you hit the "OK" button. Makes sure that you really want
        to exit the GUI
        """

        # Check if you really want to exit the GUI, and gives a firm reminder
        #     to set fit parameters before continuing
        ret = messagebox.askquestion('Question', 'Are you sure of this fit? '
            + 'If not, then click "No", adjust the fit parameters, and '
            + 'REMEMBER TO HIT THE "Set fit range and knots (SPM)" BUTTON!!!')

        # Destroy GUI if yes, return to GUI if no
        if ret == 'yes':
            self.quit()
            self.parent.destroy()
        else:
            return

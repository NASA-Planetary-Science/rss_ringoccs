"""
    Submodule:
        diffraction_correction
    Purpose:
        Provide tools for executing bash scripts
        from within python, and creating summary
        documents using LaTeX.
    Dependencies:
        #. subprocess
        #. os
        #. time
"""
import subprocess
import os
import time

from . import error_check

def shell_execute(script):
    """
        Purpose:
            Execute a shell script from within Python.
        Variables:
            :script (*str*):
                A list containing the path to the shell
                script and variables. For example:
                    script = ['path/to/script','v1',...,'vn']
                Elements must be strings.
        Outputs:
            :Process (*object*):
                An instance of the Popen class from the
                subprocess module. This contains attributes
                such as the arguments passed to it, and other
                system details.
        Notes:
            This routine has only been tested using scripts
            written in Bash, and on standard Unix commands.
        Examples:
            Suppose we have the following shell script test.sh:
                #!/bin/bash
                printf "Hello, World! My name is %s!" "$1"
            Run this shell script inside of Python:
                In [1]: from rss_ringoccs.tools import sys_tools
                In [2]: sys_tools.shell_execute(['./test.sh','Bob'])
                        Hello World! My name is Bob!
            We can also execute simple Unix commands.
                In [1]: import diffcorr as dc
                In [2]: a = dc.shell_execute(["echo","Bob"])
                        Bob
                In [3]: a.args
                Out[3]: [' echo Bob']
    """
    fname = "tools.sys_tools.shell_execute"
    error_check.check_type(script, list, "script", fname)

    string = ""
    for x in script:
        error_check.check_type(x, str, "script[i]", fname)
        string = "%s %s" % (string, x)

    return subprocess.Popen([string],shell=True)

def date_string():
    """
        Purpose:
            Create string of the form "yyyy_mm_dd_HH_MM_SS_"
        Outputs:
            :date (*str*):
                Current date "year/month/day/hour/minute/second"
    """
    strings = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = strings.split(',')
    date=""
    for x in t:
        date = "%s%s_" % (date,x)
    return date

def make_executable(path):
    error_check.check_type(path, str, "path", "tools.sys_tools.make_executable")

    mode = os.stat(path).st_mode

    # Copy R bits to X.
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)

def latex_summary_doc(pdffil, res_km, outfilename):
    fname = "tools.sys_tools.latex_summary_doc"
    error_check.check_type(pdffil, str, "pdffil", fname)
    error_check.check_type(outfilename, str, "outfilename", fname)
    res = error_check.check_type_and_convert(res_km, float, "res_km", fname)

    res = "%sM" % str(int(res_km*1000.0))

    var = pdffil.split("/")[-1]
    var = var.split("_")
    rev = var[0][3:6]
    doy = var[3]
    occ = var[5]
    if occ == 'I':
        profdir = 'Ingress'
    if occ == 'E':
        profdir = 'Egress'
    year = var[2]
    band = var[4][0:3]
    dsn = var[4][1:]

    geo = "RSS\_%s\_%s\_%s\_%s\_GEO.TAB" % (year, doy, band, occ)
    cal = "RSS\_%s\_%s\_%s\_%s\_CAL.TAB" % (year, doy, band, occ)
    tau = "RSS\_%s\_%s\_%s\_%s\_TAU\_%s.TAB" % (year, doy, band, occ, res)

    LaTeXFile = r"""
        \documentclass{article}
        \usepackage{geometry}
        \geometry{a4paper, margin = 1.0in}
        \usepackage{graphicx, float}
        \usepackage{lscape}
        \usepackage[english]{babel}
        \usepackage[dvipsnames]{xcolor}
        \usepackage[font={normalsize}, labelsep=colon]{caption}
        \addto\captionsenglish{\renewcommand{\figurename}{Fig.}}
        \newcommand{\thePDF}{%s}
        \newcommand{\theREV}{%s}
        \newcommand{\theDOY}{%s}
        \newcommand{\theRES}{%s}
        \newcommand{\theOCC}{%s}
        \newcommand{\theGEO}{%s}
        \newcommand{\theCAL}{%s}
        \newcommand{\theTAU}{%s}
        \newcommand{\theYEAR}{%s}
        \newcommand{\theBAND}{%s}
        \newcommand{\thePROFDIR}{%s}
        \newcommand{\theDSN}{%s}
        \setlength{\parindent}{0em}
        \setlength{\parskip}{0em}
        \begin{document}
            \pagenumbering{gobble}
            \begin{center}
                \LARGE{\texttt{
                    RSS\textunderscore\theYEAR%%
                    \textunderscore\theDOY\textunderscore\theBAND%%
                    \textunderscore\theOCC}\\[2.0ex]
                    Rev\theREV\
                    Cassini Radio Science Ring Occultation:\\[1.0ex]
                    Geometry, Data Calibration,
                    and Reconstructed\\[1.0ex]
                    Optical Depth and
                    Phase Shift Profiles\\[1.0ex]
                    at \theRES\
                    Resolution\\[2.5ex]}
                \large{\today}
            \end{center}
            \vspace{2ex}
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[
                        page=1,
                        trim={0.0cm 0.0cm 0.0cm 0.9cm},
                        clip,
                        width=\textwidth
                    ]{\thePDF}
                }
                \caption*{
                    The radio occultation track as seen
                    looking down on the ring plane.
                    The solid
                    \textcolor{red}{red} track is relative to
                    a reference direction defined by the
                    direction to Earth, The solid
                    \textcolor{blue}{blue}
                    track is relative to an initial reference
                    direction defined by the ascending
                    node of J2000 on the ring plane.
                }
            \end{figure}
            \newpage
            \pagenumbering{arabic}
            \begin{table}[H]
                \centering
                \begin{tabular}{ll}
                    \hline
                    Symbol&Parameter Name\\
                    \hline
                    $t_{\scriptsize{\textrm{ERT}}}$&
                    OBSERVED EVENT TIME (Earth Receiving Time)\\
                    $t_{\scriptsize{\textrm{RET}}}$&
                    RING EVENT TIME\\
                    $t_{\scriptsize{\textrm{SCET}}}$&
                    SPACECRAFT EVENT TIME\\
                    $\rho$&RING RADIUS\\
                    $\phi_{\scriptsize{\textrm{J2K}}}$&
                    RING LONGITUDE\\
                    $\phi_{\scriptsize{\textrm{E}}}$&
                    OBSERVED RING AZIMUTH\\
                    $B$&OBSERVED RING ELEVATION\\
                    $D$&SPACECRAFT TO RING INTERCEPT DISTANCE\\
                    $V_{\scriptsize{\textrm{rad}}}$&
                    RING INTERCEPT RADIAL VELOCITY\\
                    $V_{\scriptsize{\textrm{az}}}$&
                    RING INTERCEPT AZIMUTHAL VELOCITY\\
                    $F$&FRESNEL SCALE\\
                    $R_{\scriptsize{\textrm{imp}}}$&
                    IMPACT RADIUS\\
                    $r_{X}$&SPACECRAFT POSITION X\\
                    $r_{y}$&SPACECRAFT POSITION Y\\
                    $r_{z}$&SPACECRAFT POSITION Z\\
                    $v_{x}$&SPACECRAFT VELOCITY X\\
                    $v_{y}$&SPACECRAFT VELOCITY Y\\
                    $v_{z}$&SPACECRAFT VELOCITY Z\\
                    $\theta_{\scriptsize{\textrm{EL}}}$&
                    OBSERVED SPACECRAFT ELEVATION\\
                    \hline
                \end{tabular}
                \caption[Glossary of Parameters from the Geo File]{
                    Glossary of parameters in file \theGEO.
                    See companion label (.LBL) file for description
                    of parameters.
                }
                \label{tab:easydata_glossary_of_geo_file}
            \end{table}
            \begin{table}[H]
                \centering
                \begin{tabular}{l l}
                    \hline
                    Symbol&Parameter Name\\
                    \hline
                    $t_{\scriptsize{\textrm{ERT}}}$&
                    OBSERVED EVENT TIME\\
                    $f_{\scriptsize{\textrm{sky}}}$&
                    SKY FREQUENCY\\
                    $f_{\scriptsize{\textrm{resid}}}$&
                    RESIDUAL FREQUENCY\\
                    $P_{\scriptsize{\textrm{free}}}$&
                    FREESPACE POWER\\
                    \hline
                \end{tabular}
                \caption[Glossary of Data from the Cal File]{
                    Glossary of calibration data in file
                    \theCAL. See companion label (.LBL)
                    file for description of the data.
                }
                \label{tab:easydata_glossary_from_cal_file}
            \end{table}
            \begin{table}[H]
                \centering
                \begin{tabular}{l l}
                    \hline
                    Symbol&Parameter Name\\
                    \hline
                    $\rho$&RING RADIUS\\
                    $\Delta\rho_{\scriptsize{\textrm{IP}}}$
                    &RADIUS CORRECTION DUE TO IMPROVED POLE\\
                    $\Delta\rho_{\scriptsize{\textrm{TO}}}$&
                    RADIUS CORRECTION DUE TO TIMING OFFSET\\
                    $\phi_{\scriptsize{\textrm{RL}}}$&
                    RING LONGITUDE\\
                    $\phi_{\scriptsize{\textrm{ORA}}}$&
                    OBSERVED RING AZIMUTH\\
                    $P$&NORMALIZED SIGNAL POWER\\
                    $\tau$&NORMAL OPTICAL DEPTH\\
                    $\phi$&PHASE SHIFT\\
                    $\tau_{\scriptsize{\textrm{TH}}}$
                    &NORMAL OPTICAL DEPTH THRESHOLD\\
                    $t_{\scriptsize{\textrm{ERT}}}$&
                    OBSERVED EVENT TIME
                    (Earth Recieving Time)\\
                    $t_{\footnotesize{\textrm{RET}}}$&
                    RING EVENT TIME\\
                    $t_{\footnotesize{\textrm{SCET}}}$&
                    SPACECRAFT EVENT TIME\\
                    $B$&OBSERVED RING ELEVATION\\
                    \hline
                \end{tabular}
                \caption[Glossary of Parameters in Tau File]{
                    Glossary of optical depth, phase shift,
                    and selected geometry and time parameters
                    in \theTAU.
                    See companion label
                    (.LBL) files for description of the data.
                }
                \label{tab:easydata_parameters_from_tau_file}
            \end{table}
            \newpage
            \begin{figure}[H]
                \centering
                \large{\textbf{View from Earth}}\par
                \includegraphics[
                    page=2,
                    trim={0.0in 1.0in 0.0in 1.0in},
                    clip,
                    width=0.75\textwidth
                ]{\thePDF}
                \caption{View of RSS ring occultation from Earth,
                         with occultation track plotted in blue
                         (dashed when blocked by planet). 30-min markers are
                         plotted as solid blue dots. Black solid line
                         represents direction to Earth.}
            \end{figure}
            \vspace{30ex}
            \begin{figure}[H]
                \centering
                \large{\textbf{View from North Pole}}\par
                \includegraphics[
                    page=3,
                    trim={0.0in 0.65in 0.0in 0.65in},
                    clip,
                    width=\textwidth
                ]{\thePDF}
                \caption{View of RSS ring occultation from north pole of planet,
                         with occultation track plotted in blue
                         (dashed when blocked by planet). 30-min markers are
                         plotted as solid blue dots. Black solid line
                         represents direction to Earth.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \includegraphics[page=4, width=\textwidth]{\thePDF}
                \caption{Rev \theREV\ - \thePROFDIR;
                         selected occultation parameters.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \includegraphics[page=5, width=\textwidth]{\thePDF}
                \caption{Calibration data in file \theCAL.
                         The frequency residuals data
                         (the smooth curve, in the second panel)
                         is used to steer the carrier signal to the
                         middle of the recording bandwidth. The
                         free-space power data (the smooth curve in
                         the third panel) is used to normalize signal
                         power measurements so that the corresponding
                         optical depth has nearly zero value in the
                         absence of rings. Least-square fitting techniques
                         to frequency and power estimates of the direct
                         signal (the green curves in the second
                         and third panels, respectively) are used to
                         compute the calibration data.}
            \end{figure}
            \newpage
            \begin{landscape}
                \begin{figure}[H]
                    \centering
                    \includegraphics[page=6]{\thePDF}
                    \caption{Observing DSN station (DSS-\theDSN) elevation
                             angle (in \textcolor{magenta}{magenta})
                             superimposed on ring profile at \theRES M
                             resolution from data file \theTAU.
                             \textcolor{blue}{Blue}: Reconstructed normal
                             optical depth; \textcolor{cyan}{Cyan}:
                             Free-space baseline;
                             \textcolor{red}{Red}:Threshold optical depth
                             (measurement SNR $\simeq$ 1). Optical depth is
                             plotted increasing downward, the same direction
                             as increasing direct signal extinction.}
                \end{figure}
            \end{landscape}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=7, width=\textwidth]{\thePDF}
                }
                \caption{Observing DSN station (DSS-\theDSN) elevation angle
                         (in \textcolor{magenta}{magenta}) superimposed on ring
                         profile at \theRES M resolution from data file \theTAU.
                         \textcolor{blue}{Blue}: Reconstructed normal optical
                         depth; \textcolor{cyan}{Cyan}: Free-space baseline;
                         \textcolor{red}{Red}:Threshold optical depth
                         (measurement SNR $\simeq$ 1). Optical depth is plotted
                         increasing downward, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=8, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical depth
                         (measurement SNR $\simeq$ 1). Optical depth is plotted
                         downword, the same direction as increasing direct
                         signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=9, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=10, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=11, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=12, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=13, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical depth
                         (measurement SNR $\simeq$ 1). Optical depth is plotted
                         downword, the same direction as increasing direct
                         signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=14, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical depth
                         (measurement SNR $\simeq$ 1). Optical depth is plotted
                         downword, the same direction as increasing direct
                         signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=15, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical depth
                         (measurement SNR $\simeq$ 1). Optical depth is plotted
                         downword, the same direction as increasing direct
                         signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=16, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=17, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=18, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=19, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=20, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=21, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=22, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=23, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ $\tau$-profile from data file
                         \theTAU. \textcolor{blue}{Blue}: Reconstructed normal
                         optical depth; \textcolor{cyan}{Cyan}: Free-space
                         baseline; \textcolor{red}{Red}: Threshold optical
                         depth (measurement SNR $\simeq$ 1). Optical depth is
                         plotted downword, the same direction as increasing
                         direct signal extinction.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=24, width=\textwidth]{\thePDF}
                }
                \caption{Rev\theREV\-\theOCC\ DSS-\theDSN\ $\phi$-profile
                         from data file \theTAU. Phase wrapping occurs at
                         the $\pm$180$^{\circ}$ boundaries.}
            \end{figure}
        \end{document}
    """ % (pdffil, rev, doy, res, occ, geo, cal, tau, year, band, profdir, dsn)
    if (len(outfilename) >= 4):
        ext = outfilename[-4:]
        ext = ext.lower()
        if ((ext == ".pdf") or (ext == ".tex")):
            outfilename = outfilename[:-4]
        else:
            pass
    else:
        pass

    texvar = outfilename.split("/")
    if (len(texvar) > 1):
        out = "%s/" % (texvar[0])
        for i in range(1, len(texvar)-1):
            out = "%s/%s/" % (out, texvar[i])
    else:
        out = "$PWD"

    TexName = texvar[-1]
    TexFileName = '%s.tex' % TexName
    TexFile = open(TexFileName,'w')
    TexFile.write(LaTeXFile)
    TexFile.close()

    TexScript = """
        #!/bin/bash
        pdflatex %s.tex > tex.out
        pdflatex %s.tex > tex.out
        rm tex.out
        mv %s.pdf %s
        rm %s.log
        rm %s.aux
        rm %s.tex
        rm %s.sh
    """ % (TexName, TexName, TexName, out, TexName, TexName, TexName, TexName)

    TexSh = open("%s.sh" % TexName,'w')
    TexSh.write(TexScript)
    TexSh.close()

    make_executable('%s.sh' % TexName)
    shell_execute(['./%s.sh' % TexName])

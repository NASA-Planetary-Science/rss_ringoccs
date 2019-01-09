import subprocess
import os
import time

def shell_execute(script):
    """
        Function:
            shell_execute
        Purpose:
            Execute a shell script from within Python.
        Variables:
            script:     A list containing the path to the shell
                        script and variables. For example:
                            script = ['path/to/script','v1',...,'vn']
                        Elements must be strings.
        Outputs:
            Process:    An instance of the Popen class from the
                        subprocess module. This contains attributes
                        such as the arguments passed to it, and other
                        system details.
        Dependencies:
            [1] subprocess
        Notes:
            This routine has only been tested using scripts written
            in Bash, and on standard Unix commands.
        References:
            [1] https://docs.python.org/3/library/subprocess.html
            [2] https://stackoverflow.com/questions/
                3777301/how-to-call-a-shell-script-from-python-code
        Examples:
            Suppose we have the following shell script test.sh:
                #!/bin/bash
                printf "Hello, World! My name is %s!" "$1"
            Run this shell script inside of Python:
                In [1]: import diffcorr as dc
                In [2]: dc.shell_execute(['./test.sh','Bob'])
                        Hello World! My name is Bob!
            We can also execute simple Unix commands.
                In [1]: import diffcorr as dc
                In [2]: a = dc.shell_execute(["echo","Bob"])
                        Bob
                In [3]: a.args
                Out[3]: [' echo Bob']
        History:
            Created: RJM - 2018/05/16 5:49 P.M.
    """
    if not isinstance(script, list):
        raise TypeError(
            "Error: rss_ringoccs.tools.sys_tools: shell_execute\n"
            "\n\tscript must be a list of strings\n"
            "\tYour input has type: %s\n"
            % (type(script).__name__)
        )

    string = ""
    for x in script:
        if (not isinstance(x, str)):
            raise TypeError(
                "Error: rss_ringoccs.tools.sys_tools: shell_execute\n"
                "\n\tscript must be a list of strings\n"
                "\tOne of your inputs has type: %s\n"
                % (type(x).__name__)
            )
        else:
            string = "%s %s" % (string,x)
    Process=subprocess.Popen([string],shell=True)
    return Process

def date_string():
    """
        Function:
            date_string
        Purpose:
            Create string of the form "yyyy_mm_dd_HH_MM_SS_"
        Variables:  There are no variables to this function.
        Outputs:
            date:   Current date "year/month/day/hour/minute/second"
        Dependencies:
            [1] time
        Notes:
            The end of the string has an underscore "_"
        Examples:
            Get the current date.
                In [1]: import diffcorr as dc
                In [2]: dc.date_string()
                Out[2]: '2018_05_27_11_28_22_'
        History:
            Created: RJM - 2018/05/27 11:27 A.M.
    """
    strings = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = strings.split(',')
    date=""
    for x in t:
        date = "%s%s_" % (date,x)
    return date

def make_executable(path):
    if not isinstance(path, str):
        raise TypeError(
            "Error: rss_ringoccs.tools.sys_tools: make_executable\n"
            "\n\tpath must be a string\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(path).__name__)
        )
    else:
        pass

    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

def latex_summary_doc(pdffil, resolution, outfilename):
    if not isinstance(pdffil, str):
        raise TypeError(
            "Error: rss_ringoccs.tools.sys_tools: latex_summary_doc\n"
            "\n\tpdffil must be a string\n"
            "\tYour input has type: %s\n"
            % (type(pdffil).__name__)
        )
    elif not isinstance(outfilename, str):
        raise TypeError(
            "Error: rss_ringoccs.tools.sys_tools: latex_summary_doc\n"
            "\n\toutfilename must be a string\n"
            "\tYour input has type: %s\n"
            % (type(outfilename).__name__)
        )
    try:
        res = str(resolution)
    except:
        raise TypeError(
            "Error: rss_ringoccs.tools.sys_tools: latex_summary_doc\n"
            "\n\tresolution must be a floating point number\n"
            "\tYour input has type: %s\n"
            % (type(resolution).__name__)
        )
    var = pdffil.split("/")
    var = var[-1]
    var = var.split("_")
    rev = var[0][3:6]
    doy = var[3]
    occ = var[5]
    year = var[2]
    band = var[4][0:3]

    geo = "RSS\_%s\_%s\_%s\_%s\_GEO.TAB" % (year, doy, band, occ)
    cal = "RSS\_%s\_%s\_%s\_%s\_CAL.TAB" % (year, doy, band, occ)
    tau = "RSS\_%s\_%s\_%s\_%s\_TAU\_%sKM.TAB" % (year, doy, band, occ, res)

    LaTeXFile = r"""
        \documentclass{article}
        \usepackage{geometry}
        \geometry{a4paper, margin = 1.0in}
        \usepackage{graphicx, float}
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
        \setlength{\parindent}{0em}
        \setlength{\parskip}{0em}
        \begin{document}
            \pagenumbering{gobble}
            \begin{center}
                        \LARGE{\texttt{RSS\textunderscore\theYEAR%%
                               \textunderscore\theDOY\textunderscore%%
                               \theBAND\textunderscore\theOCC}\\[2.0ex]
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
                    direction to Earth, The dashed
                    \textcolor{blue}{blue}
                    track is relative to a reference
                    direction defined by the ascending
                    node of J2000 on the ring plane.
                }
            \end{figure}
            \newpage
            \pagenumbering{arabic}
            \begin{table}[H]
                \centering
                \begin{tabular}{l l}
                    \hline
                    Symbol&Parameter Name\\
                    \hline
                    $t_{OET}$&OBSERVED EVENT TIME\\
                    $t_{RET}$&RING EVENT TIME\\
                    $t_{SET}$&SPACECRAFT EVENT TIME\\
                    $\rho$&RING RADIUS\\
                    $\phi_{RL}$&RING LONGITUDE\\
                    $\phi_{ORA}$&OBSERVED RING AZIMUTH\\
                    $B$&OBSERVED RING ELEVATION\\
                    $D$&SPACECRAFT TO RING INTERCEPT DISTANCE\\
                    $\partial\rho/\partial{t}$&
                    RING INTERCEPT RADIAL VELOCITY\\
                    $\partial\theta/\partial t$&
                    RING INTERCEPT AZIMUTHAL VELOCITY\\
                    $F$&FRESNEL SCALE\\
                    $R_{impact}$&IMPACT RADIUS\\
                    $r_x$&SPACECRAFT POSITION X\\
                    $r_y$&SPACECRAFT POSITION Y\\
                    $r_z$&SPACECRAFT POSITION Z\\
                    $v_x$&SPACECRAFT VELOCITY X\\
                    $v_y$&SPACECRAFT VELOCITY Y\\
                    $v_z$&SPACECRAFT VELOCITY Z\\
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
                    $t_{OET}$&OBSERVED EVENT TIME\\
                    $f_{sky}$&SKY FREQUENCY\\
                    $f_{resid}$&RESIDUAL FREQUENCY\\
                    $P_{free}$&FREESPACE POWER\\
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
                    $\Delta\rho$&RADIUS CORRECTION\\
                    $\phi_{RL}$&RING LONGITUDE\\
                    $\phi_{ORA}$&OBSERVED RING AZIMUTH\\
                    $\tau$&NORMAL OPTICAL DEPTH\\
                    $\phi$&PHASE SHIFT\\
                    $\tau_{TH}$&NORMAL OPTICAL DEPTH THRESHOLD\\
                    $t_{OET}$&OBSERVED EVENT TIME\\
                    $t_{RET}$&RING EVENT TIME\\
                    $t_{SET}$&SPACECRAFT EVENT TIME\\
                    $B$&OBSERVED RING ELEVATION\\
                    \hline
                \end{tabular}
                \caption[Glossary of Parameters in Tau File]{
                    Glossary of optical depth, phase shift,
                    and selected geometry parameters
                    contained in files \theTAU.
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
                    trim={0.5in 1.5in 0.3in 1.5in},
                    clip,
                    width=\textwidth
                ]{\thePDF}
                \caption{Earth view of the occultation geometry
                         parameters in \theGEO.}
            \end{figure}
            \vspace{32ex}
            \begin{figure}[H]
                \centering
                \large{\textbf{View from North Pole}}\par
                \includegraphics[
                    page=3,
                    trim={0.5in 0.8in 0.3in 0.8in},
                    clip,
                    width=\textwidth
                ]{\thePDF}
                \caption{See caption of Figure 1a. $\oplus$}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \includegraphics[page=4, width=\textwidth]{\thePDF}
                \caption[Calibration Data from Cal File]{
                    Calibration data in file \theCAL.
                    The frequency residuals data
                    (the smooth curve, in the second panel)
                    is used to steer the carrier signal to the middle
                    of the recording bandwidth. The free-space power
                    data (the smooth curve in the third panel) is
                    used to normalize signal power measurements so that
                    the corresponding optical depth has nearly zero
                    value in the absence of rings. Least-square fitting
                    techniques to frequency and power estimates of
                    the direct signal (the green curves in the second
                    and third panels, respectively) are used to
                    compute the calibration data.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \includegraphics[page=5, width=\textwidth]{\thePDF}
                \caption[Ring Radius Correction from
                         Selected Occultation Geometry]{
                    Ring radius correction and
                    selected occultation geometry parameters
                    contained in the file \theTAU\
                    (solid green).
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=6, width=\textwidth]{\thePDF}
                }
                \caption[Normal Optical Depth Profiles 70000-85000km]{
                    Rev7-E normal optical depth profiles reconstructed
                    to remove diffraction effects at 1 km resolution
                    contained in the file
                    \theTAU
                    The 1 km resolution profile is plotted in green.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=7, width=\textwidth]{\thePDF}
                }
                \caption[Normal Optical Depth Profiles 85000-100000km]{
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=8, width=\textwidth]{\thePDF}
                }
                \caption[Normal Optical Depth Profiles 100000-115000km]{
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=9, width=\textwidth]{\thePDF}
                }
                \caption[Normal Optical Depth Profiles 115000-130000km]{
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution
                    (file \theTAU).
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=10, width=\textwidth]{\thePDF}
                }
                \caption[Normal Optical Depth Profiles 130000-145000km]{
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction
                    effects at 1 km resolution contained
                    in the file
                    \theTAU.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=11, width=\textwidth]{\thePDF}
                }
                \caption[Phase Shift Profile]{
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=12, width=\textwidth]{\thePDF}
                }
                \caption[Phase Shift Profile]{
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=13, width=\textwidth]{\thePDF}
                }
                \caption[Phase Shift Profile]{
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=14, width=\textwidth]{\thePDF}
                }
                \caption[Phase Shift Profile]{
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{
                    \includegraphics[page=15, width=\textwidth]{\thePDF}
                }
                \caption[Phase Shift Profile]{
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    \theTAU.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
        \end{document}
    """ % (pdffil, rev, doy, res, occ, geo, cal, tau, year, band)
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
        pdflatex %s.tex
        pdflatex %s.tex
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

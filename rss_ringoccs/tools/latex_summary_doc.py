def latex_summary_doc(pdfdir, pdffil, resolution, outfilename):

    var = 'Rev133E_RSS_2010_170_X43_E_Summary'


    LaTeXFile = r"""
        \documentclass{article}
        %---------------------------Packages----------------------------%
        \usepackage{geometry}
        \geometry{a4paper, margin = 1.0in}
        \usepackage{graphicx, float}            % Graphics/Images.
        \usepackage[english]{babel}             % Language typesetting.
        \usepackage[dvipsnames]{xcolor}         % Color names.
        \usepackage{listings, lstlinebgrd}      % Verbatim-Like Tools.
        \usepackage{mathtools, esint, mathrsfs} % amsmath and integrals.
        \usepackage{amsthm, amsfonts}           % Fonts and theorems.
        \usepackage{wrapfig}                    % Wrap text around figs.
        \usepackage{hyperref}                   % Allows for hyperlinks.
        \hypersetup{
            colorlinks=true,
            linkcolor=blue,
            filecolor=magenta,
            urlcolor=Cerulean,
            citecolor=SkyBlue
        }                                       % Colors for hyperref.
        \usepackage[
            font={normalsize},
            hypcap=true,
            labelsep=colon
        ]{caption}                              % Figure captions.
        \usepackage[%
            font=normalsize,
            labelformat=simple,
            labelsep=colon%
        ]{subcaption}                           % Subfigure captions.
        \usepackage[
            toc,
            acronym,
            nogroupskip
        ]{glossaries}                           %Glossaries and acronyms.
        %------------------------New Commands---------------------------%
        \renewcommand\labelitemii{$\circ$}
        \renewcommand\thesubfigure{%
            \arabic{section}.\arabic{figure}.\arabic{subfigure}%
        }
        \addto\captionsenglish{\renewcommand{\figurename}{Fig.}}
        %--------------------------LENGTHS------------------------------%
        % Indent and paragraph spacing.
        \setlength{\parindent}{0em}
        \setlength{\parskip}{0em}

        % CARL IS RSS\_2005\_123\_X43\_E
        % BOB IS THE PDF FILE (Or /path/to/bob.pdf)
        % Ex: Bob = /path/to/Rev133E_RSS_2010_170_X43_E_Summary.pdf

        \begin{document}
            \pagenumbering{gobble}
            \begin{center}
                \LARGE{%
                    \texttt{%
                        {RSS}\textunderscore{YEAR}%
                        \textunderscore{DOY}\textunderscore{BAND}%
                        \textunderscore{OCC}
                    }\\[2.0ex]
                    REVNUMBER Cassini Radio Science Ring
                    Occultation:\\[1.0ex]
                    Geometry, Data Calibration,
                    and Reconstructed\\[1.0ex]
                    Optical Depth and
                    Phase Shift Profiles\\[1.0ex]
                    at RESOLUTION\\[2.5ex]
                }
                \large{\today}
            \end{center}
            \vspace{2ex}
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=1,
                        trim={0.0cm 0.0cm 0.0cm 0.9cm},
                        clip,
                        width=\textwidth
                    ]{Bob}
                }
                \caption*{%
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
                    Symbol      & Parameter Name\\
                    \hline
                    $t_{OET}$ & OBSERVED EVENT TIME\\
                    $t_{RET}$ & RING EVENT TIME\\
                    $t_{SET}$ & SPACECRAFT EVENT TIME\\
                    $\rho$ & RING RADIUS\\
                    $\phi_{RL}$ & RING LONGITUDE\\
                    $\phi_{ORA}$ & OBSERVED RING AZIMUTH\\
                    $B$ & OBSERVED RING ELEVATION\\
                    $D$ & SPACECRAFT TO RING INTERCEPT DISTANCE\\
                    $\partial\rho/\partial t$
                    & RING INTERCEPT RADIAL VELOCITY\\
                    $\partial\theta/\partial t$
                    & RING INTERCEPT AZIMUTHAL VELOCITY\\
                    $F$ & FRESNEL SCALE\\
                    $R_{impact}$ & IMPACT RADIUS\\
                    $r_x$ & SPACECRAFT POSITION X\\
                    $r_y$ & SPACECRAFT POSITION Y\\
                    $r_z$ & SPACECRAFT POSITION Z\\
                    $v_x$ & SPACECRAFT VELOCITY X\\
                    $v_y$ & SPACECRAFT VELOCITY Y\\
                    $v_z$ & SPACECRAFT VELOCITY Z\\
                    \hline
                \end{tabular}
                \caption[Glossary of Parameters from the Geo File]{%
                    Glossary of parameters in file CARL\_GEO.TAB.
                    See companion label (.LBL) file for description
                    of parameters.
                }
                \label{tab:easydata_glossary_of_geo_file}
            \end{table}
            \begin{table}[H]
                \centering
                \begin{tabular}{l l}
                    \hline
                    Symbol & Parameter Name\\
                    \hline
                    $t_{OET}$ & OBSERVED EVENT TIME\\
                    $f_{sky}$ & SKY FREQUENCY\\
                    $f_{resid}$ & RESIDUAL FREQUENCY\\
                    $P_{free}$ & FREESPACE POWER\\
                    \hline
                \end{tabular}
                \caption[Glossary of Data from the Cal File]{%
                    Glossary of calibration data in file
                    CARL\_CAL.TAB. See companion label (.LBL)
                    file for description of the data.
                }
                \label{tab:easydata_glossary_from_cal_file}
            \end{table}
            \begin{table}[H]
                \centering
                \begin{tabular}{l l}
                    \hline
                    Symbol & Parameter Name\\
                    \hline
                    $\rho$ & RING RADIUS\\
                    $\Delta\rho$ & RADIUS CORRECTION\\
                    $\phi_{RL}$ & RING LONGITUDE\\
                    $\phi_{ORA}$ & OBSERVED RING AZIMUTH\\
                    $\tau$ & NORMAL OPTICAL DEPTH\\
                    $\phi$ & PHASE SHIFT\\
                    $\tau_{TH}$ & NORMAL OPTICAL DEPTH THRESHOLD\\
                    $t_{OET}$ & OBSERVED EVENT TIME\\
                    $t_{RET}$ & RING EVENT TIME\\
                    $t_{SET}$ & SPACECRAFT EVENT TIME\\
                    $B$ & OBSERVED RING ELEVATION\\
                    \hline
                \end{tabular}
                \caption[Glossary of Parameters in Tau File]{%
                    Glossary of optical depth, phase shift,
                    and selected geometry parameters
                    contained in files CARL\_TAU\_01KM.TAB
                    and RSS\_2005\_123\_X43\_E\_TAU\_10KM.TAB.
                    See companion label
                    (.LBL) files for description of the data.
                }
                \label{tab:easydata_parameters_from_tau_file}
            \end{table}
            \newpage
            \begin{figure}[H]
                \centering
                \large{\textbf{View from Earth}}\par
                \includegraphics[%
                    page=2,
                    trim={0.5in 1.5in 0.3in 1.5in},
                    clip,
                    width=\textwidth
                ]{Bob}
                \caption{%
                    Earth view of the occultation geometry
                    parameters in file CARL\_E\_GEO.tab
                }
            \end{figure}
            \vspace{32ex}
            \begin{figure}[H]
                \centering
                \large{\textbf{View from North Pole}}\par
                \includegraphics[%
                    page=3,
                    trim={0.5in 0.8in 0.3in 0.8in},
                    clip,
                    width=\textwidth
                ]{Bob}
                \caption{See caption of Figure 1a.}
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \includegraphics[%
                    page=4,
                    width=\textwidth
                ]{Bob}
                \caption[Calibration Data from Cal File]{%
                    Calibration data in file Rev007\_E\_X43\_CAL.TAB.
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
                \includegraphics[%
                    page=5,
                    width=\textwidth
                ]{Bob}
                \caption[%
                    Ring Radius Correction from
                    Selected Occultation Geometry
                ]{%
                    Ring radius correction and
                    selected occultation geometry parameters
                    contained in the file CARL\_TAU\_RESOLUTION.TAB
                    (solid green).
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=6,
                        width=\textwidth,
                    ]{Bob}
                }
                \caption[Normal Optical Depth Profiles 70000-85000km]{%
                    Rev7-E normal optical depth profiles reconstructed to
                    remove diffraction effects at 1 km resolution contained
                    in the file RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=7,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[%
                    Normal Optical Depth Profiles 85000-100000km%
                ]{%
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=8,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[%
                    Normal Optical Depth Profiles 100000-115000km
                ]{%
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=9,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[%
                    Normal Optical Depth Profiles 115000-130000km
                ]{%
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction effects
                    at 1 km resolution
                    (file RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab).
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=10,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[%
                    Normal Optical Depth Profiles 130000-145000km
                ]{%
                    Rev7-E normal optical depth profiles
                    reconstructed to remove diffraction
                    effects at 1 km resolution contained
                    in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab.
                    The 1 km resolution profile is plotted in green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=11,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[Phase Shift Profile]{%
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=12,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[Phase Shift Profile]{%
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=13,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[Phase Shift Profile]{%
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=14,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[Phase Shift Profile]{%
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
            \newpage
            \begin{figure}[H]
                \centering
                \resizebox{\textwidth}{!}{%
                    \includegraphics[%
                        page=15,
                        width=\textwidth
                    ]{Bob}
                }
                \caption[Phase Shift Profile]{%
                    Rev7-E Phase shift profile reconstructed
                    to remove diffraction effects at 1 km
                    resolution contained in the file
                    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB.
                    The 1 km resolution profile is plotted
                    in solid green.
                }
            \end{figure}
        \end{document}
    """ % (pdfdir, pdffil, pdffil, pdffil, pdffil,
           pdffil, pdffil, pdffil, pdffil, pdffil,
           pdffil, pdffil)
    TexName = outfilename
    TexFileName = '%s.tex' % TexName
    TexFile = open(TexFileName,'w')
    TexFile.write(LaTeXFile)
    TexFile.close()    

    TexScript = """
        #!/bin/bash
        pdflatex %s.tex
        pdflatex %s.tex
        rm -f %s.log
        rm -f %s.aux
        rm -f %s.tex
        rm -f %s.sh
    """ % (TexName,TexName,TexName,TexName,TexName,TexName)

    TexSh   = open("%s.sh" % TexName,'w')
    TexSh.write(TexScript)
    TexSh.close()

    make_executable('%s.sh' % TexName)
    shell_execute(['./%s.sh' % TexName])


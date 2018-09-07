def create_summary_doc(pdfdir,pdffil,outfilename):
    LaTeXFile       = r'''
    \documentclass{article}
    \usepackage{geometry}
    \geometry{a4paper, margin = 1.0in}
    \usepackage[T1]{fontenc}
    \usepackage{graphicx,float}
    \graphicspath{{%s}}
    \usepackage[dvipsnames]{xcolor}
    \usepackage{mathtools,esint,mathrsfs}
    \usepackage{amsthm,amsfonts,upgreek}
    \usepackage{wrapfig}
    \usepackage[font=scriptsize]{subcaption}
    \usepackage[font={scriptsize,it}]{caption}
    \usepackage{hyperref}
    \hypersetup{
    colorlinks=true,linkcolor=blue,filecolor=magenta,
    urlcolor=Cerulean,citecolor=SkyBlue}                           
    \title{test}
    \author{Ryan Maguire}
    \date{June 2018}
    \begin{document}
    \begin{center}
        \LARGE{RSS\_2005\_123\_X43\_E \par
        Rev7-E Cassini Radio Science Ring Occultation: Geometry, Data 
Calibration, and Reconstructed Optical Depth and Phase Shift Profiles at 
1 and 10km Resolution \par
        February 9, 2018\par}
    \end{center}
    \begin{figure}[H]
        \centering
        \includegraphics[page=1,trim = {0.67in 0.5in 0.5in 
3.1in},clip,width=\textwidth]{%s}
        \caption[Radio Occultation Track]{The radio occultation track as 
seen looking down on the ring plane. The solid red track is relative to 
a reference direction defined by the direction to Earth, The dashed blue 
track is relative to a reference direction defined by the ascending node 
of J2000 on the ring plane.}
    \end{figure}
        \begin{table}
        \centering
        \begin{tabular}{l l}
            \hline
            Symbol                      & Parameter Name                        
\\
            \hline
            $t_{OET}$                   & OBSERVED EVENT TIME                   
\\
            $t_{RET}$                   & RING EVENT TIME                       
\\
            $t_{SET}$                   & SPACECRAFT EVENT TIME                 
\\
            $\rho$                      & RING RADIUS                           
\\
            $\phi_{RL}$                 & RING LONGITUDE                        
\\
            $\phi_{ORA}$                & OBSERVED RING AZIMUTH                 
\\
            $B$                         & OBSERVED RING ELEVATION               
\\
            $D$                         & SPACECRAFT TO RING INTERCEPT 
DISTANCE \\
            $\partial\rho/\partial t$   & RING INTERCEPT RADIAL VELOCITY        
\\
            $\partial\theta/\partial t$ & RING INTERCEPT AZIMUTHAL 
VELOCITY     \\
            $F$                         & FRESNEL SCALE                         
\\
            $R_{impact}$                & IMPACT RADIUS                         
\\
            $r_x$                       & SPACECRAFT POSITION X                 
\\
            $r_y$                       & SPACECRAFT POSITION Y                 
\\
            $r_z$                       & SPACECRAFT POSITIION Z                
\\
            $v_x$                       & SPACECRAFT VELOCITY X                 
\\
            $v_y$                       & SPACECRAFT VELOCITY Y                 
\\
            $v_z$                       & SPACECRAFT VELOCITY Z                 
\\
            \hline
        \end{tabular}
        \caption[Glossary of Parameters from the Geo File]{Glossary of 
parameters in file
        RSS\_2005\_123\_X43\_E\_GEO.TAB. See companion label (.LBL) file 
for description 
        of parameters.}
        \label{tab:easydata_glossary_of_geo_file}
    \end{table}
    \begin{table}
        \centering
        \begin{tabular}{l l}
            \hline
            Symbol      & Parameter Name        \\
            \hline
            $t_{OET}$   & OBSERVED EVENT TIME   \\
            $f_{sky}$   & SKY FREQUENCY         \\
            $f_{resid}$ & RESIDUAL FREQUENCY    \\
            $P_{free}$  & FREESPACE POWER       \\
            \hline
        \end{tabular}
        \caption[Glossary of Data from the Cal File]{Glossary of 
calibration data in file
        RSS\_2005\_123\_X43\_E\_CAL.TAB. See companion label (.LBL) file 
for description of the data.}
        \label{tab:easydata_glossary_from_cal_file}
    \end{table}
    \begin{table}
        \centering
        \begin{tabular}{l l}
            \hline
            Symbol          & Parameter Name                    \\
            \hline
            $\rho$          & RING RADIUS                       \\
            $\Delta\rho$    & RADIUS CORRECTION                 \\
            $\phi_{RL}$     & RING LONGITUDE                    \\
            $\phi_{ORA}$    & OBSERVED RING AZIMUTH             \\
            $\tau$          & NORMAL OPTICAL DEPTH              \\
            $\phi$          & PHASE SHIFT                       \\
            $\tau_{TH}$     & NORMAL OPTICAL DEPTH THRESHOLD    \\
            $t_{OET}$       & OBSERVED EVENT TIME               \\
            $t_{RET}$       & RING EVENT TIME                   \\
            $t_{SET}$       & SPACECRAFT EVENT TIME             \\
            $B$             & OBSERVED RING ELEVATION           \\
            \hline
        \end{tabular}
        \caption[Glossary of Parameters in Tau File]{Glossary of
        optical depth, phase shift, and selected geometry parameters
        contained in files RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB
        and RSS\_2005\_123\_X43\_E\_TAU\_10KM.TAB. See companion label
        (.LBL) files for description of the data.}
        \label{tab:easydata_parameters_from_tau_file}
    \end{table}
    \clearpage
    \begin{table}[H]
        \centering
        \begin{tabular}{l l} 
            \hline
            Symbol 			& Parameter Name \\
            \hline
            $t_{OET}$		& OBSERVED EVENT TIME \\ 
            $t_{RET}$ 		& RING EVENT TIME \\
            $t_{SET}$ 		& SPACECRAFT EVENT TIME \\
            $\rho$	 		& RING RADIUS \\
            $\phi_{RL}$		& RING LONGITUDE \\
            $\phi_{ORA}$		& OBSERVED RING AZIMUTH \\
            $B$				& OBSERVED RING ELEVATION \\
            $D$				& SPACECRAFT TO RING INTERCEPT 
DISTANCE \\
            $\partial\rho/\partial t$	& RING INTERCEPT RADIAL VELOCITY 
\\
            $\partial\theta/\partial t$	& RING INTERCEPT AZIMUTHAL 
VELOCITY \\
            $F$				& FRESNEL SCALE \\
            $R_{impact}$		& IMPACT RADIUS \\
            $r_x$			& SPACECRAFT POSITION X \\
            $r_y$			& SPACECRAFT POSITION Y \\
            $r_z$			& SPACECRAFT POSITIION Z \\
            $v_x$			& SPACECRAFT VELOCITY X \\
            $v_y$			& SPACECRAFT VELOCITY Y \\
            $v_z$			& SPACECRAFT VELOCITY Z \\
            \hline
        \end{tabular}
        \caption[Glossary of Parameters from the Geo File]{Glossary of 
parameters in file RSS\_2005\_123\_X43\_E\_GEO.TAB. See companion label 
(.LBL) file for description of parameters.}
        \label{tab:easydata_glossary_of_geo_file}
    \end{table}
    \begin{table}[H]
        \centering
        \begin{tabular}{l l}
            \hline
            Symbol			& Parameter Name \\
            \hline
            $t_{OET}$			& OBSERVED EVENT TIME \\
            $f_{sky}$			& SKY FREQUENCY \\
            $f_{resid}$		& RESIDUAL FREQUENCY \\
            $P_{free}$		& FREESPACE POWER \\
            \hline
        \end{tabular}
        \caption[Glossary of Data from the Cal File]{Glossary of 
calibration data in file RSS\_2005\_123\_X43\_E\_CAL.TAB. See companion 
label (.LBL) file for description of the data.}
        \label{tab:easydata_glossary_from_cal_file}
    \end{table}
    \begin{table}[H]
        \centering
        \begin{tabular}{l l}
            \hline
            Symbol			& Parameter Name \\
            \hline
            $\rho$			& RING RADIUS \\
            $\Delta\rho$		& RADIUS CORRECTION \\
            $\phi_{RL}$		& RING LONGITUDE \\
            $\phi_{ORA}$		& OBSERVED RING AZIMUTH \\
            $\tau$			& NORMAL OPTICAL DEPTH \\
            $\phi$			& PHASE SHIFT \\
            $\tau_{TH}$		& NORMAL OPTICAL DEPTH THRESHOLD \\
            $t_{OET}$			& OBSERVED EVENT TIME \\
            $t_{RET}$			& RING EVENT TIME \\
            $t_{SET}$			& SPACECRAFT EVENT TIME \\
            $B$				& OBSERVED RING ELEVATION \\
            \hline
        \end{tabular}
        \caption[Glossary of Parameters in Tau File]{Glossary of optical 
depth, phase shift, and selected geometry parameters contained in files 
RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB and 
RSS\_2005\_123\_X43\_E\_TAU\_10KM.TAB. See companion label (.LBL) files 
for description of the data.}
        \label{tab:easydata_parameters_from_tau_file}
    \end{table}
    \begin{figure}[H]
        \centering
        \includegraphics[page=2,trim = {0.8in 0.5in 0.21in 
0.45in},clip,width=\textwidth]{%s}
        \caption{Occultation geometry parameters in file 
RSS\_2005\_123\_X43\_E\_GEO.tab}
    \end{figure}
    \begin{figure}[H]
        \centering
        \includegraphics[page=3,trim = {0.8in 0.5in 0.21in 
0.45in},clip,width=\textwidth]{%s}
        \caption{See caption of Figure 1a.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \includegraphics[page=4,trim = {0.67in 0.5in 0.45in 
0.5in},clip,width=\textwidth]{%s}
        \caption[Calibration Data from Cal File]{Calibration data in 
file Rev007\_E\_X43\_CAL.TAB. The frequency residuals data (the smooth 
curve, in the second panel) is used to steer the carrier signal to the 
middle of the recording bandwidth. The free-space power data (the smooth 
curve in the third panel) is used to normalize signal power measurements 
so that the corresponding optical depth has nearly zero value in the 
absence of rings. Least-square fitting techniques to frequency and power 
estimates of the direct signal (the green curves in the second and third 
panels, respectively) are used to compute the calibration data.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \includegraphics[page=5,trim = {0.8in 0.5in 0.21in 
0.45in},clip,width=\textwidth]{%s}
        \caption[Ring Radius Correction from Selected Occultation 
Geometry]{Ring radius correction and selected occultation geometry 
parameters contained in the file
    RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB (solid green).}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=1,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Normal Optical Depth Profiles 70000-85000km]{Rev7-E 
normal optical depth profiles reconstructed to remove diffraction 
effects at 1 km resolution contained in the file 
RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab. The 1 km resolution profile is 
plotted in green.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=2,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Normal Optical Depth Profiles 85000-100000km]{Rev7-E 
normal optical depth profiles reconstructed to remove diffraction 
effects at 1 km resolution contained in the file 
RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab. The 1 km resolution profile is 
plotted in green.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=3,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Normal Optical Depth Profiles 100000-115000km]{Rev7-E 
normal optical depth profiles reconstructed to remove diffraction 
effects at 1 km resolution contained in the file 
RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab. The 1 km resolution profile is 
plotted in green.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=4,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Normal Optical Depth Profiles 115000-130000km]{Rev7-E 
normal optical depth profiles reconstructed to remove diffraction 
effects at 1 km resolution (file RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab). 
The 1 km resolution profile is plotted in green.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=5,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Normal Optical Depth Profiles 130000-145000km]{Rev7-E 
normal optical depth profiles reconstructed to remove diffraction 
effects at 1 km resolution contained in the file 
RSS\_2005\_123\_X43\_E\_TAU\_01KM.tab. The 1 km resolution profile is 
plotted in green.}
    \end{figure}
    \begin{figure}[H]
        \centering
        \resizebox{\textwidth}{0.4\textheight}{
        \includegraphics[page=6,width=\textwidth,trim = {1.1in 0.85in 
0.8in 0.0in},clip]{%s}}
        \caption[Phase Shift Profile]{Rev7-E Phase shift profile 
reconstructed to remove diffraction effects at 1 km resolution contained 
in the file RSS\_2005\_123\_X43\_E\_TAU\_01KM.TAB. The 1 km resolution 
profile is plotted in solid green.}
    \end{figure}
    \end{document}
    ''' % 
(pdfdir,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil,pdffil)
    TexName         = outfilename
    TexFileName     = '%s.tex' % TexName
    TexFile         = open(TexFileName,'w')
    TexFile.write(LaTeXFile)
    TexFile.close()    

    TexScript = '''#!/bin/bash
    pdflatex %s.tex
    pdflatex %s.tex
    rm %s.log
    rm %s.aux
    rm %s.tex
    rm %s.sh''' % (TexName,TexName,TexName,TexName,TexName,TexName)

    TexSh   = open("%s.sh" % TexName,'w')
    TexSh.write(TexScript)
    TexSh.close()

    make_executable('%s.sh' % TexName)
    shell_execute(['./%s.sh' % TexName])


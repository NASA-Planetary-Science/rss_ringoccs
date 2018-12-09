#! /bin/bash

for dir in $( ls ../../output/ ); do

    if [ -d ../output/${dir} ]
    then

        for subdir in $( ls ../output/$dir ); do

            echo "Merging calibration FORFIT PDFs in ../../output/${dir}/${subdir}/"
            "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ../output/${dir}/${subdir}_FORFITS.PDF ../output/${dir}/${subdir}/*/*FORFIT*.PDF
            echo "Merging calibration FORFIT PDFs in ../../output/${dir}/${subdir}/"
            "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ../output/${dir}/${subdir}_FSPFITS.PDF ../output/${dir}/${subdir}/*/*FSPFIT*.PDF

        done

    fi

done

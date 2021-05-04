#! /bin/bash

echo -e "\n\nGeo\n\n"
python geocal.py

echo -e "\n\nDLP\n\n"
python dlp.py

echo -e "\n\nReconstruction\n\n"
python reconstruction.py

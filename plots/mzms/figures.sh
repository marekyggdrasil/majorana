#!/bin/sh
convert -background none psi1_45.png psi1_135.png +append figure_psi1.png
convert -background none psi2_45.png psi2_135.png +append figure_psi2.png
convert -background none psi3_45.png psi3_135.png +append figure_psi3.png
convert figure_psi1.png figure_psi2.png figure_psi3.png -append figure_mzms.png
rm figure_p*

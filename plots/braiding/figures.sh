#!/bin/sh
convert -background none psi0L_45.png psi0L_135.png +append figure_0L.png
convert -background none psi1L_45.png psi1L_135.png +append figure_1L.png
rm figure_p*

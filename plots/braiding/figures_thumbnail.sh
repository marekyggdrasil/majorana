#!/bin/sh
convert -background none psi0L_45.png psi0L_135.png +append figure_psi0L.png
convert -background none psi1L_45.png psi1L_135.png +append figure_psi1L.png
convert -background none figure_psi0L.png blog_mzms_pairings_logical_zero.png -gravity center -append figure_mzms_logical_zero.png
convert -background none figure_psi1L.png blog_mzms_pairings_logical_one.png -gravity center -append figure_mzms_logical_one.png
convert -background none figure_mzms_logical_zero.png figure_mzms_logical_one.png +append figure_thumbnail.png
rm figure_p*
rm figure_mzms*

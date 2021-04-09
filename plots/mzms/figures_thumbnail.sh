#!/bin/sh
convert -background none psi1_45.png psi1_135.png +append figure_psi1.png
convert -background none psi3_45.png psi3_135.png +append figure_psi3.png
convert -background none figure_psi1.png blog_mzms_pairings_trivial.png -gravity center -append figure_mzms_trivial.png
convert -background none figure_psi3.png blog_mzms_pairings_topological.png -gravity center -append figure_mzms_topological.png
convert -background none figure_mzms_trivial.png figure_mzms_topological.png +append figure_thumbnail.png
rm figure_p*
rm figure_mzms*

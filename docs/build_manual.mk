# Build the CompLB3D User Guide.
#
# Saved as build_manual.mk rather than Makefile because some sync tools refuse
# to write a bare "Makefile". Rename it, or run:  make -f build_manual.mk
#
#   make          -> CompLB3D_User_Guide.pdf
#   make clean    -> remove intermediates
#
# Needs a LaTeX distribution with tikz, listings, tcolorbox, booktabs and
# lmodern. On Debian or Ubuntu:
#   apt install texlive-latex-recommended texlive-latex-extra \
#               texlive-fonts-recommended lmodern
#
# Two passes are required: the first writes the table of contents, the second
# resolves the cross-references against it.

TARGET = CompLB3D_User_Guide

$(TARGET).pdf: main.tex part1.tex part2.tex part3.tex part4.tex
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex
	mv main.pdf $(TARGET).pdf

clean:
	rm -f main.aux main.log main.out main.toc

.PHONY: clean

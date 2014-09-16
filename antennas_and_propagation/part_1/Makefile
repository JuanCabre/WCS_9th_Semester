##################################################################
# Makefile for LaTeX
##################################################################
# Use:
# make
# make clean
# options for ps2pdf: http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm

TEX:=$(shell ls *.tex)
OTHER = *~ *.acr *.acn *.alg *.aux *.gls *.glo *.glg *.toc *.tdo *.bbl *.blg *.gz *.out *.thm *.idx *.ilg *.ind *.xdy

latex: master.tex
	latex --synctex=1 master.tex
	#bibtex master
	#makeglossaries master
	latex --synctex=1 master.tex
	latex --synctex=1 master.tex
	dvips master.dvi
	ps2pdf master.ps
	rm -f $(OTHER) $(PS)
clean:
	rm -f $(OTHER) $(PS)

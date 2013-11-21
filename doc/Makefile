VB: vbSPT_main.tex vbSPT_derivation.tex  vbSPT_statistics.tex  vbSPT_userManual.tex vbSPT_tables.tex
	pdflatex vbSPT_main
	bibtex vbSPT_main
	pdflatex vbSPT_main
	pdflatex vbSPT_main
	mv vbSPT_main.pdf ../vbSPT_userguide.pdf
HTML:    vbSPT_main.tex vbSPT_derivation.tex  vbSPT_statistics.tex  vbSPT_userManual.tex vbSPT_tables.tex
	rm -rf ../html/* *aux *bbl *blg
	htlatex vbSPT_main '' '' -d../html/




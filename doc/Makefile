VB: vbSPT_main.tex vbSPT_derivation.tex  vbSPT_statistics.tex  vbSPT_userManual.tex vbSPT_tables.tex
	pdflatex vbSPT_main
	bibtex vbSPT_main
	pdflatex vbSPT_main
	pdflatex vbSPT_main
	mv vbSPT_main.pdf ../vbSPT_userguide.pdf
html:    vbSPT_main.tex vbSPT_derivation.tex  vbSPT_statistics.tex  vbSPT_userManual.tex vbSPT_tables.tex
	rm -rf ../html/* 
	htlatex vbSPT_main '' '' -d../html/
clear: 
	rm *aux *bbl *blg *html *css *tmp *4ct *4tc *dvi *idv *lg *log *out vbSPT_main.pdf





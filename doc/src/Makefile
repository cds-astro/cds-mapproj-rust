name=proj

all : pdf view clean
update: pdfupdate clean

pdfupdate:
	pdflatex ${name}
pdf :
	pdflatex ${name}
	bibtex   ${name}
	pdflatex ${name}
	pdflatex ${name}
	@mv ${name}.pdf ..

clean :
	@rm *.aux *.log *.out *.bbl *.blg *.toc

view:
	evince ../${name}.pdf


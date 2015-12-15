FILE=programming

default:
	- latex $(FILE)
	- bibtex $(FILE)
	- latex $(FILE)
	- latex $(FILE)
	- dvips -o $(FILE).ps  $(FILE)

pdf:
	- latex $(FILE)
#	- dvips -Ppdf -updftex.map -G0 -tA4 -o $(FILE).ps  $(FILE)
	- dvips -Ppdf -G0 -tA4 -o $(FILE).ps  $(FILE)
	- ps2pdf $(FILE).ps $(FILE).pdf
	- evince $(FILE).pdf &

clean:
	- rm -f *.aux */*.aux
	- rm -f *.log */*.log
	- rm -f *.blg */*.blg
	- rm -f *.bbl */*.bbl
	- rm -f *.bak */*.bak
	- rm -f $(FILE).dvi $(FILE).ps $(FILE).pdf

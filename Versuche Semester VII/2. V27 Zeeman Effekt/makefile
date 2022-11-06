all: build/Protokoll_V27.pdf



#hier Python-Skripte


build/plots/dreieck.pdf:	|build
	python python/V27.py



# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/Protokoll_V27.pdf:  latex/auswertung.tex latex/theorie.tex latex/durchf.tex latex/diskussion.tex build/plots/dreieck.pdf| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V27.tex
	biber build/Protokoll_V27.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V27.tex


build:
	mkdir -p build

clean:
	rm -rf build


#Zusätzliche Ordner/Dateien löschen
aux2:
	rm -rf Protokoll_V27.aux Protokoll_V27.bcf Protokoll_V27.log Protokoll_V27.pdf Protokoll_V27.run.xml Protokoll_V27.toc

auxv:
	rm -rf python/variables

auxt:
	rm -rf build/Tabellen	
FORCE:

.PHONY: all clean
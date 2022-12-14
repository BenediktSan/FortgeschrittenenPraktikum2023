all: build/Protokoll_V46.pdf



#hier Python-Skripte


build/plots/probe1.pdf: python/V46.py	header-matplotlib.tex  matplotlibrc |build
	TEXINPUTS=$$(pwd): python python/V46.py



# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/Protokoll_V46.pdf:  latex/auswertung.tex latex/theorie.tex latex/durchf.tex latex/diskussion.tex build/plots/probe1.pdf| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V46.tex
	biber build/Protokoll_V46.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V46.tex


build:
	mkdir -p build

clean:
	rm -rf build


#Zusätzliche Ordner/Dateien löschen
aux2:
	rm -rf Protokoll_V46.aux Protokoll_V46.bcf Protokoll_V46.log Protokoll_V46.pdf Protokoll_V46.run.xml Protokoll_V46.toc

auxv:
	rm -rf python/variables

auxt:
	rm -rf build/Tabellen	
FORCE:

.PHONY: all clean
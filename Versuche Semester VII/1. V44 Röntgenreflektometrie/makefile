all: build/Protokoll_V44.pdf



#hier Python-Skripte


build/plots/plot1.pdf:	|build
	python python/V44.py



# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/Protokoll_V44.pdf:  latex/auswertung.tex latex/theorie.tex latex/durchf.tex latex/diskussion.tex build/plots/plot1.pdf| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V44.tex
	biber build/Protokoll_V44.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V44.tex


build:
	mkdir -p build

clean:
	rm -rf build


#Zusätzliche Ordner/Dateien löschen
aux2:
	rm -rf Protokoll_V44.aux Protokoll_V44.bcf Protokoll_V44.log Protokoll_V44.pdf Protokoll_V44.run.xml Protokoll_V44.toc

auxv:
	rm -rf python/variables

auxt:
	rm -rf build/Tabellen	
FORCE:

.PHONY: all clean
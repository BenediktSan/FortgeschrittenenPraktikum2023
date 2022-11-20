all: build/Protokoll_V49.pdf



#hier Python-Skripte


build/plots/plot1.pdf:	|build
	python python/V49.py



# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/Protokoll_V49.pdf:  latex/auswertung.tex latex/theorie.tex latex/durchf.tex latex/diskussion.tex build/plots/plot1.pdf| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V49.tex
	biber build/Protokoll_V49.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V49.tex


build:
	mkdir -p build

clean:
	rm -rf build


#Zusätzliche Ordner/Dateien löschen
aux2:
	rm -rf Protokoll_V49.aux Protokoll_V49.bcf Protokoll_V49.log Protokoll_V49.pdf Protokoll_V49.run.xml Protokoll_V49.toc

auxv:
	rm -rf python/variables

auxt:
	rm -rf build/Tabellen	
FORCE:

.PHONY: all clean
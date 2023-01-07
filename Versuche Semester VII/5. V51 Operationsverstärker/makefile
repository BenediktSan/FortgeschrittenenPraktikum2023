all: build/Protokoll_V51.pdf



#hier Python-Skripte


build/plots/plot1.pdf:	|build
	python python/V51.py



# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/Protokoll_V51.pdf:  latex/auswertung.tex latex/theorie.tex latex/durchf.tex latex/diskussion.tex build/plots/plot1.pdf| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V51.tex
	biber build/Protokoll_V51.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error latex/Protokoll_V51.tex


build:
	mkdir -p build

clean:
	rm -rf build


#Zusätzliche Ordner/Dateien löschen
aux2:
	rm -rf Protokoll_V51.aux Protokoll_V51.bcf Protokoll_V51.log Protokoll_V51.pdf Protokoll_V51.run.xml Protokoll_V51.toc

auxv:
	rm -rf python/variables

auxt:
	rm -rf build/Tabellen	
FORCE:

.PHONY: all clean
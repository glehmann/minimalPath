
# you should use the title of the article instead
ARTICLE_TITLE = Minimal Path Image Filter

# spaces are not allowed
empty:=
space:= $(empty) $(empty)
PDF_NAME:= $(subst $(space),_,$(ARTICLE_TITLE))

NAME = Article

BIBS=$(wildcard *.bib)

PNGS=$(notdir $(wildcard ../ImageResults/*.png))
JPGS=$(notdir $(wildcard ../ImageResults/*.jpg))
DIAS=$(notdir $(wildcard ../ImageResults/*.dia))
FIGS=$(notdir $(wildcard ../ImageResults/*.fig))
EPSS=$(notdir $(wildcard ../ImageResults/*.eps))
TIFS=$(notdir $(wildcard ../ImageResults/*.tif))
PLOTS=$(notdir $(wildcard ../ImageResults/*.gnuplot))

OUTPUT_EPS=$(PNGS:png=eps) $(JPGS:jpg=eps) $(DIAS:dia=eps) $(FIGS:fig=eps) $(TIFS:tif=eps) $(PLOTS:gnuplot=eps) $(EPSS)

all: pdf


ps: ${NAME}.ps
pdf: ${PDF_NAME}.pdf
dvi: ${NAME}.dvi
bbl: ${NAME}.bbl
eps: $(OUTPUT_EPS)


%.eps: ../ImageResults/%.png
	convert $< $@

%.eps: ../ImageResults/%.jpg
	convert $< $@

%.eps: ../ImageResults/%.dia
	dia --nosplash -e $@ $<

%.eps: ../ImageResults/%.fig
	fig2dev -L eps $< $@

%.eps: ../ImageResults/%.gnuplot
	gnuplot $< > $@


${NAME}.ps: ${NAME}.dvi
	dvips -R0 -o ${NAME}.ps ${NAME}.dvi

%.eps: ../ImageResults/%.tif
	convert $< $@


${PDF_NAME}.pdf: ${NAME}.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.3 \
         -dSubsetFonts=true -dEmbedAllFonts=true \
         -dAutoFilterColorImages=false \
         -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode \
	 -dMonoImageFilter=/FlateEncode \
	 ${NAME}.ps ${PDF_NAME}.pdf


${NAME}.bbl: ${BIBS} ${NAME}.aux
	bibtex ${NAME}


${NAME}.dvi: ${NAME}.tex  ${BIBS}
	$(MAKE) eps
	latex -interaction=nonstopmode ${NAME}.tex
	bibtex ${NAME}
	latex -interaction=nonstopmode ${NAME}.tex
	latex -interaction=nonstopmode ${NAME}.tex


${NAME}.aux: ${NAME}.tex  ${BIBS}
	latex ${NAME}.tex


clean:
	rm -f ${NAME}.dvi ${NAME}.ps ${PDF_NAME}.pdf ${NAME}.aux ${NAME}.log ${NAME}.bbl ${NAME}.blg ${NAME}.toc ${NAME}.out  ${NAME}.brf $(DIAS:dia=eps) $(PNGS:png=eps) $(FIGS:fig=eps)


INDEX=$(shell ls */README.md | sed 's/README.md/index.html/')
TARGET=$(INDEX)
PANDOCOPT=--mathjax -t html --template=template
PANDOC_TEXOPT=--highlight-style tango --latex-engine=lualatex -V documentclass=ltjarticle -V geometry:margin=1in
README=$(shell ls */README.md)
PDF=$(README:.md=.pdf)

all: $(TARGET) index.html

pdf: $(PDF)

index.md: README.md
	sed 's/README.md/index.html/' $< > $@

index.html: index.md
	pandoc -s $< -o $@ --mathjax -t html --template=template
	rm -f index.md

%/index.md: %/README.md
	sed '2a [[Up]](../index.html)' $< > $@
	sed -i '3a [[Repository]](https://github.com/kaityo256/md2019)\n' $@

%/index.html: %/index.md
	pandoc -s $< -o $@ $(PANDOCOPT)

%/README.pdf: %/README.md
	cd $(dir $@);pandoc $(notdir $<) -s -o $(notdir $@) $(PANDOC_TEXOPT)


clean:
	rm -f $(TARGET) index.html

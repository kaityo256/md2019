INDEX=$(shell ls */README.md | sed 's/README.md/index.html/')
TARGET=$(INDEX)
PANDOCOPT=--mathjax -t html --template=template
all: $(TARGET) index.html

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


clean:
	rm -f $(TARGET) index.html

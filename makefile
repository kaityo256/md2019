all: index.html

index.html: README.md
	pandoc -s $< -o $@ --mathjax -t html --template=template

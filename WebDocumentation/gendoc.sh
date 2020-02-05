#!/bin/sh
/opt/local/bin/pandoc --toc --toc-depth=2 --shift-heading-level-by=-1 --css=igmdoc.css --mathjax --standalone --to=html5 -fmarkdown-implicit_figures --metadata-file=igraphm.yaml --lua-filter=local_toc.lua --syntax-definition=mathematica.xml --highlight-style=mathematica.theme igraphm.md -o index.html

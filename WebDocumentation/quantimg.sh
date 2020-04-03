#!/bin/sh

/opt/local/bin/pngquant --quality=50-80 --ext .png -f --skip-if-larger --strip img/*.png

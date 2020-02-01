#!/bin/sh

find img -iname '*.png' -print0 | xargs -0 -n 1 -P4 /opt/local/bin/optipng -o4

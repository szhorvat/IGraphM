#!/bin/sh
/opt/local/bin/mogrify -filter Lanczos -resize 50% -unsharp 0x0.5+0.4 img/*.png

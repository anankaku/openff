#!/bin/bash

for f in *.com; do
    cp "$f" "$f.bak"
    sed -i '34c\
D 2 4 5 6 F\
D 4 5 6 8 S 37 10.000000' "$f"
done
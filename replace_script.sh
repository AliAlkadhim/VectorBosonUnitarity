#!/bin/sh

sed -i 's/\[/(/g' m4
sed -i 's/\]/)/g' m4

sed -i 's/\^/\**/g' m4

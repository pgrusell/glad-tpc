#!/bin/bash
find  -type f | grep -e '.(\.cpp\|\.cxx\|\.h)$'  | grep -v ^./macros | grep -v ^./params | grep -v ^./fitter | xargs -L 1 clang-format-15 -i 
echo "Use git add -A ; git commit -m \"clang-format all files\" --author=\"white space <whitespace@example.com>\" to commit changes."

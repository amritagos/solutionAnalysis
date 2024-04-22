#!/bin/sh
find ./subprojects/solvlib -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec clang-format -i {} \;
find ./python_bindings -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec clang-format -i {} \;
find ./soluanalysis -regex '.*py' -exec black {} \;

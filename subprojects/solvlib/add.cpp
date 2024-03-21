#include "add.hpp"
#include <fmt/core.h>
#include <fmt/format.h>

int my_add(int i, int j) {

    fmt::print( "{} + {} = {}\n",i, j, i+j );
    return i + j;
}
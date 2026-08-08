#pragma once

#include <cstring>


static float reference_from_matrix(const char* matname)
{
    if (!strcmp(matname, "t99") || !strcmp(matname, "t04")) {
        return 559.0;
    }
    if (!strcmp(matname, "su95") || !strcmp(matname, "su95_noGU")) {
        return 249.0;
    }
    return 363.0;
}

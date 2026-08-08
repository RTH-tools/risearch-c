#pragma once

#include <string.h>

/* can save first traversal as length is known before! */
static int reverse_inplace(char *str_beg, int j)
{
    char tmp;
    char *str_end = &str_beg[j];

    while (str_end > str_beg) {
        tmp = *str_beg;
        *str_beg++ = *str_end;
        *str_end-- = tmp;
    }

    return 0;
}

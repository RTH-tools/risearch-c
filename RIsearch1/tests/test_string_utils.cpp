#include <gtest/gtest.h>

#include "string_util.h"


TEST(StringUtils, TestReverseSanity)
{
    char str[] = "abcde";
    reverse_inplace(str, 2);
    EXPECT_STREQ(str, "cbade");
    reverse_inplace(str, 2);
    EXPECT_STREQ(str, "abcde");
}


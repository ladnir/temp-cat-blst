#include "gtest/gtest.h"

#include "runtime_constants_validation_support.h"

TEST(AvxRuntimeConstants, Validate52BitBakedConstantsAgainstGmp)
{
    EXPECT_NO_THROW(validate_runtime_constants_52());
}

TEST(AvxRuntimeConstants, Validate50BitBakedConstantsAgainstGmp)
{
    EXPECT_NO_THROW(validate_runtime_constants_50());
}

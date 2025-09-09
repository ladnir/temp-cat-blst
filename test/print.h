

extern "C" {
    #include "src/fields.h"
}


void print_vec384(const vec384 v) {
    for (int i = 0; i < 6; ++i)
        printf("%016llx%s", (unsigned long long)v[i], i < 5 ? " " : "\n");
}
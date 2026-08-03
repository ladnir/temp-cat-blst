/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#if defined(_WIN64) && (defined(_M_ARM64) || defined(_M_ARM64EC) || defined(BLST_DLL_MAIN))
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#endif

#if defined(_MSC_VER) && (defined(_M_X64) && !defined(_M_ARM64EC))
# include <intrin.h>
#endif

extern "C"
{
#if (defined(__GNUC__) || defined(__clang__) || defined(__SUNPRO_C)) && !defined(_WIN32)
__attribute__((visibility("hidden")))
#endif
int blst_platform_cap = 0;
}

#if defined(__x86_64__) || defined(__x86_64) || (defined(_M_X64) && !defined(_M_ARM64EC))

static void blst_cpuid_query(int info[4], int func, int sub)
{
#if (defined(__GNUC__) || defined(__clang__) || defined(__SUNPRO_C)) && !defined(_MSC_VER)
    int eax, ebx, ecx, edx;

    __asm__("cpuid" : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                    : "a"(func), "c"(sub));

    info[0] = eax;
    info[1] = ebx;
    info[2] = ecx;
    info[3] = edx;
#else
    __cpuidex(info, func, sub);
#endif
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((constructor))
#endif
static int blst_cpuid_init(void)
{
    int info[4], cap = 0;

    blst_cpuid_query(info, 0, 0);
    if (info[0] > 6) {
        blst_cpuid_query(info, 7, 0);
        cap |= (info[1] >> 19) & 1; /* ADX */
        cap |= (info[1] >> 28) & 2; /* SHA */
    }

    blst_platform_cap = cap;

    return 0;
}

#if defined(_MSC_VER) && !defined(__clang__) && !defined(BLST_DLL_MAIN)
# pragma section(".CRT$XCU", read)
__declspec(allocate(".CRT$XCU")) static int (*blst_cpuid_initializer)(void) = blst_cpuid_init;
#elif defined(__SUNPRO_C)
# pragma init(blst_cpuid_init)
#endif

#elif defined(__aarch64__) || defined(__aarch64) || defined(_M_ARM64) || defined(_M_ARM64EC)

#if defined(__linux__) && (defined(__GNUC__) || defined(__clang__))
extern "C" unsigned long getauxval(unsigned long type) __attribute__((weak));

__attribute__((constructor))
static int blst_cpuid_init(void)
{
    int cap = 0;

    if (getauxval) {
        unsigned long hwcap_ce = getauxval(16);
        cap = (hwcap_ce >> 6) & 1; /* SHA256 */
    }

    blst_platform_cap = cap;

    return 0;
}
#elif defined(__APPLE__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((constructor))
static int blst_cpuid_init(void)
{
    blst_platform_cap = 1; /* SHA256 */
    return 0;
}
#elif defined(__FreeBSD__) && __FreeBSD__ >= 12
# include <sys/auxv.h>
__attribute__((constructor))
static int blst_cpuid_init(void)
{
    unsigned long cap;

    if (elf_aux_info(AT_HWCAP, &cap, sizeof(cap)) == 0)
        blst_platform_cap = (cap & HWCAP_SHA2) != 0;

    return 0;
}
#elif defined(_WIN64)
static int blst_cpuid_init(void)
{
    blst_platform_cap = IsProcessorFeaturePresent(PF_ARM_V8_CRYPTO_INSTRUCTIONS_AVAILABLE);

    return 0;
}

# if defined(__GNUC__) || defined(__clang__)
__attribute__((constructor))
static int blst_cpuid_constructor(void)
{
    return blst_cpuid_init();
}
# elif defined(_MSC_VER) && !defined(BLST_DLL_MAIN)
#  pragma section(".CRT$XCU", read)
__declspec(allocate(".CRT$XCU")) static int (*blst_cpuid_initializer)(void) = blst_cpuid_init;
# endif
#endif

#endif

#if defined(_WIN64) && defined(BLST_DLL_MAIN)
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD dwReason, LPVOID lpvReserved)
{
    (void)lpvReserved;

    if (dwReason == DLL_PROCESS_ATTACH) {
        DisableThreadLibraryCalls(hinstDLL);
        blst_cpuid_init();
    }

    return TRUE;
}

#if defined(_MSC_VER)
/*
 * Even though we don't have memcpy/memset anywhere, MSVC compiler
 * generates calls to them as it recognizes corresponding patterns.
 */
# pragma function(memcpy)
extern "C" void* __cdecl memcpy(void* dst_void, const void* src_void, size_t n)
{
    unsigned char* dst = static_cast<unsigned char*>(dst_void);
    const unsigned char* src = static_cast<const unsigned char*>(src_void);
    void* ret = dst_void;

    while (n--)
        *dst++ = *src++;

    return ret;
}

# pragma function(memset)
extern "C" void* __cdecl memset(void* dst_void, int c, size_t n)
{
    unsigned char* dst = static_cast<unsigned char*>(dst_void);
    void* ret = dst_void;

    while (n--)
        *dst++ = static_cast<unsigned char>(c);

    return ret;
}
#endif
#endif

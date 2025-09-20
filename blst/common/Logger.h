#pragma once
#include "blst/common/defines.h"

namespace blst {

    struct IndexLogger {

        // Raw storage passed in from host/device
        u8* begin;
        u8* end;
        u8* iter;

        int watch_idx;       // the only index allowed to log

        BLST_HD IndexLogger() {}

        BLST_HD IndexLogger(void* buf, uint32_t bytes, int watch)
            : begin((u8*)buf),
            end((u8*)buf + bytes),
            iter((u8*)buf),
            watch_idx(watch)
        {
        }

        // Is this logger active for this thread?
        BLST_HD bool active() const { 
#ifdef __CUDA_ARCH__
			int my_idx = blockIdx.x * blockDim.x + threadIdx.x;
#else
			int my_idx = 0; // single-threaded on host
#endif
            return begin && my_idx == watch_idx; 
        
        }

        // Append a C-string (without requiring strlen from <cstring>)
        BLST_HD void append(const char* s) {
            if (!active())
                return;
            uint32_t n = 0;
            while (s[n] != 0 && iter != end)
                *iter++ = s[n++];

        }

        // Append a hex string
        BLST_HD void append_hex(const void* value, uint64_t length) {

            if (!active() || length == 0) return;

            const char* hex_chars = "0123456789abcdef";
            for (uint64_t i = 0; i < length && iter != end; i++) {
                u8 byte = ((const u8*)value)[i];
                if (iter < end)
                    *iter++ = hex_chars[(byte >> 4) & 0xF];
                if (iter < end)
                    *iter++ = hex_chars[byte & 0xF];
            }
        }



		// Append a hash. hashes the input string using a simple xor-based hash
        BLST_HD void append_hash(const void* value, uint64_t length) {

            if (!active() || length == 0) return;


			u64 hash = 0xcbf29ce484222325; // FNV offset basis
			const u64 prime = 0x100000001b3; // FNV prime
			const u8* data = (const u8*)value;
            for (uint64_t i = 0; i < length; i++) {
                hash ^= data[i];
                hash *= prime;
            }
			append_hex(&hash, sizeof(hash));
        }



        BLST_HD uint32_t size()      const { return iter - begin; }
        BLST_HD uint32_t capacity()  const { return end - begin; }
        BLST_HD bool     overflowed()const { return size() == capacity(); }
    };


    // Host global
	extern IndexLogger host_logger;

    // Device global
#if defined(__CUDACC__)
    extern __device__ IndexLogger device_logger;
#endif


    // Unified accessor
    BLST_HD inline IndexLogger& get_logger() {
#ifdef __CUDA_ARCH__
        return device_logger;
#else
        return host_logger;
#endif
    }

    inline void set_device_logger(IndexLogger l)
    {
#if defined(__CUDACC__)
        cudaMemcpyToSymbol(device_logger, &l, sizeof(l));
#endif
    }


} // namespace blst

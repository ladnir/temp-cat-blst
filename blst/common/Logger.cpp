// logger.cu
#include "logger.h"


namespace blst {
	IndexLogger host_logger;
	__device__ IndexLogger device_logger;

	int init_host_logger = []() {
		host_logger.begin = nullptr;
		host_logger.end = nullptr;
		host_logger.iter = nullptr;
		host_logger.watch_idx = 0;
		return 0;
		}();

	int init_device_logger = []() {
		// use cudaMemcpy to set up device_logger from host if needed
		cudaMemset(&device_logger, 0, sizeof(device_logger));
		return 0;
		}();
}
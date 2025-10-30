
#include <benchmark/benchmark.h>

//BENCHMARK_MAIN();
void FpMain();

void user()
{
	FpMain();
}

int main(int argc, char** argv) {

	if (argc > 1 && argv[1] && strcmp(argv[1], "--user") == 0) {
		user();
		return 0;
	}

	char arg0_default[] = "benchmark";
	char* args_default = arg0_default;
	if (!argv) {

		argc = 1;
		argv = &args_default;
	}
	::benchmark::Initialize(&argc, argv);
	if (::benchmark::ReportUnrecognizedArguments(argc, argv)) 
		return 1;
	::benchmark::RunSpecifiedBenchmarks();
	::benchmark::Shutdown();
	return 0;
}


#include "find_multibase_primes.hpp"
#include "io/io.hpp"

int main(int argc, [[maybe_unused]] char* argv[])
{
	mbp::print_preamble(argc == 2);
	// mbp::print_config();

	mbp::find_multibase_primes mbp;
	mbp.run(argc == 2);
}

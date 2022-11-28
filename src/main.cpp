
#include "find_multibase_primes.hpp"
#include "io/io.hpp"

int main()
{
	mbp::print_preamble();
	// mbp::print_config();

	mbp::find_multibase_primes mbp;
	mbp.run();
}

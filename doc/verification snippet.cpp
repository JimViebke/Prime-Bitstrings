#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <iostream>

using namespace boost::multiprecision;

cpp_int to_base(size_t n, size_t base)
{
	cpp_int number{0};
	for (size_t i = 0; i < 64; ++i)
		number += pow(cpp_int{base} * ((n >> i) & 1), i);
	return number;
}

int main()
{
	const size_t bitstring = 0b1000000010000011110100010001000101001010110111001;

	for (size_t base = 2;; ++base)
	{
		const cpp_int number = to_base(bitstring, base);
		
		if (miller_rabin_test(number, 25))
		{
			std::cout << "Prime in base " << base << ": " << number << '\n';
		}
		else
		{
			std::cout << "Not prime in base " << base << '\n';
			break;
		}
	}
}

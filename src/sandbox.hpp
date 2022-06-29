#pragma once



#include <iostream>

namespace mbp
{
	namespace detail
	{
		void compare_implementations()
		{
			size_t num = 282607273285049; // p11 - too large for 32-bit BPSW
			// size_t num = 113;
			std::cout << "Naive + native implementation:\n";
			std::cout << num << " is " << (mpir_is_prime(num) ? "prime" : "not prime") << std::endl;

			std::cout << "BPSW + native implementation:\n";
			std::cout << num << " is " << (bpsw_1_native::prime((int)num, 50) ? "prime" : "not prime") << std::endl;
		}

		void test_bpsw_1_native()
		{
			size_t primes_found = 0;

			for (size_t num = 2; ; ++num)
			{
				if (bpsw_1_native::prime((int)num, 50))
				{
					std::cout << num << ' ';
					++primes_found;
				}
			}
		}

		void test_binary_to_decimal()
		{
			for (size_t i = 0; i < 50; ++i)
			{
				std::cout << i << '\t' << binary_to_base((size_t)i, 10) << '\n';
			}
		}

		void find_p2_8()
		{
			auto start = current_time_in_ms();

			for (size_t binary = 3; ; binary += 2)
			{
				if (!mpir_is_prime(binary)) continue;

				const size_t b3 = binary_to_base(binary, 3);
				if (!mpir_is_prime(b3)) continue;

				const size_t b4 = binary_to_base(binary, 4);
				if (!mpir_is_prime(b4)) continue;

				const size_t b5 = binary_to_base(binary, 5);
				if (!mpir_is_prime(b5)) continue;

				const size_t b6 = binary_to_base(binary, 6);
				if (!mpir_is_prime(b6)) continue;

				const size_t b7 = binary_to_base(binary, 7);
				if (!mpir_is_prime(b7)) continue;

				const size_t b8 = binary_to_base(binary, 8);
				if (!mpir_is_prime(b8)) continue;

				const size_t b10 = binary_to_base(binary, 10);

				std::cout << b10 << " is prime in bases 2-8 (" <<
					binary << ", " <<
					b3 << ", " <<
					b4 << ", " <<
					b5 << ", " <<
					b6 << ", " <<
					b7 << ", " <<
					b8 << ")\n";
				break;
			}

			std::cout << "Found in " << current_time_in_ms() - start << " ms.\n";
		}

		void mpir_testing()
		{
			gmp_randclass r{ gmp_randinit_mt };
			r.seed(rand());

			const mpz_class million{ 1000 * 1000 };
			const mpz_class billion{ million * 1000 };
			const mpz_class trillion{ billion * 1000 };
			const mpz_class quadrillion{ trillion * 1000 };
			const mpz_class quintillion{ quadrillion * 1000 };
			const mpz_class sextillion{ quintillion * 1000 };
			const mpz_class septillion{ sextillion * 1000 };
			const mpz_class octillion{ septillion * 1000 };
			const mpz_class nonillion{ octillion * 1000 };
			const mpz_class decillion{ nonillion * 1000 };

			const size_t target = 10'000;

			std::cout << mpz_class{ "1000000000000000035000061" } << std::endl;

			const mpz_class num{ "1000000000000000035000061" };
			std::cout << num << std::endl;

			size_t count = 0;
			auto start = current_time_in_ms();
			for (mpz_class i = num; count < target; ++i)
				if (mpz_likely_prime_p(i.get_mpz_t(), r.get_randstate_t(), 0))
					++count;

			std::cout << "A linear search found the next " << count << " primes in " << current_time_in_ms() - start << "ms\n";

			start = current_time_in_ms();
			count = 0;
			for (mpz_class i = num; count < target; mpz_next_prime_candidate(i.get_mpz_t(), i.get_mpz_t(), r.get_randstate_t()))
				if (mpz_likely_prime_p(i.get_mpz_t(), r.get_randstate_t(), 0))
					++count;

			std::cout << "A next-candidate search found the next " << count << " primes in " << current_time_in_ms() - start << "ms\n";
		}

		void mpir_testing_2()
		{
			gmp_randclass r{ gmp_randinit_mt };
			r.seed(rand());

			const mpz_class num{ "413" };
			std::cout << num;

			auto start = current_time_in_ms();

			if (mpz_likely_prime_p(num.get_mpz_t(), r.get_randstate_t(), 0))
			{
				std::cout << " is a prime number\n";
			}
			else
			{
				std::cout << " is not a prime number, but\n";

				mpz_class next_prime = num;

				do
				{
					mpz_next_prime_candidate(next_prime.get_mpz_t(), next_prime.get_mpz_t(), r.get_randstate_t());
				} while (!mpz_likely_prime_p(next_prime.get_mpz_t(), r.get_randstate_t(), 0));

				std::cout << next_prime << " is the next prime up." << std::endl;
			}

			std::cout << "Found in " << current_time_in_ms() - start << " ms\n";
		}

		void pk_testing()
		{
			// print n largest primes, counting down from 2^64
			for (size_t i = -1, n = 0; n < 1000; i--)
			{
				if (pk::is_prime(i))
				{
					n++;
					std::cout << i << '\n';
				}
			}
		}

		void multibase_gap_tests()
		{
			const size_t p11 = 0b1000000010000011110100010001000101001010110111001; // large starting point

			size_t previous_p4 = 0;
			size_t found = 0;

			mpz_class b3;

			for (size_t b2 = p11; found < 100; b2 += 2)
			{
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], b2, 2);
				*result.ptr = '\0';

				b3.set_str(bin_str, 3);
				// b4.set_str(bin_str, 4);
				// b5.set_str(bin_str, 5);

				if (mpir_is_prime(b2) && mpir_is_prime(b3)) // && mpir_is_prime(b4) && mpir_is_prime(b5))
				{
					std::cout << std::bitset<32>(b2) << '\t' << b2;

					if (previous_p4 != 0)
					{
						std::cout << "\t+" << b2 - previous_p4;
						++found;
					}
					previous_p4 = b2;

					std::cout << '\n';
				}

				continue;

				for (int i = 2; i < 10; ++i)
				{
					size_t v = binary_to_base(b2, i);
					std::cout << '\t' << v << (mpir_is_prime(v) ? '\'' : ' ');
				}

				std::cout << '\n';
			}
		}

		void calculate_static_sieve_sizes()
		{
			size_t sieve_size = 1;
			for (auto p : small_primes_lookup)
			{
				if (p == 2) continue;
				sieve_size *= p;
				std::cout << "Static sieve size would be size " << sieve_size << " using product of primes up to " << p << '\n';
				if (sieve_size > 1'000'000'000) return;
			}
			std::cout << "(no suitable sieve size found)\n";
		}
	}
}

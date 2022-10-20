## Prime Bitstrings
 
Search for binary sequences representing primes in multiple bases.

Used to find the smallest multibase 2-12 bitstring at 53 bits:

```
10100001011000101000110101011011011101111110100101011

Prime in base 2: 5678228814363947
Prime in base 3: 7182373369076078560696363
Prime in base 4: 21551394938288350689430921610309
Prime in base 5: 2309293677550186597171330299375393881
Prime in base 6: 29906514580824071121844506525456400608319
Prime in base 7: 89923361635195477013169537997440980347877207
Prime in base 8: 92771144347686777516300137625571104334947582473
Prime in base 9: 42260965547715496390905094752471441113979492989269
Prime in base 10: 10100001011000101000110101011011011101111110100101011
Prime in base 11: 1432168478679934708683628145065883323888809796182285059
Prime in base 12: 131956356859807386250777405024493135451219537385261159117

Not prime in base 13: 8464796925606551530766977421873448346565029813972573631633 =
139277326203487467720073 * 60776561098245690157424928651883721
```

## License

MIT

## Acknowledgments

A big shout-out to the users on MersenneForum.org for invaluable discussion and analysis of the problem.

Libraries used:

* [wbhart/mpir](https://github.com/wbhart/mpir) for prime-testing functions
* [going-digital/Prime](https://github.com/going-digital/Prime) for modular exponentiation utility
* [vectorclass/version2](https://github.com/vectorclass/version2) for utilities and reference use of SIMD intrinsics

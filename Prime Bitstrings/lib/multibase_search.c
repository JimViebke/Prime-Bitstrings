#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <strings.h>
#include <gmp.h>
#include <time.h>

int report=0;
uint64_t replast=0;

clock_t   currentTime, startTime;
double    executionTime;

uint32_t starttime;
uint64_t startp;

mpz_t mpp;
mpz_t n;

void signal_alarm( int dummy )
{
	report=1;
}

int maxbase=0;
int isprime( mpz_t x )
{
	if( mpz_probab_prime_p( x, 5 ) == 0 ) {
		return( 0 );
	}
	return(1);
}

uint32_t bgcd (uint32_t a, uint32_t b)
/* Stolen from proth_sieve.c - Originally by Paul Jobling */
{
	uint32_t k = 0;

	while (!(a&1) && !(b&1))
	{
		a >>= 1;
		b >>= 1;
		k++;
	}

	if (!(a&1))
	{
		a >>= 1;
		while (!(a&1)) a >>= 1;
	}
	else if (!(b&1))
	{
		b >>= 1;
		while (!(b&1)) b >>= 1;
	}

	int32_t t = a-b;

	for (;;)
	{
		if (!t)
		{
			return a<<k;
		}

		t >>= 1;
		while (!(t&1)) t >>= 1;
		if (t < 0) 
		{
			b = -t; 
			t += a;
		}
		else 
		{
			a = t;
			t -= b;
		}
	}
}

int gcd[1155];

void initgcd( void )
{
	int ct=1;
	gcd[0]=1155;
	while( ct <= 1155 ) {
		gcd[ct]=bgcd( ct, 1155 );
		ct++;
	} 
}

int pcok[]={ 0,0,0,0,0,0,0,0,0,0,0, /* 10 */
	     1,0,1,0,0,0,1,0,1,0,0,0, /* 22 */
	     1,0,0,0,0,0,1,0,1,0,0,0, /* 34 */
	     0,0,1,0,0,0,1,0,1,0,0,0, /* 46 */
	     1,0,0,0,0,0,1,0,0,0,0,0, /* 58 */
	     1,0,1,0,0, /* 63 */ };

uint32_t popcount( uint64_t b )
{
	b = (b & 0x5555555555555555LU) + (b >> 1 & 0x5555555555555555LU);
	b = (b & 0x3333333333333333LU) + (b >> 2 & 0x3333333333333333LU);
	b = b + (b >> 4) & 0x0F0F0F0F0F0F0F0FLU;
	b = b + (b >> 8);
	b = b + (b >> 16);
	b = b + (b >> 32) & 0x0000007F;

	return (uint32_t) b;
}

void doit( uint64_t p )
{
	uint64_t a,b;
	uint32_t pc,pca,pcb;
	char nbuf[512];
	int ct;

	if( ! pcok[popcount(p)] ) {
		return;
	}
	a = p& 0xAAAAAAAAAAAAAAAA;
	b = p& 0x5555555555555555;
	pca=popcount( a );
	pcb=popcount( b );
	pc=abs(pca-pcb);
	if( gcd[pc] != 1 ) {
		return;
	}
	mpz_set_ui( mpp, p );
        mpz_get_str( nbuf, 2, mpp );
        ct=3;
	while( 1 ) {
		mpz_set_str( n, nbuf, ct );
		if( ! isprime( n ) ) {
			break;
		}
		ct++;
	}
	if( ct > maxbase ) {
		gmp_printf( "\n%s (%Zd) is a prime bases 2 to %d (%Zd)\n", mpz_get_str( NULL, 2, mpp ), mpp,  ct-1, n );
		maxbase=ct;
	}
	if( report ) {
		if( replast == 0 ) {
			gmp_printf( "%s %llu @ ???kp/s\r", nbuf, p ); fflush( stdout );
		} else {
			uint64_t diff=(uint64_t)(p-startp);
			currentTime = clock();
			executionTime = (double) (currentTime - startTime) / CLOCKS_PER_SEC;
			diff/=executionTime;
			diff/=1000;
			gmp_printf( "%s %llu @ %ukp/s\r", nbuf, p, diff ); fflush( stdout );
		}
		report=0;
		replast=p;
		alarm(3);
	}
}

// This is the number of primes in the sieve.  This will allow us to find
// all primes up to MAX_PRIME^2.
#define MAX_PRIME              70000000

// This is an approximate count of all primes (except 2) less than MAX_PRIME
// for the initial allocation of primeTable.  If the count of primes less
// than MAX_PRIME is less than this value, then some memory will be wasted.
// If the count is greater than, this will limit the amount of memory used
// and the upper limit of primes used in sieving for the search.  The
// amount of memory allocated is approximately 12x this value.
#define MAX_PRIMES_IN_TABLE    4000000
uint32_t *primeTable;
uint32_t  primesInPrimeTable;

// This table is the same size of primeTable.  The values in this table
// correspond to the smallest composite greater than the low end of the
// range being sieved.
uint64_t *compositeTable;

// RANGE_SIZE contains the current range of candidates that has been
// sieved.  The size of this variable should not exceed the size
// of level 2 cache as that will dramatically hurt performance.
#define RANGE_SIZE             500000

// This is the current number of primes in primeTable that are being
// used to sieve the current range.
uint32_t  primesUsedInRange;

uint64_t tests;

/* function prototypes */
void  InitializePrimeSieve();
void  SetupSieve(uint64_t lowEndOfRange);
void  Sieve(uint64_t lowPrime, uint64_t highPrime, void (*fPtr)(uint64_t));

uint64_t get_clocks_per_sec(void) {
        return 33333333;
}

uint64_t get_clocks(void) {
        unsigned long retval;
        asm volatile("mftb %0": "=r"(retval));
        return retval;
}

uint64_t startc, stopc, cc, mulc;

void    InitializePrimeSieve()
{
   uint64_t  maxPrime;
   uint32_t  sqrtMax, i, composite;
   uint32_t  p, minp, *lowPrimes, lowPrimeCount;
   uint8_t  *sieve, *sievePtr;

   // Find all primes less than sqrt(MAX_PRIME)
   sqrtMax = (uint32_t) sqrt((double) MAX_PRIME) + 10;

   // Make sure this is large enough to hold all of the low
   // primes we need.
   lowPrimes = malloc(1000000 * sizeof(uint32_t));
   lowPrimes[0] = 3;
   lowPrimeCount = 1;
   for (p=5; p<sqrtMax; p+=2)
   {
      for (minp=0; minp<=lowPrimeCount; minp++)
      {
         if (lowPrimes[minp] * lowPrimes[minp] > p)
         {
            lowPrimes[lowPrimeCount] = p;
            lowPrimeCount++;
            break;
         }
         if (p % lowPrimes[minp] == 0)
             break;
      }
   }

   // Divide MAX_PRIME by 2 to save memory, also because already know
   // that all even numbers in the sieve are composite
   sieve = malloc((MAX_PRIME >> 1) * sizeof(uint8_t));
   memset(sieve, 1, (MAX_PRIME >> 1));

   for (i=0; i<lowPrimeCount; i++)
   {
      // Get the current low prime.  Start sieving at 3x that prime
      // since 1x is prime and 2x is divisible by 2.
      // sieve[1] = 3, sieve[2] = 5, etc.
      composite = lowPrimes[i] * 3;
      sievePtr = &sieve[(composite - 1) >> 1];

      while (composite < MAX_PRIME)
      { 
         // composite will always be odd, so add 2*lowPrimes[i]
         *sievePtr = 0;
         sievePtr += lowPrimes[i];
         composite += (lowPrimes[i] << 1);
      } 
   }

   primeTable = malloc(MAX_PRIMES_IN_TABLE * sizeof(uint32_t));
   primesInPrimeTable = 0;
   for (i=1; i<(MAX_PRIME >> 1); i++)
   {
      if (sieve[i])
      {
         // Convert the value back to an actual prime
         primeTable[primesInPrimeTable] = (i << 1) + 1;
         primesInPrimeTable++;
         if (primesInPrimeTable == MAX_PRIMES_IN_TABLE)
            break;
      }
   }

   free(sieve);
   free(lowPrimes);

   maxPrime = primeTable[primesInPrimeTable-1];
   maxPrime *= maxPrime;
   printf("The largest prime in the sieve is %d (using %d primes).\n", primeTable[primesInPrimeTable-1], primesInPrimeTable);
   printf("The largest guaranteed prime is no larger than %lld.  Numbers larger than that might be composite\n", maxPrime);
}

void    Sieve(uint64_t lowPrime, uint64_t highPrime, void (*fPtr)(uint64_t))
{
   uint64_t  composite, prime, lowEndOfRange, candidate;
   uint32_t  i;
   uint8_t  *sieve, *sievePtr;

   if (lowPrime < primeTable[primesInPrimeTable-1])
      for (i=0; i<primesInPrimeTable; i++)
      {
         if (primeTable[i] > highPrime)
            return;
         if (primeTable[i] > lowPrime)
            fPtr(primeTable[i]);
      }

   compositeTable = (uint64_t *) malloc((primesInPrimeTable + 1) * sizeof(uint64_t));
   sieve = (uint8_t *) malloc((RANGE_SIZE + 1) * sizeof(uint8_t));
   primesUsedInRange = 0;

   lowEndOfRange = lowPrime - (lowPrime % RANGE_SIZE);

   while (lowEndOfRange < highPrime)
   {
      SetupSieve(lowEndOfRange);

      memset(sieve, 1, RANGE_SIZE);

      for (i=0; i<primesUsedInRange; i++)
      {
         prime = primeTable[i];
         composite = compositeTable[i];
         sievePtr = &sieve[(composite - lowEndOfRange) >> 1];

         while (composite < lowEndOfRange + (RANGE_SIZE << 1))
         {
            *sievePtr = 0;
            sievePtr += prime;
            composite += (prime << 1);
         }

         compositeTable[i] = composite;
     }

     candidate = lowEndOfRange + 1;

     /* printf("\rSearched to %llu  ", candidate); fflush(stdout); */

     for (i=0; i<RANGE_SIZE; i++, candidate+=2)
     {
        if (candidate > highPrime)
           break;
        if (candidate > lowPrime && sieve[i])
           fPtr(candidate);
     }

      lowEndOfRange += (RANGE_SIZE << 1);
   }
}

void    SetupSieve(uint64_t lowEndOfRange)
{
   uint32_t  i, saveUsedInRange;
   uint64_t  maxPrime, lastComposite;
   int64_t   i64;

   if (primesUsedInRange >= primesInPrimeTable)
      return;

   saveUsedInRange = primesUsedInRange;

   if (primesUsedInRange == 0)
      primesUsedInRange = 100000;

   while (primesUsedInRange < primesInPrimeTable)
   {
      maxPrime = (uint64_t) primeTable[primesUsedInRange - 1];
      maxPrime *= maxPrime;
      // The sieve range does not include even numbers, so it
      // contains twice as many candidates
      if (maxPrime > lowEndOfRange + (RANGE_SIZE << 1))
         break;
      primesUsedInRange += 50000;

      if (primesUsedInRange >= primeTable[primesInPrimeTable-1])
      {
         printf("The sieving limit has been reached.  Some non-prime candidates might be evaluated.\n");
         primesUsedInRange = primesInPrimeTable;
         break;
      }
   }

   if (saveUsedInRange == primesUsedInRange)
      return;

   // Find the largest composite greater than lowEndOfRange 
   for (i=saveUsedInRange; i<primesUsedInRange; i++)
   {
      maxPrime = primeTable[i];
      lastComposite = (lowEndOfRange / maxPrime) * maxPrime;
      compositeTable[i] = lastComposite + maxPrime;

      // We only care about odd composites since the 
      // sieve range only refers to odd values
      if (!(compositeTable[i] & 1))
         compositeTable[i] += maxPrime;
   }
}

int main(int argc, char **argv)
{
	uint64_t x;
	int nb,ct;
	uint64_t pc;
	uint64_t pull;

   long      i;
   uint64_t  lowPrime, highPrime;

   lowPrime *= 1000000000;
   highPrime *= 1000000000;

   InitializePrimeSieve();

	initgcd();

	mpz_init( n );
	mpz_init( mpp );

	if( argc == 2 ) {
		nb=atoi( argv[1] );

		x=1; x<<=nb;
		lowPrime=x;

		x<<=1;
		highPrime=lowPrime<<1;
		x>>=1;
	} else if( argc == 3 ) {
		lowPrime=atol( argv[1] );
		highPrime=atol( argv[2] );
	}

	printf( "Doing primes (%llu to %llu)\n", lowPrime, highPrime);
	report=0;
	signal( SIGALRM, signal_alarm );
	alarm(3);

	mpz_set_ui( mpp, lowPrime );

   tests = 0;
   startTime = clock();
   stopc = startc = 0;

   printf("The search will use a hand-written ASM routine that combines p-adic arithmetic with Montgomery multiplication.\n");
   printf("This algorithm has an upper limit of 2^52.  Due to overflows, numbers over that value can result in invalid results.\n");

   startp=lowPrime;
   Sieve(lowPrime, highPrime, doit);

   currentTime = clock();
   executionTime = (double) (currentTime - startTime) / CLOCKS_PER_SEC;
   printf("The search completed in %0.2f seconds (%d primes tested) (%lf clocks per expmod)\n", executionTime, tests,
              (double) startc * 2e9 / tests / get_clocks_per_sec());

   free(primeTable);
   free(compositeTable);
   return 0;
}



/* 
   This is an implementation of the Mersenne Twister which was written by  Michael
   Brundage. The original source code, which was copied and pasted into this file,
   can be found at http://qbrundage.com/michaelb/pubs/essays/random_number_generation,
   and a description of the Mersenne Twister random number generator can be found here:
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   It has a claimed periodicity of 2^{19937}-1, which is pretty incredible.  The algorithm
   as-is is mainly intended for Monte Carlo realizations, and is not intended to be used for
   cryptography.  See the website for more information.

   This code was placed into the public domain by Michael Brundage, and was put into Enzo
   by Brian W. O'Shea on 11 December 2007.  It uses the system random number generator as
   an initial seed, and the user must specify a seed for the system random number generator.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>

#define MT_LEN          624
#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

int mt_index;
unsigned long int mt_buffer[MT_LEN];

void mt_init(unsigned int seed) {

  srand(seed);

  for (int i = 0; i < MT_LEN; i++)
    mt_buffer[i] = ((unsigned long int) rand() );

  mt_index = 0;
}

void mt_read(std::ifstream& input)
{
  input >> mt_index;

  for (int i = 0; i < MT_LEN; i++)
    input >> mt_buffer[i];
}    

void mt_write(std::ofstream& output)
{
  output << mt_index << '\n';

  for (int i = 0; i < MT_LEN; i++)
    output << mt_buffer[i] << '\n';
}
 
unsigned long int mt_random() {
  unsigned long int * b = mt_buffer;
  int idx = mt_index;
  unsigned long int s;
  int i;

  if (idx == MT_LEN*sizeof(unsigned long int))
    {
      idx = 0;
      i = 0;
      for (; i < MT_IB; i++) {
	s = TWIST(b, i, i+1);
	b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
      }
      for (; i < MT_LEN-1; i++) {
	s = TWIST(b, i, i+1);
	b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
      }
      
      s = TWIST(b, MT_LEN-1, 0);
      b[MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
  mt_index = idx + sizeof(unsigned long int);
  return *(unsigned long int *)((unsigned char *)b + idx);
  /*
    Matsumoto and Nishimura additionally confound the bits returned to the caller
    but this doesn't increase the randomness, and slows down the generator by
    as much as 25%.  So I omit these operations here.
    
    r ^= (r >> 11);
    r ^= (r << 7) & 0x9D2C5680;
    r ^= (r << 15) & 0xEFC60000;
    r ^= (r >> 18);
  */
}

/****************************************************************
 *
 * random_dsfmt.h: Contains the dSFMT, a prng
 *	adapted from Mutsuo Saito and Makoto Matsumoto
 *
 *****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************/

#ifndef RANDOM_DSFMT_H_INCLUDED
#define RANDOM_DSFMT_H_INCLUDED

/*
Copyright (c) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima
University.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of the Hiroshima University nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * @file dSFMT.h
 *
 * @brief double precision SIMD oriented Fast Mersenne Twister(dSFMT)
 * pseudorandom number generator based on IEEE 754 format.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University. All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#ifndef DSFMT_H
#define DSFMT_H

#include <stdio.h>
#include <assert.h>

#define DSFMT_MEXP 19937

/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^DSFMT_MEXP-1.
 * #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define DSFMT_N32 (DSFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#if !defined(DSFMT_BIG_ENDIAN)
#if defined(__BYTE_ORDER) && defined(__BIG_ENDIAN)
#if __BYTE_ORDER == __BIG_ENDIAN
#define DSFMT_BIG_ENDIAN 1
#endif
#elif defined(_BYTE_ORDER) && defined(_BIG_ENDIAN)
#if _BYTE_ORDER == _BIG_ENDIAN
#define DSFMT_BIG_ENDIAN 1
#endif
#elif defined(__BYTE_ORDER__) && defined(__BIG_ENDIAN__)
#if __BYTE_ORDER__ == __BIG_ENDIAN__
#define DSFMT_BIG_ENDIAN 1
#endif
#elif defined(BYTE_ORDER) && defined(BIG_ENDIAN)
#if BYTE_ORDER == BIG_ENDIAN
#define DSFMT_BIG_ENDIAN 1
#endif
#elif defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN) || \
    defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)
#define DSFMT_BIG_ENDIAN 1
#endif
#endif

#if defined(DSFMT_BIG_ENDIAN) && defined(__amd64)
#undef DSFMT_BIG_ENDIAN
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
#if !defined(DSFMT_UINT32_DEFINED) && !defined(SFMT_UINT32_DEFINED)
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#define UINT64_C(v) (v##ui64)
#define DSFMT_UINT32_DEFINED
#if !defined(inline)
#define inline __inline
#endif
#endif
#else
#include <inttypes.h>
#if !defined(inline)
#if defined(__GNUC__)
#define inline __inline__
#else
#define inline
#endif
#endif
#endif

#ifndef PRIu64
#if defined(_MSC_VER) || defined(__BORLANDC__)
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#else
#define PRIu64 "llu"
#define PRIx64 "llx"
#endif
#endif

#ifndef UINT64_C
#define UINT64_C(v) (v##ULL)
#endif

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#if !defined(__APPLE__)
#include <altivec.h>
#endif
/** 128-bit data structure */
union W128_T {
  vector unsigned int s;
  uint64_t u[2];
  uint32_t u32[4];
  double d[2];
};

#elif defined(HAVE_SSE2)
#include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
  __m128i si;
  __m128d sd;
  uint64_t u[2];
  uint32_t u32[4];
  double d[2];
};
#else /* standard C */
/** 128-bit data structure */
union W128_T {
  uint64_t u[2];
  uint32_t u32[4];
  double d[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/** the 128-bit internal state array */
struct DSFMT_T {
  w128_t status[DSFMT_N + 1];
  int idx;
};
typedef struct DSFMT_T dsfmt_t;

/** dsfmt internal state vector */
extern dsfmt_t dsfmt_global_data;
/** dsfmt mexp for check */
extern const int dsfmt_global_mexp;

void dsfmt_gen_rand_all(dsfmt_t* dsfmt);
void dsfmt_fill_array_open_close(dsfmt_t* dsfmt, double array[], int size);
void dsfmt_fill_array_close_open(dsfmt_t* dsfmt, double array[], int size);
void dsfmt_fill_array_open_open(dsfmt_t* dsfmt, double array[], int size);
void dsfmt_fill_array_close1_open2(dsfmt_t* dsfmt, double array[], int size);
void dsfmt_chk_init_gen_rand(dsfmt_t* dsfmt, uint32_t seed, int mexp);
void dsfmt_chk_init_by_array(dsfmt_t* dsfmt, uint32_t init_key[],
                             int key_length, int mexp);
const char* dsfmt_get_idstring(void);
int dsfmt_get_min_array_size(void);

#if defined(__GNUC__)
#define DSFMT_PRE_INLINE inline static
#define DSFMT_PST_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#define DSFMT_PRE_INLINE __forceinline static
#define DSFMT_PST_INLINE
#else
#define DSFMT_PRE_INLINE inline static
#define DSFMT_PST_INLINE
#endif
DSFMT_PRE_INLINE uint32_t dsfmt_genrand_uint32(dsfmt_t* dsfmt) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close1_open2(dsfmt_t* dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close_open(dsfmt_t* dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_close(dsfmt_t* dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_open(dsfmt_t* dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE uint32_t dsfmt_gv_genrand_uint32(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_close1_open2(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_close_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_open_close(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_open_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_open_close(double array[],
                                                     int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_close_open(double array[],
                                                     int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_open_open(double array[],
                                                    int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_close1_open2(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_init_gen_rand(uint32_t seed) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_init_by_array(uint32_t init_key[],
                                             int key_length) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_init_gen_rand(dsfmt_t* dsfmt,
                                          uint32_t seed) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_init_by_array(dsfmt_t* dsfmt, uint32_t init_key[],
                                          int key_length) DSFMT_PST_INLINE;

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static uint32_t dsfmt_genrand_uint32(dsfmt_t* dsfmt)
{
  uint32_t r;
  uint64_t* psfmt64 = &dsfmt->status[0].u[0];

  if (dsfmt->idx >= DSFMT_N64) {
    dsfmt_gen_rand_all(dsfmt);
    dsfmt->idx = 0;
  }
  r = psfmt64[dsfmt->idx++] & 0xffffffffU;
  return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).  This is
 * the primitive and faster than generating numbers in other ranges.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close1_open2(dsfmt_t* dsfmt)
{
  double r;
  double* psfmt64 = &dsfmt->status[0].d[0];

  if (dsfmt->idx >= DSFMT_N64) {
    dsfmt_gen_rand_all(dsfmt);
    dsfmt->idx = 0;
  }
  r = psfmt64[dsfmt->idx++];
  return r;
}

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function.  This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static uint32_t dsfmt_gv_genrand_uint32(void)
{
  return dsfmt_genrand_uint32(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_close1_open2(void)
{
  return dsfmt_genrand_close1_open2(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close_open(dsfmt_t* dsfmt)
{
  return dsfmt_genrand_close1_open2(dsfmt) - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_close_open(void)
{
  return dsfmt_gv_genrand_close1_open2() - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_close(dsfmt_t* dsfmt)
{
  return 2.0 - dsfmt_genrand_close1_open2(dsfmt);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_open_close(void)
{
  return 2.0 - dsfmt_gv_genrand_close1_open2();
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_open(dsfmt_t* dsfmt)
{
  double* dsfmt64 = &dsfmt->status[0].d[0];
  union {
    double d;
    uint64_t u;
  } r;

  if (dsfmt->idx >= DSFMT_N64) {
    dsfmt_gen_rand_all(dsfmt);
    dsfmt->idx = 0;
  }
  r.d = dsfmt64[dsfmt->idx++];
  r.u |= 1;
  return r.d - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_open_open(void)
{
  return dsfmt_genrand_open_open(&dsfmt_global_data);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [1, 2) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_fill_array_close1_open2() except that this function uses
 * \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_close1_open2(double array[], int size)
{
  dsfmt_fill_array_close1_open2(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1] to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() and \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_open_close(double array[], int size)
{
  dsfmt_fill_array_open_close(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_close_open(double array[], int size)
{
  dsfmt_fill_array_close_open(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_open_open(double array[], int size)
{
  dsfmt_fill_array_open_open(&dsfmt_global_data, array, size);
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 */
inline static void dsfmt_init_gen_rand(dsfmt_t* dsfmt, uint32_t seed)
{
  dsfmt_chk_init_gen_rand(dsfmt, seed, DSFMT_MEXP);
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed. This function uses \b global variables.
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_init_gen_rand()
 */
inline static void dsfmt_gv_init_gen_rand(uint32_t seed)
{
  dsfmt_init_gen_rand(&dsfmt_global_data, seed);
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * @param dsfmt dsfmt state vector
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 */
inline static void dsfmt_init_by_array(dsfmt_t* dsfmt, uint32_t init_key[],
                                       int key_length)
{
  dsfmt_chk_init_by_array(dsfmt, init_key, key_length, DSFMT_MEXP);
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * This function uses \b global variables.
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_init_by_array()
 */
inline static void dsfmt_gv_init_by_array(uint32_t init_key[], int key_length)
{
  dsfmt_init_by_array(&dsfmt_global_data, init_key, key_length);
}

#if !defined(DSFMT_DO_NOT_USE_OLD_NAMES)
DSFMT_PRE_INLINE const char* get_idstring(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE int get_min_array_size(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void init_gen_rand(uint32_t seed) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void init_by_array(uint32_t init_key[],
                                    int key_length) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_close1_open2(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_close_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_open_close(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_open_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_open_close(double array[],
                                            int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_close_open(double array[],
                                            int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_open_open(double array[],
                                           int size) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_close1_open2(double array[],
                                              int size) DSFMT_PST_INLINE;

/**
 * This function is just the same as dsfmt_get_idstring().
 * @return id string.
 * see also \sa dsfmt_get_idstring()
 */
inline static const char* get_idstring(void) { return dsfmt_get_idstring(); }
/**
 * This function is just the same as dsfmt_get_min_array_size().
 * @return minimum size of array used for fill_array functions.
 * see also \sa dsfmt_get_min_array_size()
 */
inline static int get_min_array_size(void)
{
  return dsfmt_get_min_array_size();
}

/**
 * This function is just the same as dsfmt_gv_init_gen_rand().
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_gv_init_gen_rand(), \sa dsfmt_init_gen_rand().
 */
inline static void init_gen_rand(uint32_t seed)
{
  dsfmt_gv_init_gen_rand(seed);
}

/**
 * This function is just the same as dsfmt_gv_init_by_array().
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_gv_init_by_array(), \sa dsfmt_init_by_array().
 */
inline static void init_by_array(uint32_t init_key[], int key_length)
{
  dsfmt_gv_init_by_array(init_key, key_length);
}

/**
 * This function is just the same as dsfmt_gv_genrand_close1_open2().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close1_open2() \sa
 * dsfmt_gv_genrand_close1_open2()
 */
inline static double genrand_close1_open2(void)
{
  return dsfmt_gv_genrand_close1_open2();
}

/**
 * This function is just the same as dsfmt_gv_genrand_close_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close_open() \sa
 * dsfmt_gv_genrand_close_open()
 */
inline static double genrand_close_open(void)
{
  return dsfmt_gv_genrand_close_open();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_close().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_close() \sa
 * dsfmt_gv_genrand_open_close()
 */
inline static double genrand_open_close(void)
{
  return dsfmt_gv_genrand_open_close();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_open() \sa
 * dsfmt_gv_genrand_open_open()
 */
inline static double genrand_open_open(void)
{
  return dsfmt_gv_genrand_open_open();
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_close().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_close(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_open_close(double array[], int size)
{
  dsfmt_gv_fill_array_open_close(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_close_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_close_open(double array[], int size)
{
  dsfmt_gv_fill_array_close_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_open_open(double array[], int size)
{
  dsfmt_gv_fill_array_open_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close1_open2().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_close1_open2(double array[], int size)
{
  dsfmt_gv_fill_array_close1_open2(array, size);
}
#endif /* DSFMT_DO_NOT_USE_OLD_NAMES */

#endif /* DSFMT_H */

/*----------------------
  the parameters of DSFMT
  following definitions are in dSFMT-paramsXXXX.h file.
  ----------------------*/
/** the pick up position of the array.
#define DSFMT_POS1 122
*/

/** the parameter of shift left as four 32-bit registers.
#define DSFMT_SL1 18
 */

/** the parameter of shift right as four 32-bit registers.
#define DSFMT_SR1 12
*/

/** A bitmask, used in the recursion.  These parameters are introduced
 * to break symmetry of SIMD.
#define DSFMT_MSK1 (uint64_t)0xdfffffefULL
#define DSFMT_MSK2 (uint64_t)0xddfecb7fULL
*/

/** These definitions are part of a 128-bit period certification vector.
#define DSFMT_PCV1	UINT64_C(0x00000001)
#define DSFMT_PCV2	UINT64_C(0x00000000)
*/

#define DSFMT_LOW_MASK UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR 12

/* for sse2 */
#if defined(HAVE_SSE2)
#define SSE2_SHUFF 0x1b
#elif defined(HAVE_ALTIVEC)
#if defined(__APPLE__) /* For OSX */
#define ALTI_SR (vector unsigned char)(4)
#define ALTI_SR_PERM \
  (vector unsigned char)(15, 0, 1, 2, 3, 4, 5, 6, 15, 8, 9, 10, 11, 12, 13, 14)
#define ALTI_SR_MSK \
  (vector unsigned int)(0x000fffffU, 0xffffffffU, 0x000fffffU, 0xffffffffU)
#define ALTI_PERM \
  (vector unsigned char)(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3)
#else
#define ALTI_SR \
  {             \
    4           \
  }
#define ALTI_SR_PERM                                      \
  {                                                       \
    15, 0, 1, 2, 3, 4, 5, 6, 15, 8, 9, 10, 11, 12, 13, 14 \
  }
#define ALTI_SR_MSK                                    \
  {                                                    \
    0x000fffffU, 0xffffffffU, 0x000fffffU, 0xffffffffU \
  }
#define ALTI_PERM                                        \
  {                                                      \
    12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 \
  }
#endif
#endif

#define DSFMT_POS1 117
#define DSFMT_SL1 19
#define DSFMT_MSK1 UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2 UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1 0x000ffaffU
#define DSFMT_MSK32_2 0xfffffb3fU
#define DSFMT_MSK32_3 0x000ffdffU
#define DSFMT_MSK32_4 0xfc90fffdU
#define DSFMT_FIX1 UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2 UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1 UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2 UINT64_C(0x0000000000000001)
#define DSFMT_IDSTR "dSFMT2-19937:117-19:ffafffffffb3f-ffdfffc90fffd"

/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__) /* For OSX */
#define ALTI_SL1 (vector unsigned int)(3, 3, 3, 3)
#define ALTI_SL1_PERM \
  (vector unsigned char)(2, 3, 4, 5, 6, 7, 30, 30, 10, 11, 12, 13, 14, 15, 0, 1)
#define ALTI_SL1_MSK \
  (vector unsigned int)(0xffffffffU, 0xfff80000U, 0xffffffffU, 0xfff80000U)
#define ALTI_MSK                                                     \
  (vector unsigned int)(DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, \
                        DSFMT_MSK32_4)
#else /* For OTHER OSs(Linux?) */
#define ALTI_SL1 \
  {              \
    3, 3, 3, 3   \
  }
#define ALTI_SL1_PERM                                      \
  {                                                        \
    2, 3, 4, 5, 6, 7, 30, 30, 10, 11, 12, 13, 14, 15, 0, 1 \
  }
#define ALTI_SL1_MSK                                   \
  {                                                    \
    0xffffffffU, 0xfff80000U, 0xffffffffU, 0xfff80000U \
  }
#define ALTI_MSK                                               \
  {                                                            \
    DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4 \
  }
#endif

#endif  // RANDOM_DSFMT_H_INCLUDED

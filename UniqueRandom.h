/**************************************************************************************************
 *            THIS CLASS IS ONLY DESIGNED FOR COMPARISON WITH OLD VERSION NON SMP                 *
 **************************************************************************************************/
#ifndef _UNIQUE_RANDOM_H_
#define _UNIQUE_RANDOM_H_
#include <cstddef>
#include <cmath>
#include "config.h"
#include "icsilog.h"
#include <cstring>
//#include "cond_var.h"
#include <stdint.h>
#include <limits.h>

extern int NumberOfThreads;

//// A routine to generate poisson deviates with mean xm
//extern double poidev(double xm, long *idum);
////returns gamma deviates
//extern double gamma_dev(double a);
//extern double igamma_dev(int ia);

// extern "C" void zigset(unsigned long jsrseed);
// extern "C" double rexp();

//typedef double my_double
typedef double my_double; // Be careful, program crashes if we change double to real

const int numRanNumToKeep = 10000;
        

/*Taylor (deg 3) implementation of the log: http://www.flipcode.com/cgi-bin/fcarticles.cgi?show=63828*/
//Not extremely precise and slower than when using LOOKUP_TABLE2
inline float fast_log (float val)
{
   register int *const     exp_ptr = ((int*)&val);
   register int            x = *exp_ptr;
   register const int      log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return ((val + log_2)* 0.69314718);
}

#ifdef _RAN3_
  const long _MBIG=1000000000L;
  const long _MBIG_1=_MBIG-1L;
  const long _MSEED=161803398L;
  const long _MZ=0L;
  const double _FAC=1.0/_MBIG;
  const double _OVER=0.999999999;
  const double _LOG_OVER=-0.000000001;
  const long _HALF_MBIG_=_MBIG/2;
  //Loro_05_11_18 Transforming double into float to prevent a bug leading to log of random numbers >0 to be generated
//  const double _LOG_FAC=log(_FAC);
  const float _LOG_FAC=log(_FAC);
  const double _ALMOST_ONE=1.0*(_MBIG-1.0)*_FAC;
#else
  
//Loro_29_04_20
// const double _inv_max_uint64_ = 1.0/_UI64_MAX;
 const double _inv_max_uint64_ = 1.0/UINT64_MAX;
 const double _inv_max_uint32_ = 1.0/UINT32_MAX;
  //Loro_05_11_18 Correction to prevent generation of logRan numbers >0
//  const double log_inv_max_uint64_ = log(_inv_max_uint64_);
 
//Loro_29_04_20
//  const float log_inv_max_uint64_ = log(_inv_max_uint64_); 
  const float log_inv_max_uint64_ = log(_inv_max_uint64_); 
//  const uint64_t _half_max_ui64_= _UI64_MAX / 2;
  const uint64_t _half_max_ui64_= UINT64_MAX / 2;
  const uint32_t _half_max_ui32_= UINT32_MAX / 2;
  
#endif
  
  

class GlobalRandom { 
  /************************************************************************************************
   *                                     Definitions                                              *
   ************************************************************************************************/
private:
//   typedef unsigned long uint64_t;
  /************************************************************************************************
   *                                  Local variables                                             *
   ************************************************************************************************/

#ifdef _XOSHIRO_
   my_double ranArray[numRanNumToKeep];
#else      
   my_double ranArray[numRanNumToKeep];
#endif  
   float logRanArray[numRanNumToKeep];
   int numRanGiven;
   int numLogRanGiven;
   #ifdef _ICSI_LOG_
   static float * LOOKUP_TABLE2;
   static float * LOOKUP_TABLE23; //Lookup table with full precision... needed to compute logs precisely with icsilog())
   #endif

private: 
#ifdef _RAN3_
   int _inext;
   int _inextp;
   long _ma[56]; 
#else
   #ifdef _XOROSHIRO_
   uint64_t  s[2]={0,0};
   #endif

   #ifdef _XOSHIRO_
   uint32_t s[4]={0,0,0,0};
   #endif
#endif
    
  long _Seed;
   
  /* Poisson variable */
  my_double sq,alxm,g, oldm;
  
  /************************************************************************************************
   *                                  Constructor/Destructor                                      *
   ************************************************************************************************/
public:
   
   static unsigned  logPrecision;
     
   // -----------------NEW RANDOM NUMBER GENERATORS ----------------------------

   //based on https://nullprogram.com/blog/2017/09/21/                                 //Comparison of random number generators

   //   xoroshiro128+ is the obvious winner in this benchmark and it seems to be the best 
   //   64-bit simulation PRNG available. If you need a fast, 
   //   quality PRNG, just drop these 11 lines into your C or C++ program:
   #ifdef _XOROSHIRO_
   uint64_t
   xoroshiro128plus()
   {
      uint64_t s0 = s[0];
      uint64_t s1 = s[1];
      uint64_t result = s0 + s1;
      s1 ^= s0;
      s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
      s[1] = (s1 << 36) | (s1 >> 28);
      return result;
   }
   double
   xoroshiro128plus_real()
   {
      uint64_t s0 = s[0];
      uint64_t s1 = s[1];
      uint64_t result = s0 + s1;
      s1 ^= s0;
      s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
      s[1] = (s1 << 36) | (s1 >> 28);
      return (result*_inv_max_uint64_);
   }
   #endif
   //If it weren’t for xoroshiro128+, then xorshift128+ would have been the 
   //winner of the benchmark and my new favorite choice.

   #ifdef _XORSHIFT_
   uint64_t
   xorshift128plus()
   {
       uint64_t x = s[0];
       uint64_t y = s[1];
       s[0] = y;
       x ^= x << 23;
       s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
       return s[1] + y;
   }
   double
   xorshift128plus()
   {
      uint64_t x = s[0];
      uint64_t y = s[1];
      s[0] = y;
      x ^= x << 23;
      s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
      return ((s[1] + y)*_inv_max_uint64_);
   }
   #endif

   
   #ifdef _XOSHIRO_

//   static uint32_t s[4]; //Loro_22_12_21 Already defined above as a private variable

   //https://prng.di.unimi.it/xoshiro128plus.c
   
   static inline uint32_t rotl(const uint32_t x, int k) {
      return (x << k) | (x >> (32 - k));
   }
   
   uint32_t xoshiro128plus(void) {
      const uint32_t result = s[0] + s[3];

      const uint32_t t = s[1] << 9;

      s[2] ^= s[0];
      s[3] ^= s[1];
      s[1] ^= s[2];
      s[0] ^= s[3];

      s[2] ^= t;

      s[3] = rotl(s[3], 11);

      return result;
   }
   
   
   double xoshiro128plus_real(void) {
      const uint32_t result = s[0] + s[3];

      const uint32_t t = s[1] << 9;

      s[2] ^= s[0];
      s[3] ^= s[1];
      s[1] ^= s[2];
      s[0] ^= s[3];

      s[2] ^= t;

      s[3] = rotl(s[3], 11);

//      if (result== UINT32_MAX) return 0;
      return (result*_inv_max_uint32_);
   }

   #endif
   
   #ifdef _RAN3_
   GlobalRandom() : oldm(-1.0) {
      _inext=0;
      _inextp=0; 
      _Seed=-1;
      sq=0;
      alxm=0;
      g=0;
   }
   GlobalRandom(const GlobalRandom & GR) {
      _inext=GR._inext;
      _inextp=GR._inextp;
      memcpy(_ma, GR._ma, 56*sizeof(GR._ma[0])); 
      _Seed=GR._Seed;
      sq=GR.sq;
      alxm=GR.alxm;
      g=GR.g;
      oldm=GR.oldm;
   }
   #else
   GlobalRandom() : oldm(-1.0) {
      sq=0;
      alxm=0;
      g=0;
   }
   GlobalRandom(const GlobalRandom & GR) {
      _Seed=GR._Seed;
      sq=GR.sq;
      alxm=GR.alxm;
      g=GR.g;
      oldm=GR.oldm;
   }
   #endif

   ~GlobalRandom() {
   #ifdef _ICSI_LOG_
   #pragma omp critical
      {
      if (GlobalRandom::LOOKUP_TABLE2) delete [] GlobalRandom::LOOKUP_TABLE2;
      GlobalRandom::LOOKUP_TABLE2=NULL;
      if (GlobalRandom::LOOKUP_TABLE23) delete [] GlobalRandom::LOOKUP_TABLE23;
      GlobalRandom::LOOKUP_TABLE23=NULL;
      }
   #endif
   }
   
  
  /************************************************************************************************
   *                                     Member functions                                         *
   ************************************************************************************************/
private:

public:
  long seed() {return _Seed;}
  void InitializeGenerator(const long Seed);
  bool YesOrNo() {
     
#ifdef _RAN3_
    if (++_inext == 56) _inext=1;
    if (++_inextp == 56) _inextp=1;
    register long mj =_ma[_inext]-_ma[_inextp];
    if (mj < _MZ) mj += _MBIG;
    _ma[_inext]=mj;
    return (mj < _HALF_MBIG_);
#endif
    
#ifdef _XOROSHIRO_ 
    return (xoroshiro128plus() < _half_max_ui64_); 
#endif
    
#ifdef _XORSHIFT_
    return (xorshift128plus() < _half_max_ui64_);
#endif
     
#ifdef _XOSHIRO_ 
    return (xoshiro128plus() < _half_max_ui32_); 
#endif
  }
  
#ifdef _XOSHIRO_ 
   inline double NormalizedFlatDistributionValue() {
      if (numRanGiven==numRanNumToKeep) updateRanNums();
      return ranArray[numRanGiven++];
   }  
#else   
   inline my_double NormalizedFlatDistributionValue() {
      if (numRanGiven==numRanNumToKeep) updateRanNums();
      return ranArray[numRanGiven++];
   }  
#endif
   inline uint64_t integerValue(uint64_t max) {
      if (numRanGiven==numRanNumToKeep) updateRanNums();
      uint64_t val=ranArray[numRanGiven++]*max;
      return (max==val? val-1: val);
   } 
   
   inline float logNormalizedFlatDistributionValue() {      
      if (numLogRanGiven==numRanNumToKeep) updateLogRanNums();
      return logRanArray[numLogRanGiven++];
   }
      
   inline float icsi_log23(my_double & val) {      
   #ifndef _ICC_
      return icsi_log_v2_23(val,GlobalRandom::LOOKUP_TABLE23);
   #else
      return log(val); 
   #endif
   }

   //Returns a geometrically distributed random variate between zero and infinity
   //Number of failures before first success
   int GeometricDistributionValue(const my_double p) {
      if (p < 1.e-10 || p>=1.0) return 0;
      my_double rand;
      //    while ((rand=NormalizedFlatDistributionValue())==0);
      rand=NormalizedFlatDistributionValue();
      if (rand) return (int) floor(log(rand)/log(p));
      return 0;
   }
   my_double PoissonDistributionValue(my_double xm); 
  
   long getSeed() const { return _Seed; } 
//  void PrintState();
 
   
   #ifdef _ICSI_LOG_
    void initializeLookUpTables();
   #endif
    void updateRanNums();
    void updateLogRanNums();

};

#endif /* _UNIQUE_RANDOM_H_ */

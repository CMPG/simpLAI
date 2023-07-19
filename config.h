#ifndef _CONFIG_
#define _CONFIG_

//#define _DEBUG_ANCIENT_LINEAGES_

//**************************  VARIOUS GENERAL SETTINGS  ************************ 
   //#define _DEBUG_  
   //#define _GCC_
   //#define _PRINT_TREE_PAUP_            //Output generated trees in PAUP format (one tree per locus)
//   #define _LOG_FILE_
//   #define _ADVANCED_FEATURES_
//   #define _POLYNOMIAL_FIT_

//    #define _DEBUG_DEAD_DEMES_ 

//*****************************  RECOMBINATION  ******************************** 
//   #define _OLD_REC_IMPL_
#define _NEW_FAST_REC_IMPL_
//******************************************************************************

//**************************  RANDOM NUMBERS AND LOGS  *************************
#define _ICSI_LOG_ 
#ifdef _ICSI_LOG_
   //For  _ICSI_LOG_ fast log generation
   const unsigned PRECISIONBITS=23; //Cannot be larger than 23
          //PRECISIONBITS is the number of bits used from the mantissa (0<=n<=23). 
          //Higher n means higher accuracy but slower execution. We found that a good value for n is 14 
#endif

//#define _RAN3_
#define _XOROSHIRO_
//Loro_22_12_21
//#define _XOSHIRO_  //Does not seem really better that xoroshiro, and does not work with floats so no real benefit from it...
//#define _XORSHIFT_

///*************************************************************************/
///*                         RANDOM GENERATOR                              */
///*************************************************************************/
//#define RANDOM_GENERATOR_USE_VECTOR
//#if defined(RANDOM_GENERATOR_USE_VECTOR)
//# define RANDOM_USE_VECTOR
//#endif
//
///* Use of vectors is not possible for unique generator */
//#define _UNIQUE_GENERATOR_
//#if defined(_UNIQUE_GENERATOR_)
//#undef RANDOM_USE_VECTOR
//#endif
//
//#define RANDOM_LOGARITHM
//#if defined(RANDOM_LOGARITHM) && defined(RANDOM_USE_VECTOR)
//# define USE_STANDARD_MATH_LIBRARY
//#endif
//
////#define USE_VECTOR_MATH_LIBRARY Needs Intel math kernel library for windows or linux
//
///********************* GNU C++ compiler *************************/
//#define ALIGNED(var, alignment) var __attribute__((aligned(alignment)))
//#define __ALWAYS_INLINE __attribute__((always_inline))
//#define TIMEVAL_T struct timeval
   
   
//Loro_27_09_18
#ifdef _LINUX_
#define _UI64_MAX 0xffffffffffffffffull
#endif

#endif //_CONFIG_

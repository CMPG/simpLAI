#include "config.h"
//#include "cond_var.h"

#include "UniqueRandom.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
//#include "simSettings.h"

#define PI 3.141592654

#ifdef _ICSI_LOG_
   float * GlobalRandom::LOOKUP_TABLE2 = NULL;
   float * GlobalRandom::LOOKUP_TABLE23 = NULL;
   unsigned  GlobalRandom::logPrecision = PRECISIONBITS;
#endif
   
using namespace std;
    
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
//------------------------------------------------------------------------------
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
//------------------------------------------------------------------------------
void GlobalRandom::InitializeGenerator(const long Seed) {
    
#ifdef _RAN3_
  long mj,mk;
  int i,ii,k;
    
  oldm=-1.0; //Just to be sure...
  
  //Initialize normal random numbers -------------------------------------------
  long lSeed = Seed;

  mj=_MSEED-(lSeed < 0 ? -lSeed : lSeed);
  mj %= _MBIG;
  _ma[55]=mj;
  mk=1;
  for (i=1;i<=54;++i) {
    ii=(21*i) % 55;
    _ma[ii]=mk;
    mk=mj-mk;
    if (mk < _MZ) mk += _MBIG;
    mj=_ma[ii];
  }
  for (k=1;k<=4;++k) {
    for (i=1;i<=55;++i) {
      _ma[i] -= _ma[1+(i+30) % 55];
      if (_ma[i] < _MZ) _ma[i] += _MBIG;
    }
  }
  
  _Seed    = Seed;
  _inext   = 0;
  _inextp  = 31;
  
  //Loro_25_12_13 Adding this here to keep syncing with ran3...
   if (++_inext == 56)  _inext=1;
   if (++_inextp == 56) _inextp=1;
   mj =_ma[_inext]-_ma[_inextp];
   if (mj < _MZ) mj += _MBIG;
   _ma[_inext]=mj;
  
#endif
   
#ifdef _XOROSHIRO_
   s[0]=Seed; s[1]=Seed+1234567;
   uint64_t first_rand=xoroshiro128plus(); 
#endif
   
#ifdef _XOSHIRO_
   s[0]=Seed; s[1]=Seed+123456;
   s[2]=Seed+2*123456; s[3]=Seed+3*123456;
   uint32_t first_rand=xoshiro128plus(); 
#endif
   
#ifdef _XORSHIFT_
   uint64_t first_rand=xorshift128plus(); 
#endif
   
   
updateRanNums();
updateLogRanNums();   
   
}

my_double GlobalRandom::PoissonDistributionValue(my_double xm)
{
  my_double em,t,y;

  if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= NormalizedFlatDistributionValue();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*NormalizedFlatDistributionValue());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (NormalizedFlatDistributionValue() > t);
	}
  return em;
}

//#ifdef _RAN_VECTORIZE_

#ifdef _ICSI_LOG_
void GlobalRandom::initializeLookUpTables() {    
   #ifndef _ICC_
   if (!GlobalRandom::LOOKUP_TABLE2) GlobalRandom::LOOKUP_TABLE2 = new float[(int)(pow(2.0,(int)GlobalRandom::logPrecision))];
   fill_icsi_log_table2(GlobalRandom::logPrecision,GlobalRandom::LOOKUP_TABLE2);
   if (!GlobalRandom::LOOKUP_TABLE23) GlobalRandom::LOOKUP_TABLE23 = new float[(int)(pow(2.0,23))];
   fill_icsi_log_table2(23,GlobalRandom::LOOKUP_TABLE23);
   #endif
}
#endif

void GlobalRandom::updateRanNums() {
#ifdef _RAN3_
   for (int i=0; i<numRanNumToKeep; ++i) {
      if (++_inext == 56)  _inext=1;
      if (++_inextp == 56) _inextp=1;
      register long mj =_ma[_inext]-_ma[_inextp];
      if (mj < _MZ) mj += _MBIG;
      _ma[_inext]=mj;
      ranArray[i]=mj;
   }
   for (int i=0; i<numRanNumToKeep; ++i) {
      ranArray[i]*=_FAC;
   }
#endif
   for (int i=0; i<numRanNumToKeep; ++i) {
#ifdef _XOROSHIRO_
      ranArray[i]= xoroshiro128plus_real();
#endif
      
#ifdef _XORSHIFT_
      ranArray[i]= xorshift128plus_real();
#endif
 
#ifdef _XOSHIRO_
//      //Debug
//      float temp=xoshiro128plus_real();
//      ranArray[i]= temp;
      ranArray[i]= xoshiro128plus_real();
#endif     
   }
      
   numRanGiven=0; 
}

//Loro_06_04-17 New version with some for loop vectorization
void GlobalRandom::updateLogRanNums() {
   
//   //Debug
//   my_string fileName("LogRandTable_");
//   fileName+=this->_Seed;
//   fileName+=(char*)".txt";
//   ofstream debugRanNumbers(fileName.c_str(), ios::app);
//   debugRanNumbers << "IN:  _log_inext " << _log_inext << "\t_log_inextp " << _log_inextp << "\t_log_ma[_log_inext] "  << _log_ma[_log_inext] << "\t_log_ma[_log_inextp] "  << _log_ma[_log_inextp] << "\n"; 
   
#ifdef _RAN3_
   #ifdef _ICSI_LOG_
      if (GlobalRandom::logPrecision==23) {
         for (int i=0; i<numRanNumToKeep; ++i) {
            if (++_inext == 56)  _inext=1;
            if (++_inextp == 56) _inextp=1;
            register long mj =_ma[_inext]-_ma[_inextp];         
            if (mj < _MZ) mj += _MBIG;
            _ma[_inext]=mj;
   //         //Loro_25_06_15 To avoid log(zero))
   //         if (mj <= _MZ) mj += _MBIG_1;
   //         if(mj<=0.0) {
   //            bool bugFound=true;
   //         }
    #ifndef _ICC_
//            logRanArray[i]=icsi_log_v2_23(mj,GlobalRandom::LOOKUP_TABLE2)+_LOG_FAC;
            logRanArray[i]=icsi_log_v2_23(mj,GlobalRandom::LOOKUP_TABLE2);
    #else
//            logRanArray[i]=log(mj)+_LOG_FAC;
            logRanArray[i]=log(mj);
    #endif  //_ICC_     

//         debugRanNumbers << logRanArray[i] << "\n";
          } 
      }
      else {
         float curlogRand; 
         for (int i=0; i<numRanNumToKeep; ++i) {
            if (++_inext == 56)  _inext=1;
            if (++_inextp == 56) _inextp=1;
            register long mj =_ma[_inext]-_ma[_inextp];         
            if (mj < _MZ) mj += _MBIG;
      //      if (mj < _MZ) {
      //         mj += _MBIG;
      //      }
   //         //Loro_25_06_15 To avoid log(zero))
   //         if (mj <= _MZ) mj += _MBIG_1;
   //         if(mj<=0.0) {
   //            bool bugFound=true;
   //         }
            _ma[_inext]=mj;
    #ifndef _ICC_
//            curlogRand=icsi_log_v2(mj,GlobalRandom::LOOKUP_TABLE2,GlobalRandom::logPrecision)+_LOG_FAC;
            curlogRand=icsi_log_v2(mj,GlobalRandom::LOOKUP_TABLE2,GlobalRandom::logPrecision);
    #else
//            curlogRand=log(mj)+_LOG_FAC;
            curlogRand=log(mj);
    #endif //_ICC_
//            logRanArray[i]=(curlogRand>0 ? 0 : curlogRand);
            logRanArray[i]=(curlogRand>-_LOG_FAC ? -_LOG_FAC : curlogRand);
            
//            debugRanNumbers << logRanArray[i] << "\n";
         } 
      }

#else //_ICSI_LOG_ not defined
      //#pragma omp simd
      for (int i=0; i<numRanNumToKeep; ++i) {
      if (++_inext == 56)  _inext=1;
      if (++_inextp == 56) _inextp=1;
      long mj =_ma[_inext]-_ma[_inextp];
   //      if (mj < _MZ) {
   //         mj += _MBIG;
   //      }
         //Loro_25_06_15 To avoid log(zero))
         if (mj <= _MZ) mj += _MBIG_1;
      _ma[_inext]=mj;
         logRanArray[i]=log(mj)+_LOG_FAC;
      }
   #endif //_ICSI_LOG_ not defined

   for (int i=0; i<numRanNumToKeep; ++i) {
      logRanArray[i]+=_LOG_FAC;
//      //Loro_05_11_18 Debug
//      if (logRanArray[i]>0) {
//         cerr << "Invalid logRan number : " << logRanArray[i] << "\n";
//      }
   }   
#else //_RAN3_ not defined
   #ifdef _ICSI_LOG_
      if (GlobalRandom::logPrecision==23) {
         for (int i=0; i<numRanNumToKeep; ++i) {
            #ifndef _ICC_
               #ifdef _XOROSHIRO_
               logRanArray[i]=icsi_log_v2_23(xoroshiro128plus(),GlobalRandom::LOOKUP_TABLE2);
               #endif
               #ifdef _XOSHIRO_
               logRanArray[i]=icsi_log_v2_23(xoshiro128plus(),GlobalRandom::LOOKUP_TABLE2);
               #endif
               #ifdef _XORSHIFT_
               logRanArray[i]=icsi_log_v2_23(xorshift128plus(),GlobalRandom::LOOKUP_TABLE2);
               #endif
            #else
               logRanArray[i]=log(xoroshiro128plus);
            #endif         
         }  
//         #pragma omp simd         
//         for (int i=0; i<numRanNumToKeep; ++i) {
//            logRanArray[i]+=log_inv_max_uint64_;
//         }
      }
      else {
         for (int i=0; i<numRanNumToKeep; ++i) {
            #ifndef _ICC_
               #ifdef _XOROSHIRO_
               logRanArray[i]=icsi_log_v2(xoroshiro128plus(),GlobalRandom::LOOKUP_TABLE2,GlobalRandom::logPrecision);
               #endif
               #ifdef _XOSHIRO_
               logRanArray[i]=icsi_log_v2(xoshiro128plus(),GlobalRandom::LOOKUP_TABLE2,GlobalRandom::logPrecision);
               #endif
               #ifdef _XORSHIFT_
               logRanArray[i]=icsi_log_v2(xorshift128plus(),GlobalRandom::LOOKUP_TABLE2,GlobalRandom::logPrecision);
               #endif
            #else
               #ifdef _XOROSHIRO_
               logRanArray[i]=log(xoroshiro128plus());
               #endif
               #ifdef _XOSHIRO_
               logRanArray[i]=log(xoshiro128plus());
               #endif
               #ifdef _XORSHIFT_
               logRanArray[i]=log(xorshift128plus());
               #endif
            #endif
         } 
      }
           
      #pragma omp simd     
      for (int i=0; i<numRanNumToKeep; ++i) {
         logRanArray[i]+=log_inv_max_uint64_;
         //Loro_05_11_18 Debug
//         if (logRanArray[i]>0) {
//            cerr << "Invalid logRan number : " << logRanArray[i] << "\n";
//         }
      }
   #else 
      logRanArray[i]=log(xoroshiro128plus_real());
   #endif
      
#endif      
   numLogRanGiven=0;
}

//#endif
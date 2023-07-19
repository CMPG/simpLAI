#ifndef chrom_painter_h              // Sentry, use file only if it's not already included.
#define chrom_painter_h

#include "arrays.h"
#include "cstring.h"
#include <iostream>

class TDef {
   public:
      int n1, n2, n3, n4, nadm, NPop1, NPop2, NPop3, NPop4, NAnc, tADm, tDiv;
      TDef() {
         n1=n2=n3=n4=nadm=NPop1=NPop2=NAnc=tADm=tDiv=0;
      }  
      TDef& operator=(const TDef & TD) {
         if (this!=&TD) {
            n1=TD.n1;
            n2=TD.n2;
            n3=TD.n3;
            n4=TD.n4;
            nadm=TD.nadm;
            NPop1=TD.NPop1;
            NPop2=TD.NPop2;
            NPop4=TD.NPop3;
            NPop2=TD.NPop4;
            NAnc=TD.NAnc;
            tADm=TD.tADm;
            tDiv=TD.tDiv;
         }
         return (*this);
      }
      
   friend istream & operator>>(istream& is, TDef& def);
};

   istream & operator>>(istream& is, TDef& def);
   
#endif
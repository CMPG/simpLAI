/*//////////////////////////////////////////////////////////////////////////////

Copyright Laurent Excoffier, CMPG, University of Berne, August 2022
E-mail: laurent.excoffier@unibe.ch

The use of the following source code in only allowed if it is kept intact
with its original banner
 * 
 * 22.08.22  I have noticed that the windows of inferred segments tended to be shifted to the right as compared to the true segments
 *           One would need to revise the way window beginnings and ends are assigned
 *           If segment changes its assignment, the overlapping part with the previous segment
 *           Reports center of the segment instead of the beginning? What about beginning and end of chromosome?
 * 
 * 18.10.22 I have added the computation of FSTs in windows of size winInc
 * 
 * 18.10.22 Compute No. of differences in winInc size windows and then add them to compute them in winSize windows: DONE
 * 
 * 19.10.22 Remove private mutations from the computation of distances as it could lead to 
 *          wrong inferences regarding the evolutionary proximity of sequences: DONE (but does not bring much)
 * 
 * 18.10.22 Trying to add assigned windows of admixed chromosomes to reference panels. 
 *          The idea is to have a first pass where chromosome segments are assigned to a given reference only on the basis of the reference panel 
 *          and then having a second pass where the assigned admixed segments are used as additional references to reassign the other admixed segments
 * 
 * 18.12.22 Plan is to introduce a recombination path among source chromosomes that best explains the observed chromosome. 
 *          The aim is to find within each window a combination of source chromosome segments that can best the target admixed chromosome.
 * 
 *          The algorithm would go like, for each target chromosome
 *          - Remove all private mutations from the target chromosome
 *          - For each remaining SNP of the target chromosome and for each source population:
 *             - Examine all source chromosome with a SNP matching that on the target chromosome, and select the source chrom. with the least mismatches with the current target chromosome. 
 *               Record this number. 
 *             - In case of ex-aequo, if possible keep the same chromosome as for the previous SNP.
 *             - Jump to the closest chromosome
 *             - Record jump history among source chromosomes.
 * 
 * 09.02.23 Plan to extend the program to test for the admixture of up to four populations
 * 10.02.23 There is a BUG, as matching polym sites are not found on pop3 and pop4 
 *          in function void computeMinNumDiffs(bool * admChrom, ...
 *          NEED TO RECONSIDER COMPUTATION OF DIFFERENCES FOR SHORTEST RECOMBINATION PATH
 * 
 *          DONE on 13.02!
 * 
 * 20.02.23 Corrected a bug when there was only 3 populations and not 4 (n4==0)
 * 
 * 09.03.23 Added the possibility to define the non-recombining segment length (in no. of snps) 
 * 

//////////////////////////////////////////////////////////////////////////////*/


#include "config.h"
#include "anyoption.h"

#include <sstream>
#include "cstring.h"
#include <ctime>
#include "UniqueRandom.h"

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include "simpLAI.h"

using namespace std;

int  RETURN_OKAY = 0;
int  RETURN_PROBLEM = 1;

char ** newArg = NULL;
char strTime[80];

typedef MY_TArrayAsVector<TDef> TSimDefinitions;
typedef MY_TArrayAsVector<int> TIntArray;

extern double gammln(double xx);
typedef bool** TGenArray;

typedef std::vector<int> TPolymPosArray;

//------------------------------------------------------------------------------
istream & operator >> (istream& is, TDef& def) {
   double pAdm;
   int NAdm;
   is >> def.n1 >> def.n2 >> def.n3 >> def.n4 >> def.nadm >> def.NPop1 >> def.NPop2 >> NAdm >> def.NAnc >> def.tADm >> pAdm >> def.tDiv;
   return is;   
}
//------------------------------------------------------------------------------
//Reading def file
void readDefFile(const my_string &fileName, TSimDefinitions & SimDefs) {  
   my_string curLine;
   ifstream ifs(fileName.c_str());
   //Reading header (SS1	SS2	SS3	NPOP1	NPOP2	NADM      NANC      tadm	padm      tdiv)
   curLine.read_line(ifs);
   while (ifs) {
      curLine.read_line(ifs);
//      curLine.remove_blanks();
//      curLine.remove_comments();
      if (curLine.length()) {         
         stringstream ss;   
         ss << curLine;
         TDef curDef;
         ss >> curDef;
         SimDefs.Add(curDef); 
      }     
   }   
}
//------------------------------------------------------------------------------
//Reading genotype file
void readGenotypeFile(const my_string &fileName, TGenArray & GA, int*  pos, int numChroms) {  
   my_string curLine;
   ifstream ifs(fileName.c_str()); 
  
   curLine.read_line(ifs);
   int curSite=0;
   int chrom;
   char anc_all, der_all;
   while (ifs) {
//      curLine.read_line(ifs);
//      curLine.remove_blanks();
//      curLine.remove_comments();
//      if (curLine.length()) {         
//         stringstream ss;   
//         ss << curLine;
//         int chrom;
//         
//         char anc_all, der_all;
//         
//         ss >> chrom >> pos[curSite] >> anc_all >> der_all;
//         for (int nc=0; nc<numChroms; ++nc) {
//            ss >> GA[nc][curSite];
//         }         
//         curSite++;
//         if (curSite%10000==0) cout << "\t\tReading line " << curSite <<  "\r";
//      }  
      ifs >> chrom >> pos[curSite] >> anc_all >> der_all;
      for (int nc=0; nc<numChroms; ++nc) {
         ifs >> GA[nc][curSite];
      } 
      curLine.read_line(ifs);
      curSite++;
      if (curSite%10000==0) cout << "\t\tReading line " << curSite <<  "\r";
      
   } 
   cout << "\t\tReading line " << curSite <<  "\n";
   return;
}
//------------------------------------------------------------------------------
int getLineNumber(char* fn) {
   ifstream myfile(fn);
   my_string curLine;
   int number_of_lines=0;
   if(myfile.is_open()){
      while(!myfile.eof()){
         curLine.read_line(myfile);
//         curLine.remove_blanks();
//         curLine.remove_comments();
         if (curLine.length()) {     
            number_of_lines++;
         }
      }
      myfile.close();
   }
   return number_of_lines;
}
//------------------------------------------------------------------------------
int computeNumDiff(bool *h1, bool *h2, int size) {
   int numDiff=0;
   for (int i=0; i<size; ++i) {
      numDiff+=(*h1++)^(*h2++);
   }
   return numDiff;
}
//------------------------------------------------------------------------------
void computeMismatchDistribution(bool * admChrom, bool ** targetChroms, int beg, int end, int targetSize, TIntArray & mmD) {
   for (int i=0; i<targetSize; ++i)  {
      int size=end-beg+1;
      int curNumDiff=computeNumDiff(&admChrom[beg], &targetChroms[i][beg], size);
      mmD.Add(curNumDiff);
   }
}
//------------------------------------------------------------------------------
double computePiX(bool ** X, int beg, int end, int nX) {
   if (nX<=1) return 0;
   double sum=0;
   int size=end-beg+1;
   for (int i=0; i<nX; ++i)  {
      for (int j=0; j<i; ++j)  {
         sum+=computeNumDiff(&X[i][beg], &X[j][beg], size);
      }
   }
   double piX=sum/(nX*(nX-1.0)/2.0);
   return piX;
}
//------------------------------------------------------------------------------
double computePiXY(bool ** X, bool ** Y, int beg, int end, int nX, int nY) {
   double sum=0;
   int size=end-beg+1;
   for (int i=0; i<nX; ++i)  {
      for (int j=0; j<nY; ++j)  {
         sum+=computeNumDiff(&X[i][beg], &Y[j][beg], size);
      }
   }
   double piXY=sum/(nX*nY);
   return piXY;
}
//------------------------------------------------------------------------------
double computeMeanDiff(TIntArray & curArray) {
   double sum=0;
   int size=curArray.GetItemsInContainer();
   if (!size) return 0;
   for (int i=0; i<size; ++i) {
      sum+=curArray[i];
   }
   return 1.0*sum/size;
}
//------------------------------------------------------------------------------
double computeMeanDiff(TIntArray & curArray, int size) {
   double sum=0;
   for (int i=0; i<size; ++i) {
      sum+=curArray[i];
   }
   return 1.0*sum/size;
}
//------------------------------------------------------------------------------
double computeMinDiff(TIntArray & curArray) {
   double min=1e50;
   int size=curArray.GetItemsInContainer();
   if (!size) return 0;
   for (int i=0; i<size; ++i) {
      if (curArray[i]<min) min=curArray[i];
   }
   return min;
}
//------------------------------------------------------------------------------
double computeMinDiff(TIntArray & curArray, int size) {
   double min=1e50;
   for (int i=0; i<size; ++i) {
      if (curArray[i]<min) min=curArray[i];
   }
   return min;
}
//------------------------------------------------------------------------------
double computeProbaSdiffIfAdmixed(int S, double & theta1, double & theta0, double & tA, double tD) {
   double B=expl(-(tD-tA)/theta1)*expl(-tD);
   double log_theta1_ratio=log(theta1/(theta1+1));
   double log_theta0_ratio=log(theta0/(theta0+1));
   double sum1=0, sum2=0;
   for (int j=0; j<=S; j++) {
      double lnfactj=gammln(j+1);
      sum1 += (expl(j*log(tA)-lnfactj-tA) -expl(j*log(tD)-lnfactj)*B)*expl((S-j)*log_theta1_ratio)/(theta1+1);
      sum2 += expl(j*log(tD)-lnfactj)*B*expl((S-j)*log_theta0_ratio)/(theta0+1);
   }
   return sum1+sum2;
}
//------------------------------------------------------------------------------
double computeProbaSdiffIfNOTAdmixed(int S, double & theta0, double tD) {
   double log_theta0_ratio=log(theta0/(theta0+1));
   double sum=0;
   for (int j=0; j<=S; j++) {
      double lnfactj=gammln(j+1);
      sum += expl(j*log(tD)-lnfactj)*expl((S-j)*log_theta0_ratio)/(theta0+1);
   }
   double prob=expl(-tD)*sum;
   return prob;
}
//------------------------------------------------------------------------------
//double logProduct(double* a, int size) {
//   double lp=0;
//   for (int i=0; i<size; ++i) lp+=log(a[i]);
//   return lp;
//}
//------------------------------------------------------------------------------
void findSitePositionsForCurrentWindows(int* positions, int numSites, int winBeg, int winEnd, int & posBeg, int & posEnd ) {
 for (int i=0; i<numSites; ++i) {
    if (positions[i]>=winBeg) {
       posBeg=i;
       for (int j=i+1; j<numSites; ++j) {
          if (positions[j]>=winEnd) {
             posEnd=j-1;
             return;
          }
       }
    }
 }  
}
//------------------------------------------------------------------------------
int assignNewCmdLine(const my_string & newCmdLine, int & argc, char ** & argv) {
   
//   //debug 
//   cout << endl << "Old argv: " << endl;
//   for (int i=0; i<argc; ++i) {
//      cout << "Arg. " << i+1 << " : " << argv[i] << "\n";
//   }
   
   #ifdef _GCC_
   stringstream ss;
   #else
   strstream ss;
   #endif
   ss << newCmdLine;
   vector<string> strVec;
   while (ss) {
      string curArg;
      ss >> curArg;
      if (curArg.length()) strVec.push_back(curArg);
   }
   int numArg=strVec.size()+1;
   newArg = new char * [numArg];
   newArg[0] = new char[strlen(argv[0])+1];
   strcpy(newArg[0], argv[0]);
   for (int i=1; i<numArg; ++i) {
//      cout << strVec[i-1].c_str() << endl;
      newArg[i] = new char[strlen(strVec[i-1].c_str())+1];
      strcpy(newArg[i], strVec[i-1].c_str());
//      cout << newArg[i] << endl;
   }  
   //if (argv) delete argv;
   argv=newArg;
   argc=numArg;
   //debug 
//   cout << "New argv: " << endl;
//   for (int i=0; i<argc; ++i) {
//      cout << "Arg. " << i+1 << " : " << argv[i] << endl;
//   }
   if (argc==1) {
      cerr << "No options specified on command line" << endl;
      return 0;
   }
   return 1;
}
//------------------------------------------------------------------------------
bool isInteger(const std::string & s) {
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}
//------------------------------------------------------------------------------
void computeNumDiffs(bool * admChrom, bool ** targetChroms, int beg, int end, int targetSize, int* numDiffs) {
   for (int i=0; i<targetSize; ++i)  {
      int size=end-beg+1;
      numDiffs[i]=computeNumDiff(&admChrom[beg], &targetChroms[i][beg], size);
   }
}
//------------------------------------------------------------------------------
void computeMinNumDiffs(bool * admChrom, bool ** targetChroms, int beg, int end, int targetSize, int & minNumDiffs, int & closestSourceChrom, bool requestMatchingPolymSite){
   minNumDiffs=INT_MAX;
   int curDiffs=0;
   bool foundMatchingPolymSite=false;
   for (int i=0; i<targetSize; ++i)  {
      if (!requestMatchingPolymSite || admChrom[end]==targetChroms[i][end]) {
         foundMatchingPolymSite=true;
         int size=end-beg+1;
         curDiffs=computeNumDiff(&admChrom[beg], &targetChroms[i][beg], size);
         if (curDiffs<minNumDiffs) {
            minNumDiffs=curDiffs;
            closestSourceChrom=i;
         }
      }
   }
   if (!foundMatchingPolymSite) {
      minNumDiffs=-1;
      closestSourceChrom=-1;
   }
}   
//------------------------------------------------------------------------------
void computeMinNumDiffs_new(bool * admChrom, bool ** targetChroms, int beg, int winSize, int targetSize, int & minNumDiffs, int & closestSourceChrom){
   minNumDiffs=INT_MAX;
   int curDiffs=0;
   for (int i=0; i<targetSize; ++i)  {
      curDiffs=computeNumDiff(&admChrom[beg], &targetChroms[i][beg], winSize);
      if (curDiffs<minNumDiffs) {
         minNumDiffs=curDiffs;
         closestSourceChrom=i;
      }
   }
}   
//------------------------------------------------------------------------------
/* Function that needs to
  1) Find source chromosome that carry the current derived allele
  2) Identify the source chromosome with the least derived allele from the last position, not counting this last position
  3) Return this number of mismatches
*/
int findShortestPathFromLastPolymPos(bool * admChrom, bool ** targetChroms, int lastPos, int curPos, int targetSize, bool ignoreCurrentPolymPos) {   
   int minDiff=INT_MAX; 
   int size=curPos-lastPos;
   bool foundSharedDerivedAllele=false;
   for (int i=0; i<targetSize; ++i) { 
      //Only look at source chromosomes sharing the derived allele with the admixed chrom at the current position
      if (ignoreCurrentPolymPos || admChrom[curPos]==targetChroms[i][curPos]) { 
         foundSharedDerivedAllele=true;
         //Compute number of differences since the last polym pos. on the admixed chromosome (not counting this last polym pos)
         int numDiffs=computeNumDiff(&admChrom[lastPos+1], &targetChroms[i][lastPos+1], size); 
         if (numDiffs<minDiff) minDiff=numDiffs;
      } 
   }
   if (foundSharedDerivedAllele) return minDiff; 
   //If derived allele is not found among any source chromosome
   else return -1;   
}

//------------------------------------------------------------------------------
int computeDiffsWithClosestSourceChrom(bool * admChrom, bool ** targetChroms, int beg, int end, int targetSize) {   
   int minDiff=INT_MAX; 
   int size=end-beg;
   for (int i=0; i<targetSize; ++i) { 
      //Compute number of differences since the last polym pos. on the admixed chromosome (not counting this last polym pos)
      int numDiffs=computeNumDiff(&admChrom[beg], &targetChroms[i][beg], size); 
      if (numDiffs<minDiff) minDiff=numDiffs;
   }
   return minDiff;   
}
//------------------------------------------------------------------------------
int findShortestRecombPath(bool * admChrom, bool ** targetChroms, TPolymPosArray & polymPositions, int beg, int end, int targetSize) {
   int numDiffs=0;
   //Go over all polym. positions of current admixed chromosome
   int lastPos=0;
   int curPos=beg;
   bool polymPosInCurrentSegment=false;
   for (int i=0; i<polymPositions.size(); ++i) {
      curPos=polymPositions[i];
      //Now compute the shortest possible path across all source chromosomes 
      if (curPos >=beg && curPos<=end) {
         polymPosInCurrentSegment=true;
         bool ignoreCurrentPolymPos=true; //Flag to explore all source chromosomes and not only those having the same allele at the polymorphic position
         int curNumDiffs=findShortestPathFromLastPolymPos(admChrom, targetChroms, lastPos, curPos, targetSize, ignoreCurrentPolymPos); 
         if (curNumDiffs>-1) { 
            numDiffs+=curNumDiffs;
            lastPos=curPos;
         } else { 
           //if there is no shared derived allele in the current population, we do not do anything and we examine the next polym position 
         }
      }
   }
   //Now handle the case where there is no polym pos in this segment and 
   //if there are still some mutations between the current polym pos and the end of the segment
   if (!polymPosInCurrentSegment) {
      curPos=beg;
   }
   numDiffs+=computeDiffsWithClosestSourceChrom(admChrom, targetChroms, curPos, end, targetSize); 
   return numDiffs;   
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

double mu=1.25e-8;
bool removeSingletons=false;
bool removePrivateAllelesOnAdmixedChroms=true;

int main(int argc, char *argv[]) {
   
   //Handling of command line --------------------------------------------------
   AnyOption opt;
   
   opt.addUsage("");
   char versionSpec[400];
   int MAJOR_VERSION=0;
   int MINOR_VERSION=9.1;
   my_string RELEASE_DATE = (char*) "21.03.23";
   sprintf(versionSpec, " genomePainter ver %i.%i - (%s)", MAJOR_VERSION, MINOR_VERSION, RELEASE_DATE.c_str()); 
   opt.addUsage(versionSpec);
   opt.addUsage("");
   opt.addUsage("Usage:");
   opt.addUsage("");
   opt.addUsage(" -h  --help               : prints this usage ");
   opt.addUsage(" -g  --genFile test.gen   : name of genotype file (phased)");
   opt.addUsage("     --ss1  10            : number of chromosome in source 1");
   opt.addUsage("     --ss2  10            : number of chromosome in source 2");
   opt.addUsage("     --ss3  10            : number of chromosome in source 3");
   opt.addUsage("     --ss4  10            : number of chromosome in source 4");
   opt.addUsage("     --ssa  20            : number of chromosome in admixed sample");
   opt.addUsage(" -l  --genomeLength  1e8  : total length of genome");
   opt.addUsage(" -s  --winSize       1e6  : size of sliding window (optional: default = 1Mb)");
   opt.addUsage(" -i  --winInc        5e5  : increment for sliding window (optional: default = 500Kb)");
   opt.addUsage(" -n  --numPolSites   1e4  : number of polymorphic sites on which to search for");
   opt.addUsage("                            the shortest recombination path (default = 10,000 sites) ");
   opt.addUsage(" -m  --incRec        1e3  : increment for sliding window when searching for");
   opt.addUsage("                            the shortest recombination path (default = 1,000 sites) ");
   opt.addUsage(" -t                  5    : no. of linked polymorphic sites between potential");
   opt.addUsage("                            recombination breakpoints (default = 5 sites) ");
//   opt.addUsage(" -p  --penRec        0    : penalty for recombination event when searching for");
//   opt.addUsage("                            the shortest recombination path (default = 0 ) ");
//   opt.addUsage(" -q  --penPol        1    : penalty for lack of matching polymorphic site when searching for");
//   opt.addUsage("                            the shortest recombination path (default = 1 ) ");
   

   opt.setFlag("help", 'h');
   opt.setOption("genFile", 'g'); 
   opt.setOption("winSize", 's'); 
   opt.setOption("incWin", 'i');  
   opt.setOption("ss1");
   opt.setOption("ss2");
   opt.setOption("ss3");
   opt.setOption("ss4");
   opt.setOption("ssa"); 
   opt.setOption("genomeLength", 'l');  
   opt.setOption("numPolSites", 'n'); 
   opt.setOption("incRec", 'm'); 
   opt.setOption('t'); 
//   opt.setOption("penRec", 'p');  
//   opt.setOption("penPol", 'q');  
   
   /* go through the command line and get the options  */
   opt.processCommandArgs(argc, argv);
   
   if (argc>1) {
      cout << "\ngenomePainter was invoked with the following command line arguments:\n";
      for (int i=0; i<argc ; ++i) {
         cout << argv[i] << " ";
      }
      cout << "\n";
   }

   if (!opt.hasOptions()) { /* print usage if no options */
      opt.printUsage();
      return 1;
   }
   
   if (opt.getFlag("help") || opt.getFlag('h'))  {
      opt.printUsage();
      return 1;
   }
      
   my_string genFileName;
   my_string genericName;
   my_string admFileName;
   
   //Get genotype file name
   if (opt.getValue('g') != NULL || opt.getValue("genFile") != NULL) {
      genFileName = (char*) opt.getValue('g');
      genericName= genFileName;
      genericName.remove_extension();
   }
   else {
      cerr <<  "Genotype file name not provided. Exiting program";      
      opt.printUsage();
      return 1;
   }
   
   int chromLength;
   //Get genome length
   if (opt.getValue('l') != NULL || opt.getValue("genomeLength") != NULL) {
      chromLength = atof(opt.getValue('l'));
   }
   else {
      cerr <<  "Chromosome length not specified. Exiting program";      
      opt.printUsage();
      return 1;
   }
   
   //Default values for window size and increment size
   int winSize=1e6, winInc=5e5;
   
   //Read optional window size
   if (opt.getValue('s') != NULL || opt.getValue("winSize") != NULL) {
      winSize = atof(opt.getValue('s'));
   }
   
   //Read optional window increment size
   if (opt.getValue('i') != NULL || opt.getValue("incWin") != NULL) {
      winInc = atof(opt.getValue('i'));
   }
   
   TDef curDef;  
   
   if (opt.getValue("ss1") != NULL) {
      if (isInteger(opt.getValue("ss1"))) {
         curDef.n1=atof(opt.getValue("ss1"));
      } else {
        cerr << "Bad value for source population 1 sample size. Need an integer value. Exiting program!";
         opt.printUsage();
         return 1;
      }
   } else {
         cerr << "Sample size for source population 1 is not defined. Exiting program!";
         opt.printUsage();
         return 1;
   }

   if (opt.getValue("ss2") != NULL) {
      if (isInteger(opt.getValue("ss2"))) {
         curDef.n2=atof(opt.getValue("ss2"));
      } else {
        cerr << "Bad value for source population 2 sample size. Need an integer value. Exiting program!";
         opt.printUsage();
         return 1;
      }
   } else {
         cerr << "Sample size for source population 2 is not defined. Exiting program!";
         opt.printUsage();
         return 1;
   }
   
   if (opt.getValue("ss3") != NULL) {
      if (isInteger(opt.getValue("ss3"))) {
         curDef.n3=atof(opt.getValue("ss3"));
      } else {
         cerr << "Bad value for source population 3 sample size. Need an integer value. Exiting program!";
         opt.printUsage();
         return 1;
      }
   } else {
      curDef.n3=0;
   }
   
   if (opt.getValue("ss4") != NULL) {
      if (isInteger(opt.getValue("ss4"))) {
         curDef.n4=atof(opt.getValue("ss4"));
      } else {
         cerr << "Bad value for source population 4 sample size. Need an integer value. Exiting program!";
         opt.printUsage();
         return 1;
      }
   } else {
      curDef.n4=0;
   }

   if (opt.getValue("ssa") != NULL) {
      if (isInteger(opt.getValue("ssa"))) {
         curDef.nadm=atof(opt.getValue("ssa"));
      } else {
        cerr << "Bad value for admixed population sample size. Need an integer value. Exiting program!";
         opt.printUsage();
         return 1;
      }
   } else {
         cerr << "Sample size for admixed population is not defined. Exiting program!";
         opt.printUsage();
         return 1;
   }
   
    int numConsecPolSitesToConsider=10000;
   
   if (opt.getValue('n') != NULL || opt.getValue("numPolSites") != NULL) {
         numConsecPolSitesToConsider=atof(opt.getValue('n'));
   }
   int incNumPolSites=1000;
   if (opt.getValue('m') != NULL || opt.getValue("incRec") != NULL) {
         incNumPolSites=atof(opt.getValue('m'));
   }
   
   int nonRecTract=5;
   if (opt.getValue('t') != NULL) {
         nonRecTract=atof(opt.getValue('t'));
   }
   
//   int recombPenalty=0; // Penalty for a recombination event to be added to number of mutations to compute total no. of mismatches
//   if (opt.getValue('p') != NULL || opt.getValue("penRec") != NULL) {
//         recombPenalty=atof(opt.getValue('p'));
//   }
//   
//   int lackOfPolymSitePenalty=1;
//   if (opt.getValue('q') != NULL || opt.getValue("penPol") != NULL) {
//         lackOfPolymSitePenalty=atof(opt.getValue('q'));
//   }

   //End handling of command line ----------------------------------------------
      
   struct tm *time_now;
   time_t secs_now;
   time(&secs_now);
   time_now = localtime(&secs_now);
   strftime(strTime, 80, "%d/%m/20%y at %H:%M:%S", time_now);
     
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   long randSeed = ((tv.tv_sec ^ tv.tv_usec) ^ getpid()) % 1000000; //Implies only  one million possible seeds!
     

   int numSites=getLineNumber(genFileName.c_str())-1;      
   int numChroms=curDef.n1+curDef.n2+curDef.n3+curDef.n4+curDef.nadm;

   TGenArray curGenotypes = new bool*[numChroms];
   for (int s=0; s<numChroms; ++s) {
      curGenotypes[s] = new bool[numSites]; 
   }  
  
//   TPolymPosArray polymPositions;

   int * positions = new int[numSites];

   int numWin=(chromLength-winSize)/winInc+1;

   int ** seq_origins_from_mean  = new int*[curDef.nadm];
   int ** seq_origins_from_min   = new int*[curDef.nadm];
   int ** seq_origins_from_rec   = new int*[curDef.nadm];
      
   //Lists of adm segments assigned to a given reference panel 
//   int ** index_adm_ref1         = new int*[numWin]; 
//   int ** index_adm_ref2         = new int*[numWin];
   //No. of admixed seq assigned to the two ref panels for each window
//   int * num_adm_in_ref1         = new int[numWin]; 
//   int * num_adm_in_ref2         = new int[numWin];

   for (int na=0; na<curDef.nadm; ++na) {         
      seq_origins_from_mean[na]  = new int[numWin]; 
      seq_origins_from_min[na]   = new int[numWin];
   }
   
//   for (int nw=0; nw<numWin; ++nw) {    
//      index_adm_ref1[nw]         = new int[curDef.nadm]; 
//      index_adm_ref2[nw]         = new int[curDef.nadm];
//      num_adm_in_ref1[nw] =0;
//      num_adm_in_ref2[nw] =0;
//   }
   
   bool ** GenotypesP1=&curGenotypes[0];
   bool ** GenotypesP2=&curGenotypes[curDef.n1];
   bool ** GenotypesP3=&curGenotypes[curDef.n1+curDef.n2];
   bool ** GenotypesP4=&curGenotypes[curDef.n1+curDef.n2+curDef.n3];
   
   int winSizeForFST=winInc;
   
   int numFSTWin=(chromLength-winSizeForFST)/winSizeForFST+1;
   
   double * pi12 = new double[numFSTWin];
   double * pi1 = new double[numFSTWin];
   double * pi2 = new double[numFSTWin];
   double * pi3 = new double[numFSTWin];
   double * pi4 = new double[numFSTWin];
   double * FSTBetweenRefs = new double[numFSTWin];
   
   int numWinInc=(chromLength-winInc)/winInc+1;
   int *** numDiffsRef1 = new int**[numWinInc];
   int *** numDiffsRef2 = new int**[numWinInc];
   int *** numDiffsRef3 = new int**[numWinInc];
   int *** numDiffsRef4 = new int**[numWinInc];   
   
   int ** recPathLengthRef1 = new int*[numWinInc];
   int ** recPathLengthRef2 = new int*[numWinInc];
   int ** recPathLengthRef3 = new int*[numWinInc];
   int ** recPathLengthRef4 = new int*[numWinInc];
   
   //Creating arrays to store number of differences between admixed and reference chromosome segments
   for (int i=0; i<numWinInc; ++i) {
      recPathLengthRef1[i]= new int[curDef.nadm];
      recPathLengthRef2[i]= new int[curDef.nadm];
      recPathLengthRef3[i]= new int[curDef.nadm];
      recPathLengthRef4[i]= new int[curDef.nadm];
      numDiffsRef1[i]= new int*[curDef.nadm];
      numDiffsRef2[i]= new int*[curDef.nadm];
      numDiffsRef3[i]= new int*[curDef.nadm];
      numDiffsRef4[i]= new int*[curDef.nadm];
      for (int j=0; j<curDef.nadm; ++j) {
         numDiffsRef1[i][j] = new int[curDef.n1];
         numDiffsRef2[i][j] = new int[curDef.n2];
         if (curDef.n3) numDiffsRef3[i][j] = new int[curDef.n3]; else numDiffsRef3[i][j]=NULL;
         if (curDef.n4) numDiffsRef4[i][j] = new int[curDef.n4]; else numDiffsRef4[i][j]=NULL;
      }
   }
   
   //Reading genotype file -----------------------------------------------------
   cout << "\tReading genotype file " << genFileName << "...\n";
   readGenotypeFile(genFileName, curGenotypes, positions, numChroms); 
   cout << " \tDONE!\n\n";
   
   //Compute FST in all incremental windows ------------------------------------  
   cout << "\tComputing FSTs in all windows ..."; 
//   int begWin=0, endWin=winSizeForFST, beg, end;  
   #pragma omp parallel for 
   for (int w=0; w<numFSTWin; ++w) { 
//      if (w%10==0 | w==numFSTWin-1) cout << "\t\twindow " << (w+1) << "/" << numFSTWin << "\r"; 
   
      int curWinSize=winSizeForFST;
      int begWin=w*winSizeForFST, endWin=begWin+winSizeForFST, beg, end;
      if (endWin>chromLength) {
         endWin=chromLength;
         curWinSize=endWin-begWin;
      }   
      
      findSitePositionsForCurrentWindows(positions, numSites, begWin, endWin, beg, end);      
      
      pi12[w]=computePiXY(GenotypesP1, GenotypesP2, beg, end, curDef.n1, curDef.n2);
      pi1[w]=computePiX(GenotypesP1, beg, end, curDef.n1);
      pi2[w]=computePiX(GenotypesP2,  beg, end, curDef.n2);
      if (curDef.n3) pi3[w]=computePiX(GenotypesP3, beg, end, curDef.n3);
      if (curDef.n4) pi4[w]=computePiX(GenotypesP4, beg, end, curDef.n4);
      
      if (pi12[w]>0) FSTBetweenRefs[w]=(pi12[w]-(pi1[w]+pi2[w])/2.0)/pi12[w];
      else FSTBetweenRefs[w]=-2;         
      
      begWin+=winSizeForFST;
      endWin+=winSizeForFST;
   }
   cout << " DONE!\n\n";
   
   //Remove singletons ---------------------------------------------------------   
   TGenArray nonSingletonGenotypes = new bool*[numChroms];       
   TGenArray nonPrivateGenotypes = new bool*[numChroms];     
   for (int s=0; s<numChroms; ++s) {
      nonSingletonGenotypes[s] = new bool[numSites]; 
      nonPrivateGenotypes[s] = new bool[numSites]; 
   }    
   int * positionsNonSingletons = new int[numSites];
   int * positionsNonPrivate = new int[numSites]; 
         
   int * curPos;
   TGenArray * curGenot;
   int curNumSites; 
   
   int numPolymSitesInSources=0; 
      
   for (int cases=0; cases<2; ++cases) {
      
      //Let's skip cases without singletons and cases where we look for the shortest recombination path
      if (cases==1 || cases==2) continue;
      
      removePrivateAllelesOnAdmixedChroms=false;
      removeSingletons=false;
      switch (cases) {
         case 1: removeSingletons=true; 
                 break;
      }  
      
      if (removeSingletons)  {
         cout << "\tInferring chromosome origins without singletons ...\n";
      } else {   
         cout << "\tInferring chromosome origins with singletons ...\n";
      }
      if (removeSingletons) {
         int numMonoSites=0;
         int numNonSingletons=0; 
         for (int i=0; i<numSites ; ++i) {
            int numDerivedAlleles=0;
            for (int j=0; j<numChroms; ++j) {
               if (curGenotypes[j][i]) ++numDerivedAlleles;
            }
            if (numDerivedAlleles>1 && numDerivedAlleles<numChroms-1) {
               positionsNonSingletons[numNonSingletons]=positions[i];
               for (int j=0; j<numChroms; ++j) {
                  nonSingletonGenotypes[j][numNonSingletons]=curGenotypes[j][i];
               }
               numNonSingletons++;
            }
            else {
               if (numDerivedAlleles==0 || numDerivedAlleles==numChroms) ++numMonoSites;
            } 
         }      
         curPos=positionsNonSingletons;
         curGenot=&nonSingletonGenotypes;
         curNumSites=numNonSingletons;
         cout << "\t\tFound "<< numNonSingletons << " non-singleton sites, " << (numSites-numNonSingletons-numMonoSites) << " singleton sites and " << numMonoSites << " monomorphic sites\n";
      } else {
         curPos=positions;
         curGenot=&curGenotypes;
         curNumSites=numSites;
      }
      GenotypesP1=&(*curGenot)[0];
      GenotypesP2=&(*curGenot)[curDef.n1];
      GenotypesP3=&(*curGenot)[curDef.n1+curDef.n2];
      GenotypesP4=&(*curGenot)[curDef.n1+curDef.n2+curDef.n3];
         
      //Computing number of pairwise differences in small windows -----------------     
      cout << "\t\tComputing number of pairwise differences in small windows ..."; 
      #pragma omp parallel for 
      for (int w=0; w<numWinInc; ++w) {

         int begWin=w*winInc, endWin=begWin+winInc, beg, end;
         if (endWin>chromLength) {
            endWin=chromLength;
         }     

         findSitePositionsForCurrentWindows(curPos, curNumSites, begWin, endWin, beg, end);

         for (int j=0; j<curDef.nadm; ++j) {
            int curAdmChrom=curDef.n1+curDef.n2+curDef.n3+curDef.n4+j;
            computeNumDiffs((*curGenot)[curAdmChrom], GenotypesP1, beg, end, curDef.n1, numDiffsRef1[w][j]);
            computeNumDiffs((*curGenot)[curAdmChrom], GenotypesP2, beg, end, curDef.n2, numDiffsRef2[w][j]);
            if (curDef.n3) computeNumDiffs((*curGenot)[curAdmChrom], GenotypesP3, beg, end, curDef.n3, numDiffsRef3[w][j]);
            if (curDef.n4) computeNumDiffs((*curGenot)[curAdmChrom], GenotypesP4, beg, end, curDef.n4, numDiffsRef4[w][j]);
         }   
      }
      cout << " DONE!\n\n";
      
      //------------------------------------------------------------------------
      
//      cout << "\tInferring chromosome origins from shortest recombination paths ... ";
//      
//      int numSourceChroms=curDef.n1+curDef.n2+curDef.n3+curDef.n4;
//      
//      //Step 1: Identify the non-private mutations in admixed chromosome
//         
//      for (int i=0; i<numSites ; ++i) {
//         int numDerivedAllelesInSources=0;
//         //Count derived allele frequency in source chromosomes
//         for (int j=0; j<numSourceChroms; ++j) {
//            if (curGenotypes[j][i]) ++numDerivedAllelesInSources;
//         } 
//         bool monoInSources=false;
//         if ((numDerivedAllelesInSources==0 || numDerivedAllelesInSources==numSourceChroms))  monoInSources=true;
//         double curFreqInSource=numDerivedAllelesInSources/numSourceChroms;
//         if (!monoInSources) {
//            //Then this position must be polymorphic in all admixed individuals
//            polymPositions.push_back(i);
//            numPolymSitesInSources++;
//         }
//      }   
//      
//      //Pretend that the last site is polymorphic in sources
//      polymPositions.push_back(numSites-1);
//      
//      int numRecSegments=(numPolymSitesInSources-numConsecPolSitesToConsider)/incNumPolSites+1;
//      for (int na=0; na<curDef.nadm; ++na) { 
//         seq_origins_from_rec[na]   = new int[numRecSegments+1];
//      }
//      
//      //Loop over all admixed chromosomes
//      #pragma omp parallel for
//      for (int j=0; j<curDef.nadm; ++j) {
//         
////         cout << "\t\tExamining individual " << j+1 << "/" << curDef.nadm <<"\r";
//      
//         int curAdmChrom=curDef.n1+curDef.n2+curDef.n3+curDef.n4+j;         
//         
//         int startSite=0;
//         int endSite=numConsecPolSitesToConsider;
//         for (int rs=0; rs<numRecSegments+1; ++rs) {
//            int lastSite=polymPositions[startSite];
//            int lastCloseChromP1=-1,lastCloseChromP2=-1, lastCloseChromP3=-1,lastCloseChromP4=-1;
//            int totalScoreP1=0, totalScoreP2=0, totalScoreP3=0, totalScoreP4=0;
//            int numRecP1=0, numRecP2=0, numRecP3=0, numRecP4=0 ;
//            //Loop over all polymorphic sites in the current segment
//            for (int p=startSite; p<endSite; ++p) {
//               int curSite=polymPositions[p];
//               //Examine each source pop in turn and find the source chromosome that is closest to the admixed chromosome
//               int curScoreP1, curScoreP2,  curScoreP3, curScoreP4;
//               int closestChromP1, closestChromP2, closestChromP3, closestChromP4;
//               bool requestMatchingPolymSite=true;
//               if (p==numPolymSitesInSources) { //Examine the last segment after the last polymorphic site
//                  curSite=numSites-1;
//                  requestMatchingPolymSite=false;
//               }
////               computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP1, lastSite, curSite, curDef.n1, curScoreP1, closestChromP1, requestMatchingPolymSite);
////               computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP2, lastSite, curSite, curDef.n2, curScoreP2, closestChromP2, requestMatchingPolymSite);
////               if (curDef.n3) computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP3, lastSite, curSite, curDef.n3, curScoreP3, closestChromP3, requestMatchingPolymSite);
////               if (curDef.n4) computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP4, lastSite, curSite, curDef.n4, curScoreP4, closestChromP4, requestMatchingPolymSite);
//               
//               
//               computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP1, lastSite, curSite, curDef.n1, curScoreP1, closestChromP1);
//               computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP2, lastSite, curSite, curDef.n2, curScoreP2, closestChromP2);
//               if (curDef.n3) computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP3, lastSite, curSite, curDef.n3, curScoreP3, closestChromP3);
//               if (curDef.n4) computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP4, lastSite, curSite, curDef.n4, curScoreP4, closestChromP4);
//               
//               if (curScoreP1>-1 && curScoreP2>-1) {                  
//                  //We introduce a penalty in case of recombination
//                  if (p && closestChromP1!=lastCloseChromP1) curScoreP1+=recombPenalty;
//                  if (p && closestChromP2!=lastCloseChromP2) curScoreP2+=recombPenalty;
//                  if (curDef.n3) {
//                     if (curScoreP3>-1 && p && closestChromP3!=lastCloseChromP3) curScoreP3+=recombPenalty; 
//                     else if (curScoreP3<0) {
//                        curScoreP3=curScoreP3+lackOfPolymSitePenalty;
//                     }
//                  }
//                  if (curDef.n4) {
//                     if (curScoreP4>-1 && p && closestChromP4!=lastCloseChromP4) curScoreP3+=recombPenalty; 
//                     else if (curScoreP3<0) {
//                        curScoreP4=curScoreP4+lackOfPolymSitePenalty;
//                     }
//                  }
//               } else { //one of the two source pop is fixed for the allele that is not found on the admixed chromosome at the polymorphic position
//                  if (curScoreP1<0) {
//                    curScoreP1=curScoreP2+lackOfPolymSitePenalty;
//                  } else {
//                     if (curScoreP2<0) {
//                        curScoreP2=curScoreP1+lackOfPolymSitePenalty;
//                     } else { //both P1 and P2 are fixed for another allele, which should not happen
//                        bool bugFound=true;
//                     }
//                  }                        
//                  //We introduce a penalty in case of recombination
//                  if (p && closestChromP1!=lastCloseChromP1) curScoreP1+=recombPenalty;
//                  if (p && closestChromP2!=lastCloseChromP2) curScoreP2+=recombPenalty;
//                  if (curDef.n3 && p && closestChromP3!=lastCloseChromP3) curScoreP3+=recombPenalty;
//                  if (curDef.n4 && p && closestChromP4!=lastCloseChromP4) curScoreP4+=recombPenalty;
//               }
//               
//               totalScoreP1+=curScoreP1;
//               totalScoreP2+=curScoreP2;
//               if (curDef.n3) totalScoreP3+=curScoreP3;
//               if (curDef.n4) totalScoreP4+=curScoreP4;
////               //Assign chromosome segment to source pop with closest segment
////               if (curScoreP1<curScoreP2 || (closestChromP1>-1 && closestChromP2<0) ) {
////                  seq_origins_from_rec[j][p]=1; 
////               }
////               else {
////                  if (curScoreP2<curScoreP1 || (closestChromP2>-1 && closestChromP1<0)) {
////                     seq_origins_from_rec[j][p]=2;
////                  }
////                  else // it means that (curScoreP1==curScoreP2) {
////                  if (closestChromP1==lastCloseChromP1) {
////                     seq_origins_from_rec[j][p]=1;
////                  }
////                  if (closestChromP2==lastCloseChromP2) {
////                     seq_origins_from_rec[j][p]=2;
////                  }
////                  else seq_origins_from_rec[j][p]=2; 
////               }
//               if (closestChromP1!=lastCloseChromP1) ++numRecP1;
//               if (closestChromP2!=lastCloseChromP2) ++numRecP2;
//               if (closestChromP3!=lastCloseChromP3) ++numRecP3;
//               if (closestChromP4!=lastCloseChromP4) ++numRecP4;
//               lastCloseChromP1=closestChromP1;
//               lastCloseChromP2=closestChromP2;
//               lastCloseChromP3=closestChromP3;
//               lastCloseChromP4=closestChromP4;
//               lastSite=curSite+1;
//            } 
//            //Assign admixed chromosome segment to source pop with closest recombining segment
//            int assignedPop=0;
//            if (curDef.n3) {
//              if (curDef.n4) {
//                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3 && totalScoreP1<=totalScoreP4) {
//                   assignedPop=1;
//                }
//                else {                   
//                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3 && totalScoreP2<=totalScoreP4) {
//                      assignedPop=2;
//                   }
//                   else {
//                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2 && totalScoreP3<=totalScoreP4) {
//                        assignedPop=3;
//                     }
//                     else {
//                        if (totalScoreP4<=totalScoreP1 && totalScoreP4<=totalScoreP2 && totalScoreP4<=totalScoreP3) {
//                           assignedPop=4;
//                        }
//                        else assignedPop=1; 
//                     }
//                   }
//                } 
//              } 
//              else { //only 3 possible source populations
//                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3) {
//                   assignedPop=1; 
//                }
//                else {
//                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3) {
//                      assignedPop=2; 
//                   }
//                   else {
//                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2) {
//                        assignedPop=3; 
//                     }
//                     else assignedPop=1; 
//                   } 
//                }
//              } 
//            } else {  //only 2 possible source populations
//               if (totalScoreP1<=totalScoreP2) {
//                  assignedPop=1; 
//               }
//               else assignedPop=2;  
//            }
//            seq_origins_from_rec[j][rs]=assignedPop;
//            
////            if (totalScoreP1<totalScoreP2) {
////               seq_origins_from_rec[j][rs]=1; 
////            }
////            else {
////               if (totalScoreP1==totalScoreP2) { //Then select the one with the fewest recombinations 
////                  if (numRecP1<numRecP2) {
////                     seq_origins_from_rec[j][rs]=1;
////                  }
////                  else {
////                     seq_origins_from_rec[j][rs]=2;
////                  }
////               }
////               else seq_origins_from_rec[j][rs]=2;
////            }
//            startSite+=incNumPolSites;
//            endSite+=incNumPolSites;
//         }
////         cout << endl;
//      }
//      cout << " DONE!\n\n";
//
//      
//      cout << "\tInferring chromosome origins from shortest recombination paths ... ";
//      
//      int numSourceChroms=curDef.n1+curDef.n2+curDef.n3+curDef.n4;
//      
//      //Step 1: Identify the non-private mutations in admixed chromosome
//         
//      for (int i=0; i<numSites ; ++i) {
//         int numDerivedAllelesInSources=0;
//         //Count derived allele frequency in source chromosomes
//         for (int j=0; j<numSourceChroms; ++j) {
//            if (curGenotypes[j][i]) ++numDerivedAllelesInSources;
//         } 
//         bool monoInSources=false;
//         if ((numDerivedAllelesInSources==0 || numDerivedAllelesInSources==numSourceChroms))  monoInSources=true;
//         double curFreqInSource=numDerivedAllelesInSources/numSourceChroms;
//         if (!monoInSources) {
//            //Then this position must be polymorphic in all admixed individuals
//            polymPositions.push_back(i);
//            numPolymSitesInSources++;
//         }
//      }   
//      
//      //Pretend that the last site is polymorphic in sources
//      polymPositions.push_back(numSites-1);
//      
//      int numRecSegments=(numPolymSitesInSources-numConsecPolSitesToConsider)/incNumPolSites+1;
//      for (int na=0; na<curDef.nadm; ++na) { 
//         seq_origins_from_rec[na]   = new int[numRecSegments+1];
//      }
//      
//      //Loop over all admixed chromosomes
//      #pragma omp parallel for
//      for (int j=0; j<curDef.nadm; ++j) {
//         
////         cout << "\t\tExamining individual " << j+1 << "/" << curDef.nadm <<"\r";
//      
//         int curAdmChrom=curDef.n1+curDef.n2+curDef.n3+curDef.n4+j;         
//         
//         int startSite=0;
//         int endSite=numConsecPolSitesToConsider;
//         for (int rs=0; rs<numRecSegments+1; ++rs) {
//            int lastSite=polymPositions[startSite];
//            int lastCloseChromP1=-1,lastCloseChromP2=-1, lastCloseChromP3=-1,lastCloseChromP4=-1;
//            int totalScoreP1=0, totalScoreP2=0, totalScoreP3=0, totalScoreP4=0;
//            int numRecP1=0, numRecP2=0, numRecP3=0, numRecP4=0 ;
//            //Loop over all polymorphic sites in the current segment
//            for (int p=startSite; p<endSite; ++p) {
//               int curSite=polymPositions[p];
//               //Examine each source pop in turn and find the source chromosome that is closest to the admixed chromosome
//               int curScoreP1, curScoreP2,  curScoreP3, curScoreP4;
//               int closestChromP1, closestChromP2, closestChromP3, closestChromP4;
//               bool requestMatchingPolymSite=true;
//               if (p==numPolymSitesInSources) { //Examine the last segment after the last polymorphic site
//                  curSite=numSites-1;
//                  requestMatchingPolymSite=false;
//               }
//               computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP1, lastSite, curSite, curDef.n1, curScoreP1, closestChromP1, requestMatchingPolymSite);
//               computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP2, lastSite, curSite, curDef.n2, curScoreP2, closestChromP2, requestMatchingPolymSite);
//               if (curDef.n3) computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP3, lastSite, curSite, curDef.n3, curScoreP3, closestChromP3, requestMatchingPolymSite);
//               if (curDef.n4) computeMinNumDiffs((*curGenot)[curAdmChrom], GenotypesP4, lastSite, curSite, curDef.n4, curScoreP4, closestChromP4, requestMatchingPolymSite);
//               
//               if (curScoreP1>-1 && curScoreP2>-1) {                  
//                  //We introduce a penalty in case of recombination
//                  if (p && closestChromP1!=lastCloseChromP1) curScoreP1+=recombPenalty;
//                  if (p && closestChromP2!=lastCloseChromP2) curScoreP2+=recombPenalty;
//                  if (curDef.n3) {
//                     if (curScoreP3>-1 && p && closestChromP3!=lastCloseChromP3) curScoreP3+=recombPenalty; 
//                     else if (curScoreP3<0) {
//                        curScoreP3=curScoreP3+lackOfPolymSitePenalty;
//                     }
//                  }
//                  if (curDef.n4) {
//                     if (curScoreP4>-1 && p && closestChromP4!=lastCloseChromP4) curScoreP4+=recombPenalty; 
//                     else if (curScoreP3<0) {
//                        curScoreP4=curScoreP4+lackOfPolymSitePenalty;
//                     }
//                  }
//               } else { //one of the two source pop is fixed for the allele that is not found on the admixed chromosome at the polymorphic position
//                  if (curScoreP1<0) {
//                    curScoreP1=curScoreP2+lackOfPolymSitePenalty;
//                  } else {
//                     if (curScoreP2<0) {
//                        curScoreP2=curScoreP1+lackOfPolymSitePenalty;
//                     } else { //both P1 and P2 are fixed for another allele, which should not happen
//                        bool bugFound=true;
//                     }
//                  }                        
//                  //We introduce a penalty in case of recombination
//                  if (p && closestChromP1!=lastCloseChromP1) curScoreP1+=recombPenalty;
//                  if (p && closestChromP2!=lastCloseChromP2) curScoreP2+=recombPenalty;
//                  if (curDef.n3 && p && closestChromP3!=lastCloseChromP3) curScoreP3+=recombPenalty;
//                  if (curDef.n4 && p && closestChromP4!=lastCloseChromP4) curScoreP4+=recombPenalty;
//               }
//               
//               totalScoreP1+=curScoreP1;
//               totalScoreP2+=curScoreP2;
//               if (curDef.n3) totalScoreP3+=curScoreP3;
//               if (curDef.n4) totalScoreP4+=curScoreP4;
////               //Assign chromosome segment to source pop with closest segment
////               if (curScoreP1<curScoreP2 || (closestChromP1>-1 && closestChromP2<0) ) {
////                  seq_origins_from_rec[j][p]=1; 
////               }
////               else {
////                  if (curScoreP2<curScoreP1 || (closestChromP2>-1 && closestChromP1<0)) {
////                     seq_origins_from_rec[j][p]=2;
////                  }
////                  else // it means that (curScoreP1==curScoreP2) {
////                  if (closestChromP1==lastCloseChromP1) {
////                     seq_origins_from_rec[j][p]=1;
////                  }
////                  if (closestChromP2==lastCloseChromP2) {
////                     seq_origins_from_rec[j][p]=2;
////                  }
////                  else seq_origins_from_rec[j][p]=2; 
////               }
//               if (closestChromP1!=lastCloseChromP1) ++numRecP1;
//               if (closestChromP2!=lastCloseChromP2) ++numRecP2;
//               if (closestChromP3!=lastCloseChromP3) ++numRecP3;
//               if (closestChromP4!=lastCloseChromP4) ++numRecP4;
//               lastCloseChromP1=closestChromP1;
//               lastCloseChromP2=closestChromP2;
//               lastCloseChromP3=closestChromP3;
//               lastCloseChromP4=closestChromP4;
//               lastSite=curSite+1;
//            } 
//            //Assign admixed chromosome segment to source pop with closest recombining segment
//            int assignedPop=0;
//            if (curDef.n3) {
//              if (curDef.n4) {
//                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3 && totalScoreP1<=totalScoreP4) {
//                   assignedPop=1;
//                }
//                else {                   
//                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3 && totalScoreP2<=totalScoreP4) {
//                      assignedPop=2;
//                   }
//                   else {
//                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2 && totalScoreP3<=totalScoreP4) {
//                        assignedPop=3;
//                     }
//                     else {
//                        if (totalScoreP4<=totalScoreP1 && totalScoreP4<=totalScoreP2 && totalScoreP4<=totalScoreP3) {
//                           assignedPop=4;
//                        }
//                        else assignedPop=1; 
//                     }
//                   }
//                } 
//              } 
//              else { //only 3 possible source populations
//                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3) {
//                   assignedPop=1; 
//                }
//                else {
//                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3) {
//                      assignedPop=2; 
//                   }
//                   else {
//                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2) {
//                        assignedPop=3; 
//                     }
//                     else assignedPop=1; 
//                   } 
//                }
//              } 
//            } else {  //only 2 possible source populations
//               if (totalScoreP1<=totalScoreP2) {
//                  assignedPop=1; 
//               }
//               else assignedPop=2;  
//            }
//            seq_origins_from_rec[j][rs]=assignedPop;
//            
////            if (totalScoreP1<totalScoreP2) {
////               seq_origins_from_rec[j][rs]=1; 
////            }
////            else {
////               if (totalScoreP1==totalScoreP2) { //Then select the one with the fewest recombinations 
////                  if (numRecP1<numRecP2) {
////                     seq_origins_from_rec[j][rs]=1;
////                  }
////                  else {
////                     seq_origins_from_rec[j][rs]=2;
////                  }
////               }
////               else seq_origins_from_rec[j][rs]=2;
////            }
//            startSite+=incNumPolSites;
//            endSite+=incNumPolSites;
//         }
////         cout << endl;
//      }      
//      cout << " DONE!\n\n";
      
      cout << "\tInferring chromosome origins from shortest recombination paths ... ";
      
//      int numSourceChroms=curDef.n1+curDef.n2+curDef.n3+curDef.n4;
      
      //Step 1: Identify the non-private mutations in admixed chromosome
//         
//      for (int i=0; i<numSites ; ++i) {
//         polymPositions.push_back(i);
//      }   
      
      int numRecSegments=(numSites-numConsecPolSitesToConsider)/incNumPolSites+1;
      for (int na=0; na<curDef.nadm; ++na) { 
         seq_origins_from_rec[na]   = new int[numRecSegments+1];
      }
      
      //Loop over all admixed chromosomes
      #pragma omp parallel for
      for (int j=0; j<curDef.nadm; ++j) {
         
//         cout << "\t\tExamining individual " << j+1 << "/" << curDef.nadm <<"\r";
      
         int curAdmChrom=curDef.n1+curDef.n2+curDef.n3+curDef.n4+j; 
         int winBeg=0;
         
         for (int rs=0; rs<numRecSegments+1; ++rs) {
            int winEnd=winBeg+numConsecPolSitesToConsider;
            if (winEnd>numSites-1) winEnd=numSites-1;
            
            int totalScoreP1=0, totalScoreP2=0, totalScoreP3=0, totalScoreP4=0;                                   
            
            //Loop over all polymorphic sites in the current segment
            int numSitesToCompare=nonRecTract;
            
            int lastClosestChromP1=-1, lastClosestChromP2=-1, lastClosestChromP3=-1, lastClosestChromP4=-1; 
            int closestChromP1, closestChromP2, closestChromP3, closestChromP4; 
            for (int p=winBeg; p<winEnd; p+=numSitesToCompare) {
               int beg=p;
               if ( beg+numSitesToCompare > winEnd) numSitesToCompare=winEnd-beg;
               //Examine each source pop in turn and find the source chromosome that is closest to the admixed chromosome
               int curScoreP1=0, curScoreP2=0,  curScoreP3=0, curScoreP4=0;
               
               computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP1, beg, numSitesToCompare, curDef.n1, curScoreP1, closestChromP1);
               computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP2, beg, numSitesToCompare, curDef.n2, curScoreP2, closestChromP2);
               if (curDef.n3) computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP3, beg, numSitesToCompare, curDef.n3, curScoreP3, closestChromP3);
               if (curDef.n4) computeMinNumDiffs_new((*curGenot)[curAdmChrom], GenotypesP4, beg, numSitesToCompare, curDef.n4, curScoreP4, closestChromP4);
               
//               if (p>winBeg) {
//                  if (lastClosestChromP1!=closestChromP1) {
//                     curScoreP1+=recombPenalty;
//                  }
//                  if (lastClosestChromP2!=closestChromP2) {
//                     curScoreP2+=recombPenalty;
//                  }
//                  if (curDef.n3 && lastClosestChromP3!=closestChromP3) {
//                     curScoreP3+=recombPenalty;
//                  }
//                  if (curDef.n4 && lastClosestChromP4!=closestChromP4) {
//                     curScoreP4+=recombPenalty;
//                  }
//               } 
               
               lastClosestChromP1=closestChromP1;
               lastClosestChromP2=closestChromP2;
               lastClosestChromP3=closestChromP3;
               lastClosestChromP4=closestChromP4;
               
               totalScoreP1+=curScoreP1;
               totalScoreP2+=curScoreP2;
               if (curDef.n3) totalScoreP3+=curScoreP3;
               if (curDef.n4) totalScoreP4+=curScoreP4;
               
            } 
            winBeg+=incNumPolSites;
            //Assign admixed chromosome segment to source pop with closest recombining segment
            int assignedPop=0;
            if (curDef.n3) {
              if (curDef.n4) {
                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3 && totalScoreP1<=totalScoreP4) {
                   assignedPop=1;
                }
                else {                   
                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3 && totalScoreP2<=totalScoreP4) {
                      assignedPop=2;
                   }
                   else {
                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2 && totalScoreP3<=totalScoreP4) {
                        assignedPop=3;
                     }
                     else {
                        if (totalScoreP4<=totalScoreP1 && totalScoreP4<=totalScoreP2 && totalScoreP4<=totalScoreP3) {
                           assignedPop=4;
                        }
                        else assignedPop=1; 
                     }
                   }
                } 
              } 
              else { //only 3 possible source populations
                if (totalScoreP1<=totalScoreP2 && totalScoreP1<=totalScoreP3) {
                   assignedPop=1; 
                }
                else {
                   if (totalScoreP2<=totalScoreP1 && totalScoreP2<=totalScoreP3) {
                      assignedPop=2; 
                   }
                   else {
                     if (totalScoreP3<=totalScoreP1 && totalScoreP3<=totalScoreP2) {
                        assignedPop=3; 
                     }
                     else assignedPop=1; 
                   } 
                }
              } 
            } else {  //only 2 possible source populations
               if (totalScoreP1<=totalScoreP2) {
                  assignedPop=1; 
               }
               else assignedPop=2;  
            }
            seq_origins_from_rec[j][rs]=assignedPop;
            
//            startSite+=incNumPolSites;
//            endSite+=incNumPolSites;
         }
      }      
      cout << " DONE!\n\n";
      
      //Inferring chromosome origins ----------------------------------------------
      cout << "\t\tInferring chromosome origins from population average or maximum similarities..."; 
      int begWin=0, endWin=winSize, beg, end;
      int numIncsInWindows=winSize/winInc;
      #pragma omp parallel for
      for (int ncAdm=0; ncAdm<curDef.nadm; ++ncAdm) { 

         TIntArray mismatchWithP1(curDef.n1, 0), mismatchWithP2(curDef.n2, 0), mismatchWithP3(curDef.n3, 0), mismatchWithP4(curDef.n4, 0);;
         for (int w=0; w<numWin; ++w) {       

            //SumNumDiff for Current window
            if (w==0) {            
               for (int k=0; k<curDef.n1; ++k) {
                  mismatchWithP1[k]=0;
                  for (int j=0; j<numIncsInWindows; ++j) {
                     mismatchWithP1[k]+=numDiffsRef1[j][ncAdm][k];
                  }
               }                       
               for (int k=0; k<curDef.n2; ++k) {
                  mismatchWithP2[k]=0;
                  for (int j=0; j<numIncsInWindows; ++j) {
                     mismatchWithP2[k]+=numDiffsRef2[j][ncAdm][k];
                  }
               }           
               for (int k=0; k<curDef.n3; ++k) {
                  mismatchWithP3[k]=0;
                  for (int j=0; j<numIncsInWindows; ++j) {
                     mismatchWithP3[k]+=numDiffsRef3[j][ncAdm][k];
                  }
               }                       
               for (int k=0; k<curDef.n4; ++k) {
                  mismatchWithP4[k]=0;
                  for (int j=0; j<numIncsInWindows; ++j) {
                     mismatchWithP4[k]+=numDiffsRef4[j][ncAdm][k];
                  }
               }
            } else {     
               for (int k=0; k<curDef.n1; ++k) {
                  mismatchWithP1[k]+=numDiffsRef1[w+numIncsInWindows-1][ncAdm][k]-numDiffsRef1[w-1][ncAdm][k];
               }            
               for (int k=0; k<curDef.n2; ++k) {
                  mismatchWithP2[k]+=numDiffsRef2[w+numIncsInWindows-1][ncAdm][k]-numDiffsRef2[w-1][ncAdm][k];
               }        
               for (int k=0; k<curDef.n3; ++k) {
                  mismatchWithP3[k]+=numDiffsRef3[w+numIncsInWindows-1][ncAdm][k]-numDiffsRef3[w-1][ncAdm][k];
               }            
               for (int k=0; k<curDef.n4; ++k) {
                  mismatchWithP4[k]+=numDiffsRef4[w+numIncsInWindows-1][ncAdm][k]-numDiffsRef4[w-1][ncAdm][k];
               }        
            }
            double meanDiffWithP1=0, meanDiffWithP2=0, meanDiffWithP3=0, meanDiffWithP4=0;
            
            meanDiffWithP1=computeMeanDiff(mismatchWithP1, curDef.n1);
            meanDiffWithP2=computeMeanDiff(mismatchWithP2, curDef.n2);
            if (curDef.n3) meanDiffWithP3=computeMeanDiff(mismatchWithP3, curDef.n3);
            if (curDef.n4) meanDiffWithP4=computeMeanDiff(mismatchWithP4, curDef.n4);          

            double minDiffWithP1=0, minDiffWithP2=0, minDiffWithP3=0, minDiffWithP4=0; 
            minDiffWithP1=computeMinDiff(mismatchWithP1, curDef.n1);
            minDiffWithP2=computeMinDiff(mismatchWithP2, curDef.n2); 
            if (curDef.n3) minDiffWithP3=computeMinDiff(mismatchWithP3, curDef.n3);
            if (curDef.n4) minDiffWithP4=computeMinDiff(mismatchWithP4, curDef.n4); 
                  
            int assignedSource=0;
            
            if (curDef.n3) {
               if (curDef.n4) { //4 sources
                  if (meanDiffWithP1<=meanDiffWithP2 && meanDiffWithP1<=meanDiffWithP3 && meanDiffWithP1<=meanDiffWithP4) assignedSource=1;
                  else {
                     if (meanDiffWithP2<=meanDiffWithP1 && meanDiffWithP2<=meanDiffWithP3 && meanDiffWithP2<=meanDiffWithP4) assignedSource=2;
                     else {
                        if (meanDiffWithP3<=meanDiffWithP1 && meanDiffWithP3<=meanDiffWithP2 && meanDiffWithP3<=meanDiffWithP4) assignedSource=3;
                        else {
                           if (meanDiffWithP4<=meanDiffWithP1 && meanDiffWithP4<=meanDiffWithP2 && meanDiffWithP4<=meanDiffWithP3) assignedSource=4;
                           else assignedSource=1;
                        }
                     }
                  }
               }
               else { //Only 3 sources                  
                  if (meanDiffWithP1<=meanDiffWithP2 && meanDiffWithP1<=meanDiffWithP3) assignedSource=1;
                  else {
                     if (meanDiffWithP2<=meanDiffWithP1 && meanDiffWithP2<=meanDiffWithP3 ) assignedSource=2;
                     else {
                        if (meanDiffWithP3<=meanDiffWithP1 && meanDiffWithP3<=meanDiffWithP2 ) assignedSource=3;
                        else assignedSource=1;
                     }
                  }
               }
            }
            else { //Only 2 sources
               if (meanDiffWithP1<=meanDiffWithP2) assignedSource=1;
               else {
                  if (meanDiffWithP2<=meanDiffWithP1) assignedSource=2;
                  else  assignedSource=1;
               }
            }
            seq_origins_from_mean[ncAdm][w]=assignedSource;
            
            assignedSource=0;
            
            if (curDef.n3) {
               if (curDef.n4) { //4 sources
                  if (minDiffWithP1<=minDiffWithP2 && minDiffWithP1<=minDiffWithP3 && minDiffWithP1<=minDiffWithP4) assignedSource=1;
                  else {
                     if (minDiffWithP2<=minDiffWithP1 && minDiffWithP2<=minDiffWithP3 && minDiffWithP2<=minDiffWithP4) assignedSource=2;
                     else {
                        if (minDiffWithP3<=minDiffWithP1 && minDiffWithP3<=minDiffWithP2 && minDiffWithP3<=minDiffWithP4) assignedSource=3;
                        else {
                           if (minDiffWithP4<=minDiffWithP1 && minDiffWithP4<=minDiffWithP2 && minDiffWithP4<=minDiffWithP3) assignedSource=4;
                           else assignedSource=1;
                        }
                     }
                  }
               }
               else { //Only 3 sources                  
                  if (minDiffWithP1<=minDiffWithP2 && minDiffWithP1<=minDiffWithP3) assignedSource=1;
                  else {
                     if (minDiffWithP2<=minDiffWithP1 && minDiffWithP2<=minDiffWithP3 ) assignedSource=2;
                     else {
                        if (minDiffWithP3<=minDiffWithP1 && minDiffWithP3<=minDiffWithP2 ) assignedSource=3;
                        else assignedSource=1;
                     }
                  }
               }
            }
            else { //Only 2 sources
               if (minDiffWithP1<=minDiffWithP2) assignedSource=1;
               else {
                  if (minDiffWithP2<=minDiffWithP1) assignedSource=2;
                  else  assignedSource=1;
               }
            }            
            seq_origins_from_min[ncAdm][w]=assignedSource;
            
//            if (meanDiffWithP1<meanDiffWithP2) seq_origins_from_mean[ncAdm][w]=1; else seq_origins_from_mean[ncAdm][w]=2;
//            if (minDiffWithP1<minDiffWithP2) seq_origins_from_min[ncAdm][w]=1; else seq_origins_from_min[ncAdm][w]=2;

            int stop=1;
         }

         begWin+=winInc;
         endWin+=winInc;

      }

      cout << " DONE!\n\n";

      char suffix[50];
      if (removeSingletons) sprintf(suffix, "%s", "noSingl");
      else sprintf(suffix, "%s", "withSingl");

      //Writing chromosomes origins -----------------------------------------------
      ofstream outputChromOriginsFromMean, outputChromOriginsFromMin, outputChromOriginsFromRec, outputFST, outputPiXY;
      char  nameOutputChromOriginsFromMean[200],nameOutputChromOriginsFromMin[200] , nameOutputChromOriginsFromRec[200] ,nameOutputFST[200], nameOutputPiXY[200];
      sprintf(nameOutputChromOriginsFromMean,  "%s_win_%i_inc_%i_fromMean_%s.adm", genericName.c_str(),winSize, winInc, suffix);
      sprintf(nameOutputChromOriginsFromMin,  "%s_win_%i_inc_%i_fromMin_%s.adm", genericName.c_str(),winSize, winInc, suffix);
//      sprintf(nameOutputChromOriginsFromRec,  "%s_win_%i_inc_%i_fromRec_%s.adm", genericName.c_str(),winSize, winInc, suffix);
//      sprintf(nameOutputChromOriginsFromRec,  "%s_n_%i_m_%i_p_%i_q_%i_fromRec_%s.adm", genericName.c_str(),numConsecPolSitesToConsider, incNumPolSites, recombPenalty, lackOfPolymSitePenalty, suffix);
      sprintf(nameOutputChromOriginsFromRec,  "%s_n_%i_m_%i_t_%i_fromRec_%s.adm", genericName.c_str(),numConsecPolSitesToConsider, incNumPolSites, nonRecTract, suffix);
      sprintf(nameOutputFST,  "%s_win_%i_inc_%i_FST.txt", genericName.c_str(),winSize, winInc);
      sprintf(nameOutputPiXY,  "%s_win_%i_inc_%i_PiXY.txt", genericName.c_str(),winSize, winInc);

      outputChromOriginsFromMean.open(nameOutputChromOriginsFromMean);
      outputChromOriginsFromMin.open(nameOutputChromOriginsFromMin);
      outputChromOriginsFromRec.open(nameOutputChromOriginsFromRec);
      outputFST.open(nameOutputFST);
      outputPiXY.open(nameOutputPiXY);

      int * last_adm_from_mean = new int[curDef.nadm];
      int * last_adm_from_min = new int[curDef.nadm];
      int * last_adm_from_rec = new int[curDef.nadm];


      //Writing FST result to output files-----------------------------------------
      cout << "\t\tWriting FST to output files ...";
      outputFST << "Pos.\tpi12\tpi1\tpi2\tFST\n";

      outputFST <<  winSizeForFST/2 << "\t" << pi12[0] << "\t" << pi1[0] << "\t"  << pi2[0] << "\t"  << FSTBetweenRefs[0] << "\n"; 

      int newBeg=winSizeForFST;
      for (int w=1; w<numFSTWin; ++w) {
         outputFST <<  newBeg + winSizeForFST/2 << "\t" << pi12[w] << "\t" << pi1[w] << "\t"  << pi2[w] << "\t"  << FSTBetweenRefs[w] << "\n" ;    
         newBeg+=winSizeForFST;   
      }
      outputFST.close();
      cout << "\t\tDONE!\n\n";
      
      //Writing PiXY result to output files-----------------------------------------
      cout << "\t\tWriting PiXY to output files ...";
      outputPiXY << "Pos.\tpi12\tpi1\tpi2\tNei'D per nucleotide\n";
      
      int numPops=2;
      if ( pi1[0]==0  | pi2[0]==0 ) numPops=1;

      outputPiXY <<  winSizeForFST/2 << "\t" << pi12[0] << "\t" << pi1[0] << "\t"  << pi2[0] << "\t"  << (pi12[0]-(pi1[0]+pi2[0])/numPops)/winSizeForFST << "\n"; 

      newBeg=winSizeForFST;
      for (int w=1; w<numFSTWin; ++w) {
         if ( pi1[w]==0  | pi2[w]==0 ) numPops=1;
         outputPiXY <<  newBeg + winSizeForFST/2 << "\t" << pi12[w] << "\t" << pi1[w] << "\t"  << pi2[w] << "\t"  << (pi12[w]-(pi1[w]+pi2[w])/numPops)/winSizeForFST << "\n" ;    
         newBeg+=winSizeForFST;   
      }
      outputPiXY.close();
      cout << "\t\tDONE!\n\n";

      //Writing chromosome origins to output files --------------------------------
      cout << "\t\tWriting chromosome origins to output files ...";
      //Writing header
      outputChromOriginsFromMean << "beg";
      outputChromOriginsFromMin << "beg";  
      outputChromOriginsFromRec << "beg";  
      for (int na=0; na<curDef.nadm; na++) {
         outputChromOriginsFromMean << "\ta" << (na+1);
         outputChromOriginsFromMin << "\ta" << (na+1);
         outputChromOriginsFromRec << "\ta" << (na+1);
      }      
      outputChromOriginsFromMean << "\n";
      outputChromOriginsFromMin << "\n";
      outputChromOriginsFromRec << "\n";

      //Writing first line
      outputChromOriginsFromMean << "0";
      outputChromOriginsFromMin << "0";  
      outputChromOriginsFromRec << "0";   
      for (int na=0; na<curDef.nadm; na++) {
         outputChromOriginsFromMean << "\t" << seq_origins_from_mean[na][0];
         outputChromOriginsFromMin << "\t" << seq_origins_from_min[na][0];
         outputChromOriginsFromRec << "\t" << seq_origins_from_rec[na][0];
         last_adm_from_mean[na]=seq_origins_from_mean[na][0];
         last_adm_from_min[na]=seq_origins_from_min[na][0];
         last_adm_from_rec[na]=seq_origins_from_rec[na][0];
      }      
      outputChromOriginsFromMean << "\n";
      outputChromOriginsFromMin << "\n";
      outputChromOriginsFromRec << "\n";

      bool writeNewLineInMeanFile=false, writeNewLineInMinFile=false, writeNewLineInRecFile=false; 
      
      //Output of chromosome segments for mean diff and min diff based approaches
      //-------------------------------------------------------------------------
      newBeg=winInc;;
      for (int w=1; w<numWin; ++w) {   

         for (int na=0; na<curDef.nadm; na++) {

            if (last_adm_from_mean[na]!=seq_origins_from_mean[na][w]) {
               last_adm_from_mean[na]=seq_origins_from_mean[na][w];
               writeNewLineInMeanFile=true;
            }   
            if (last_adm_from_min[na]!=seq_origins_from_min[na][w]) {
               last_adm_from_min[na]=seq_origins_from_min[na][w];
               writeNewLineInMinFile=true;
            }     
         }
         if (writeNewLineInMeanFile) {
            outputChromOriginsFromMean << newBeg + winSize/2;
            for (int na=0; na<curDef.nadm; na++) {
               outputChromOriginsFromMean << "\t" << seq_origins_from_mean[na][w];
            }
            outputChromOriginsFromMean << "\n";
            writeNewLineInMeanFile=false;
         }
         if (writeNewLineInMinFile) {
            outputChromOriginsFromMin << newBeg + winSize/2;
            for (int na=0; na<curDef.nadm; na++) {
               outputChromOriginsFromMin << "\t" << seq_origins_from_min[na][w];
            }
            outputChromOriginsFromMin << "\n";
            writeNewLineInMinFile=false;
         }
         
         newBeg+=winInc;
      }
      
      //Output of chromosome segments for min rec path based approach
      //-------------------------------------------------------------------------
      
//      int numRecSegments=(numPolymSitesInSources-numConsecPolSitesToConsider)/incNumPolSites+1;
      int curPolSiteBeg=0;
      for (int rs=0; rs<numRecSegments+1; ++rs) { 
//         newBeg=positions[polymPositions[curPolSiteBeg]];
         //Loro_16_01_23 Using the SNP in the middle the window rather than SNP at the beginning
//         int midSNP=curPolSiteBeg+incNumPolSites/2;
         int midSNP=curPolSiteBeg+numConsecPolSitesToConsider/2;
         
         if (midSNP>numSites) newBeg=(chromLength+positions[curPolSiteBeg])/2;
         else newBeg=positions[midSNP];
         for (int na=0; na<curDef.nadm; na++) {
            if (last_adm_from_rec[na]!=seq_origins_from_rec[na][rs]) {
               last_adm_from_rec[na]=seq_origins_from_rec[na][rs];
               writeNewLineInRecFile=true;
            }       
         }
         
         if (writeNewLineInRecFile) {
            outputChromOriginsFromRec << newBeg;
            for (int na=0; na<curDef.nadm; na++) {
               outputChromOriginsFromRec << "\t" << seq_origins_from_rec[na][rs];
            }
            outputChromOriginsFromRec << "\n";
            writeNewLineInRecFile=false;
         }
         curPolSiteBeg+=incNumPolSites;
      }

      cout << "\t\tDONE!\n";
   
      outputChromOriginsFromMean << endl;
      outputChromOriginsFromMin << endl;
      outputChromOriginsFromRec << endl;
      outputFST <<  endl; 


      outputChromOriginsFromMean.close();
      outputChromOriginsFromMin.close();
      outputChromOriginsFromRec.close();
      outputFST.close();
   
      cout << "\tDONE!\n\n";
   }
   //=======================
   //Deleting dynamic arrays
   //=======================

   for (int s=0; s<numChroms; ++s) {
      if (curGenotypes[s]) delete [] curGenotypes[s];
      if (nonSingletonGenotypes[s]) delete [] nonSingletonGenotypes[s]; 
      
   }
   delete [] curGenotypes;
   delete [] nonSingletonGenotypes;
   delete [] nonPrivateGenotypes;

   delete[] positions;
   delete[] positionsNonSingletons;
   delete[] positionsNonPrivate;

   for (int na=0; na<curDef.nadm; ++na) {         
      delete [] seq_origins_from_mean[na];        
      delete [] seq_origins_from_min[na];         
      delete [] seq_origins_from_rec[na];
      
   }

   delete []  seq_origins_from_mean;
   delete []  seq_origins_from_min;
   delete []  seq_origins_from_rec;
   delete []  FSTBetweenRefs;
   delete [] pi12;
   delete [] pi1;
   delete [] pi2;
   delete [] pi3;
   delete [] pi4;
   
    for (int i=0; i<numWinInc; ++i) {
      for (int j=0; j<curDef.nadm; ++j) {
         delete [] numDiffsRef1[i][j];
         delete [] numDiffsRef2[i][j];
         if (numDiffsRef3[i][j]) delete [] numDiffsRef3[i][j];
         if (numDiffsRef4[i][j]) delete [] numDiffsRef4[i][j];
         
      }
      delete [] numDiffsRef1[i];
      delete [] numDiffsRef2[i];
      delete [] numDiffsRef3[i];
      delete [] numDiffsRef4[i];
      
      delete []  recPathLengthRef1[i];
      delete []  recPathLengthRef2[i];
      delete []  recPathLengthRef3[i];
      delete []  recPathLengthRef4[i];
   }
   delete [] numDiffsRef1;
   delete [] numDiffsRef2;   
   delete [] numDiffsRef3;
   delete [] numDiffsRef4;
   
   
   delete []  recPathLengthRef1;
   delete []  recPathLengthRef2;
   delete []  recPathLengthRef3;
   delete []  recPathLengthRef4;
 
   return 0;
};
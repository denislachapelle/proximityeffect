//                                TwoDWiresGenerator
//   
/*

*/

#ifndef TWODWIRESGENERATOR
#define TWODWIRESGENERATOR

#include "mytools.hpp"
#include <iostream>
#include <math.h>
#include <filesystem>
#include "gmsh.h"

using namespace std;

class WireInfo
{
   public:
   enum WireTypes {roundwire, rectangular, awg};
   WireTypes type;
   double dimensions[10];
   double center[2];
   double current[2];
};

class TwoDWiresGenerator
{

   protected:
      WireInfo *wiresInfo;

   private:
      // 1. Parse command-line options.
      const char *configFile = "fourwires.txt";
      const char *meshFile = "fourwires.msh";
      int nbrwires = 2;
            
      real_t domainRadius = -1.0;
      int refineTimes = 0;

   public:
      //parse the options.
      int Parser(int argc, char *argv[]);

      int ReadConfigFile();
      int CreateMeshFile();
      int getNbrWires();
      WireInfo *getWireInfo();

};
#endif //2DWIRESGENERATOR

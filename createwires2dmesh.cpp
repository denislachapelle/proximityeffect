//                                createwires2dmesh
//                                
// compile with: make createwires2dmesh, need gmsh.
//
// Sample runs:  ./createwires2dmesh -cf wiresconfig.txt
//
/*
Description:  

*/

#include "twodwiresgenerator.hpp"

int main(int argc, char *argv[])
{

   TwoDWiresGenerator WG;
 
   WG.Parser(argc, argv);

   WG.ReadConfigFile();

   WG.CreateMeshFile();
   
   return 0;
}

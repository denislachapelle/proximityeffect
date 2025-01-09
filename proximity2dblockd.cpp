//                                proximity2dblockc
//                                inspired from proximity2dblock and se2dblock
//                                inspired from proximity2d
//                                based on MFEM Example 22 prob 1 (case 0), ex5...
//
// Compile with: make proximity2dblockc, need MFEM version 4.7 and GLVIS-4.3.
//
// Sample runs:  ./proximity2dblockd
//
/*
Description:  

Implementation suggested by paul hilsher in https://github.com/mfem/mfem/issues/4584
from paper "Specialized conductor models for finite element eddy current simulation".
https://www.iem.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaabcfymi

The equation is 
- div(grad(Az)) + i w u s Az - u s V = 0

Az : magnetic vector potential z-component.
V : voltage potential.
u : permeability.
s : conductivity.
w : pulsation, 2 pi f.
Assume 1m long, so s V is electric field.

then split in real and imag parts
- div(grad(Azr)) - w u s Azi - u s Vr = 0
- div(grad(Azi)) + w u s Azr - u s Vi = 0

Let define matrix A1, A2, A3 and A4 as implementing ....
A1 -div(grad((.))
A2 - u w s (.)
A3 - u s (.)
A4   w u s (.)

We also want to enforce the current in each wires..
intS(- i s1 w A) + V1/R1 = I1
same for wire 2.

then split real and imaginary equation...
intS(s1 w Ai) + V1r/R1 = I1r
intS(- s1 w Ar) + V1i/R1 = I1i
same for wire 2.

let define A5_1 and A6_1 as implementing ...
A5 intS(s1 w (.))
A6 1/R1 (.)
A7 intS(-s1 w (.))

intS means surface integral over the wire n.

Then we can write the assembled matrix...

[A1   A2   A3_1 0    A3_2 0   ] [Azr] = [0]
[A4   A1   0    A3_1 0    A3_2] [Azi] = [0]
[0    A5_1 A6_1 0    0    0   ] [V1r]  = [I1r]
[A7_1 0    0    A6_1 0    0   ] [V1i]  = [I1i]
[0    A5_2 0    0    A6_2    0] [V2r]  = [I2r]
[A7_2 0    0    0    0    A6_2] [V2i]  = [I2i]

Sub block dimensions:
[dofxdof dofxdof dofx1 dofx1 ...] [dofx1] = [dofx1]
[dofxdof dofxdof dofx1 dofx1 ...] [dofx1] = [dofx1]
[1xdof   1xdof   1x1   1x1 ...]   [1x1] =   [1x1]
[1xdof   1xdof   1x1   1x1 ...]   [1x1] =   [1x1]
[1xdof   1xdof   1x1   1x1 ...]   [1x1] =   [1x1]
[1xdof   1xdof   1x1   1x1 ...]   [1x1] =   [1x1]
I1r being total real current in wire 1.
I1i being total imaginary current in wire 1.

Once solved the current density can be computed...
J = - i w s A - s V

Jr =   w s Ai + s Vr
Ji = - w s Ar + s Vi

||J|| = sqrt(Jr^2+Ji^2)

u: permeability.
e: permitivity.
s: conductivity.
i: sqrt(-1)
*/

#include <mfem.hpp>
#include <linalg/hypre.hpp>
#include "mytools.hpp"
#include <iostream>
#include <math.h>
#include <filesystem>
#include "/usr/local/include/gmsh.h"

//Value expected from mesh file.
#define AIR        (nbrwires+1)
#define AIRCONTOUR (nbrwires+2)

using namespace std;
using namespace mfem;


class WireInfo
{
   public:
   enum WireTypes {roundwire, rectangular};
   WireTypes type;
   double dimensions[10];
   double center[2];
   double current[2];
};


class ProximityEffect
{

   protected:
      WireInfo *wiresInfo;


   private:
      // 1. Parse command-line options.
      const char *configFile = "";
      const char *meshFile = "tworoundwires2d.msh";
      int order = 1;
      double freq = -1.0;
      int nbrwires = 2;
            
      real_t mu_ = 1.256637061E-6;
      real_t epsilon_ = 8.8541878188e-12;
      real_t sigma_ = 59.59E6; // at 20C.
      real_t omega_ = 2.0*M_PI*60;
      real_t domainRadius = -1.0;

      Mesh *mesh;
      int dim;
      int nbrel;  //number of element.

      FiniteElementCollection *fec;
      FiniteElementSpace *fespace;
      int nbrdof;   //number of degree of freedom.

      Array<int> *ess_tdof_list_block; //essential dof list.

      //matrix pointer and matrix array pointer
      //for holding the matrix forming the block operator.
      SparseMatrix *A1, *A2, **A3, *A4, **A5, **A6, **A7;
      Vector *rhs, *x;

      BlockOperator *A;
      Array<int> *blockOffset;
      BlockOperator *ProxOp;

      Operator *A_ptr;
      BlockVector *B, *X;
      
      BlockDiagonalPreconditioner *block_prec;

      GridFunction *AzrGF, *AziGF;  // magnetic vector potential z-axis, real and imaginary.

      //Space and gridfunction for the current density.
      FiniteElementCollection *JFec;
      FiniteElementSpace *JFESpace;
      GridFunction *JrGF, *JiGF, *JGF;

   public:
      //delete all files in out dir.
      int CleanOutDir();
      //parse the options.
      int Parser(int argc, char *argv[]);

      int ReadConfigFile();
      int CreateMeshFile();

      int LoadMeshFile();
      int CreateFESpace();
      int CreateEssentialBoundary();
      int CreateOperatorA1();
      int CreateOperatorA2();
      int CreateOperatorA3();
      int CreateOperatorA4();
      int CreateOperatorA5();
      int CreateOperatorA6();
      int CreateOperatorA7();
      int CreaterhsVector();
      int CreatexVector();
      int CreateBlockOperator();
      int CreatePreconditionner();
      int Solver();
      int PostPrecessing();
      int DisplayResults();
//      bool IsCreatedMeshFile();
};

/*
bool ProximityEffect::IsCreatedMeshFile()
{
   if(meshFile == "") return true;
   else return false;
}
*/
int ProximityEffect::Parser(int argc, char *argv[])
{

   OptionsParser args(argc, argv);
   args.AddOption(&meshFile, "-mf", "--meshfile",
                  "file to use as mesh file.");
   args.AddOption(&configFile, "-cf", "--configfile",
                  "file to use as wires config file.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&mu_, "-mu", "--permeability",
                  "Permeability of free space (or 1/(spring constant)).");
   args.AddOption(&epsilon_, "-eps", "--permittivity",
                  "Permittivity of free space (or mass constant).");
   args.AddOption(&sigma_, "-sigma", "--conductivity",
                  "Conductivity (or damping constant).");
   args.AddOption(&freq, "-f", "--frequency",
                  "Frequency (in Hz).");
                
   args.Parse();

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, occ::, RAJA and OpenMP based on command line options.
   Device device("cpu");
   device.Print();

   if ( !(freq < 0.0) ) omega_ = 2.0 * M_PI * freq;

   if (args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);
   return 0;
}


int ProximityEffect::ReadConfigFile()
{
   std::string s;
   std::stringstream ss;
   nbrwires = HowManyWire(configFile);
   cout << nbrwires << " nbrWire\n";
   assert(nbrwires>0);
   getchar();

   if(nbrwires>0) wiresInfo = new WireInfo[nbrwires];
 
   for(int wc=0; wc<nbrwires; wc++)
   {
      ss.str("");
      ss << "wire" << wc+1 << "shape";
      s = findWordN(configFile, ss.str(), 2);
      if(s =="round")
      {
         wiresInfo[wc].type = WireInfo::WireTypes::roundwire;

         ss.str("");
         ss << "wire" << wc+1 << "dimensions";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].dimensions[0] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "center";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].center[0] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "center";
         s = findWordN(configFile, ss.str(), 3);
         wiresInfo[wc].center[1] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "current";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].current[0] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "current";
         s = findWordN(configFile, ss.str(), 3);
         wiresInfo[wc].current[1] = stof(s);
      }
      else if(s =="rectangular")
      {
         wiresInfo[wc].type = WireInfo::WireTypes::rectangular;

         ss.str("");
         ss << "wire" << wc+1 << "dimensions";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].dimensions[0] = stof(s);
         
         ss.str("");
         ss << "wire" << wc+1 << "dimensions";
         s = findWordN(configFile, ss.str(), 3);
         wiresInfo[wc].dimensions[1] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "center";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].center[0] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "center";
         s = findWordN(configFile, ss.str(), 3);
         wiresInfo[wc].center[1] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "current";
         s = findWordN(configFile, ss.str(), 2);
         wiresInfo[wc].current[0] = stof(s);

         ss.str("");
         ss << "wire" << wc+1 << "current";
         s = findWordN(configFile, ss.str(), 3);
         wiresInfo[wc].current[1] = stof(s);
      }
   }
   
   ss.str("");
   ss << "domainradius";
   s = findWordN(configFile, ss.str(), 2);
   domainRadius = stof(s);
   
   return 1;
}


int ProximityEffect::CreateMeshFile()
{
   // Before using any functions in the C++ API, gmsh::must be initialized:
   gmsh::initialize();

   gmsh::model::add("proxmesh");
   int wc;
   int *wtag;
   wtag = new int[nbrwires];

   for(wc=0; wc<nbrwires; wc++)
   {
      if(wiresInfo[wc].type==WireInfo::WireTypes::roundwire)
      {
         int circleTag =gmsh::model::occ::addCircle(wiresInfo[wc].center[0], wiresInfo[wc].center[1], 0, wiresInfo[wc].dimensions[0], -1);
         int curveLoopTag = gmsh::model::occ::addCurveLoop({circleTag}, -1);
         wtag[wc] = gmsh::model::occ::addPlaneSurface({curveLoopTag}, -1);
      }
      else if(wiresInfo[wc].type==WireInfo::WireTypes::rectangular)
      {
         wtag[wc] = gmsh::model::occ::addRectangle(wiresInfo[wc].center[0] - wiresInfo[wc].dimensions[0]/2.0,
                                              wiresInfo[wc].center[1] - wiresInfo[wc].dimensions[1]/2.0,
                                              0,
                                              wiresInfo[wc].dimensions[0],
                                              wiresInfo[wc].dimensions[1],
                                              -1,
                                              0);
      }

   }
   gmsh::model::occ::synchronize();

   //domain limit

   int circleTag = gmsh::model::occ::addCircle(0, 0, 0, domainRadius, -1);

   int curveLoopTag = gmsh::model::occ::addCurveLoop({circleTag}, -1);

   std::vector<int> vec(nbrwires+1);
   vec.at(0)=curveLoopTag;
   for(wc=0; wc<nbrwires; wc++ )
   {
      vec.at(wc+1)=-(wtag[wc]);
   /*
      if(wiresInfo[wc].type==WireInfo::WireTypes::roundwire)
      {
         array.at(wc+1)=-(wc+100);
      }
      else if(wiresInfo[wc].type==WireInfo::WireTypes::rectangular)
      {
         // std::vector<std::pair<int, int>> outDimTag;
         // Get the boundary of the surface (dimension 2 -> 1, surface to edges)
         // gmsh::model::getBoundary({{2, wc+100}}, outDimTag, true, true, false);
         std::vector<int> curveLoopTags;
         std::vector<std::vector<int> > curveTags;
         gmsh::model::occ::getCurveLoops(wc+100,
                                  curveLoopTags,
                                  curveTags);

         array.at(wc+1)=-(curveLoopTags.at(0));
      }
   */   
   }

   
   int domainTag = gmsh::model::occ::addPlaneSurface(vec, -1);

   //
   //synchronize prior to add physical group.
   //
   gmsh::model::occ::synchronize();
   //
   //add physical groups.
   //
   for(wc=0; wc<nbrwires; wc++)
   {
      char s[10];
      sprintf(s, "wire_%d",  wc+1);
      gmsh::model::addPhysicalGroup(2, {wtag[wc]}, wc+1, s);
   }

   gmsh::model::addPhysicalGroup(2, {domainTag}, nbrwires+1 , "air");
   gmsh::model::addPhysicalGroup(1, {circleTag}, nbrwires+2, "aircontour");

   // We can then generate a 2D mesh...
 //  gmsh::option::setNumber("Mesh.Algorithm", 6);

   gmsh::model::mesh::generate(2);
   gmsh::model::mesh::refine();
   gmsh::model::mesh::refine();

   // glvis can read mesh version 2.2
   gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

   // ... and save it to disk
   gmsh::write("mesh.msh");

   // start gmsh
   gmsh::fltk::run();

   //before leaving.
   gmsh::finalize();

   return 1;
}


/*
int ProximityEffect::CreateMeshFile3()
{
   // Before using any functions in the C++ API, gmsh::must be initialized:
   gmsh::initialize();

   gmsh::model::add("proxmesh");
   int wc;
   
   //assemble the domain entities.
   gmsh::model::occ::addDisk(0, 0, 0, domainRadius, domainRadius, 10);
   
   // prepare the cutting tools which are the wires.

   for(wc=0; wc<nbrwires; wc++)
   {
      if(wiresInfo[wc].type==WireInfo::WireTypes::roundwire)
      {
         gmsh::model::occ::addDisk(wiresInfo[wc].center[0], wiresInfo[wc].center[1], 0, wiresInfo[wc].dimensions[0], wiresInfo[wc].dimensions[0], wc+100);
      }
      else if(wiresInfo[wc].type==WireInfo::WireTypes::rectangular)
      {
         gmsh::model::occ::addRectangle(wiresInfo[wc].center[0] - wiresInfo[wc].dimensions[0]/2.0, wiresInfo[wc].center[1] - wiresInfo[wc].dimensions[1]/2.0, 0, wiresInfo[wc].dimensions[0], wiresInfo[wc].dimensions[1], wc+100, 0);
      }
   }

   gmsh::model::occ::synchronize();  

   std::vector<std::pair<int, int>> cuttingSurface;
   for(wc=0; wc<nbrwires; wc++ )
   {
      cuttingSurface.emplace_back(2, wc+100);
   }
      
   std::vector<std::pair<int, int>> outSurfaces;
   std::vector<gmsh::vectorpair> outDimTagsMap;
   gmsh::model::occ::cut({{2, 10}}, cuttingSurface, outSurfaces, outDimTagsMap);
   gmsh::model::occ::synchronize();

// assemble the wires entities.

   for(wc=0; wc<nbrwires; wc++)
   {
      if(wiresInfo[wc].type==WireInfo::WireTypes::roundwire)
      {
         gmsh::model::occ::addDisk(wiresInfo[wc].center[0], wiresInfo[wc].center[1], 0, wiresInfo[wc].dimensions[0], wiresInfo[wc].dimensions[0], wc+200);
      }
      else if(wiresInfo[wc].type==WireInfo::WireTypes::rectangular)
      {
         gmsh::model::occ::addRectangle(wiresInfo[wc].center[0] - wiresInfo[wc].dimensions[0]/2.0, wiresInfo[wc].center[1] - wiresInfo[wc].dimensions[1]/2.0, 0, wiresInfo[wc].dimensions[0], wiresInfo[wc].dimensions[1], wc+200, 0);
      }
   }

   //
   //synchronize prior to add physical group.
   //
   gmsh::model::occ::synchronize();
   //
   //add physical groups.
   //
   for(wc=0; wc<nbrwires; wc++)
   {
      char s[10];
      sprintf(s, "wire_%d",  wc+1);
      gmsh::model::addPhysicalGroup(2, {wc+200}, wc+1, s);
   }


   gmsh::model::addPhysicalGroup(2, {outSurfaces[0].second}, nbrwires+1 , "air");
   gmsh::model::addPhysicalGroup(1, {nbrwires+2}, nbrwires+2, "aircontour");

// Get all entities
    std::vector<std::pair<int, int>> entities2;
    gmsh::model::getEntities(entities2, -1);

    // Print all entities
    std::cout << "Entities in the model:\n";
    for (const auto &entity : entities2) {
        std::cout << "Dimension: " << entity.first << ", Tag: " << entity.second << "\n";
    }

   // We can then generate a 2D mesh...
   gmsh::option::setNumber("Mesh.Algorithm", 6);

   gmsh::model::mesh::generate(2);
   gmsh::model::mesh::refine();
   gmsh::model::mesh::refine();

   // glvis can read mesh version 2.2
   gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

   // ... and save it to disk
   gmsh::write("mesh.msh");

   // start gmsh
   gmsh::fltk::run();

   //before leaving.
   gmsh::finalize();

   return 1;
}
*/


int ProximityEffect::LoadMeshFile()
{

   // 3. Read the mesh from the given mesh file.
//   if(IsCreatedMeshFile()) mesh = new Mesh("mesh.msh", 1, 1);
//   else 
   mesh = new Mesh(meshFile, 1, 1);
   
   dim = mesh->Dimension();

   cout << mesh->bdr_attributes.Max() << " bdr attr max\n"
        << mesh->Dimension() << " dimensions\n"
        << mesh->GetNV() << " vertices\n"
        << mesh->GetNE() << " elements\n"
        << mesh->GetNBE() << " boundary elements\n"
        << mesh->GetNEdges() << " edges\n";

   nbrel = mesh->GetNE();

   return 1;
}

int ProximityEffect::CreateFESpace()
{
   fec = new H1_FECollection(order, dim);
   fespace = new FiniteElementSpace(mesh, fec);
   nbrdof = fespace->GetNDofs(); 
   cout << fespace->GetNDofs() << " degree of freedom\n"
        << fespace->GetVDim() << " vectors dimension\n\n";   

return 1;
}

int ProximityEffect::CreateEssentialBoundary()
{

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   // real and imag are the same because they refer to mesh nodes.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   assert(mesh->bdr_attributes.Max()==nbrwires+2);
   ess_bdr.SetSize(nbrwires+2);
   ess_bdr = 0;
   ess_bdr[AIRCONTOUR-1]=1;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   {
      std::ofstream out("out/ess_tdof_list.txt");
      ess_tdof_list.Print(out, 10);
   }

//duplicate the essential boundary condition for the imaginary part of Az.
//it shall be + nbrdof to align with the second block.
   ess_tdof_list_block = new Array<int>(2*ess_tdof_list.Size());
   for(int i=0, size=ess_tdof_list.Size(); i<size; i++)
   {
      (*ess_tdof_list_block)[i]=ess_tdof_list[i];
      (*ess_tdof_list_block)[i+size]=ess_tdof_list[i]+nbrdof;
   }

   {
      std::ofstream out("out/ess_tdof_list_block.txt");
      ess_tdof_list_block->Print(out, 10);
   }
   return 1;
}

int ProximityEffect::CreateOperatorA1()
{
   ConstantCoefficient One(1.0);
// note DiffusionIntegrator is "- div(grad(.))".
   BilinearForm BLFA1(fespace);
   BLFA1.AddDomainIntegrator(new DiffusionIntegrator(One));
   BLFA1.Assemble();
   BLFA1.Finalize();
   
   A1 = new SparseMatrix(BLFA1.SpMat());

   std::ofstream out("out/A1.txt");
   A1->Print(out, 10);
  
   cout << A1->Height() << " A1 Height()\n " 
        << A1->Width()  << " A1 Width()\n\n ";

   return 1;
}

int ProximityEffect::CreateOperatorA2()
{
   
   PWConstCoefficient K(2+nbrwires);
   Vector CoeffVector(2+nbrwires);
   CoeffVector = 0.0;
   for(int wc = 0; wc < nbrwires; wc++)
   {
      CoeffVector[wc] = -mu_*omega_*sigma_;
   }   
   K.UpdateConstants(CoeffVector);

   BilinearForm BLFA2(fespace);
   BLFA2.AddDomainIntegrator(new MassIntegrator(K));
   BLFA2.Assemble();
   BLFA2.Finalize();
   
   A2 = new SparseMatrix(BLFA2.SpMat());

   std::ofstream out("out/A2.txt");
   A2->Print(out, 10);
  
   cout << A2->Height() << " A2 Height()\n " 
        << A2->Width()  << " A2 Width()\n\n ";

   return 1;
}

int ProximityEffect::CreateOperatorA3()
{  
   A3 = new SparseMatrix*[nbrwires];
   stringstream ss;
   char fn[256];
   
   for(int wc = 0; wc < nbrwires; wc++)
   {
      Vector CoeffVector(2+nbrwires);
      CoeffVector = 0.0;
      CoeffVector[wc] = -mu_*sigma_;
      PWConstCoefficient K(CoeffVector);

      BilinearForm BLFA3(fespace);
      BLFA3.AddDomainIntegrator(new MassIntegrator(K));
      BLFA3.Assemble();
      BLFA3.Finalize();

      
      sprintf(fn,"out/BLFA3_%d.txt", wc);
      std::ofstream out4(fn);
      BLFA3.SpMat().Print(out4, 10);
      
      SparseMatrix TempSM(BLFA3.SpMat());
      TempSM.Finalize();
      
      Vector TempVEC(nbrdof);
      TempSM.GetRowSums(TempVEC);

      A3[wc] = new SparseMatrix(nbrdof, 1);
      *(A3[wc]) = 0.0;
      
      for(int i=0; i<nbrdof; i++)
      {
         A3[wc]->Add(i, 0, TempVEC[i]);
      }

      A3[wc]->Finalize();

      sprintf(fn,"out/A3_%d.txt", wc );
      std::ofstream out3(fn);
      A3[wc]->Print(out3, 10);
      
      cout << A3[wc]->Height() << " A3 Height()\n " 
         << A3[wc]->Width()  << " A3 Width()\n\n ";
   }
   return 1;
}


int ProximityEffect::CreateOperatorA4()
{
   A4 = new SparseMatrix(*A2);
   *A4 *= -1.0;
   A4->Finalize();
      
   std::ofstream out("out/A4.txt");
   A4->Print(out, 10);

   cout << A4->Height() << " A4 Height()\n " 
        << A4->Width()  << " A4 Width()\n\n ";

   return 1;
}

int ProximityEffect::CreateOperatorA5()
{
   /*
This section of code compute the operator performing the 
integration.

For each element compute the current which is the integral of J x s.
Note s, the conductivity, is a PWCoefficient
*/
   A5 = new SparseMatrix*[nbrwires];
   stringstream ss; 
   char fn[256];
   for(int wc = 0; wc < nbrwires; wc++)
   {
      Vector CoeffVector(2+nbrwires);
      CoeffVector = 0.0;
      CoeffVector[wc] = omega_*sigma_;
      PWConstCoefficient K(CoeffVector);
         
      //surface integral.
      LinearForm LFA5(fespace);
      LFA5.AddDomainIntegrator(new DomainLFIntegrator(K));
      LFA5.Assemble();
      
      sprintf(fn,"out/LFA5_%d.txt", wc );
      std::ofstream out1(fn);
      LFA5.Print(out1, 10);

      A5[wc] = new SparseMatrix(1, nbrdof);
      for(int k=0; k<nbrdof; k++)
      {
         A5[wc]->Set(0, k, LFA5[k]);
      }

      A5[wc]->Finalize();

      sprintf(fn,"out/A5_%d.txt", wc );
      std::ofstream out2(fn);
      A5[wc]->Print(out2, 10);

      cout << A5[wc]->Height() << " A5 Height()\n " 
         << A5[wc]->Width()  << " A5 Width()\n\n ";
   }

   return 1;
}


int ProximityEffect::CreateOperatorA6()
{
   A6 = new SparseMatrix*[nbrwires];
   stringstream ss; 
   char fn[256];
   for(int wc = 0; wc < nbrwires; wc++)
   {
    
      ConstantCoefficient One(1.0);
      real_t wireArea = IntegrateScalar(*fespace, One, wc+1);
      A6[wc] = new SparseMatrix(1, 1);
      A6[wc]->Set(0, 0, sigma_ * wireArea);
      A6[wc]->Finalize();

      sprintf(fn,"out/A6_%d.txt", wc );
      std::ofstream out(fn);
      A6[wc]->Print(out, 10);

      cout << A6[wc]->Height() << " A6 Height()\n " 
         << A6[wc]->Width()  << " A6 Width()\n\n ";
   }
   return 1;
}


int ProximityEffect::CreateOperatorA7()
{
   A7 = new SparseMatrix*[nbrwires];

   stringstream ss; 
   char fn[256];
   for(int wc = 0; wc < nbrwires; wc++)
   {
      A7[wc] = new SparseMatrix(*(A5[wc]));
      *(A7[wc]) *= -1.0;
      A7[wc]->Finalize();

      sprintf(fn,"out/A7_%d.txt", wc );        
      std::ofstream out(ss.str());
      A7[wc]->Print(out, 10);

      cout << A7[wc]->Height() << " A7 Height()\n " 
         << A7[wc]->Width()  << " A7 Width()\n\n ";
   }
   return 1;
}


int ProximityEffect::CreaterhsVector()
{
// computing rhs 
   rhs = new Vector(2*nbrdof+2*nbrwires);
   *rhs = 0.0;
   for(int wc=0; wc < nbrwires; wc++)
   {
      (*rhs)[2*nbrdof + 2*wc + 0] = wiresInfo[wc].current[0] * cos(2*M_PI*wiresInfo[wc].current[1]/360.0);
      (*rhs)[2*nbrdof + 2*wc + 1] = wiresInfo[wc].current[0] * sin(2*M_PI*wiresInfo[wc].current[1]/360.0);;     
   }

   std::ofstream out("out/rhs.txt");
   rhs->Print(out, 10);

   cout << rhs->Size() << " rhs size\n\n ";
   return 1;

}

int ProximityEffect::CreatexVector()
{
   x = new Vector(2*nbrdof+2*nbrwires);
   *x = 0.0; 
   std::ofstream out("out/x.txt");
   x->Print(out, 10);
   cout << x->Size() << " x size\n\n ";
   return 1;
}

int ProximityEffect::CreateBlockOperator()
{

// 6. Define the BlockStructure of the problem
   int blockoffsetsize = 2 + 2 * nbrwires + 1;
   blockOffset = new Array<int>(blockoffsetsize);  //n+1
   (*blockOffset)[0]=0;
   (*blockOffset)[1]=nbrdof; 
   (*blockOffset)[2]=nbrdof;
   for(int i = 3; i<blockoffsetsize; i++) (*blockOffset)[i]=1;
   blockOffset->PartialSum();
   {
      std::ofstream out("out/blockOffset.txt");
      blockOffset->Print(out, 10);
   }
  
   // 8. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
  
  Device device("cpu");
  MemoryType mt = device.GetMemoryType();

   ProxOp = new BlockOperator(*blockOffset);

   
// Build the operator
// row 0 ...
      ProxOp->SetBlock(0, 0, A1);
      ProxOp->SetBlock(0, 1, A2);
      for(int i = 0; i<nbrwires; i++) ProxOp->SetBlock(0, 2+2*i, A3[i]);
// row 1
      ProxOp->SetBlock(1, 0, A4);
      ProxOp->SetBlock(1, 1, A1);
      for(int i = 0; i<nbrwires; i++) ProxOp->SetBlock(1, 3+2*i, A3[i]);
// col 0
      for(int i = 0; i<nbrwires; i++) ProxOp->SetBlock(3+2*i, 0, A7[i]);
// col 1
      for(int i = 0; i<nbrwires; i++) ProxOp->SetBlock(2+2*i, 1, A5[i]);
// diagonal (2, 2)
      for(int i = 0; i<nbrwires; i++)
      {
         ProxOp->SetBlock(2+2*i, 2+2*i, A6[i]);
         ProxOp->SetBlock(3+2*i, 3+2*i, A6[i]);
      }

      {
      std::ofstream out("out/ProxOp.txt");
      ProxOp->PrintMatlab(out);
      }

      assert(ProxOp->Height() == 2*nbrdof+2*nbrwires);
      assert(ProxOp->Width() == 2*nbrdof+2*nbrwires);
      
//DL241125: I check ProxOp it contains all the BLFA1 to 4 in proper order.

      assert(2*A1->NumRows()+2*nbrwires==ProxOp->NumRows());
      assert(2*A1->NumCols()+2*nbrwires==ProxOp->NumCols());


   A = new BlockOperator(*blockOffset);
   A_ptr = A;

   B = new BlockVector(*blockOffset, mt);
   X = new BlockVector(*blockOffset, mt);

   // note the A_ptr do not point at the same operator after !!!
   ProxOp->FormLinearSystem(*ess_tdof_list_block, *x, *rhs, A_ptr, *X, *B);

   {
      std::ofstream out("out/X.txt");
      X->Print(out, 10);
   }

   {
      std::ofstream out("out/B.txt");
      B->Print(out, 10);
   }

   {
      std::ofstream out("out/A_ptr.txt");
      A_ptr->PrintMatlab(out);
   }

   cout << A_ptr->Height() << " Operator Height()\n " 
         << A_ptr->Width()  << " Operator Width()\n "
         << X->Size() << " X.Size()\n " 
         << B->Size() << " B.Size()\n\n ";

   return 1;
}
   
int ProximityEffect::CreatePreconditionner()
{

   // 10. Construct the operators for preconditioner

   // ************************************
   // Here I have no idea how to build a preconditionner.
   // ex5 is not appropriated as per my understanding.
   // ************************************

   // Create smoothers for diagonal blocks
   int wc;
   GSSmoother gs1(*A1); // Gauss-Seidel smoother for A11
   GSSmoother gs2(*A1); // Gauss-Seidel smoother for A22
   GSSmoother* gs[20];
   for(wc=0;wc<nbrwires; wc++)
   {
      gs[2*wc+0] = new GSSmoother(*A6[wc]);
      gs[2*wc+1] = new GSSmoother(*A6[wc]);
   }

   block_prec = new BlockDiagonalPreconditioner(*blockOffset);
   block_prec->SetDiagonalBlock(0, &gs1); // Set smoother for A11
   block_prec->SetDiagonalBlock(1, &gs2); // Set smoother for A22
   for(wc=0;wc<nbrwires; wc++)
   {
      block_prec->SetDiagonalBlock(2+2*wc+0, gs[wc]);
      block_prec->SetDiagonalBlock(2+2*wc+1, gs[wc]);
   }
   return 1;
}

int ProximityEffect::Solver()
{
   // Solve system Ax = b
   GMRESSolver solver1;
   solver1.SetOperator(*A_ptr);
   //solver1.SetPreconditioner(*block_prec);
   solver1.SetRelTol(1e-16);
   //   solver.SetAbsTol(1e-8);
   solver1.SetMaxIter(5000);
   solver1.SetPrintLevel(1);

   //x = 0.0;       // Initial guess
   solver1.Mult(*B, *X);

   

   A->RecoverFEMSolution(*X, *rhs, *x);
   {
      std::ofstream out("out/xsol.txt");
      x->Print(out, 10);
   }

   return 1;
   
}


int ProximityEffect::PostPrecessing()
{
   AzrGF = new GridFunction(fespace);
   AziGF = new GridFunction(fespace);
  
   // rebuild GFR and GFI from x.
   AzrGF->MakeRef(fespace, *x, 0);
   AziGF->MakeRef(fespace, *x, nbrdof);

   return 1;

   GridFunctionCoefficient AzrGFCoeff(AzrGF);
   GridFunctionCoefficient AziGFCoeff(AziGF);

   // compute Jr 

   Vector CoeffVector(2+nbrwires);
   CoeffVector = 0.0;
   for(int wc=0; wc<nbrwires; wc++) CoeffVector[wc] = omega_*sigma_;
   PWConstCoefficient K1(CoeffVector);
   
   ProductCoefficient Jr1Coeff(K1, AziGFCoeff);

   PWConstCoefficient K2;
   
   PWConstCoefficient Jr2Coeff;
   Vector CoeffVector2(2+nbrwires);
   CoeffVector2 = 0.0;
   for(int wc=0; wc<nbrwires; wc++) CoeffVector2[wc] = sigma_*(*x)[2*nbrdof+wc*2];
   Jr2Coeff.UpdateConstants(CoeffVector2);
  
   SumCoefficient JrCoeff(Jr1Coeff, Jr2Coeff);

   // Compute Ji
   PWConstCoefficient K3;   
   Vector CoeffVector3(2+nbrwires);
   CoeffVector3 = 0.0;
   for(int wc=0; wc<nbrwires; wc++) CoeffVector3[wc] = -omega_*sigma_;
   K3.UpdateConstants(CoeffVector3);
   
   ProductCoefficient Ji1Coeff(K3, AzrGFCoeff);

   PWConstCoefficient Ji2Coeff;
   Vector CoeffVector4(2+nbrwires);
   CoeffVector4 = 0.0;
   for(int wc=0; wc<nbrwires; wc++) CoeffVector4[wc] = sigma_*(*x)[2*nbrdof+2*wc+1];
   Ji2Coeff.UpdateConstants(CoeffVector4);
   
   SumCoefficient JiCoeff(Ji1Coeff, Ji2Coeff);

   //Compute J
   PowerCoefficient jrSquareCoeff(JrCoeff, 2.0);
   PowerCoefficient jiSquareCoeff(JiCoeff, 2.0);
   SumCoefficient JSquare(jrSquareCoeff, jiSquareCoeff);
   PowerCoefficient JCoeff(JSquare, 0.5);

   JFec = new DG_FECollection(order, dim);
   JFESpace = new FiniteElementSpace(mesh, JFec);

   JiGF = new GridFunction(JFESpace);
   JiGF->ProjectCoefficient(JiCoeff);
   JrGF = new GridFunction(JFESpace);
   JrGF->ProjectCoefficient(JrCoeff);
   JGF = new GridFunction(JFESpace);
   JGF->ProjectCoefficient(JCoeff);

   return 1;
}

int ProximityEffect::DisplayResults()
{

   Glvis(mesh, AzrGF, "A-Real" );
   Glvis(mesh, AziGF, "A-Imag" );

   return 1;

   Glvis(mesh, JrGF, "J-Real" );
   Glvis(mesh, JiGF, "J-Imag" );
   Glvis(mesh, JGF, "J" );

   ConstantCoefficient One(1.0);
   
   cout << "\n";
   for(int wc=0; wc<nbrwires; wc++)
   {
      cout << "wire " << wc << ": "<< (*x)[2*nbrdof+2*wc] << " Vr, "
           << (*x)[2*nbrdof+2*wc+1] << " Vi, "
           << sqrt(pow((*x)[2*nbrdof+2*wc], 2)+pow((*x)[2*nbrdof+2*wc+1], 2)) << " V\n";
   
   real_t WireArea = IntegrateScalar(*fespace, One, wc+1);
   real_t Rdc = 1.0/(WireArea * sigma_);
   real_t Rac = (*x)[2*nbrdof+2*wc] / (*rhs)[2*nbrdof+2*wc] ;
   real_t RacdcRatio = Rac/Rdc;
   real_t Lw = (*x)[2*nbrdof+2*wc+1] / ((*rhs)[2*nbrdof+2*wc]  * omega_);

   cout << Rdc << " Rdc\n"
        << Rac << " Rac\n"
        << RacdcRatio << " AC DC Ratio at " << omega_/2.0/M_PI << "Hz\n"
        << 1E6*Lw << " L uH\n\n";
   }
   return 1;
}

int ProximityEffect::CleanOutDir()
{
    system("rm -f out/*");
    return 1;
}



int main(int argc, char *argv[])
{

   ProximityEffect PE;

   PE.CleanOutDir();
   
   PE.Parser(argc, argv);

   PE.ReadConfigFile();

//   if(PE.IsCreatedMeshFile())
     PE.CreateMeshFile();

   PE.LoadMeshFile();
  
   PE.CreateFESpace();

   PE.CreateEssentialBoundary();

   PE.CreateOperatorA1();
   PE.CreateOperatorA2();
   PE.CreateOperatorA3();
   PE.CreateOperatorA4();
   PE.CreateOperatorA5();
   PE.CreateOperatorA6();
   PE.CreateOperatorA7();

   PE.CreaterhsVector();
   PE.CreatexVector();

   PE.CreateBlockOperator();

 //  PE.CreatePreconditionner();

   PE.Solver();

   PE.PostPrecessing();

   PE.DisplayResults();

   cout << "return to end\n";
   
   return 0;
}

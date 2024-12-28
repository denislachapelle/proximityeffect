//                                proximity2dblockint
//                                inspired from proximity2dblock
//                                inspired from proximity2d
//                                based on MFEM Example 22 prob 1 (case 0), ex5...
//
// Compile with: make proximity2dblock, need MFEM version 4.7 and GLVIS-4.3.
//
// Sample runs:  ./proximity2dblock -m ../mesh/ProxRoundRoundWires2d.msh
//
/*
Description:  

Assuming no displacement current, the equation for magnetic vector potential is 
Div(grad(Az)) + i w u s Az = 0

div(grad((Az)) + i c Az = 0 with c = w u s

and then split in real and imag parts
div(grad(Azr)) + i div(grad((Azi)) + i c Azr - c Azi = 0
and in two equations
div(grad(Azr)) - c Azi = 0
div(grad((Azi)) + c Azr = 0

Now computing Jz.
Jz = 1/u curl(curl(A))

Assembly will be ...

Let define matrix A1, A2, A3 and A4 as implementing ....
A1 div(grad((Azr)) + b Azr
A2 div(grad(Azi)) + b Azi
A3 -c Azi
A4 c Azr

Then we can write the assembled matrix...

[A1 A3] [Azr] = [0]
[A4 A2] [Azi] = [0]

A fix A-Field is set at the wire countour.

The result I get is an A-field which respect the equation.
But computing J is again quite complicated.
J = -1/u div(grad(A))
To compute J another system need to be solved.
This is not done yet.

I learn a lot with block structure.

This is not a good approach because we do not know the values 
of the BC to get the required current.

Starting from proximity2dblock I add a bilinear to compute J and another block to compute 
the integral of J over a given domain section which is the current a a given wire.

[A1 A3 0  0]   [Azr] = [0]
[A4 A2 0  0]   [Azi] = [0]
[A5 A7 0  0]   [Jzr] = [0]
[A8 A6 0  0]   [Jzi] = [0]
[0  0  A9 A10]       = [Ir]
[0  0 A11 A12]       = [Ii]

If the system as N dof, A1 to A9 will be NxN, A9 to A12 will be 1xN.

I cannot for Ir and Ii since I don't know theiur relation. I should constraint
sqrt(I1r^2 + I1i^2) = I1. I have no idea how to do this.
Wire many wires I will have to constraint many different current... such as 
sqrt(I2r^2 + I2i^2) = I2.

When more than one wires extra row will be added instead of A9 and A10 for wire 1
it will be A11 and A12 for wire 2.

Here i am stuck ....

I think NewtonSolver is the way to go.

u: permeability.
e: permitivity.
s: conductivity.
i: sqrt(-1)
*/

#include <mfem.hpp>
#include <linalg/hypre.hpp>
#include "mytools.hpp"
#include <fstream>
#include <iostream>
#include <math.h>

#define wire_1 1
#define wire_2 2
#define air 3
#define aircontour 4
#define wire1contour 5
#define wirecenter 6

using namespace std;
using namespace mfem;

static double mu_ = 1.257E-6;
static double epsilon_ = 8.854E-12;
static double sigma_ = 58E6;
static double omega_ = 2.0*M_PI*60;


int main(int argc, char *argv[])
{
   StopWatch chrono;
   tic();
   // 1. Parse command-line options.
   const char *mesh_file = "../mesh/proxroundroundwires2d.msh";
   int ref_levels = 0;
   int order = 1;
   int prob = 0;
   double freq = -1.0;
   double a_coef = 0.0;
   bool visualization = 1;
   bool herm_conv = true;
   bool pa = false;
   bool noDispl = true;
   const char *device_config = "cpu";
   double alpha = 1e-10; // penalty scaling factor.

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
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
   args.AddOption(&herm_conv, "-herm", "--hermitian", "-no-herm",
                  "--no-hermitian", "Use convention for Hermitian operators.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&noDispl, "-nd", "--no-displ", "-dincl", "--displ-include",
                  "Do not include displacement current in the computation.");
   args.AddOption(&alpha, "-al", "--alpha", 
                  "Penalty scaling factor.");
                  
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   int nbrWire = HowManyWire("../mesh/proxmeshconfig.txt");
   cout << "\n nbrWire = " << nbrWire <<endl;

   if ( freq >= 0.0 )
   {
      omega_ = 2.0 * M_PI * freq;
   }

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
   //    with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1); // 
   int dim = mesh->Dimension();

   // 4. Refine the mesh to increase resolution. In this example we do
   //    'ref_levels' of uniform refinement where the user specifies
   //    the number of levels with the '-r' option.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   FiniteElementCollection *fecR = NULL;
   fecR = new H1_FECollection(order, dim);
   FiniteElementSpace *fespaceR = new FiniteElementSpace(mesh, fecR);

   FiniteElementCollection *fecI = NULL;
   fecI = new H1_FECollection(order, dim);
   FiniteElementSpace *fespaceI = new FiniteElementSpace(mesh, fecI);  

   cout << "\n" 
        << fespaceR->GetTrueVSize() << " fespaceR->GetTrueVSize()\n " 
        << fespaceI->GetTrueVSize() << " fespaceI->GetTrueVSize()\n ";
        
   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined based on the type
   //    of mesh and the problem type.

   // real and imag are the same because they refr to mesh nodes.
   Array<int> ess_tdof_list_R;
   Array<int> ess_bdr_R;
   assert(mesh->bdr_attributes.Max()==5);
   ess_bdr_R.SetSize(mesh->bdr_attributes.Max());
   ess_bdr_R = 0;
   ess_bdr_R[wire1contour-1]=1;
   fespaceR->GetEssentialTrueDofs(ess_bdr_R, ess_tdof_list_R);

   {
      std::ofstream out("out/ess_tdof_list_R.txt");
      ess_tdof_list_R.Print(out, 10);
   }

   // the diffusion integrator is negative, -div(grad)).
   ConstantCoefficient a_(-1.0);
   ConstantCoefficient b_(-omega_ * omega_ * epsilon_*mu_);

   PWConstCoefficient *c_ = new PWConstCoefficient;
   {
      double CoeffArray[]={omega_*mu_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      c_->UpdateConstants(CoeffVector);
   }
   
   PWConstCoefficient *negc_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*mu_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      negc_->UpdateConstants(CoeffVector);
   }

   PWConstCoefficient *d_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      d_->UpdateConstants(CoeffVector);
   }

/* current density in conductor left and conductor right.
   PWConstCoefficient *CurrentDensityCoeff = new PWConstCoefficient;
   {
      double CoeffArray[]={1.0, 0.0, 0.0, 0.0, 0.0};  // wire_1 1A/m2, wire_2 -1A/m2.
      Vector CoeffVector(CoeffArray, 5);
      CurrentDensityCoeff->UpdateConstants(CoeffVector);
   }
*/

   // 9. Assemble the finite element matrices for the operator
   
   BilinearForm BLFA1(fespaceR);
   BLFA1.AddDomainIntegrator(new DiffusionIntegrator(a_));
   if (!noDispl) BLFA1.AddDomainIntegrator(new MassIntegrator(b_));
   BLFA1.Assemble();
   BLFA1.Finalize();

   {
      std::ofstream out("out/BLFA1.txt");
      const SparseMatrix &SpPtr = BLFA1.SpMat();
      SpPtr.Print(out, 10);
   }
   
   BilinearForm BLFA2(fespaceI);
   BLFA2.AddDomainIntegrator(new DiffusionIntegrator(a_));
   if (!noDispl) BLFA2.AddDomainIntegrator(new MassIntegrator(b_));
   BLFA2.Assemble();
   BLFA2.Finalize();

   {
      std::ofstream out("out/BLFA2.txt");
      const SparseMatrix &SpPtr = BLFA2.SpMat();
      SpPtr.Print(out, 10);
   }

   BilinearForm BLFA3(fespaceI);
   BLFA3.AddDomainIntegrator(new MassIntegrator(*negc_)); 
   BLFA3.Assemble();
   BLFA3.Finalize();

   {
      std::ofstream out("out/BLFA3.txt");
      const SparseMatrix &SpPtr = BLFA3.SpMat();
      SpPtr.Print(out, 10);
   }

   BilinearForm BLFA4(fespaceR);
   BLFA4.AddDomainIntegrator(new MassIntegrator(*c_)); 
   BLFA4.Assemble();
   BLFA4.Finalize();
   
   {
      std::ofstream out("out/BLFA4.txt");
      const SparseMatrix &SpPtr = BLFA4.SpMat();
      SpPtr.Print(out, 10);
   }

   LinearForm LFR(fespaceR);
   //LFR.AddDomainIntegrator(new DomainLFIntegrator(*CurrentDensityCoeff));
   LFR.Assemble();

   {
      std::ofstream out("out/LFR.txt");
      LFR.Print(out, 10);
   }
   
   
   LinearForm LFI(fespaceI);
   LFI.Assemble();

   {
      std::ofstream out("out/LFI.txt");
      LFI.Print(out, 10);
   }
      cout << BLFA1.Height() << " BLFA1.Height()\n " 
         << BLFA1.Width()  << " BLFA1.Width()\n "

         << BLFA2.Height() << " BLFA2.Height()\n " 
         << BLFA2.Width()  << " BLFA2.Width()\n "
         
         << BLFA3.Height() << " BLFA3.Height()\n " 
         << BLFA3.Width()  << " BLFA3.Width()\n "
         
         << BLFA4.Height() << " BLFA4.Height()\n " 
         << BLFA4.Width()  << " BLFA4.Width()\n "

         << LFR.Size() << " LFR.Size()\n " 
         << LFI.Size() << " LFI.Size()\n ";
        
// 12. Create the grid functions u and p. Compute the L2 error norms.
   GridFunction GFR(fespaceR), GFI(fespaceI);
   GFR=0.0;
   GFI=0.0;

   ConstantCoefficient One(1.0);
   GFR.ProjectBdrCoefficient(One, ess_bdr_R);
   //GFI.ProjectBdrCoefficient(Zero, ess_bdr_R);

//   SparseMatrix A1, A2, A3, A4;


// 6. Define the BlockStructure of the problem.
   Array<int> blockRowOffset(3);
   blockRowOffset[0]=0;
   blockRowOffset[1]=BLFA1.NumRows(); 
   blockRowOffset[2]=2*BLFA1.NumRows();

   assert(BLFA1.NumRows()==BLFA1.Height());

   Array<int> blockColOffset(3);
   blockColOffset[0]=0;
   blockColOffset[1]=BLFA1.NumCols();
   blockColOffset[2]=2*BLFA1.NumCols();

   assert(BLFA1.NumCols()==BLFA1.Width());
   
   // 8. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
   MemoryType mt = device.GetMemoryType();
   BlockVector x(blockRowOffset, mt), rhs(blockRowOffset, mt);

   BlockOperator ProxOp(blockRowOffset, blockColOffset);

      SparseMatrix SpBLFA1(BLFA1.SpMat());
      SparseMatrix SpBLFA2(BLFA2.SpMat());
      SparseMatrix SpBLFA3(BLFA3.SpMat());
      SparseMatrix SpBLFA4(BLFA4.SpMat());

      {
      std::ofstream out("out/SpBLFA1.txt");
      SpBLFA1.Print(out, 10);
      }

      ProxOp.SetBlock(0,0, &SpBLFA1);
      ProxOp.SetBlock(0,1, &SpBLFA3);
      ProxOp.SetBlock(1,0, &SpBLFA4);
      ProxOp.SetBlock(1,1, &SpBLFA2);

      {
      std::ofstream out("out/ProxOp.txt");
      ProxOp.PrintMatlab(out);
      }

//DL241125: I check ProxOp it contains all the BLFA1 to 4 in proper order.

      assert(2*BLFA1.NumRows()==ProxOp.NumRows());
      assert(2*BLFA1.NumCols()==ProxOp.NumCols());


 //  GFR.GetTrueDofs(x.GetBlock(0)); // Copy x1 into the first block of X
 //  GFI.GetTrueDofs(x.GetBlock(1)); // Copy x2 into the second block of X
   
   for(int i=0; i<BLFA1.Height(); i++)
   {
      rhs[i]=LFR[i];
      rhs[i+BLFA1.Height()]=LFI[i];

      x[i]=GFR[i];
      x[i+BLFA1.Height()]=GFI[i];
   }

   {
      std::ofstream out("out/rhs.txt");
      rhs.Print(out, 10);
   }

   {
      std::ofstream out("out/x.txt");
      x.Print(out, 10);
   }

   Array<int> ess_tdof_list_block(2*ess_tdof_list_R.Size());


   for(int i=0, size=ess_tdof_list_R.Size(); i<size; i++)
   {
      ess_tdof_list_block[i]=ess_tdof_list_R[i];
      ess_tdof_list_block[i+size]=ess_tdof_list_R[i]+BLFA1.Height();
   }

   {
      std::ofstream out("out/ess_tdof_list_block.txt");
      ess_tdof_list_block.Print(out, 10);
   }



   BlockOperator A(blockRowOffset, blockColOffset);
   mfem::Operator *A_ptr = &A;
   BlockVector B(blockRowOffset, mt), X(blockRowOffset, mt);

   // note the A_ptr do not point at the same operator after !!!
   ProxOp.FormLinearSystem(ess_tdof_list_block, x, rhs, A_ptr, X, B);

   {
      std::ofstream out("out/X1.txt");
      X.Print(out, 10);
   }

   {
      std::ofstream out("out/B.txt");
      B.Print(out, 10);
   }

   {
      std::ofstream out("out/A_ptr.txt");
      A_ptr->PrintMatlab(out);
   }

      cout << A_ptr->Height() << " A_ptr->Height()\n " 
           << A_ptr->Width()  << " _ptr->Width()\n "

         << X.Size() << " X.Size()\n " 
         << B.Size() << " B.Size()\n ";


   // 10. Construct the operators for preconditioner

// ************************************
// Here I have no idea how to build a preconditionner.
// ex5 is not appropriated as per my understanding.
// ************************************

SparseMatrix A11_mat = BLFA1.SpMat();
SparseMatrix A22_mat = BLFA2.SpMat();

// Create smoothers for diagonal blocks
GSSmoother gs1(A11_mat); // Gauss-Seidel smoother for A11
GSSmoother gs2(A22_mat); // Gauss-Seidel smoother for A22

DSmoother ds1(A11_mat); // Diagonal smoother for A11
DSmoother ds2(A22_mat); // Diagonal smoother for A22

BlockDiagonalPreconditioner block_prec(blockRowOffset);
block_prec.SetDiagonalBlock(0, &gs1); // Set smoother for A11
block_prec.SetDiagonalBlock(1, &gs2); // Set smoother for A22

   chrono.Clear();
   chrono.Start();



   if(1)
   {
     // Solve system Ax = b
      GMRESSolver solver;
      solver.SetOperator(*A_ptr);
     // solver.SetPreconditioner(block_prec);
      solver.SetRelTol(1e-6);
      solver.SetAbsTol(1e-8);
      solver.SetMaxIter(10000);
      solver.SetPrintLevel(1);

      //x = 0.0;       // Initial guess
      solver.Mult(B, X);

   }
   chrono.Stop();

A.RecoverFEMSolution(X, rhs, x);

   PWConstCoefficient *f_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*sigma_, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 4);
      f_->UpdateConstants(CoeffVector);
   }
   
   PWConstCoefficient *negf_ = new PWConstCoefficient;
   {
      double CoeffArray[]={omega_*sigma_, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 4);
      negf_->UpdateConstants(CoeffVector);
   }

// rebuild GFR and GFI from x.
GFR.MakeRef(fespaceR, x.GetBlock(0), 0);
GFI.MakeRef(fespaceI, x.GetBlock(1), 0);

Glvis(mesh, &GFR, "A-field: Real Part" );
Glvis(mesh, &GFI, "A-field: Imag Part" );
 
   // 15. Free the used memory.
delete fecR;
delete fecI;
delete fespaceR;
delete fespaceI;
delete mesh;
delete c_;
delete negc_;
delete d_;

   cout << "\n time elapsed = " << toc() << endl;

   return 0;
}

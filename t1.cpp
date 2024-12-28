//                                proximity2dblockb
//                                inspired from proximity2dblock
//                                inspired from proximity2d
//                                based on MFEM Example 22 prob 1 (case 0), ex5...
//
// Compile with: make proximity2dblockb, need MFEM version 4.7 and GLVIS-4.3.
//
// Sample runs:  ./proximity2dblock -m ../mesh/ProxRoundRoundWires2d.msh
//
/*
Description:  

Implementation suggested by paul hilsher in https://github.com/mfem/mfem/issues/4584
from paper "Specialized conductor models for finite element eddy current simulation".
https://www.iem.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaabcfymi

The equation is 
- div(grad(Az)) + i w u s Az + u s V = 0

Az : magnetic vector potential z-component.
V : voltage potential.
u : permeability.
s : conductivity.
w : pulsation, 2 pi f.
Assume 1m long, so s V is electric field.

then split in real and imag parts
- div(grad(Azr)) - w u s Azi + u s Vr = 0
- div(grad(Azi)) + w u s Azr + u s Vi = 0

Let define matrix A1, A2, A3 and A4 as implementing ....
A1 -div(grad((.))
A2 - u w s (.)
A3 u s (.)
A4 w u s (.)

We also want to enforce the current ..
intS(- i s w A) + V/R = I

then split real and imaginary equation...
intS(s w Ai) + Vr/R = Ir
intS(- s w Ar) + Vi/R = Ii

let define A5 and A6 as implementing ...
A5 intS(s w (.))
A6 1/R (.)
A7 intS(-s w (.))

Then we can write the assembled matrix...

[A1 A2 A3 0 ] [Azr] = [0]
[A4 A1 0  A3] [Azi] = [0]
[0  A5 A6 0 ] [Vr]  = [Ir]
[A7 0  0  A6] [Vi]  = [Ii]

[dofxdof dofxdof dofx1 dofx1] [dofx1] = [dofx1]
[dofxdof dofxdof dofx1 dofx1] [dofx1] = [dofx1]
[1xdof   1xdof   1x1   1x1]   [dofx1] = [dofx1]
[1xdof   1xdof   1x1   1x1]   [dofx1] = [dofx1]

Ir being total real current in wire 1.

Ii being total imaginary current in wire 1.

Computing A5 and A7...
 
Bilinear, MassIntegrator, sum column to form a row vector
which one represent the surface of each element.

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

void SMSumColumns(const SparseMatrix &mat, Vector &vec)
{
   int num_rows = mat.NumRows();
   for(int i = 0; i<num_rows; i++) 
   {
      // Get the data for the specified row
      const int *columns = mat.GetRowColumns(i);
      const double *values = mat.GetRowEntries(i);
      int row_size = mat.RowSize(i);
      // Search for the column index
      for (int k = 0; k < row_size; k++)
      {
         vec(columns[k]) += mat.Elem(i, k);
      }
   }
   return;    
}

double ComputeElementArea(FiniteElementSpace *fes, int element_index)
{
    ElementTransformation *trans = fes->GetMesh()->GetElementTransformation(element_index);
    const IntegrationRule &ir = IntRules.Get(fes->GetMesh()->GetElementBaseGeometry(element_index),
                                             fes->GetElementOrder(element_index) * 2);

    double area = 0.0;
    for (int i = 0; i < ir.GetNPoints(); i++) {
        const IntegrationPoint &ip = ir.IntPoint(i);
        trans->SetIntPoint(&ip);
        double detJ = trans->Jacobian().Det();
        area += ip.weight * detJ;
    }
    return area;
}

// Custom Coefficient Class
class ElementCoefficient : public Coefficient
{
private:
    int target_element; // The ID of the target element

public:
    // Constructor
    ElementCoefficient(int elem_id) : target_element(elem_id) {}

    // Override Eval to return 1.0 in the target element and 0.0 elsewhere
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
    {
        return (T.ElementNo == target_element) ? 1.0 : 0.0;
    }
};

class ProximityEffect
{
   private:
   // 1. Parse command-line options.
   const char *mesh_file = "../mesh/proxroundroundwires2d.msh";
   int ref_levels = 0;
   int order = 1;
   double freq = -1.0;

   Mesh *mesh;
   int dim;
   int nbrel;

   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   int nbrdof;

   Array<int> *ess_tdof_list_block;

   SparseMatrix *A1, *A2, *A3, *A4, *A5, *A6, *A7;


   public:
      int Parser(int argc, char *argv[]);
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

      
      int CreateRhsVector();
      int CreateXVector();
   
      


};

int ProximityEffect::Parser(int argc, char *argv[])
{

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
   
                 
   
 
   args.Parse();

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device("cpu");
   device.Print();

   if ( freq >= 0.0 ) omega_ = 2.0 * M_PI * freq;

   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);
   return 0;
}

int ProximityEffect::LoadMeshFile()
{

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
   //    with the same code.
   mesh = new Mesh(mesh_file, 1, 1); // 
   dim = mesh->Dimension();

   // 4. Refine the mesh to increase resolution. In this example we do
   //    'ref_levels' of uniform refinement where the user specifies
   //    the number of levels with the '-r' option.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   cout << mesh->Dimension() << " dimensions\n"
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
}

int ProximityEffect::CreateEssentialBoundary()
{

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined based on the type
   //    of mesh and the problem type.

   // real and imag are the same because they refer to mesh nodes.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   assert(mesh->bdr_attributes.Max()==5);
   ess_bdr.SetSize(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[aircontour-1]=1;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   {
      std::ofstream out("out/ess_tdof_list_R.txt");
      ess_tdof_list.Print(out, 10);
   }

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
}

int ProximityEffect::CreateOperatorA1()
{
   ConstantCoefficient One(1.0);

   BilinearForm BLFA1(fespace);
   BLFA1.AddDomainIntegrator(new DiffusionIntegrator(One));
   BLFA1.Assemble();
   BLFA1.Finalize();
   
   std::ofstream out("out/BLFA1.txt");
   const SparseMatrix &SpPtr = BLFA1.SpMat();
   SpPtr.Print(out, 10);
  
   A1 = new SparseMatrix(BLFA1.SpMat());
   return 1;
}

int ProximityEffect::CreateOperatorA2()
{
   
   double CoeffArray[]={-mu_*omega_*sigma_, 0.0, 0.0, 0.0, 0.0};
   Vector CoeffVector(CoeffArray, 5);
   PWConstCoefficient K(CoeffVector);
   
   BilinearForm BLFA2(fespace);
   BLFA2.AddDomainIntegrator(new MassIntegrator(K));
   BLFA2.Assemble();
   BLFA2.Finalize();
   
   std::ofstream out("out/BLFA2.txt");
   const SparseMatrix &SpPtr = BLFA2.SpMat();
   SpPtr.Print(out, 10);

   A2 = new SparseMatrix(BLFA2.SpMat());

   return 1;
}

int ProximityEffect::CreateOperatorA3()
{
   //computing SMA3
   // for each element
   //   if attr = wire 1 all dof = K5 otherwise 0
   A3 = new SparseMatrix(nbrdof, 1);
   for(int el=0; el<nbrel; el++)
   {
      if(fespace->GetAttribute(el) == wire_1)
      {
         const FiniteElement* elem = fespace->GetFE(el);
         Array<int> ElGlobalDofs(elem->GetDof());
         fespace->GetElementDofs(el, ElGlobalDofs);
         for(int k=0; k<ElGlobalDofs.Size(); k++)
            A3->Set(ElGlobalDofs[k], 0, mu_ * sigma_);
      }
   }
   
   std::ofstream out("out/SMA3.txt");
   A3->Print(out, 10);

}


int ProximityEffect::CreateOperatorA4()
{
   A4 = new SparseMatrix(*A2);
   *A4= -1.0;
      
   std::ofstream out("out/A4.txt");
   A4->Print(out, 10);
}

int ProximityEffect::CreateOperatorA5()
{
   /*
This section of code compute the operator performing the 
integration.

For each element compute the current which is the integral of J x s.
Note s, the conductivity, is a PWCoefficient
*/

 
   double CoeffArray[]={omega_*sigma_, 0.0, 0.0, 0.0, 0.0};
   Vector CoeffVector(CoeffArray, 5);
   PWConstCoefficient K(CoeffVector);
 
   //A5 method #2, linearform.
   //surface integral.
   LinearForm LFA5(fespace);
   LFA5.AddDomainIntegrator(new DomainLFIntegrator(K));
   LFA5.Assemble();
   
   std::ofstream out("out/LFA5.txt");
   LFA5.Print(out, 10);
      
   A5 = new SparseMatrix(1, nbrdof);
   for(int k=0; k<nbrdof; k++)
   {
      A5->Set(0, k, LFA5[k]);
   }

   std::ofstream out("out/A5.txt");
   A5->Print(out, 10);
}


int ProximityEffect::CreateOperatorA6()
{

   ConstantCoefficient One(1.0);
   real_t WireArea = IntegrateScalar(*fespace, One, wire_1);
   A6 = new SparseMatrix(1, 1);
   A6->Set(0, 0, sigma_ * WireArea);

   std::ofstream out("out/A6.txt");
   A6->Print(out, 10);

}


int ProximityEffect::CreateOperatorA7()
{
   A7 = new SparseMatrix(*A5);
   *A7= -1.0;
      
   std::ofstream out("out/A4.txt");
   A7->Print(out, 10);
}


int ProximityEffect::CreateRhsVector()
{

}
int ProximityEffect::CreateXVector()
{

   
}
   


static double mu_ = 1.257E-6;
static double epsilon_ = 8.854E-12;
static double sigma_ = 58E6;
static double omega_ = 2.0*M_PI*60;

int main(int argc, char *argv[])
{

   StopWatch chrono;
   tic();
   
   ProximityEffect PE;
   
   PE.Parser(argc, argv);

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

   PE.CreateRhsVector();
   PE.CreateXVector();
   







   ConstantCoefficient One(1.0);
   ConstantCoefficient MinusOne(-1.0);
   ConstantCoefficient Zero(0.0);
   
   PWConstCoefficient *K1_ = new PWConstCoefficient;
   {
      double CoeffArray[]={omega_*mu_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K1_->UpdateConstants(CoeffVector);
   }
   PWConstCoefficient *K1Neg_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*mu_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K1Neg_->UpdateConstants(CoeffVector);
   }
   
   PWConstCoefficient *K2_ = new PWConstCoefficient;
   {
      double CoeffArray[]={sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K2_->UpdateConstants(CoeffVector);
   }

   PWConstCoefficient *K3_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*sigma_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K3_->UpdateConstants(CoeffVector);
   }
 
   

   PWConstCoefficient *K4_ = new PWConstCoefficient;
   {
      real_t val1 = IntegrateScalar(*fespace, One, wire_1);
      real_t val = sigma_ * IntegrateScalar(*fespace, One, wire_1);
      double CoeffArray[]={val, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K4_->UpdateConstants(CoeffVector);
   }

   PWConstCoefficient *K5_ = new PWConstCoefficient;
   {
      double CoeffArray[]={sigma_*mu_, 0.0, 0.0, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 5);
      K5_->UpdateConstants(CoeffVector);
   }


   // 9. Assemble the finite element matrices for the operator
  



// the x vector.
   Vector x(2*nbrdof+2);
   x = 0.0;
   {
      std::ofstream out("out/x.txt");
      x.Print(out, 10);
   }

      cout << BLFA1.Height() << " BLFA1.Height()\n " 
            << BLFA1.Width()  << " BLFA1.Width()\n "

            << BLFA2.Height() << " BLFA2.Height()\n " 
            << BLFA2.Width()  << " BLFA2.Width()\n "
            
            << SMA3.Height() << " SMA3.Height()\n " 
            << SMA3.Width()  << " SMA3.Width()\n "
            
            << BLFA4.Height() << " BLFA4.Height()\n " 
            << BLFA4.Width()  << " BLFA4.Width()\n "

            << SMA5.Height() << " SMA5.Height()\n " 
            << SMA5.Width()  << " SMA5.Width()\n "

            << SMA7.Height() << " SMA7.Height()\n " 
            << SMA7.Width()  << " SMA7.Width()\n ";

// 6. Define the BlockStructure of the problem.
   Array<int> blockRowOffset(5);
   blockRowOffset[0]=0;
   blockRowOffset[1]=BLFA1.NumRows(); 
   blockRowOffset[2]=BLFA4.NumRows();
   blockRowOffset[3]=SMA5.NumRows();
   blockRowOffset[4]=SMA5.NumRows();
   blockRowOffset.PartialSum();
   {
      std::ofstream out("out/blockRowOffset.txt");
      blockRowOffset.Print(out, 10);
   }

   assert(BLFA1.NumRows()==BLFA1.Height());

   Array<int> blockColOffset(5);
   blockColOffset[0]=0;
   blockColOffset[1]=BLFA1.NumCols();
   blockColOffset[2]=BLFA2.NumCols();
   blockColOffset[3]=SMA3.NumCols();
   blockColOffset[4]=SMA3.NumCols();
   blockColOffset.PartialSum();
   {
      std::ofstream out("out/blockColOffset.txt");
      blockColOffset.Print(out, 10);
   }
   assert(BLFA1.NumCols()==BLFA1.Width());
   
   // 8. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
   MemoryType mt = device.GetMemoryType();
//   BlockVector x(blockRowOffset, mt), rhs(blockRowOffset, mt);

   BlockOperator ProxOp(blockRowOffset, blockColOffset);

      
      SparseMatrix SpBLFA2(BLFA2.SpMat());
      SparseMatrix SpBLFA4(BLFA4.SpMat());

      {
      std::ofstream out("out/SpBLFA1.txt");
      SpBLFA1.Print(out, 10);
      }

// Build the operator row by row...
      ProxOp.SetBlock(0, 0, &SpBLFA1);
      ProxOp.SetBlock(0, 1, &SpBLFA2);
      ProxOp.SetBlock(0, 2, &SMA3);

      ProxOp.SetBlock(1, 0, &SpBLFA4);
      ProxOp.SetBlock(1, 1, &SpBLFA1);
      ProxOp.SetBlock(1, 3, &SMA3);
      
      ProxOp.SetBlock(2, 1, &SMA5);
      ProxOp.SetBlock(2, 2, &DMA6);

      ProxOp.SetBlock(3, 0, &SMA7);
      ProxOp.SetBlock(3, 3, &DMA6);

      {
      std::ofstream out("out/ProxOp.txt");
      ProxOp.PrintMatlab(out);
      }

      assert(ProxOp.Height() == 2*nbrdof+2);
      assert(ProxOp.Width() == 2*nbrdof+2);
      
//DL241125: I check ProxOp it contains all the BLFA1 to 4 in proper order.

      assert(2*BLFA1.NumRows()+2==ProxOp.NumRows());
      assert(2*BLFA1.NumCols()+2==ProxOp.NumCols());


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
SparseMatrix A22_mat = BLFA1.SpMat();
SparseMatrix A33_mat(1, 1); A33_mat.Set(0, 0, sigma_);
SparseMatrix A44_mat(1, 1); A44_mat.Set(0, 0, sigma_); 


// Create smoothers for diagonal blocks
GSSmoother gs1(A11_mat); // Gauss-Seidel smoother for A11
GSSmoother gs2(A22_mat); // Gauss-Seidel smoother for A22
GSSmoother gs3(A33_mat);
GSSmoother gs4(A44_mat);

DSmoother ds1(A11_mat); // Diagonal smoother for A11
DSmoother ds2(A22_mat); // Diagonal smoother for A22

BlockDiagonalPreconditioner block_prec(blockRowOffset);
block_prec.SetDiagonalBlock(0, &gs1); // Set smoother for A11
block_prec.SetDiagonalBlock(1, &gs2); // Set smoother for A22
block_prec.SetDiagonalBlock(2, &gs3);
block_prec.SetDiagonalBlock(3, &gs4);

   chrono.Clear();
   chrono.Start();

   if(1)
   {
     // Solve system Ax = b
      GMRESSolver solver;
      solver.SetOperator(*A_ptr);
      solver.SetPreconditioner(block_prec);
      solver.SetRelTol(1e-16);
   //   solver.SetAbsTol(1e-8);
      solver.SetMaxIter(10000);
      solver.SetPrintLevel(1);

      //x = 0.0;       // Initial guess
      solver.Mult(B, X);

   }
   chrono.Stop();

A.RecoverFEMSolution(X, rhs, x);
   {
      std::ofstream out("out/xsol.txt");
      x.Print(out, 10);
   }
   // 12. Create the grid functions.
   GridFunction GFAzr(fespace), GFAzi(fespace);
   GFAzr=0.0;
   GFAzi=0.0;
//   GFAzr.ProjectBdrCoefficient(Zero, ess_bdr_R);
//   GFAzi.ProjectBdrCoefficient(Zero, ess_bdr_R);

// rebuild GFR and GFI from x.
GFAzr.MakeRef(fespace, x, 0);
GFAzi.MakeRef(fespace, x, nbrdof);

Glvis(mesh, &GFAzr, "A-field: Real Part" );
Glvis(mesh, &GFAzi, "A-field: Imag Part" );


// see https://github.com/mfem/mfem/issues/4543

   
   GridFunctionCoefficient ARealCoeff(&GFAzr), AImagCoeff(&GFAzi);

   ProductCoefficient JRealACoeff(AImagCoeff, *K3Neg_);
   ProductCoefficient JRealVCoeff(x[2*nbrdof+0], *K2_);
   SumCoefficient JRealCoeff(JRealACoeff, JRealVCoeff);
   
   ProductCoefficient JImagACoeff(ARealCoeff, *K3_);
   ProductCoefficient JImagVCoeff(x[2*nbrdof+1], *K2_);
   SumCoefficient JImagCoeff(JImagACoeff, JImagVCoeff);
  
   PowerCoefficient JRealSquareCoeff(JRealCoeff, 2.0), JImagSquareCoeff(JImagCoeff, 2.0);
   SumCoefficient JSquareCoeff(JRealSquareCoeff, JImagSquareCoeff);
   PowerCoefficient JCoeff(JSquareCoeff, 0.5);

   FiniteElementCollection *JIndFec = new DG_FECollection(order, dim);
   FiniteElementSpace *JIndFESpace = new FiniteElementSpace(mesh, JIndFec);

   GridFunction JImagGridFunction(JIndFESpace),
                JRealGridFunction(JIndFESpace),
                JGridFunction(JIndFESpace);

   JImagGridFunction.ProjectCoefficient(JRealCoeff);
   JRealGridFunction.ProjectCoefficient(JImagCoeff);
   JGridFunction.ProjectCoefficient(JCoeff);
   
  
   Glvis(mesh, &JRealGridFunction, "JRealGridFunction" );
   Glvis(mesh, &JImagGridFunction, "JImagGridFunction" );
   Glvis(mesh, &JGridFunction, "JGridFunction" );   

/*
cout << "\n ****************************** \n";
cout << "\n IntegrateScalar JImagCoeff = " << IntegrateScalar(*fespace, JImagCoeff, wire_1) << endl;
cout << "\n IntegrateScalar JRealCoeff = " << IntegrateScalar(*fespace, JRealCoeff, wire_1) << endl;
cout << "\n IntegrateScalar JCoeff = " << IntegrateScalar(*fespace, JCoeff, wire_1) << endl;

*/

   // 15. Free the used memory.
delete fec;
delete fespace;
delete mesh;
//delete c_;
//delete negc_;
//delete d_;

   cout << "\n time elapsed = " << toc() << endl;

   return 0;
}

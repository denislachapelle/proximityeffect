// rien de bon :(

if(0){
   int ndofs = fespaceR->GetNDofs();
   cout << "ndofs = " << ndofs << endl;
   SparseMatrix SMA5(ndofs, ndofs);
   Array<int> Cols(ndofs); for(int i=0; i<ndofs; i++) Cols[i]=i;
   for(int i=0; i<ndofs; i++)
   {
      ElementCoefficient *ECptr = new ElementCoefficient(i); 
      LinearForm *LFptr = new LinearForm(fespaceR);
      LFptr->AddDomainIntegrator(new DomainLFIntegrator(*ECptr));
      LFptr->Assemble();
      SMA5.AddRow(i, Cols, *LFptr);
      delete ECptr;
      delete LFptr; 
   }
   
   SMA5.Finalize();
   {
      std::ofstream out("out/SPA5_LF.txt");
      SMA5.Print(out, 10);
   }
}



   
   int ndofs = fespaceR->GetNDofs();
   cout << "ndofs = " << ndofs << endl;
   SparseMatrix SMA5(ndofs, ndofs);
   Array<int> Cols(ndofs); for(int i=0; i<ndofs; i++) Cols[i]=i;
   for(int i=0; i<ndofs; i++)
   {
      ElementCoefficient *ECptr = new ElementCoefficient(i); 
      BilinearForm *BLFptr = new BilinearForm(fespaceR);
      BLFptr->AddDomainIntegrator(new MassIntegrator(*ECptr));
      LFptr->Assemble();

      SMA5.AddRow(i, Cols, *LFptr);
      delete ECptr;
      delete LFptr; 
   }
   
   SMA5.Finalize();
   {
      std::ofstream out("out/SPA5_LF.txt");
      SMA5.Print(out, 10);
   }
}


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


/*
  Based on the note "Separating Components of Error Matrices" by M.H.Shaevitz
  systCOV = COV_shape + COV_mixed + COV_norm
*/

class TSeparate
{
public:
  TSeparate(int user_nbins) {
    nbins = user_nbins;
    matrix_pred.ResizeTo(1, nbins);
    matrix_data.ResizeTo(1, nbins);
    matrix_syst.ResizeTo(nbins, nbins);
    matrix_shape.ResizeTo(nbins, nbins);
    matrix_mixed.ResizeTo(nbins, nbins);
    matrix_norm.ResizeTo(nbins, nbins);
  }

  int nbins;
  
  TMatrixD matrix_pred;
  TMatrixD matrix_data;
  TMatrixD matrix_syst;
  
  TMatrixD matrix_shape;  
  TMatrixD matrix_mixed;
  TMatrixD matrix_norm;

  double val_chi2;
  int val_ndf;
  
  void Set_matrix_shape_mixed_norm();  
  int Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst_user, int index);
};

int TSeparate::Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst_user, int index)
{
  if( num_Y+num_X!=matrix_syst_user.GetNrows() ) {cout<<" ERROR"<<endl; exit(10);}
  
  ///////////////////////////////////////////////////////////////////////////////////////////// for no-systematics
  
  TMatrixD matrix_cov_total = matrix_syst_user;
  for(int idx=1; idx<=num_Y+num_X; idx++)
    if( matrix_cov_total(idx-1, idx-1)==0 ) matrix_cov_total(idx-1, idx-1) = 1e-6;// case inverse
  
  ///////////////////////////////////////////////////////////////////////////////////////////// noConstraint
  
  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_Y = matrix_pred.GetSub(0, num_Y-1, 0, 0);
  TMatrixD matrix_data_Y = matrix_data.GetSub(0, num_Y-1, 0, 0);
  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_YY = matrix_cov_total.GetSub(0, num_Y-1, 0, num_Y-1);
 
  ///////////////////////////// Goodness of fit, Pearson's format

  TMatrixD matrix_Goodness_cov_total_noConstraint = matrix_YY;
  for( int i=0; i<num_Y; i++ ) matrix_Goodness_cov_total_noConstraint(i,i) += matrix_pred_Y(i, 0);
  
  matrix_pred_Y.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_noConstraint = matrix_pred_Y - matrix_data_Y;
  matrix_pred_Y.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_noConstraint_T = matrix_pred_Y - matrix_data_Y;

  TMatrixD matrix_cov_noConstraint_inv = matrix_Goodness_cov_total_noConstraint;
  matrix_cov_noConstraint_inv.Invert();
  TMatrixD matrix_chi2_noConstraint = matrix_delta_noConstraint * matrix_cov_noConstraint_inv * matrix_delta_noConstraint_T;
  double val_chi2_noConstraint = matrix_chi2_noConstraint(0,0);
  double p_value_noConstraint = TMath::Prob( val_chi2_noConstraint, num_Y );  
  cout<<endl<<TString::Format(" ---> GOF noConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.6f",
			      val_chi2_noConstraint, num_Y, val_chi2_noConstraint/num_Y, p_value_noConstraint)<<endl;  
  if( num_X==0 ) {
    cout<<endl;
    return 1;
  }
 
  ///////////////////////////////////////////////////////////////////////////////////////////// wiConstraint
    
  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_X = matrix_pred.GetSub(num_Y, num_Y+num_X-1, 0, 0);
  TMatrixD matrix_data_X = matrix_data.GetSub(num_Y, num_Y+num_X-1, 0, 0);
  matrix_pred.T(); matrix_data.T();

  TMatrixD matrix_XX = matrix_cov_total.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1);
  for(int ibin=1; ibin<=num_X; ibin++) matrix_XX(ibin-1, ibin-1) += matrix_pred_X(ibin-1, 0);// Pearson's term for statistics
  TMatrixD matrix_XX_inv = matrix_XX;
  matrix_XX_inv.Invert();

  TMatrixD matrix_YX = matrix_cov_total.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1);
  TMatrixD matrix_XY(num_X, num_Y); matrix_XY.Transpose(matrix_YX);

  TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_data_X - matrix_pred_X);
  TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
  
  TMatrixD matrix_Goodness_cov_total_wiConstraint = matrix_YY_under_XX;
  for( int i=0; i<num_Y; i++ ) matrix_Goodness_cov_total_wiConstraint(i,i) += matrix_Y_under_X(i, 0);

  matrix_Y_under_X.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_wiConstraint = matrix_Y_under_X - matrix_data_Y;
  matrix_Y_under_X.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_wiConstraint_T = matrix_Y_under_X - matrix_data_Y;  

  TMatrixD matrix_cov_wiConstraint_inv = matrix_Goodness_cov_total_wiConstraint;
  matrix_cov_wiConstraint_inv.Invert();
  TMatrixD matrix_chi2_wiConstraint = matrix_delta_wiConstraint * matrix_cov_wiConstraint_inv * matrix_delta_wiConstraint_T;
  double val_chi2_wiConstraint = matrix_chi2_wiConstraint(0,0);
  double p_value_wiConstraint = TMath::Prob( val_chi2_wiConstraint, num_Y );  
  cout<<TString::Format(" ---> GOF wiConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.6f",
			val_chi2_wiConstraint, num_Y, val_chi2_wiConstraint/num_Y, p_value_wiConstraint)<<endl<<endl;
  
  return 1;
}

void TSeparate::Set_matrix_shape_mixed_norm()
{
  ///
  double N_T = 0;
  for(int idx=0; idx<nbins; idx++) N_T += matrix_pred(0, idx);

  ///
  double M_kl = 0;
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      M_kl += matrix_syst(i,j);
    }
  }

  ///
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {      
      double N_i = matrix_pred(0, i);
      double N_j = matrix_pred(0, j);

      double M_ij = matrix_syst(i,j);      
      double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
      double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);

      matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
      matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;	
      matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
    }
  }
  
}

/////////////////////////////////////////////////////////////
//////////////////////// MAIN ///////////////////////////////
/////////////////////////////////////////////////////////////

/*
  To run the script:
         root -l read_obj.cc
*/

void read_obj()
{
  TFile *file_input = new TFile("file_spectra_systCOV.root", "read");
  TMatrixD *matrix_pred = (TMatrixD*)file_input->Get("matrix_pred");
  TMatrixD *matrix_data = (TMatrixD*)file_input->Get("matrix_data");
  TMatrixD *matrix_systCOV = (TMatrixD*)file_input->Get("matrix_systCOV_FluxXs");
  int nbins = matrix_pred->GetNcols();
  
  ///////////////////

  TSeparate *sep_test = new TSeparate( nbins );  
  sep_test->matrix_pred = (*matrix_pred);
  sep_test->matrix_data = (*matrix_data);
  sep_test->matrix_syst = (*matrix_systCOV);
  sep_test->Set_matrix_shape_mixed_norm();

  TMatrixD matrix_cov_for_GOF = sep_test->matrix_syst;
  //TMatrixD matrix_cov_for_GOF = sep_test->matrix_shape + sep_test->matrix_mixed + sep_test->matrix_norm;
  //TMatrixD matrix_cov_for_GOF = sep_test->matrix_shape;
  //TMatrixD matrix_cov_for_GOF = sep_test->matrix_mixed;
  //TMatrixD matrix_cov_for_GOF = sep_test->matrix_norm;
  
  if( 1 ) {// last 11 bins constrained by first 15 bins
    TMatrixD matrix_gof_trans( sep_test->nbins, 11+15 );// oldworld, newworld
    for( int ibin=1; ibin<=11; ibin++) matrix_gof_trans(15+ibin -1, ibin -1) = 1;
    for( int ibin=1; ibin<=15; ibin++) matrix_gof_trans(ibin -1, 11+ibin -1) = 1;
        
    TMatrixD matrix_gof_trans_T = matrix_gof_trans.T();
    matrix_gof_trans.T();
    
    TMatrixD matrix_gof_pred = sep_test->matrix_pred * matrix_gof_trans;
    TMatrixD matrix_gof_data = sep_test->matrix_data * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * matrix_cov_for_GOF * matrix_gof_trans;
    
    sep_test->Exe_Goodness_of_fit( 11, 15, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 1);    
  }
  
  if( 1 ) {// (last 11 bins --> 1 bin) constrained by first 15 bins
    TMatrixD matrix_gof_trans( sep_test->nbins, 1+15 );// oldworld, newworld
    for( int ibin=1; ibin<=11; ibin++) matrix_gof_trans(15+ibin -1, 0) = 1;
    for( int ibin=1; ibin<=15; ibin++) matrix_gof_trans(ibin -1, 1+ibin -1) = 1;
        
    TMatrixD matrix_gof_trans_T = matrix_gof_trans.T();
    matrix_gof_trans.T();
    
    TMatrixD matrix_gof_pred = sep_test->matrix_pred * matrix_gof_trans;
    TMatrixD matrix_gof_data = sep_test->matrix_data * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * matrix_cov_for_GOF * matrix_gof_trans;
    
    sep_test->Exe_Goodness_of_fit( 1, 15, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 1);    
  }
  
  
}

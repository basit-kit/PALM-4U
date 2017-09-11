
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  Matlab Gateway for the Function Hessian
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 USE kchem_kpp_Model

      INTEGER nlhs, nrhs
      INTEGER plhs(*), prhs(*)
      INTEGER mxGetPr, mxCreateFull, mxGetM, mxgetN
      INTEGER VPtr, FPtr, RPtr, HESSPtr
      REAL(kind=dp) V(13), F(3), RCT(12)
      REAL(kind=dp) HESS(28)

! Check for the right number of input arguments
      IF ( nrhs .ne. 3 ) THEN
         CALL mexErrMsgTxt('Hessian requires 3 input vectors: &
     &V(13), F(3), RCT(12)')
      END IF 
! Check for the right number of output arguments
      IF ( nlhs .ne. 1 ) THEN
         CALL mexErrMsgTxt('Hessian requires 1 output vector: &
     &HESS(28)')
      END IF 

      plhs(1) = mxCreateDoubleMatrix(28,1,0)

      VPtr = mxGetPr(prhs(1));
      CALL mxCopyPtrToReal8(VPtr,V,13)
      
      FPtr = mxGetPr(prhs(2));
      CALL mxCopyPtrToReal8(FPtr,F,3)
      
      RPtr = mxGetPr(prhs(3));
      CALL mxCopyPtrToReal8(RPtr,RCT,12)

      HESSPtr = mxGetPr(plhs(1))

      CALL Hessian( V, F, RCT, HESS )

      CALL mxCopyReal8ToPtr(HESS, HESSPtr, 28)

 END SUBROUTINE mexFunction

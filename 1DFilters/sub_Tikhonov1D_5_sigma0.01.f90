!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!
!Redistribution and use in source and binary forms, with or
!without modification, are permitted provided that the following
!conditions are met:
!
!   Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
!   Redistributions in binary form must reproduce the above
!   copyright notice, this list of conditions and the following
!   disclaimer in the documentation and/or other materials provided
!   with the distribution.
!   Neither the name of the SEISCOPE project nor the names of
!   its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!Warranty Disclaimer:
!THIS SOFTWARE IS PROVIDED BY THE SEISCOPE PROJECT AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!SEISCOPE PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!POSSIBILITY OF SUCH DAMAGE.

!
! Modified Tikhonov regularization routine for edge-preserving FWI
! replace this file with "sub_Tikhonov.f90"
! and recompile TOY2DAC
!
! For 1D windows, directional Tikhonov weightings 
! lambda_x and lambda_z need to be set separately in fwi_input file
! For 2D window types, there are no particular directional
! penalties: either set lambda_x or lambda_z greater than zero. 
! If both set above zero it would cause double penalization.
!

subroutine sub_Tikhonov_fcost(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  
  real :: fcost_tikho
  integer :: i,i1,i2

  !STORE FCOST DATA
  inv%fcost_data=inv%fcost/inv%scalingfactor
  
  !COMPUTE REGUL IN Z DIRECTION
  fcost_tikho=0.
  do i=1,inv%npar
     do i2=1,pbdir%n2     
        do i1=inv%ibathy(i2),pbdir%n1-1        
           fcost_tikho=fcost_tikho+&
                inv%lambda(i)*(inv%model(i1+1,i2,i)-inv%model(i1,i2,i))**2
        enddo
     enddo
  enddo  
  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho  

  !STORE FCOST REG Z
  inv%fcost_reg=1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho

  !COMPUTE REGUL IN X DIRECTION
  fcost_tikho=0.
  do i=1,inv%npar
     do i2=1,pbdir%n2-1     
        do i1=inv%ibathy(i2),pbdir%n1-1        
           fcost_tikho=fcost_tikho+&
                inv%lambda(i)*(inv%model(i1,i2+1,i)-inv%model(i1,i2,i))**2
        enddo
     enddo
  enddo
  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  
  
  !STORE FCOST REG X
  inv%fcost_reg=inv%fcost_reg+&
       1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  

  !NORMALIZE FCOST REG
  inv%fcost_reg=inv%fcost_reg/inv%scalingfactor
  
end subroutine sub_Tikhonov_fcost

subroutine sub_Tikhonov_fgrad(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  real,allocatable,dimension(:,:,:) :: fgrad_tikho_x,fgrad_tikho_z
  integer :: i,i1,i2, i1s, i1e, ii, n, nn, vnn, vrep, cs
  real :: csum, temp
  ! coefficients sum up to 2
  ! sigma 0.001
  integer, parameter :: nwindow = 5
  real, parameter :: cf(nwindow) = [ 0.000000, 0.000000, 2.000000, 0.000000, 0.000000 ]
  real, dimension(nwindow) :: values, values_sorted
  ! half window size (total wondow size is (2*nwh+1) including the center )
  integer, parameter :: nwh = int( (nwindow-1)/2 )
  
  !Vertical term
  allocate(fgrad_tikho_z(pbdir%n1,pbdir%n2,inv%npar))
  fgrad_tikho_z(:,:,:)=0.
  
  do i=1, inv%npar
	 do i2=1, pbdir%n2
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			do ii = (i1-nwh),(i1+nwh)
				if ( ii>=i1s .AND. ii<=i1e .AND. &
				     ii>=inv%ibathy(i2) .AND. ii<=i1e) THEN 
				   nn = nn + 1
				   values(nn) = inv%model(ii,i2,i)
				endif
			enddo

			!sorting
			do vrep = 1, nn
				n = MINLOC(values,1)
				values_sorted(vrep) = values(n)
				values(n) = 1.0e32
			enddo
			! if the number of values are even
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = 1, nn
				   if ( vnn < (nn/2) ) THEN
					   values_sorted(vnn) = values_sorted(vnn);
				   else if ( vnn == (nn/2) ) THEN
					   values_sorted(vnn) = &
					   0.5*(values_sorted(vnn)+values_sorted(vnn+1))
				   else
					   values_sorted(vnn)= values_sorted(vnn+1)
				   endif
			   enddo
			   nn = nn-1
			endif
			
			csum = 0
			do vrep= INT( (nwh+1) - (nn-1)/2 ), &
										INT( (nwh+1) + (nn-1)/2 )
				csum = csum + cf(vrep)
			enddo
			
			fgrad_tikho_z(i1,i2,i) = 2 * inv%model(i1, i2, i)
			do vrep = 1, nn
			   cs = INT( ((nwh+1) - (nn-1)/2) - 1 + vrep )
			   fgrad_tikho_z(i1,i2,i) = fgrad_tikho_z(i1,i2,i) - &
				 (2/csum) * ( values_sorted(vrep) * cf(cs) )
			enddo
			
		enddo
	 enddo
  enddo
  
  
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_z*fgrad_tikho_z(:,:,i)       
  enddo
  deallocate(fgrad_tikho_z)
  
  !Horizontal term
  allocate(fgrad_tikho_x(pbdir%n1,pbdir%n2,inv%npar))
  fgrad_tikho_x(:,:,:)=0.
  do i=1, inv%npar
	 do i2=1, pbdir%n2
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			do ii = (i2-nwh),(i2+nwh)
				if ( ii>=1 .AND. ii<=pbdir%n2 .AND. &
				     ii>=1 ) THEN 
					if 	(i1>=inv%ibathy(ii) .AND. i1<=i1e) THEN
				    	nn = nn + 1
						values(nn) = inv%model(i1,ii,i)
					endif
				endif
			enddo

			!sorting
			do vrep = 1, nn
				n = MINLOC(values,1)
				values_sorted(vrep) = values(n)
				values(n) = 1.0e32
			enddo
			! if the number of values are even
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = 1, nn
				   if ( vnn < (nn/2) ) THEN
					   values_sorted(vnn) = values_sorted(vnn);
				   else if ( vnn == (nn/2) ) THEN
					   values_sorted(vnn) = &
					   0.5*(values_sorted(vnn)+values_sorted(vnn+1))
				   else
					   values_sorted(vnn)= values_sorted(vnn+1)
				   endif
			   enddo
			   nn = nn-1;
			endif
			
			csum = 0
			do vrep= INT( (nwh+1) - (nn-1)/2 ), &
										INT( (nwh+1) + (nn-1)/2 )
				csum = csum + cf(vrep)
			enddo
			
			fgrad_tikho_x(i1,i2,i) = 2 * inv%model(i1, i2, i)
			do vrep = 1, nn
			   cs = INT( ((nwh+1) - (nn-1)/2) - 1 + vrep )
			   fgrad_tikho_x(i1,i2,i) = fgrad_tikho_x(i1,i2,i) - &
				 (2/csum) * ( values_sorted(vrep) * cf(cs) )
			enddo
			
		enddo
	 enddo
  enddo
  
  
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_x*fgrad_tikho_x(:,:,i)
  enddo
  deallocate(fgrad_tikho_x)
  
end subroutine sub_Tikhonov_fgrad

subroutine sub_Tikhonov_Hv(pbdir,inv,Hv,v)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v,Hv
  real,allocatable,dimension(:) :: vect_x,vect_z
  integer :: i,i1,i2, i1s, i2s, i2e, i1e, ii, n, nm, nn, vnn, vrep, current, neighbor, cs
  real :: csum, temp
  ! coefficients sum up to 2
  ! sigma 0.001
  integer, parameter :: nwindow = 5
  real, parameter :: cf(nwindow) = [ 0.000000, 0.000000, 2.000000, 0.000000, 0.000000 ]
  real, dimension(nwindow) :: values, values_sorted
  ! half window size (total wondow size is (2*nwh+1) including the center )
  integer, parameter :: nwh = int( (nwindow-1)/2 )
    
   !Vertical term
  allocate(vect_z(pbdir%n1*pbdir%n2*inv%npar))
  vect_z(:)=0.

  ! total number of parameters    
  nm = pbdir%n1*pbdir%n2*inv%npar
  do i=1, inv%npar 
     i2s=1
	 i2e=pbdir%n2
     do i2 = i2s,i2e  
         i1s = inv%ibathy(i2)
         i1e = pbdir%n1
         do i1 = i1s , i1e
             nn = 0
             values(:) = 1.0e32
             values_sorted(:) = 0.0
             
             current = INT( i1+(i2-1)*pbdir%n1 + &
                                               (i-1)*pbdir%n1*pbdir%n2)
             do ii = (-1*nwh),nwh
                 neighbor = INT( (i1+ii) + (i2-1)*pbdir%n1 + &
                                              (i-1)*pbdir%n1*pbdir%n2 )
                 if ( (i1+ii)>=i1s .AND. (i1+ii)<=i1e .AND. &			  
                      (neighbor)>=1 .AND. (neighbor)<=nm .AND. &
					  (i1+ii)>=inv%ibathy(i2) .AND. (i1+ii)<=i1e) THEN
                    nn = nn + 1
                    values(nn) = v(neighbor)
                 endif
             enddo
             
             !sorting
             do vrep = 1,nn
                 n = MINLOC(values,1)
                 values_sorted(vrep) = values(n)
                 values(n) = 1.0e32
             enddo
             
             temp = MOD(nn,2)
             if ( temp == 0.0 ) THEN
                do vnn=1,nn
                    if ( vnn<(nn/2) ) THEN
                        values_sorted(vnn) = values_sorted(vnn);
                    else if ( vnn == (nn/2) ) THEN
                        values_sorted(vnn)= 0.5*(values_sorted(vnn) &
                                            + values_sorted(vnn+1) )
                    else
                        values_sorted(vnn)= values_sorted(vnn+1)
                    endif
                enddo
                nn = nn-1
             endif
               
             csum = 0;
             do vrep = ( (nwh+1) - (nn-1)/2 ), ( (nwh+1) + (nn-1)/2 )
                 csum = csum + cf(vrep)
             enddo
             vect_z(current)= 2*v(current)
             do vrep=1, nn
                cs = INT( ((nwh+1) - (nn-1)/2) - 1 + vrep )
                vect_z(current) = vect_z(current) - &
                           ( 2/csum ) * ( values_sorted(vrep)*cf(cs) ) 
             enddo
             vect_z(current) = inv%lambda(i)*( vect_z(current) )
         enddo
     enddo
      
     
  enddo
    
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_z*vect_z(:)
  deallocate(vect_z)
  
  !Horizontal term  
  ! total number of parameters    
  nm = pbdir%n1*pbdir%n2*inv%npar
  allocate(vect_x(pbdir%n1*pbdir%n2*inv%npar))
  vect_x(:)=0.
  do i=1,inv%npar  
     i2s=1
	 i2e=pbdir%n2
	 do i2 = i2s,i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1  
		do i1 = i1s , i1e
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			
			current = INT( i1+(i2-1)*pbdir%n1 + (i-1)*pbdir%n1*pbdir%n2)
			do ii = (-1*nwh),nwh
				neighbor = INT( (i1) + (i2-1 + ii)*pbdir%n1 + &
				                               (i-1)*pbdir%n1*pbdir%n2 )
                if ( (i2+ii)>=i2s .AND. (i2+ii)<=i2e .AND. &
				     (neighbor)>=1 .AND. (neighbor)<=nm .AND. &
					 (i2+ii)>=1 ) THEN
					if (i1>=inv%ibathy(i2+ii) .AND. i1<=i1e)  THEN
			        	nn = nn + 1
			        	values(nn) = v(neighbor)
					endif
				endif
			enddo
			
		   !sorting
			do vrep = 1,nn
				n = MINLOC(values,1)
				values_sorted(vrep) = values(n)
				values(n) = 1.0e32
			enddo

			temp = MOD(nn,2)
			if ( temp == 0.0 ) THEN
			   do vnn=1,nn
				   if ( vnn<(nn/2) ) THEN
					   values_sorted(vnn) = values_sorted(vnn);
				   else if ( vnn == (nn/2) ) THEN
					   values_sorted(vnn)= 0.5*(values_sorted(vnn) &
										   + values_sorted(vnn+1) )
				   else
					   values_sorted(vnn)= values_sorted(vnn+1)
				   endif
			   enddo
			   nn = nn-1
			endif
			
			csum = 0;
			do vrep = ( (nwh+1) - (nn-1)/2 ), ( (nwh+1) + (nn-1)/2 )
				csum = csum + cf(vrep)
			enddo
			vect_x(current)= 2*v(current)
			do vrep=1, nn
			   cs = INT( ((nwh+1) - (nn-1)/2) - 1 + vrep )
			   vect_x(current) =  vect_x(current) - &
						( 2/csum ) * ( values_sorted(vrep)*cf(cs) ) 
			enddo
			vect_x(current) =  inv%lambda(i)*( vect_x(current) )  
		enddo
	 enddo

  enddo
  
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_x*vect_x(:)  
  deallocate(vect_x)
     
end subroutine sub_Tikhonov_Hv


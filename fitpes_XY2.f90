module xy2_globfit

 implicit none

 private 
 public morbid,morse,fitting,wifin_242,morse_242,morse_tyuterev,morse_tyuterev_II,taylor_I,morbid_MEP,R_r12_dr_morse_pol,morbid_mep2
 public R_r12_morbid_switch,morbid_mep3,R_r12_dr_morse_pol_x,V_rep_disp,morbid_MEP4,potential_xy2_tyuterev_damp,dipol_xy2_p,dipol_xy2_q,dipol_xy2_p_dump
!    

! Here we go:

contains 

subroutine fitting(potential)


!
! Parameters:  
 ! 
 ! number of parameters and maximum number of ab initio points   
 !
 integer,parameter          ::  parmax=500,enermax =100000   
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10
 ! 
 ! Dinite difference differentiation with 2 or 3 points
 !
 integer,parameter          ::  findif = 3
 !
 ! parameter to control the correlation when we through out the parameters with large st. error (0.0..1.0)
 !
 double precision   :: fact_corr=0.99
 !
 ! the iteration will be stoped when the best standard deviation reached 
 !
 double precision,parameter :: stadev_target = 0.00001d0   
 !
 integer,parameter          ::  verbose = 5
 ! 
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
 !
 ! Universal constants 
 ! 
 double precision,parameter :: bohr = 0.52917715e+00,hartre=2.194746354d+05,pi=3.141592653589793d0
 ! Avareging region for the robustfit = 3 
 integer,parameter          :: ishift = 10

 ! Variables:  
 double precision  :: maxdeviat  ! st. deviation at the previos iteration 
 !
 character(len=40)  :: fit_type='dgelss'
 !
 integer           :: i,npts,nused,itmax,numpar,ncol,iter,nrow,l,itertypemax,irow,icolumn ! just integers  
 ! 
 character(len=80) :: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=11) :: parnam_(parmax) ! Names of the parameters, as they appear in the input file  
 character(len=1)  :: mark           ! used at the parameters output to mark the dependent or independent ones 
 character(len=11) :: label          ! Temporary label 
 double precision  :: wtsum,ssq      ! 
 !
 ! some variables 
 !
 double precision  :: stadev,stadev_best,rms_best,epsil,stadev_old,v,potential,dd
 double precision  :: tempx,deltax,potright,potleft,stability
 double precision  :: r12,r32,alpha1
 double precision  :: corr,ve
 integer           :: ierror,i1,i2,itertype,last_i,robustfit,nparams
 double precision  :: stab_best,acc_rms,conf_int
 integer           :: ivar_(parmax),maxivar,minivar,tivar
 double precision  :: last_stadev,factor_back,factor_out
 double precision  :: rms,s,da1,da2,a_wats,a_wats_best,da,ZPE,local(3)
 integer           :: last_ivar,ndigits,npoints_delta,rand_itmin,alloc
 integer           :: max_pnts,ipar
 logical           :: yes,still_run,ifopen
 integer           :: flag_borh, flag_hartre,ivartmp,it0,imin,imax,ilarge_sterr
 double precision  :: conv_lenght,conv_energy,paramtmp,relat_sterr
 double precision  :: param_(parmax),relat_sterr_min
 character(len=68)  :: fmt1,fmt2
 character(len=2)   :: fmt0
 integer            :: Nout,ilow_sterr,l0,i0
 
 !
 ! some matrices, to be allocated  
 !
 character(len=11),allocatable :: parnam(:)
 integer, allocatable :: ivar(:),icorr(:),ivar0(:)
 double precision, allocatable :: rjacob(:,:),last_param(:)
 double precision, allocatable :: crr(:),blen1(:),blen2(:)
 double precision, allocatable :: bang(:),energy(:)
 double precision, allocatable :: eps(:),wt(:),parold(:),sigma(:),wt_tmp(:)
 double precision, allocatable :: al(:,:),bl(:),dx(:)
 double precision, allocatable :: ai(:,:),sterr(:)
 double precision, allocatable :: sens(:),param(:)
 double precision, allocatable :: Tsing(:,:)
 integer                       :: rank0,info,j,jlistmax,neval_t(0:1)

 integer,parameter             :: lwork = 8*PARMAX
 double precision,allocatable  :: wspace(:)
 !
 ! Common block is to reduce the number of assignments within the potential routine
 ! that repeat every iteration  
 !common /potparam/ param
 !
 ! Here we go!
 !
 ! Array's allocation:
 !

! allocate (param(parmax),stat=alloc)
! if (alloc/=0) then 
!  write(6,"('parmax - out of memory')") 
! stop
! endif 

 allocate (blen1(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('blen1 - out of memory')")  
  stop
 endif 

 allocate (blen2(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('blen2 - out of memory')")  
  stop
 endif 


 allocate (bang(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('bang - out of memory')")  
  stop
 endif 

 allocate (energy(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('energy - out of memory')")  
  stop
 endif 

 allocate (eps(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('eps - out of memory')")  
  stop
 endif 

 allocate (wt(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('wt - out of memory')")  
  stop
 endif 

 allocate (wt_tmp(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('wt_tmp - out of memory')")  
  stop
 endif 


 allocate (sigma(enermax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sigma - out of memory')")  
  stop
 endif 

 allocate (Wspace(lwork),Tsing(parmax,parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('Wspace - out of memory')")  
  stop
 endif


 !
 !  initial constants 
 !
 last_stadev = 1000000.
 maxivar = 8
 maxdeviat = 200.0
 !
 ! Input data from file: 
 !
 !call skiplines(f_inp,1)
 !
 !  input the job title 
 !
 !do i=1,4
 !  read  (f_inp,"(a80)") title(i)
 !enddo
 !
 !  output the job title 
 !
 !write (f_out,"(3(132('*')/' '),16('*'),100x,16('*'))") 
 !write (f_out,"(4(' ',16('*'),10x,a80,10x,16('*')/))") (title(i), i=1,4)
 !write (f_out,"(' ',16('*'),100x,16('*')/' ',3(132('*')/' '))") 
 !
 !call skiplines(f_inp,4)
 !
 read (f_inp,"(a40)") fit_type
 read (f_inp,*) itmax
 read (f_inp,*) stab_best,fact_corr
 read (f_inp,*) conv_lenght
 read (f_inp,*) conv_energy
 read (f_inp,*) factor_out
 read (f_inp,*) factor_back
 read (f_inp,*) robustfit
 ! 
 write (f_out,"(i10,' <= number of iterations in the fit')")         itmax
 write (f_out,"(f10.6,' <= The fit is over if the stability is reached ')") stab_best
 write (f_out,"(f10.6,' <= sterr parameter out ')")                  factor_out
 write (f_out,"(f10.6,' <= stdev parameter out ')")                  factor_back
 write (f_out,"(i10,  ' <= watson robust on (1) or of (0) ')")       robustfit
 !
 ! Converiosn factors:
 ! bond length 
 !
 !if (flag_borh ==1) then 
 !  conv_lenght = bohr 
 !  write (f_out,"(' The bond lengths are in bohr')") 
 !else 
 !  conv_lenght = 1.0
 !  write (f_out,"(' The bond lengths are in Angstrom')") 
 !endif   
 !
 ! energies 
 !if (flag_hartre==1) then 
 !  conv_energy = hartre
 !  write (f_out,"(' The ab initio energies are in hartre')") 
 !else 
 !  conv_energy = 1.0 
 !  write (f_out,"(' The ab initio energies are in 1/cm')") 
 !endif   
 !
 ! input potential parameters 
 !
 nparams = parmax
 !
 call skiplines(f_inp,3)
 do i = 1,parmax
   read (f_inp,"(a80)") longlabel
   read (longlabel(1:11),"(a11)") label
   read (longlabel(12:80),*) ivartmp,paramtmp
   if (ivartmp/=-1) then 
     parnam_(i) = label 
     ivar_(i)   = ivartmp
     param_(i)  = paramtmp
   else
     !
     nparams = i-1
     exit 
     !
     !write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i,parmax
     !stop 'Too few parameters'
   endif 
 enddo 


 allocate (al(nparams,nparams),bl(nparams),crr(nparams),sterr(nparams),last_param(nparams),parold(nparams),stat=alloc)
 if (alloc/=0) then 
  write(6,"('al - out of memory')")  
  stop
 endif 

 allocate (param(nparams),parnam(nparams),ivar(nparams),rjacob(enermax,nparams),&
           sens(nparams),ai(nparams,nparams),dx(nparams),icorr(nparams),ivar0(nparams),stat=alloc)
 if (alloc/=0) then 
  write(6,"('rjacob - out of memory')")  
  stop
 endif 

 do i = 1,nparams
   parnam(i) = parnam_(i) 
   ivar(i)   = ivar_(i)
   param(i)  = param_(i)
 enddo 


 !read (f_inp,"(a80)") longlabel
 ! write(6,"(a80)") longlabel
 !read (longlabel,"(a10,i4,d18.8)") label,ivartmp,paramtmp
 if (ivartmp/=-1) then 
   write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i,nparams
   stop 'Too many parameters'
 endif 
 !
 ! Output of the potential parameters 
 ! 
 write (f_out,"('0',5x,12('*'),' potential energy expansion coefficients (cm-1) ',12('*')//)")
 do i = 1,int(nparams/3)*3,3 
 write (f_out,"(3('  ',a11,'=',f16.3))") parnam(i),param(i),     &
                                        parnam(i+1),param(i+1), &
                                        parnam(i+2),param(i+2)
 enddo 
 write (f_out,"(3('  ',a11,'=',f16.3))") (parnam(int(nparams/3)*3+i), &
                                         param(int(nparams/3)*3+i),i=1,mod(nparams,3))


!---input energy values------!
 npts = 0
 yes = .true. 
 Nout = 1
 last_i = 0
 do while( yes )
    npts=npts+1
    if (npts>enermax) then
      write(f_out,"('Too many ab initio points, increase enermax:',I6,' vs ',I6)") npts,enermax
      endif
    read (f_inp,*) blen1(npts),blen2(npts),bang(npts),energy(npts),wt(npts)
    if ( blen1(npts).lt.0.0 ) yes = .false.
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
 !
 ! For the Robust fit we need to dedine accuracy of the abinitio data "sigma"
 ! 
 if (RobustFit/=0) then 
   do nrow=1,npts
     if (wt(nrow)>0.0) then 
       sigma(nrow) = 1.0d0/wt(nrow)*0.1
       !wt(nrow) = 1.0
       !sigma(nrow) = 1.0
     else 
       sigma(nrow) = 1000.0 
     endif 
   enddo
 endif 
 ! 
 ! "normalising" the weight factors
 nused=0
 wtsum=0.0d0
 do i=1,npts
   if (wt(i) > 0.0) nused=nused+1
   wtsum=wtsum+wt(i)
 enddo 
 wtsum=wtsum/nused
 wt(:)=wt(:)/wtsum
 !
 ! We introduce a "zero point energy" shift to avoid loosin accuracy 
 ! because of too too large numbers
 ZPE = minval(energy(1:npts))
 ! Convert the energies to the internal dimension, e.g. 1/cm 
 energy(:) = ( energy(:) )*conv_energy
 !
 ve = ZPE*conv_energy
 !

 a_wats = 0.33d0
 stadev_best = huge(1.0d0)
 rms_best = huge(1.0d0)
 !
 ! Sometimes pot. parameters are dependent. They are nor to be in the fit 
 ! If we meet a dependent parameter, we get rid of it and start the fit again 
 ! This is the outer loop for
 !

 still_run = .true.
 outer_loop: do while (still_run)  
   !
   !  make a table header 
   !
   if (itmax.ne.0) then 
     write(f_out,"(/3x,59('-'))") 
     write(f_out,"('   |   iter | points |   deviat    |    rms      | stability |')")
     write(f_out,"(3x,59('-'))") 
   endif 
   !
   ! Parameters to control current and previous st.deviations and stability
   !
   stadev_old = 1.e10
   stability =  1.e10
   stadev    =  1.e10
   ! 
   ! Parameters to control excluding dependent parameters   
   ! We need to know the maximum and minimum values of the pot. parameters weights "ivar" 
   !
   ! numpar is to count the number fitted varying pot. parameters  
   !
   numpar  = 0
   maxivar = 0
   minivar = 1000000
   do i=1,nparams
     if (ivar(i) .gt. 0) numpar=numpar+1
     if (ivar(i) .gt. maxivar) maxivar = ivar(i)
     if (ivar(i) .lt. minivar .and. ivar(i).ne.0) minivar = ivar(i)
   enddo 
   !
   if (nused<numpar) then 
       write(f_out,"('worning! the number of parameter <  the number of data points',2I6)") nused,numpar
       stop 'worning! the number of parameter <  the number of data points'
    endif
   if (nused==numpar) then 
      write(f_out,"('worning! the same number of parameters and data points')")
   endif 
   !
   !
   !##################################################c
   !############ start of the fitting ################c
   !##################################################c
   rjacob = 0 
   iter = 0
   sens = 0.d0 
   do while( iter<=itmax .and. stadev>stadev_target .and. (  stability.ge.stab_best) )   
     iter = iter + 1
     ssq=0.0d+00
     rms=0.0d+00
     eps = 0
     do nrow=1,npts
       !
       !print*,nrow
       !
       ! Running bond lenght 
       !
       local(1)=blen1(nrow)*conv_lenght
       local(2)=blen2(nrow)*conv_lenght
       !
       ! Running interbond angles
       !
       local(3)=bang(nrow)*pi/1.8d+02
       ! 
       ! Value of the potential energy fucntion at the current geometry 
       !
       !call wifin(local(1),local(2),local(3),v)
       !
       v = potential(local,ZPE*conv_energy,nparams,param)
       !
       !
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (itmax/=0 .and. wt(nrow)>0.0) then
         ncol=0
         do  i=1,nparams
           if (ivar(i) .ne. 0) then
             ncol=ncol+1
             tempx=param(i)
             deltax=factordeltax*abs(tempx)
             if (deltax .le. 1e-15) deltax=1e-6
             !
             param(i)=tempx+deltax
             potright=potential(local,ZPE*conv_energy,nparams,param)
             !
             select case ( findif )
             case default
               write (f_out,"(' Bad finite diff. parameter',i)") findif
               stop 'Bad finite diff. parameter'
             case (2)
               dd=(potright-v)/(deltax)
             case (3)
               param(i)=tempx-deltax
               potleft =potential(local,ZPE*conv_energy,nparams,param)
               dd=(potright-potleft)/(2.d0*deltax)
             end select
             param(i)=tempx
             !
             rjacob(nrow,ncol)=dd
           endif
         enddo ! --- ncol
       endif
       !
       eps(nrow) = energy(nrow)-v
       !
       ! Weighted st. square deviation 
       !
       ssq=ssq+eps(nrow)*eps(nrow)*wt(nrow)
       !
       ! mean square 
       !
       rms=rms+eps(nrow)*eps(nrow)

     enddo  ! ---  nrow
     !
     ! We constract a set of linear equations A x = B
     !
     if (itmax.ne.0) then
        ! form A matrix 
        do irow=1,numpar       !==== row-...... ====!
          do icolumn=1,irow    !==== column-....====!
            al(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow)*wt(:))
            al(icolumn,irow)=al(irow,icolumn)
          enddo
        enddo
        ! form B matrix 
        do irow=1,numpar       !==== row-...... ====!
          bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wt(1:npts))
        enddo   
        !
        ! Solve the set of linear equations 
        !
        select case (trim(fit_type)) 
        case default
          write (6,"('fit_type ',a,' unknown')") trim(fit_type)
          stop 'fit_type unknown'
        case('linur') 
          !
          call linur(numpar,nparams,al,bl,dx,ierror)
          ! If ierror /=0, ierror gives the number of the dependent parameter
          ! we remove this parameter from the fit 
          if (ierror.ne.0) then 
            ncol=0
            do i=1,nparams
              if (ivar(i) .ne. 0) then
                ncol=ncol+1
                if  ( ncol.eq.ierror ) then 
                    ivar(i) = 0
                    param(i) = 0.d0
                    write(f_out,"(i,'-th is out - ',a11)") i,parnam(i)
                endif 
              endif 
            enddo 
            cycle outer_loop    
          endif            !
        case ('dgelss')
          !
          stop 'dgelss is not implemented'
          !
          ai = al 
          !call DGELSS(numpar,numpar,1,Ai(1:numpar,1:numpar),numpar,BL(1:numpar),numpar,Tsing,1.D-12,RANK0,wspace,lwork,info)
          !
          if (info/=0) then
            write(6,"('DGELSS:error',i)") info
            stop 'dgelss'
          endif
          !
          dx = bl
          !
        end select 
        !
        ! We can define the stand. deviation sigma if nused / =numpar
        ! otherwise stadev is defined as a weighted RMS 
        !
        if (nused.eq.numpar) then 
          stadev=dsqrt(ssq/float(nused))
        else 
          stadev=dsqrt(ssq/float(nused-numpar))
        endif 
        !
        ! We keep the previous values of pot. parameters 
        !
        parold=param
        !   
        ! Update the pot. parameters to the new values 
        !
        ncol=0
        do i=1,nparams
          if (ivar(i) /= 0) then
             ncol=ncol+1
             param(i)=param(i)+dx(ncol)
          endif
        enddo
        !
        ! Calcualte the inverse to A matrix, needed for the standard error of the fitted pot. parameters
        !
        call invmat(al,ai,numpar,nparams)
        !
        !
        ncol = 0 
        do i=1,nparams
          if (ivar(i) /= 0) then
             ncol=ncol+1
            !  
            ! If nused = numpar, the problem is degenerated, Ax = B is a system of ordinary equations
            ! st. error doesn't make sense  
            !
            if (nused==numpar) then  
              sterr(ncol)=0
            else
              !
              ! Standard definition of the st. error 
              !
              sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
              !
            endif
          endif
        enddo    
        !
        ! Here we define stability to see if the fit is converged 
        !      
        stability=abs( (stadev-stadev_old)/stadev )
        stadev_old=stadev
        !
        if (RobustFit/=0) then 
          !
          !Watson alpha-parameter
          !
          !if (iter==itmax) then 
           !  !
          !  a_wats = a_wats_best
          !  !
          !else
          !
          da = 1e8 
          !
          it0 = 0
          !
          do while (it0<=100.and.abs(da)>=1e-8)
            !
            it0 = it0+1
            da1 = 0
            da2 = 0
            do nrow=1,npts
              if (wt(nrow) .gt. 0.0d0) then 
                da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
                da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
              endif 
            enddo 
            !
            da =( da1 -  float(nused-numpar) )/da2
            a_wats = a_wats + da
            !
            if (a_wats<0.0000001) a_wats = 0.001+it0*0.01
          enddo
          !
          a_wats = 0.9
          !
          !endif 
          ! 
          !  adjusting the weights  
          ! 
          do nrow=1,npts
             if (wt(nrow) .gt. 0.0d0) then 
               wt(nrow) = 1.d0/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
             endif 
          enddo 
          ! 
          ! "re-normalising" the weight factors
          !
          nused=0
          wtsum=0.0d0
          do i=1,npts
            if (wt(i) > 0.0) nused=nused+1
            wtsum=wtsum+wt(i)
          enddo 
          wtsum=wtsum/nused
          wt(:)=wt(:)/wtsum

          ssq = sum(eps(1:npts)**2*wt(1:npts))
          !
          ! We can recalculated the stand. deviation with new weigh factors
          ! otherwise stadev is defined as a weighted RMS 
          !
          if (nused.eq.numpar) then 
            stadev=dsqrt(ssq/float(nused))
          else 
            stadev=dsqrt(ssq/float(nused-numpar))
          endif 
        endif ! --- robust fit
            
     else   ! Itermax = 0, so there is no fit
            ! only straightforward calculations 
        !
        stadev_old=stadev
        stadev=dsqrt(ssq/float(nused))
     endif 

     !
     ! Here we define the root mean square 
     !
     rms = dsqrt(rms/npts)
     !
     if (stadev<stadev_best.and.rms_best<rms.and.iter<itmax) then 
       !
       stadev_best = stadev
       !
       rms_best = rms
       !
       a_wats_best = a_wats
       !
     endif 
     !
     ! Do some output 
     !
     if (itmax/=0 .and. stadev < maxdeviat ) then 
         write (f_out,"('   |  ',i3,'   | ',i5,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
                iter,nused,stadev,rms,stability
     endif 
     !
    !inquire(f_res,opened=ifopen)
    !if ( ifopen ) then
    !  rewind(f_res)
    !else
    !  open(f_res,file='00.res',status='replace' )
    !endif
    !!
    !do i=1,nparams
    !   write (f_res,"(a11,i4,2x,e22.12)") parnam(i),ivar(i),param(i)
    !enddo 
    !do nrow=1,npts
    !!
    !r12=blen1(nrow)*conv_lenght
    !r32=blen2(nrow)*conv_lenght
    !!
    !alpha1=bang(nrow)
    !
    ! v = energy(nrow)-eps(nrow)
    ! write (f_res,"(2f14.8,f8.3,3f12.2,2x,d8.2)")            & 
    !       r12,r32,                                  &
    !       alpha1,                            &
    !       energy(nrow)-ve,v-ve,eps(nrow),                   &
    !       wt(nrow)
    ! !
    !enddo
     !
   enddo  ! --- iter
   !################ end of iterations  #######################
    if (RobustFit/=0) write(f_out,"(/'Watson Alpha parameter',f8.4,' fitted with rms2 = ',f8.4)") &
                                     a_wats,da1/dfloat((nused-numpar))
   !
   ! if RobustFit==2,3 after first iteration it converts to the standard fit RobustFit=0, but with robust weights
   !
   if ( RobustFit==2 ) then
     !
     write(f_out,"('watson-param = ',f18.8)") a_wats
     !
     RobustFit = 0 
   else if ( RobustFit==3 ) then ! Averaging after the first iteration 
     do i=1,npts
       imin = max(   1,i-npts/ishift) 
       imax = min(npts,i+npts/ishift) 
       ssq = sum(wt(imin:imax))/(imax-imin)
       wt_tmp(i) = ssq
     end do
     wtsum=sum(wt_tmp(:))
     wtsum=wtsum/nused
     wt(1:npts)=wt_tmp(1:npts)/wtsum
     RobustFit = 0 
     cycle outer_loop    
   end if 
   !
   ! We need to know the equilibrium value of the pot. function to substract it from the 
   ! output energies 
   !
   !ve = 0 ! param(1) !-ZPE*conv_energy

   ! Output some staqtistics and results 
   !
   !  only if we are fitting:  
   !
   if (itmax.ne.0) then
     !
     ! If any of the pot. parameters (last_i) has been removed from the fit at the previous iteration step
     ! at this moment we check whether the st. deviation is worse than it was before we got rid of the pot. parameter 
     ! in this case we take the last_i-pot. parameter back
     ! this very tricky process is controled by factor_back 
     ! 
     if (last_i/=0.and. (stadev-last_stadev)/stadev > Nout*factor_back ) then 
       ivar = ivar0     
       ivar(last_i)  =  ivar(last_i)+1
       param = last_param
       last_stadev = 10000000.
       write(f_out,"(I6,'-th is back - ',A11)") last_i,parnam(last_i)
       cycle outer_loop    
     endif 
     ! 
     ! It's time to do the output of pot. parameters 
     ! We print them out rounded with their st.errors
     ! 
      write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
      l = 0 
      do i=1,nparams
        if (ivar(i) .ne. 0) then
           l=l+1
           ndigits = 0
           conf_int = sterr(l)
           do while (conf_int.le.10.0)
             ndigits = ndigits +1 
             conf_int = conf_int*10
           enddo
           write(fmt0,"(i2)") ndigits+2
           fmt0 = adjustl(fmt0)
           fmt1 = '(a11,I4,2x,f24.'//fmt0//'    ''('',i8,'')'')'
           write (f_out,FMT=fmt1) parnam(i),ivar(i),param(i),nint(conf_int)
!          write (f_out,7007) parnam(i),ivar(i),param(i),nint(conf_int)
        else 
           ndigits =2
           if (param(i).ne.0.0) ndigits = 10

           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           fmt1 = '(a11,I4,2x,f24.'//fmt0//')'
           write (f_out,FMT=fmt1) parnam(i),ivar(i),param(i)
!          write (f_out,7008) parnam(i),ivar(i),param(i)
        endif
      enddo  ! --- i

!      7007  format(a10,I4,2x,f14.<ndigits>,'    (',i8,')')
!      7008  format(a10,I4,2x,f14.<ndigits>)
      7007  format(a11,I4,2x,f14.4,'    (',i8,')')
      7008  format(a11,I4,2x,f14.4)

      !
      ! Another statistics of the fit: 
      ! correlations between pot.parameters 
      !
      write (f_out,"(/'Correlation coefficients larger than 0.9: ')")
      do i1=1,numpar
        do i2=i1+1,numpar
          corr=ai(i1,i2)/sqrt(abs(ai(i1,i1)*ai(i2,i2)))
           if (abs(corr) .gt. 0.9) write (f_out,"(10x,'corr(',i2,',',i2,') = ',f10.7)") i1,i2,abs(corr)
        enddo
      enddo
      ! 
      ! sensitivity
      !
      ncol = 0 
      do i=1,nparams
        if (ivar(i) /=  0) then
          ncol=ncol+1
          if (nused /= numpar) then  
            ssq = sum(rjacob(:,ncol)**2*wt(:))
            sens(ncol)=stadev*0.1d0/dfloat(numpar)/dsqrt( ssq/dfloat(nused) )
          endif
        endif
      enddo    
      ! 
      ! Robust fit secion 
      !
      !
      ! Pot. parameters output with their st. errors 
      !
      write (f_out,"(/15x,80('-')/15x,':     old parm   :     new parm   :  delta parm  :   std. err.  : sensitility'/15x,80('-'))")
      l = 0
      do i=1,nparams
        if (ivar(i) /= 0) then
          l=l+1
          !
          epsil=param(i)-parold(i)
          crr(l)=(sterr(l)*al(l,l)/stadev)**2
          !
          mark  = ' '
          if (abs(sens(l)) >= abs(epsil) ) mark = '!'
          !
          if (abs(param(i)) > abs(sterr(l))) then
             write (f_out,"(i3,2x,a11,1x,2(':',f18.4,2x),':',3(f18.4,':'),a1)") &
                    l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          else
              write (f_out,"(i3,2('#'),a11,1x,2(':',f18.4,2x),':',3(f18.4,':'),a1,'#')") & 
                      l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          endif
        endif
      enddo
      !
      !  Here we check if there is a pot. parameter with too big st. error  
      !  We consider a parameter to be unneccesary and remove from the fit 
      !  keeping in mind that we might have to take it back in case the st.deviation will get worse
      !
      !
       !
      ! Search for the largest relative sterr
      !
      if(trim(fit_type)/='dgelss') then 
        !
        tivar = minivar
        ivar0 = ivar
        icorr = 0
        last_param = param
        Nout = 0
        do while ( tivar.lt.maxivar )
          if (verbose>=5) write (f_out,"(' tivar, factor_out =',i5,f18.8)") tivar,factor_out
          relat_sterr_min = 1.0e20
          ilow_sterr = -1

          relat_sterr = 0.0
          ilarge_sterr = -1

          l = numpar+1
          do i=nparams,1,-1
            if (ivar(i) /= 0) then
              l=l-1
              !
              ! Get read of parameters with large rel.st.error and their dependencies 
              !
              if (verbose>=5) write (f_out,"(' param, number, ivar, st.error ',f22.12,i5,i5,f18.7)") param(i),i,ivar(i),abs(sterr(l)/param(i))
              if (param(i)/=0.0d0.and.ivar(i) == tivar) then
                if (icorr(i)==0.and.factor_out<abs(sterr(l)/param(i))) then
                  if (verbose>=5) write (f_out,"(' param, number, st.err., ivar',f18.7,i5,f18.7,i5)") param(i),i,abs(sterr(l)/param(i)),ivar(i)
                  ! relat_sterr=abs(sterr(l)/param(i))
                  ! ilarge_sterr = i
                  icorr(i) = 1
                  l0=0
                  relat_sterr = abs(sterr(l)/param(i))
                  ilarge_sterr = i
                  do i0=1,nparams,1
                    if (ivar(i0) /= 0.and.l/=l0) then
                      l0=l0+1
                      corr=abs(ai(l,l0))/sqrt(abs(ai(l,l)*ai(l0,l0)))
                      !   if (verbose>=4) then 
                      !     write (f_out,"('corr(',i,',',i,') = ',f14.7)") l,l0,corr
                      !   endif 
                      if (i0/=i.and.param(i0)/=0.0d0.and.ivar(i0) == tivar) then
                        if (corr>fact_corr.and.icorr(i0)==0) then
                          if (verbose>=4) write (f_out,"('corr, param, (st. err.)  number',2f18.7,f18.7,i4)") corr,param(i0),sterr(l0),i0
                          icorr(i0) = 1
                          !
                          ! Here one finds the (lowest) largest st. error and which par. it belongs to
                          if ( relat_sterr<abs(sterr(l0)/param(i0))) then
                              relat_sterr=abs(sterr(l0)/param(i0))
                              ilarge_sterr = i0
                          endif
                        endif
                      endif
                    endif
                  enddo
                  ivar(ilarge_sterr)  = 0
                  param(ilarge_sterr) = 0.0d0
                  if (verbose>2) then 
                    write(6,"('relat_sterr,ilarge_sterr, tivar:',f18.8,2i5)") relat_sterr,ilarge_sterr, tivar
                  endif 
                  Nout = Nout+1
                  !
                  ! Here one finds the (lowest) largest st. error and which par. it belongs to
                  if ( relat_sterr_min>relat_sterr) then
                      relat_sterr_min=relat_sterr
                     ilow_sterr = ilarge_sterr
                  endif
                endif
              endif
            endif
          enddo
           if (verbose>2) then 
             write(6,"('relat_sterr_min,ilow_sterr, tivar:',d18.8,2i5)") relat_sterr_min,ilow_sterr, tivar
          endif 
          if (verbose>=4) then 
            do i=1,nparams
               write (f_out,"(a11,1x,i4,2x,g22.12)") parnam(i),ivar(i),param(i)
            enddo 
          endif 
          !
          if (relat_sterr_min > factor_out .and. ilow_sterr/=-1) then
            last_ivar = ilow_sterr
            ! last_param = param
            last_stadev = stadev
            last_i = ilow_sterr
            ! ivar(ilarge_sterr) = 0
            !
            ! The parameter gets zero value  
            !
            ! param(ilarge_sterr) = 0.d0
            write(f_out,"(I6,' parameters are out - ',A10)") Nout
            cycle outer_loop    
          endif 
          tivar = tivar + 1
        enddo

      endif 
      !
      still_run = .false.
    endif 
    still_run = .false.
    !


 ! Here we can stop the fit 
 enddo outer_loop   


!--  printing  out the resulted information 

   inquire(f_res,opened=ifopen)
   if ( ifopen ) then
     rewind(f_res)
   else
     open  (f_res,file='00.res',status='replace' )
   endif

  if (itmax.ne.0) then 
   !write (f_res,"(a80)") (title(i), i=1,4)
   do i=1,nparams
      write (f_res,"(a11,i4,2x,e22.12)") parnam(i),ivar(i),param(i)
   enddo 
   !
   !write (f_res,"(a1)") '!'
   !write (f_res,"(a11,1x,i4,2x,f18.8)") 'alpha     ',ivar(2),180.d0-param(2)
   !write (f_res,"(a11,1x,i4,2x,f18.8)") parnam(11),ivar(11),param(11)
   !write (f_res,"(a11,1x,i4,2x,f18.8)") parnam(12),ivar(12),param(12)
   !write (f_res,"(a11,1x,i4,2x,f18.8)") parnam(1),ivar(1),0.0
   !do i=3,10
   !   write (f_res,"(a11,1x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   !enddo 
   !do i=13,nparams
   !   write (f_res,"(a11,1x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
   !enddo 

  endif  ! ---- itmax/=0

 do nrow=1,npts
   !
   r12=blen1(nrow)*conv_lenght
   r32=blen2(nrow)*conv_lenght
   !
   alpha1=bang(nrow)

    v = energy(nrow)-eps(nrow)
    write (f_res,"(2f14.4,f11.4,3f16.9,2x,d9.2)")            & 
          r12,r32,                                  &
          alpha1,                            &
          energy(nrow)-ve,v-ve,eps(nrow),                   &
          wt(nrow)

 enddo



!--  printing  out the fitted PES in the form of the original ab initio
!--  for refitting with the lineariyed coordinates later 

 open  (f_new,file='new.en',status='replace')
 !
 !open  (f_new,file='new.en',access='append',status='old')
 !


 do nrow=1,npts
   !
   if (wt(nrow)>0) then 
      !
      alpha1=bang(nrow)
      !
      write (f_new,"(f10.2,1x,f18.8,2x,f14.8,2x,f14.8,2x,f20.12,2x,f18.6)") alpha1,param(1),param(5),param(8) ! ,param(17),param(38)
      !
      !r12=blen1(nrow)*conv_lenght
      !r32=blen2(nrow)*conv_lenght

      !v = (energy(nrow)-eps(nrow))*conv_energy
      !write (f_new,"(2f14.8,f8.3,2x,3f20.10,2x,d8.2)")            & 
      !      r12,r32,                                  &
      !      alpha1,                            &
      !      energy(nrow),v,eps(nrow),                   &
      !      wt(nrow)*1e-8
     !
     exit
     !
   endif 
   !
 enddo



 !do nrow=1,npts
 !  !
 !  r12=blen1(nrow)
 !  r32=blen2(nrow)
 !  !
 !  alpha1=bang(nrow)
 !   v = (energy(nrow)-eps(nrow)) !/conv_energy+ZPE
 !   write (f_new,"(2f10.4,f10.4,' ',f22.4,2x,f9.4)") & 
 !        r12,r32,                          &
 !        alpha1,                           &
 !        v, wt(nrow)
 ! enddo


 write (f_res,"(/3x,59('-')/'   | param  | points |   deviat    |    rms      | stability |')")
 write (f_res,"(3x,59('-'))")
 write (f_res,"('   |  ',i3,'   | ',i5,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_res,"(/92('#'))")



 write (f_out,"(/3x,59('-')/'   | param  | points |   deviat    |    rms      | stability |')")
 write (f_out,"(3x,59('-'))")
 write (f_out,"('   |  ',i3,'   | ',i5,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_out,"(/92('#'))")



 end subroutine fitting


  subroutine wifin_113(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r1,r2,th,reoh,thetae,b1,roh,alphaoh
      double precision            :: phh2,t0,ut,x0,ut2,x02,xs1,xs2,xst,rs,rm,rr1,rr2,xep1,xep2,xep3
      double precision            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,vp9,vp10,vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2,v
      integer(4) :: i 
      !
      r1 = local(1)
      r2 = local(2)
      th = local(3)

      !

      !open(unit=32,status='old',form='formatted',file='fit.pot')
      !rewind 32
      !read(32,*) npropin
      !
      reoh=0.958649d0
      thetae=104.3475d0
      b1=2.0d0
      roh=0.951961d0
      alphaoh=2.587949757553683d0
      phh2=6.70164303995d0
      t0=0.01d0
      ut=20.d0
      x0=2.5d0
      ut2=20.d0
      x02=2.5d0
           
      thetae=thetae*.314159265358979312d01*.00555555555555555555d0

      xs1=(r1+r2)*0.5d0-reoh
      xs2=(r1-r2)*0.5d0
      xst=dcos(th)-dcos(thetae)


      rs=dsqrt(0.5d0)*(r1+r2)
      rm=dsqrt(0.5d0)*(r1-r2)

      rr1=r1-roh
      rr2=r2-roh

      xep1=dexp(-2.d0*alphaoh*rr1)-2.d0*dexp(-alphaoh*rr1)+1.d0
      xep2=dexp(-2.d0*alphaoh*rr2)-2.d0*dexp(-alphaoh*rr2)+1.d0
      xep3=dexp(-b1*(rr1**2+rr2**2))
      rhh=dsqrt(r1**2+r2**2-2.d0*r1*r2*dcos(th))
      vhh=0.895714083584240987d6*dexp(-phh2*rhh)
      vpb1=0.368181683997432046d5*xep1
      vpb2=0.368181683997432046d5*xep2

      v0=-(xp(1))
      vp1=-(xp(2))*xst                  &
       -(xp(3))*xs1                     &
       +(xp(4))*xst**2                  &
       -(xp(5))*xs2**2                  &
       -(xp(6))*xs1*xst                 &
       -(xp(7))*xs1**2                  &
       +(xp(8))*xst**3                  &
       -(xp(9))*xs2**2*xst              &
       -(xp(10))*xs1*xst**2             &
       +(xp(11))*xs1*xs2**2             &
       +(xp(12))*xs1**2*xst
     vp2=(xp(13))*xs1**3                &
       +(xp(14))*xst**4                 &
       +(xp(15))*xs2**2*xst**2          &
       -(xp(16))*xs2**4                 &
       -(xp(17))*xs1*xst**3             &
       -(xp(18))*xs1*xs2**2*xst         &
       +(xp(19))*xs1**2*xst**2          &
       -(xp(20))*xs1**2*xs2**2          &
       -(xp(21))*xs1**3*xst             &
       -(xp(22))*xs1**4                 &
       -(xp(23))*xst**5                 &
       -(xp(24))*xs2**2*xst**3          &
       -(xp(25))*xs2**4*xst
     vp3=(xp(26))*xs1*xst**4            &
       -(xp(27))*xs1*xs2**2*xst**2      &
       +(xp(28))*xs1*xs2**4             &
       -(xp(29))*xs1**2*xst**3          &
       -(xp(30))*xs1**2*xs2**2*xst      &
       -(xp(31))*xs1**3*xst**2          &
       +(xp(32))*xs1**3*xs2**2          &
       +(xp(33))*xs1**4*xst             &
       +(xp(34))*xs1**5                 &
       +(xp(35))*xst**6
     vp4=(xp(36))*xs2**2*xst**4         &
       -(xp(37))*xs2**6                 &
       -(xp(38))*xs1*xst**5             &
       -(xp(39))*xs1*xs2**2*xst**3      &
       +(xp(40))*xs1**2*xst**4          &
       +(xp(41))*xs1**2*xs2**2*xst**2   &
       +(xp(42))*xs1**4*xst**2          &
       -(xp(43))*xs1**6                 &
       +(xp(44))*xs2**2*xst**5          &
       -(xp(45))*xs2**6*xst                  &
       -(xp(46))*xs1*xst**6                  &
       -(xp(47))*xs1*xs2**2*xst**4                  &
       -(xp(48))*xs1*xs2**6
     vp5=(xp(49))*xs1**2*xst**5                  &
       -(xp(50))*xs1**3*xs2**4                  &
       -(xp(51))*xs1**4*xs2**2*xst                  &
       +(xp(52))*xs1**5*xst**2                  &
       -(xp(53))*xs1**5*xs2**2                  &
       +(xp(54))*xs1**7                  &
       -(xp(55))*xs2**2*xst**6                  &
       +(xp(56))*xs2**4*xst**4                  &
       +(xp(57))*xs2**6*xst**2                  &
       +(xp(58))*xs2**8                  &
       +(xp(59))*xs1*xs2**2*xst**5                  &
       +(xp(60))*xs1*xs2**4*xst**3                  &
       -(xp(61))*xs1**2*xst**6
     vp6=-(xp(62))*xs1**2*xs2**2*xst**4                  &
       +(xp(63))*xs1**2*xs2**6                  &
       +(xp(64))*xs1**3*xst**5                  &
       +(xp(65))*xs1**3*xs2**2*xst**3                  &
       -(xp(66))*xs1**4*xst**4                  &
       +(xp(67))*xs1**4*xs2**4                  &
       +(xp(68))*xs1**5*xst**3                  &
       +(xp(69))*xs1**5*xs2**2*xst                  &
       -(xp(70))*xs1**6*xst**2
     vp7=(xp(71))*xs1**6*xs2**2                  &
       +(xp(72))*xs1**7*xst                  &
       -(xp(73))*xs2**4*xst**5                  &
       +(xp(74))*xs1*xst**8                  &
       +(xp(75))*xs1*xs2**2*xst**6                  &
       -(xp(76))*xs1*xs2**6*xst**2                  &
       -(xp(77))*xs1*xs2**8                  &
       -(xp(78))*xs1**2*xst**7                  &
       -(xp(79))*xs1**2*xs2**2*xst**5                  &
       +(xp(80))*xs1**3*xst**6                  &
       -(xp(81))*xs1**3*xs2**4*xst**2
     vp8=-(xp(82))*xs1**3*xs2**6                  &
       -(xp(83))*xs1**4*xst**5                  &
       +(xp(84))*xs1**4*xs2**2*xst**3                  &
       -(xp(85))*xs1**5*xs2**2*xst**2                  &
       -(xp(86))*xs1**5*xs2**4                  &
       -(xp(87))*xs1**6*xs2**2*xst                  &
       -(xp(88))*xs1**7*xs2**2                  &
       -(xp(89))*xs1**8*xst                  &
       +(xp(90))*xs2**2*xst**8                  &
       -(xp(91))*xs2**6*xst**4                  &
       -(xp(92))*xs2**8*xst**2
     vp9=-(xp(93))*xs1*xst**9                  &
       -(xp(94))*xs1*xs2**2*xst**7                  &
       -(xp(95))*xs1*xs2**4*xst**5                  &
       +(xp(96))*xs1**2*xst**8                  &
       +(xp(97))*xs1**2*xs2**2*xst**6                  &
       +(xp(98))*xs1**2*xs2**4*xst**4                  &
       +(xp(99))*xs1**2*xs2**6*xst**2                  &
       +(xp(100))*xs1**2*xs2**8                  &
       -(xp(101))*xs1**3*xst**7                  &
       -(xp(102))*xs1**3*xs2**2*xst**5                  &
       -(xp(103))*xs1**3*xs2**4*xst**3                  &
       +(xp(104))*xs1**4*xst**6
     vp10=(xp(105))*xs1**4*xs2**2*xst**4                  &
       +(xp(106))*xs1**4*xs2**6                  &
       -(xp(107))*xs1**5*xs2**2*xst**3                  &
       +(xp(108))*xs1**6*xst**4                  &
       +(xp(109))*xs1**6*xs2**2*xst**2                  &
       +(xp(110))*xs1**6*xs2**4                  &
       +(xp(111))*xs1**9*xst                  &
       -(xp(112))*xst**11                  &
       -(xp(113))*xs2**2*xst**9


       vp=vp1+vp2+vp3+vp4+vp5+vp6+vp7+vp8+vp10+vp9
       !
       vps1=42395.535333d0*xep1
       vps2=42395.535333d0*xep2

       y1=1.d0/(1.d0+dexp(ut*(x0-r1)))
       y2=1.d0/(1.d0+dexp(ut*(x0-r2)))
       y12=1.d0/(1.d0+dexp(ut2*(x02-r1)))
       y22=1.d0/(1.d0+dexp(ut2*(x02-r2)))

       vp=vp*xep3*(1-y12)*(1-y22)
       voh1=vpb1*(1-y1)+y1*vps1
       voh2=vpb2*(1-y2)+y2*vps2

       v=v0+vp+voh1+voh2+vhh

        return
        end subroutine wifin_113





  double precision function wifin_242(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r1,r2,th,reoh,thetae,b1,roh,alphaoh
      double precision            :: phh2,t0,ut,x0,ut2,x02,xs1,xs2,xst,rs,rm,rr1,rr2,xep1,xep2,xep3
      double precision            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2,v
      integer(4) :: i 
      !
      r1 = local(1)
      r2 = local(2)
      th = local(3)

      !
      
      !open(unit=32,status='old',form='formatted',file='fit.pot')
      !rewind 32
      !read(32,*) npropin
      !
      reoh  = xp(2) !  0.958649d0
      thetae= xp(3) ! 104.3475d0
      b1=2.0d0
      roh = xp(4) ! 0.951961d0
      alphaoh = xp(5) !2.587949757553683d0
      phh2= xp(6) !6.70164303995d0
      t0=0.01d0
      ut=20.d0
      x0=2.5d0
      ut2=20.d0
      x02=2.5d0
        
      thetae=thetae*.314159265358979312d01*.00555555555555555555d0

      xs1=(r1+r2)*0.5d0-reoh
      xs2=(r1-r2)*0.5d0
      xst=dcos(th)-dcos(thetae)

      rs=dsqrt(0.5d0)*(r1+r2)
      rm=dsqrt(0.5d0)*(r1-r2)

      rr1=r1-roh
      rr2=r2-roh

      xep1=dexp(-2.d0*alphaoh*rr1)-2.d0*dexp(-alphaoh*rr1)+1.d0
      xep2=dexp(-2.d0*alphaoh*rr2)-2.d0*dexp(-alphaoh*rr2)+1.d0
      xep3=dexp(-b1*(rr1**2+rr2**2))
      rhh=dsqrt(r1**2+r2**2-2.d0*r1*r2*dcos(th))
      vhh=xp(7)*dexp(-phh2*rhh) ! 0.900642240911285975d6
      vpb1=xp(8)*xep1  ! 0.518622556959170834d5
      vpb2=xp(8)*xep2 ! 0.518622556959170834d5


  v0= xp(  1)*xs1**0*xs2**0*xst**0
 vp1= xp(  9)*xs1**0*xs2**0*xst**1&
     +xp( 10)*xs1**1*xs2**0*xst**0&
     +xp( 11)*xs1**0*xs2**0*xst**2&
     +xp( 12)*xs1**0*xs2**2*xst**0&
     +xp( 13)*xs1**1*xs2**0*xst**1&
     +xp( 14)*xs1**2*xs2**0*xst**0&
     +xp( 15)*xs1**0*xs2**0*xst**3&
     +xp( 16)*xs1**0*xs2**2*xst**1&
     +xp( 17)*xs1**1*xs2**0*xst**2&
     +xp( 18)*xs1**1*xs2**2*xst**0&
     +xp( 19)*xs1**2*xs2**0*xst**1&
     +xp( 20)*xs1**3*xs2**0*xst**0&
     +xp( 21)*xs1**0*xs2**0*xst**4&
     +xp( 22)*xs1**0*xs2**2*xst**2&
     +xp( 23)*xs1**0*xs2**4*xst**0&
     +xp( 24)*xs1**1*xs2**0*xst**3&
     +xp( 25)*xs1**1*xs2**2*xst**1&
     +xp( 26)*xs1**2*xs2**0*xst**2&
     +xp( 27)*xs1**2*xs2**2*xst**0&
     +xp( 28)*xs1**3*xs2**0*xst**1&
     +xp( 29)*xs1**4*xs2**0*xst**0&
     +xp( 30)*xs1**0*xs2**0*xst**5&
     +xp( 31)*xs1**0*xs2**2*xst**3&
     +xp( 32)*xs1**0*xs2**4*xst**1&
     +xp( 33)*xs1**1*xs2**0*xst**4&
     +xp( 34)*xs1**1*xs2**2*xst**2&
     +xp( 35)*xs1**1*xs2**4*xst**0&
     +xp( 36)*xs1**2*xs2**0*xst**3&
     +xp( 37)*xs1**2*xs2**2*xst**1&
     +xp( 38)*xs1**3*xs2**0*xst**2&
     +xp( 39)*xs1**3*xs2**2*xst**0&
     +xp( 40)*xs1**4*xs2**0*xst**1&
     +xp( 41)*xs1**5*xs2**0*xst**0&
     +xp( 42)*xs1**0*xs2**0*xst**6&
     +xp( 43)*xs1**0*xs2**2*xst**4&
     +xp( 44)*xs1**0*xs2**4*xst**2&
     +xp( 45)*xs1**0*xs2**6*xst**0&
     +xp( 46)*xs1**1*xs2**0*xst**5&
     +xp( 47)*xs1**1*xs2**2*xst**3&
     +xp( 48)*xs1**1*xs2**4*xst**1&
     +xp( 49)*xs1**2*xs2**0*xst**4&
     +xp( 50)*xs1**2*xs2**2*xst**2&
     +xp( 51)*xs1**2*xs2**4*xst**0&
     +xp( 52)*xs1**3*xs2**0*xst**3&
     +xp( 53)*xs1**3*xs2**2*xst**1&
     +xp( 54)*xs1**4*xs2**0*xst**2&
     +xp( 55)*xs1**4*xs2**2*xst**0&
     +xp( 56)*xs1**5*xs2**0*xst**1&
     +xp( 57)*xs1**6*xs2**0*xst**0&
     +xp( 58)*xs1**0*xs2**0*xst**7&
     +xp( 59)*xs1**0*xs2**2*xst**5&
     +xp( 60)*xs1**0*xs2**4*xst**3&
     +xp( 61)*xs1**0*xs2**6*xst**1&
     +xp( 62)*xs1**1*xs2**0*xst**6&
     +xp( 63)*xs1**1*xs2**2*xst**4&
     +xp( 64)*xs1**1*xs2**4*xst**2&
     +xp( 65)*xs1**1*xs2**6*xst**0&
     +xp( 66)*xs1**2*xs2**0*xst**5&
     +xp( 67)*xs1**2*xs2**2*xst**3&
     +xp( 68)*xs1**2*xs2**4*xst**1&
     +xp( 69)*xs1**3*xs2**0*xst**4&
     +xp( 70)*xs1**3*xs2**2*xst**2&
     +xp( 71)*xs1**3*xs2**4*xst**0&
     +xp( 72)*xs1**4*xs2**0*xst**3&
     +xp( 73)*xs1**4*xs2**2*xst**1&
     +xp( 74)*xs1**5*xs2**0*xst**2&
     +xp( 75)*xs1**5*xs2**2*xst**0&
     +xp( 76)*xs1**6*xs2**0*xst**1&
     +xp( 77)*xs1**7*xs2**0*xst**0&
     +xp( 78)*xs1**0*xs2**0*xst**8&
     +xp( 79)*xs1**0*xs2**2*xst**6&
     +xp( 80)*xs1**0*xs2**4*xst**4&
     +xp( 81)*xs1**0*xs2**6*xst**2&
     +xp( 82)*xs1**0*xs2**8*xst**0&
     +xp( 83)*xs1**1*xs2**0*xst**7&
     +xp( 84)*xs1**1*xs2**2*xst**5&
     +xp( 85)*xs1**1*xs2**4*xst**3&
     +xp( 86)*xs1**1*xs2**6*xst**1&
     +xp( 87)*xs1**2*xs2**0*xst**6&
     +xp( 88)*xs1**2*xs2**2*xst**4&
     +xp( 89)*xs1**2*xs2**4*xst**2&
     +xp( 90)*xs1**2*xs2**6*xst**0&
     +xp( 91)*xs1**3*xs2**0*xst**5&
     +xp( 92)*xs1**3*xs2**2*xst**3&
     +xp( 93)*xs1**3*xs2**4*xst**1&
     +xp( 94)*xs1**4*xs2**0*xst**4&
     +xp( 95)*xs1**4*xs2**2*xst**2&
     +xp( 96)*xs1**4*xs2**4*xst**0&
     +xp( 97)*xs1**5*xs2**0*xst**3&
     +xp( 98)*xs1**5*xs2**2*xst**1&
     +xp( 99)*xs1**6*xs2**0*xst**2
 vp2= xp(100)*xs1**6*xs2**2*xst**0&
     +xp(101)*xs1**7*xs2**0*xst**1&
     +xp(102)*xs1**8*xs2**0*xst**0&
     +xp(103)*xs1**0*xs2**0*xst**9&
     +xp(104)*xs1**0*xs2**2*xst**7&
     +xp(105)*xs1**0*xs2**4*xst**5&
     +xp(106)*xs1**0*xs2**6*xst**3&
     +xp(107)*xs1**0*xs2**8*xst**1&
     +xp(108)*xs1**1*xs2**0*xst**8&
     +xp(109)*xs1**1*xs2**2*xst**6&
     +xp(110)*xs1**1*xs2**4*xst**4&
     +xp(111)*xs1**1*xs2**6*xst**2&
     +xp(112)*xs1**1*xs2**8*xst**0&
     +xp(113)*xs1**2*xs2**0*xst**7&
     +xp(114)*xs1**2*xs2**2*xst**5&
     +xp(115)*xs1**2*xs2**4*xst**3&
     +xp(116)*xs1**2*xs2**6*xst**1&
     +xp(117)*xs1**3*xs2**0*xst**6&
     +xp(118)*xs1**3*xs2**2*xst**4&
     +xp(119)*xs1**3*xs2**4*xst**2&
     +xp(120)*xs1**3*xs2**6*xst**0&
     +xp(121)*xs1**4*xs2**0*xst**5&
     +xp(122)*xs1**4*xs2**2*xst**3&
     +xp(123)*xs1**4*xs2**4*xst**1&
     +xp(124)*xs1**5*xs2**0*xst**4&
     +xp(125)*xs1**5*xs2**2*xst**2&
     +xp(126)*xs1**5*xs2**4*xst**0&
     +xp(127)*xs1**6*xs2**0*xst**3&
     +xp(128)*xs1**6*xs2**2*xst**1&
     +xp(129)*xs1**7*xs2**0*xst**2&
     +xp(130)*xs1**7*xs2**2*xst**0&
     +xp(131)*xs1**8*xs2**0*xst**1&
     +xp(132)*xs1**9*xs2**0*xst**0&
     +xp(133)*xs1**0*xs2**0*xst**10&
     +xp(134)*xs1**0*xs2**2*xst**8&
     +xp(135)*xs1**0*xs2**4*xst**6&
     +xp(136)*xs1**0*xs2**6*xst**4&
     +xp(137)*xs1**0*xs2**8*xst**2&
     +xp(138)*xs1**0*xs2**10*xst**0&
     +xp(139)*xs1**1*xs2**0*xst**9&
     +xp(140)*xs1**1*xs2**2*xst**7&
     +xp(141)*xs1**1*xs2**4*xst**5&
     +xp(142)*xs1**1*xs2**6*xst**3&
     +xp(143)*xs1**1*xs2**8*xst**1&
     +xp(144)*xs1**2*xs2**0*xst**8&
     +xp(145)*xs1**2*xs2**2*xst**6&
     +xp(146)*xs1**2*xs2**4*xst**4&
     +xp(147)*xs1**2*xs2**6*xst**2&
     +xp(148)*xs1**2*xs2**8*xst**0&
     +xp(149)*xs1**3*xs2**0*xst**7&
     +xp(150)*xs1**3*xs2**2*xst**5&
     +xp(151)*xs1**3*xs2**4*xst**3&
     +xp(152)*xs1**3*xs2**6*xst**1&
     +xp(153)*xs1**4*xs2**0*xst**6&
     +xp(154)*xs1**4*xs2**2*xst**4&
     +xp(155)*xs1**4*xs2**4*xst**2&
     +xp(156)*xs1**4*xs2**6*xst**0&
     +xp(157)*xs1**5*xs2**0*xst**5&
     +xp(158)*xs1**5*xs2**2*xst**3&
     +xp(159)*xs1**5*xs2**4*xst**1&
     +xp(160)*xs1**6*xs2**0*xst**4&
     +xp(161)*xs1**6*xs2**2*xst**2&
     +xp(162)*xs1**6*xs2**4*xst**0&
     +xp(163)*xs1**7*xs2**0*xst**3&
     +xp(164)*xs1**7*xs2**2*xst**1&
     +xp(165)*xs1**8*xs2**0*xst**2&
     +xp(166)*xs1**8*xs2**2*xst**0&
     +xp(167)*xs1**9*xs2**0*xst**1&
     +xp(168)*xs1**10*xs2**0*xst**0&
     +xp(169)*xs1**0*xs2**0*xst**11&
     +xp(170)*xs1**0*xs2**2*xst**9&
     +xp(171)*xs1**0*xs2**4*xst**7&
     +xp(172)*xs1**0*xs2**6*xst**5&
     +xp(173)*xs1**0*xs2**8*xst**3&
     +xp(174)*xs1**0*xs2**10*xst**1&
     +xp(175)*xs1**1*xs2**0*xst**10&
     +xp(176)*xs1**1*xs2**2*xst**8&
     +xp(177)*xs1**1*xs2**4*xst**6&
     +xp(178)*xs1**1*xs2**6*xst**4&
     +xp(179)*xs1**1*xs2**8*xst**2&
     +xp(180)*xs1**1*xs2**10*xst**0&
     +xp(181)*xs1**2*xs2**0*xst**9
 vp3= xp(182)*xs1**2*xs2**2*xst**7&
     +xp(183)*xs1**2*xs2**4*xst**5&
     +xp(184)*xs1**2*xs2**6*xst**3&
     +xp(185)*xs1**2*xs2**8*xst**1&
     +xp(186)*xs1**3*xs2**0*xst**8&
     +xp(187)*xs1**3*xs2**2*xst**6&
     +xp(188)*xs1**3*xs2**4*xst**4&
     +xp(189)*xs1**3*xs2**6*xst**2&
     +xp(190)*xs1**3*xs2**8*xst**0&
     +xp(191)*xs1**4*xs2**0*xst**7&
     +xp(192)*xs1**4*xs2**2*xst**5&
     +xp(193)*xs1**4*xs2**4*xst**3&
     +xp(194)*xs1**4*xs2**6*xst**1&
     +xp(195)*xs1**5*xs2**0*xst**6&
     +xp(196)*xs1**5*xs2**2*xst**4&
     +xp(197)*xs1**5*xs2**4*xst**2&
     +xp(198)*xs1**5*xs2**6*xst**0&
     +xp(199)*xs1**6*xs2**0*xst**5&
     +xp(200)*xs1**6*xs2**2*xst**3&
     +xp(201)*xs1**6*xs2**4*xst**1&
     +xp(202)*xs1**7*xs2**0*xst**4&
     +xp(203)*xs1**7*xs2**2*xst**2&
     +xp(204)*xs1**7*xs2**4*xst**0&
     +xp(205)*xs1**8*xs2**0*xst**3&
     +xp(206)*xs1**8*xs2**2*xst**1&
     +xp(207)*xs1**9*xs2**0*xst**2&
     +xp(208)*xs1**9*xs2**2*xst**0&
     +xp(209)*xs1**0*xs2**0*xst**12&
     +xp(210)*xs1**0*xs2**2*xst**10&
     +xp(211)*xs1**0*xs2**4*xst**8&
     +xp(212)*xs1**0*xs2**6*xst**6&
     +xp(213)*xs1**0*xs2**8*xst**4&
     +xp(214)*xs1**0*xs2**10*xst**2&
     +xp(215)*xs1**0*xs2**12*xst**0&
     +xp(216)*xs1**1*xs2**0*xst**11&    
     +xp(217)*xs1**1*xs2**2*xst**9&
     +xp(218)*xs1**1*xs2**4*xst**7&
     +xp(219)*xs1**1*xs2**6*xst**5&
     +xp(220)*xs1**1*xs2**8*xst**3&
     +xp(221)*xs1**1*xs2**10*xst**1&
     +xp(222)*xs1**2*xs2**0*xst**10&
     +xp(223)*xs1**2*xs2**2*xst**8&
     +xp(224)*xs1**2*xs2**4*xst**6&
     +xp(225)*xs1**2*xs2**6*xst**4&
     +xp(226)*xs1**2*xs2**8*xst**2&
     +xp(227)*xs1**2*xs2**10*xst**0&
     +xp(228)*xs1**3*xs2**0*xst**9&
     +xp(229)*xs1**3*xs2**2*xst**7&
     +xp(230)*xs1**3*xs2**4*xst**5&
     +xp(231)*xs1**3*xs2**6*xst**3&
     +xp(232)*xs1**3*xs2**8*xst**1&
     +xp(233)*xs1**4*xs2**0*xst**8&
     +xp(234)*xs1**4*xs2**2*xst**6&
     +xp(235)*xs1**4*xs2**4*xst**4&
     +xp(236)*xs1**4*xs2**6*xst**2&
     +xp(237)*xs1**4*xs2**8*xst**0&
     +xp(238)*xs1**5*xs2**0*xst**7&
     +xp(239)*xs1**5*xs2**2*xst**5&
     +xp(240)*xs1**5*xs2**4*xst**3&
     +xp(241)*xs1**5*xs2**6*xst**1&
     +xp(242)*xs1**6*xs2**0*xst**6&
     +xp(243)*xs1**6*xs2**2*xst**4&
     +xp(244)*xs1**6*xs2**4*xst**2&
     +xp(245)*xs1**6*xs2**6*xst**0&
     +xp(246)*xs1**7*xs2**0*xst**5&
     +xp(247)*xs1**7*xs2**2*xst**3



       vp=vp1+vp2+vp3
       !
       vps1=42395.535333d0*xep1
       vps2=42395.535333d0*xep2

       y1=1.d0/(1.d0+dexp(ut*(x0-r1)))
       y2=1.d0/(1.d0+dexp(ut*(x0-r2)))
       y12=1.d0/(1.d0+dexp(ut2*(x02-r1)))
       y22=1.d0/(1.d0+dexp(ut2*(x02-r2)))

       !vp=vp*xep3*(1-y12)*(1-y22)
       !voh1=vpb1*(1-y1)+y1*vps1
       !voh2=vpb2*(1-y2)+y2*vps2

       vp=vp*xep3
       voh1=vpb1 !*(1-y1)+y1*vps1
       voh2=vpb2 !*(1-y2)+y2*vps2

       wifin_242=v0+vp+voh1+voh2+vhh

        return
        end function wifin_242




  double precision function morse_tyuterev_II(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,aa1,pi,re12,alphae,xst,y1,y2,xs1,xs2,v0,vp1,vp2,vp3
      double precision            :: g1,g2,b1,b2,rhh,vhh
      integer(4) :: i 
      !
      r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)

      pi = 4.0d0 * atan2(1.0d0,1.0d0)

      re12      = xp(2)
      alphae    = xp(3)*pi/180.0d0
      !
      aa1  = xp(4)
      b1   = xp(5)
      b2   = xp(6)
      g1   = xp(7)
      g2   = xp(8)
      !
      rhh=dsqrt(r12**2+r32**2-2.d0*r12*r32*dcos(alpha))
      vhh=b1*dexp(-g1*rhh)+b2*dexp(-g2*rhh**2)
      !
! calculate potential energy function values
!
      xst=cos(alphae)-cos(alpha)
      !st=sin(alphae)-sin(alpha)
!     coro=sin(rhoe)-sin(alpha)

      y1=1.0d+00-exp(-aa1*(r12-re12))
      y2=1.0d+00-exp(-aa1*(r32-re12))
      !xst=cos(alphae)-cos(alpha)
      xst=cotan(alpha*0.5d0)
      !
      xs1 = (y1+y2)*0.5d0
      xs2 = (y1-y2)*0.5d0
      !

  v0= xp(  1)*xs1**0*xs2**0*xst**0
 vp1= xp(  9)*xs1**0*xs2**0*xst**1&
     +xp( 10)*xs1**1*xs2**0*xst**0&
     +xp( 11)*xs1**0*xs2**0*xst**2&
     +xp( 12)*xs1**0*xs2**2*xst**0&
     +xp( 13)*xs1**1*xs2**0*xst**1&
     +xp( 14)*xs1**2*xs2**0*xst**0&
     +xp( 15)*xs1**0*xs2**0*xst**3&
     +xp( 16)*xs1**0*xs2**2*xst**1&
     +xp( 17)*xs1**1*xs2**0*xst**2&
     +xp( 18)*xs1**1*xs2**2*xst**0&
     +xp( 19)*xs1**2*xs2**0*xst**1&
     +xp( 20)*xs1**3*xs2**0*xst**0&
     +xp( 21)*xs1**0*xs2**0*xst**4&
     +xp( 22)*xs1**0*xs2**2*xst**2&
     +xp( 23)*xs1**0*xs2**4*xst**0&
     +xp( 24)*xs1**1*xs2**0*xst**3&
     +xp( 25)*xs1**1*xs2**2*xst**1&
     +xp( 26)*xs1**2*xs2**0*xst**2&
     +xp( 27)*xs1**2*xs2**2*xst**0&
     +xp( 28)*xs1**3*xs2**0*xst**1&
     +xp( 29)*xs1**4*xs2**0*xst**0&
     +xp( 30)*xs1**0*xs2**0*xst**5&
     +xp( 31)*xs1**0*xs2**2*xst**3&
     +xp( 32)*xs1**0*xs2**4*xst**1&
     +xp( 33)*xs1**1*xs2**0*xst**4&
     +xp( 34)*xs1**1*xs2**2*xst**2&
     +xp( 35)*xs1**1*xs2**4*xst**0&
     +xp( 36)*xs1**2*xs2**0*xst**3&
     +xp( 37)*xs1**2*xs2**2*xst**1&
     +xp( 38)*xs1**3*xs2**0*xst**2&
     +xp( 39)*xs1**3*xs2**2*xst**0&
     +xp( 40)*xs1**4*xs2**0*xst**1&
     +xp( 41)*xs1**5*xs2**0*xst**0&
     +xp( 42)*xs1**0*xs2**0*xst**6&
     +xp( 43)*xs1**0*xs2**2*xst**4&
     +xp( 44)*xs1**0*xs2**4*xst**2&
     +xp( 45)*xs1**0*xs2**6*xst**0&
     +xp( 46)*xs1**1*xs2**0*xst**5&
     +xp( 47)*xs1**1*xs2**2*xst**3&
     +xp( 48)*xs1**1*xs2**4*xst**1&
     +xp( 49)*xs1**2*xs2**0*xst**4&
     +xp( 50)*xs1**2*xs2**2*xst**2&
     +xp( 51)*xs1**2*xs2**4*xst**0&
     +xp( 52)*xs1**3*xs2**0*xst**3&
     +xp( 53)*xs1**3*xs2**2*xst**1&
     +xp( 54)*xs1**4*xs2**0*xst**2&
     +xp( 55)*xs1**4*xs2**2*xst**0&
     +xp( 56)*xs1**5*xs2**0*xst**1&
     +xp( 57)*xs1**6*xs2**0*xst**0&
     +xp( 58)*xs1**0*xs2**0*xst**7&
     +xp( 59)*xs1**0*xs2**2*xst**5&
     +xp( 60)*xs1**0*xs2**4*xst**3&
     +xp( 61)*xs1**0*xs2**6*xst**1&
     +xp( 62)*xs1**1*xs2**0*xst**6&
     +xp( 63)*xs1**1*xs2**2*xst**4&
     +xp( 64)*xs1**1*xs2**4*xst**2&
     +xp( 65)*xs1**1*xs2**6*xst**0&
     +xp( 66)*xs1**2*xs2**0*xst**5&
     +xp( 67)*xs1**2*xs2**2*xst**3&
     +xp( 68)*xs1**2*xs2**4*xst**1&
     +xp( 69)*xs1**3*xs2**0*xst**4&
     +xp( 70)*xs1**3*xs2**2*xst**2&
     +xp( 71)*xs1**3*xs2**4*xst**0&
     +xp( 72)*xs1**4*xs2**0*xst**3&
     +xp( 73)*xs1**4*xs2**2*xst**1&
     +xp( 74)*xs1**5*xs2**0*xst**2&
     +xp( 75)*xs1**5*xs2**2*xst**0&
     +xp( 76)*xs1**6*xs2**0*xst**1&
     +xp( 77)*xs1**7*xs2**0*xst**0&
     +xp( 78)*xs1**0*xs2**0*xst**8&
     +xp( 79)*xs1**0*xs2**2*xst**6&
     +xp( 80)*xs1**0*xs2**4*xst**4&
     +xp( 81)*xs1**0*xs2**6*xst**2&
     +xp( 82)*xs1**0*xs2**8*xst**0&
     +xp( 83)*xs1**1*xs2**0*xst**7&
     +xp( 84)*xs1**1*xs2**2*xst**5&
     +xp( 85)*xs1**1*xs2**4*xst**3&
     +xp( 86)*xs1**1*xs2**6*xst**1&
     +xp( 87)*xs1**2*xs2**0*xst**6&
     +xp( 88)*xs1**2*xs2**2*xst**4&
     +xp( 89)*xs1**2*xs2**4*xst**2&
     +xp( 90)*xs1**2*xs2**6*xst**0&
     +xp( 91)*xs1**3*xs2**0*xst**5&
     +xp( 92)*xs1**3*xs2**2*xst**3&
     +xp( 93)*xs1**3*xs2**4*xst**1&
     +xp( 94)*xs1**4*xs2**0*xst**4&
     +xp( 95)*xs1**4*xs2**2*xst**2&
     +xp( 96)*xs1**4*xs2**4*xst**0&
     +xp( 97)*xs1**5*xs2**0*xst**3&
     +xp( 98)*xs1**5*xs2**2*xst**1&
     +xp( 99)*xs1**6*xs2**0*xst**2
 vp2= xp(100)*xs1**6*xs2**2*xst**0&
     +xp(101)*xs1**7*xs2**0*xst**1&
     +xp(102)*xs1**8*xs2**0*xst**0&
     +xp(103)*xs1**0*xs2**0*xst**9&
     +xp(104)*xs1**0*xs2**2*xst**7&
     +xp(105)*xs1**0*xs2**4*xst**5&
     +xp(106)*xs1**0*xs2**6*xst**3&
     +xp(107)*xs1**0*xs2**8*xst**1&
     +xp(108)*xs1**1*xs2**0*xst**8&
     +xp(109)*xs1**1*xs2**2*xst**6&
     +xp(110)*xs1**1*xs2**4*xst**4&
     +xp(111)*xs1**1*xs2**6*xst**2&
     +xp(112)*xs1**1*xs2**8*xst**0&
     +xp(113)*xs1**2*xs2**0*xst**7&
     +xp(114)*xs1**2*xs2**2*xst**5&
     +xp(115)*xs1**2*xs2**4*xst**3&
     +xp(116)*xs1**2*xs2**6*xst**1&
     +xp(117)*xs1**3*xs2**0*xst**6&
     +xp(118)*xs1**3*xs2**2*xst**4&
     +xp(119)*xs1**3*xs2**4*xst**2&
     +xp(120)*xs1**3*xs2**6*xst**0&
     +xp(121)*xs1**4*xs2**0*xst**5&
     +xp(122)*xs1**4*xs2**2*xst**3&
     +xp(123)*xs1**4*xs2**4*xst**1&
     +xp(124)*xs1**5*xs2**0*xst**4&
     +xp(125)*xs1**5*xs2**2*xst**2&
     +xp(126)*xs1**5*xs2**4*xst**0&
     +xp(127)*xs1**6*xs2**0*xst**3&
     +xp(128)*xs1**6*xs2**2*xst**1&
     +xp(129)*xs1**7*xs2**0*xst**2&
     +xp(130)*xs1**7*xs2**2*xst**0&
     +xp(131)*xs1**8*xs2**0*xst**1&
     +xp(132)*xs1**9*xs2**0*xst**0&
     +xp(133)*xs1**0*xs2**0*xst**10&
     +xp(134)*xs1**0*xs2**2*xst**8&
     +xp(135)*xs1**0*xs2**4*xst**6&
     +xp(136)*xs1**0*xs2**6*xst**4&
     +xp(137)*xs1**0*xs2**8*xst**2&
     +xp(138)*xs1**0*xs2**10*xst**0&
     +xp(139)*xs1**1*xs2**0*xst**9&
     +xp(140)*xs1**1*xs2**2*xst**7&
     +xp(141)*xs1**1*xs2**4*xst**5&
     +xp(142)*xs1**1*xs2**6*xst**3&
     +xp(143)*xs1**1*xs2**8*xst**1&
     +xp(144)*xs1**2*xs2**0*xst**8&
     +xp(145)*xs1**2*xs2**2*xst**6&
     +xp(146)*xs1**2*xs2**4*xst**4&
     +xp(147)*xs1**2*xs2**6*xst**2&
     +xp(148)*xs1**2*xs2**8*xst**0&
     +xp(149)*xs1**3*xs2**0*xst**7&
     +xp(150)*xs1**3*xs2**2*xst**5&
     +xp(151)*xs1**3*xs2**4*xst**3&
     +xp(152)*xs1**3*xs2**6*xst**1&
     +xp(153)*xs1**4*xs2**0*xst**6&
     +xp(154)*xs1**4*xs2**2*xst**4&
     +xp(155)*xs1**4*xs2**4*xst**2&
     +xp(156)*xs1**4*xs2**6*xst**0&
     +xp(157)*xs1**5*xs2**0*xst**5&
     +xp(158)*xs1**5*xs2**2*xst**3&
     +xp(159)*xs1**5*xs2**4*xst**1&
     +xp(160)*xs1**6*xs2**0*xst**4&
     +xp(161)*xs1**6*xs2**2*xst**2&
     +xp(162)*xs1**6*xs2**4*xst**0&
     +xp(163)*xs1**7*xs2**0*xst**3&
     +xp(164)*xs1**7*xs2**2*xst**1&
     +xp(165)*xs1**8*xs2**0*xst**2&
     +xp(166)*xs1**8*xs2**2*xst**0&
     +xp(167)*xs1**9*xs2**0*xst**1&
     +xp(168)*xs1**10*xs2**0*xst**0&
     +xp(169)*xs1**0*xs2**0*xst**11&
     +xp(170)*xs1**0*xs2**2*xst**9&
     +xp(171)*xs1**0*xs2**4*xst**7&
     +xp(172)*xs1**0*xs2**6*xst**5&
     +xp(173)*xs1**0*xs2**8*xst**3&
     +xp(174)*xs1**0*xs2**10*xst**1&
     +xp(175)*xs1**1*xs2**0*xst**10&
     +xp(176)*xs1**1*xs2**2*xst**8&
     +xp(177)*xs1**1*xs2**4*xst**6&
     +xp(178)*xs1**1*xs2**6*xst**4&
     +xp(179)*xs1**1*xs2**8*xst**2&
     +xp(180)*xs1**1*xs2**10*xst**0&
     +xp(181)*xs1**2*xs2**0*xst**9
 vp3= xp(182)*xs1**2*xs2**2*xst**7&
     +xp(183)*xs1**2*xs2**4*xst**5&
     +xp(184)*xs1**2*xs2**6*xst**3&
     +xp(185)*xs1**2*xs2**8*xst**1&
     +xp(186)*xs1**3*xs2**0*xst**8&
     +xp(187)*xs1**3*xs2**2*xst**6&
     +xp(188)*xs1**3*xs2**4*xst**4&
     +xp(189)*xs1**3*xs2**6*xst**2&
     +xp(190)*xs1**3*xs2**8*xst**0&
     +xp(191)*xs1**4*xs2**0*xst**7&
     +xp(192)*xs1**4*xs2**2*xst**5&
     +xp(193)*xs1**4*xs2**4*xst**3&
     +xp(194)*xs1**4*xs2**6*xst**1&
     +xp(195)*xs1**5*xs2**0*xst**6&
     +xp(196)*xs1**5*xs2**2*xst**4&
     +xp(197)*xs1**5*xs2**4*xst**2&
     +xp(198)*xs1**5*xs2**6*xst**0&
     +xp(199)*xs1**6*xs2**0*xst**5&
     +xp(200)*xs1**6*xs2**2*xst**3&
     +xp(201)*xs1**6*xs2**4*xst**1&
     +xp(202)*xs1**7*xs2**0*xst**4&
     +xp(203)*xs1**7*xs2**2*xst**2&
     +xp(204)*xs1**7*xs2**4*xst**0&
     +xp(205)*xs1**8*xs2**0*xst**3&
     +xp(206)*xs1**8*xs2**2*xst**1&
     +xp(207)*xs1**9*xs2**0*xst**2&
     +xp(208)*xs1**9*xs2**2*xst**0&
     +xp(209)*xs1**0*xs2**0*xst**12&
     +xp(210)*xs1**0*xs2**2*xst**10&
     +xp(211)*xs1**0*xs2**4*xst**8&
     +xp(212)*xs1**0*xs2**6*xst**6&
     +xp(213)*xs1**0*xs2**8*xst**4&
     +xp(214)*xs1**0*xs2**10*xst**2&
     +xp(215)*xs1**0*xs2**12*xst**0&
     +xp(216)*xs1**1*xs2**0*xst**11&    
     +xp(217)*xs1**1*xs2**2*xst**9&
     +xp(218)*xs1**1*xs2**4*xst**7&
     +xp(219)*xs1**1*xs2**6*xst**5&
     +xp(220)*xs1**1*xs2**8*xst**3&
     +xp(221)*xs1**1*xs2**10*xst**1&
     +xp(222)*xs1**2*xs2**0*xst**10&
     +xp(223)*xs1**2*xs2**2*xst**8&
     +xp(224)*xs1**2*xs2**4*xst**6&
     +xp(225)*xs1**2*xs2**6*xst**4&
     +xp(226)*xs1**2*xs2**8*xst**2&
     +xp(227)*xs1**2*xs2**10*xst**0&
     +xp(228)*xs1**3*xs2**0*xst**9&
     +xp(229)*xs1**3*xs2**2*xst**7&
     +xp(230)*xs1**3*xs2**4*xst**5&
     +xp(231)*xs1**3*xs2**6*xst**3&
     +xp(232)*xs1**3*xs2**8*xst**1&
     +xp(233)*xs1**4*xs2**0*xst**8&
     +xp(234)*xs1**4*xs2**2*xst**6&
     +xp(235)*xs1**4*xs2**4*xst**4&
     +xp(236)*xs1**4*xs2**6*xst**2&
     +xp(237)*xs1**4*xs2**8*xst**0&
     +xp(238)*xs1**5*xs2**0*xst**7&
     +xp(239)*xs1**5*xs2**2*xst**5&
     +xp(240)*xs1**5*xs2**4*xst**3&
     +xp(241)*xs1**5*xs2**6*xst**1&
     +xp(242)*xs1**6*xs2**0*xst**6&
     +xp(243)*xs1**6*xs2**2*xst**4&
     +xp(244)*xs1**6*xs2**4*xst**2&
     +xp(245)*xs1**6*xs2**6*xst**0&
     +xp(246)*xs1**7*xs2**0*xst**5&
     +xp(247)*xs1**7*xs2**2*xst**3



       morse_tyuterev_II=v0+vp1+vp2+vp3+vhh

        return
        end function morse_tyuterev_II





  double precision function potential_xy2_tyuterev_damp(local,ZPE,N,force) result (f)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  N
      double precision,intent(in) ::  force(N)
      double precision            :: r12,r32,alpha,aa1,pi,re12,alphae,xst,y1,y2,y3,y(3),beta,damp,vlong,vdamp,de
      double precision            :: v4,v5,v6,v7,v8,v0,v1,v2,v3
      !
      r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
      !
      pi = 4.0d0 * atan2(1.0d0,1.0d0)
      !
      re12      = force(2)
      alphae    = force(3)*pi/180.0d0
   !
   De     = force(4)
   beta   = 2.1d0
   damp   = force(5)
   !
   y(1) = r12-re12
   y(2) = r32-re12
   y(3) = alpha-alphae
   !
   vlong = De*(1.0d0-exp(-beta*y(1)))**2 + De*(1.0d0-exp(-beta*y(2)))**2
   vdamp = exp( -damp*( y(1)**2+y(2)**2)-0.0*damp*( y(1)**4+y(2)**4 )-0.*damp*y(3)**2 )
   !
   ! calculate potential energy function values
   !
   y1=r12-re12 ! 1.0d+00-exp(-beta*(r12-re12))
   y2=r32-re12 ! 1.0d+00-exp(-beta*(r32-re12))
   !
   y3=cos(alpha)-cos(alphae)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
 v0 = force( 1)*y1**0*y2**0*y3**0
 v1 = force( 6)*y1**0*y2**0*y3**1& 
    + force( 7)*y1**1*y2**0*y3**0& 
    + force( 7)*y1**0*y2**1*y3**0
 v2 = force( 8)*y1**0*y2**0*y3**2& 
    + force( 9)*y1**1*y2**0*y3**1& 
    + force( 9)*y1**0*y2**1*y3**1& 
    + force(10)*y1**1*y2**1*y3**0& 
    + force(11)*y1**2*y2**0*y3**0& 
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3& 
    + force(13)*y1**1*y2**0*y3**2& 
    + force(13)*y1**0*y2**1*y3**2& 
    + force(14)*y1**1*y2**1*y3**1& 
    + force(15)*y1**2*y2**0*y3**1& 
    + force(15)*y1**0*y2**2*y3**1& 
    + force(16)*y1**2*y2**1*y3**0& 
    + force(16)*y1**1*y2**2*y3**0& 
    + force(17)*y1**3*y2**0*y3**0& 
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then 
  v4 =force(18)*y1**0*y2**0*y3**4& 
    + force(19)*y1**1*y2**0*y3**3& 
    + force(19)*y1**0*y2**1*y3**3& 
    + force(20)*y1**1*y2**1*y3**2& 
    + force(21)*y1**2*y2**0*y3**2& 
    + force(21)*y1**0*y2**2*y3**2& 
    + force(22)*y1**2*y2**1*y3**1& 
    + force(22)*y1**1*y2**2*y3**1& 
    + force(23)*y1**2*y2**2*y3**0& 
    + force(24)*y1**3*y2**0*y3**1& 
    + force(24)*y1**0*y2**3*y3**1& 
    + force(25)*y1**3*y2**1*y3**0& 
    + force(25)*y1**1*y2**3*y3**0& 
    + force(26)*y1**4*y2**0*y3**0& 
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then 
  v5 =force(27)*y1**0*y2**0*y3**5& 
    + force(28)*y1**1*y2**0*y3**4& 
    + force(28)*y1**0*y2**1*y3**4& 
    + force(29)*y1**1*y2**1*y3**3& 
    + force(30)*y1**2*y2**0*y3**3& 
    + force(30)*y1**0*y2**2*y3**3& 
    + force(31)*y1**2*y2**1*y3**2& 
    + force(31)*y1**1*y2**2*y3**2& 
    + force(32)*y1**2*y2**2*y3**1& 
    + force(33)*y1**3*y2**0*y3**2& 
    + force(33)*y1**0*y2**3*y3**2& 
    + force(34)*y1**3*y2**1*y3**1& 
    + force(34)*y1**1*y2**3*y3**1& 
    + force(35)*y1**3*y2**2*y3**0& 
    + force(35)*y1**2*y2**3*y3**0& 
    + force(36)*y1**4*y2**0*y3**1& 
    + force(36)*y1**0*y2**4*y3**1& 
    + force(37)*y1**4*y2**1*y3**0& 
    + force(37)*y1**1*y2**4*y3**0& 
    + force(38)*y1**5*y2**0*y3**0& 
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then 
  v6 = force(39)*y1**0*y2**0*y3**6& 
    + force(40)*y1**1*y2**0*y3**5& 
    + force(40)*y1**0*y2**1*y3**5& 
    + force(41)*y1**1*y2**1*y3**4& 
    + force(42)*y1**2*y2**0*y3**4& 
    + force(42)*y1**0*y2**2*y3**4& 
    + force(43)*y1**2*y2**1*y3**3& 
    + force(43)*y1**1*y2**2*y3**3& 
    + force(44)*y1**2*y2**2*y3**2& 
    + force(45)*y1**3*y2**0*y3**3& 
    + force(45)*y1**0*y2**3*y3**3& 
    + force(46)*y1**3*y2**1*y3**2& 
    + force(46)*y1**1*y2**3*y3**2& 
    + force(47)*y1**3*y2**2*y3**1& 
    + force(47)*y1**2*y2**3*y3**1& 
    + force(48)*y1**3*y2**3*y3**0& 
    + force(49)*y1**4*y2**0*y3**2& 
    + force(49)*y1**0*y2**4*y3**2& 
    + force(50)*y1**4*y2**1*y3**1& 
    + force(50)*y1**1*y2**4*y3**1& 
    + force(51)*y1**4*y2**2*y3**0& 
    + force(51)*y1**2*y2**4*y3**0& 
    + force(52)*y1**5*y2**0*y3**1& 
    + force(52)*y1**0*y2**5*y3**1& 
    + force(53)*y1**5*y2**1*y3**0& 
    + force(53)*y1**1*y2**5*y3**0& 
    + force(54)*y1**6*y2**0*y3**0& 
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then 
 v7 = force(55)*y1**0*y2**0*y3**7& 
    + force(56)*y1**1*y2**0*y3**6& 
    + force(56)*y1**0*y2**1*y3**6& 
    + force(57)*y1**1*y2**1*y3**5& 
    + force(58)*y1**2*y2**0*y3**5& 
    + force(58)*y1**0*y2**2*y3**5& 
    + force(59)*y1**2*y2**1*y3**4& 
    + force(59)*y1**1*y2**2*y3**4& 
    + force(60)*y1**2*y2**2*y3**3& 
    + force(61)*y1**3*y2**0*y3**4& 
    + force(61)*y1**0*y2**3*y3**4& 
    + force(62)*y1**3*y2**1*y3**3& 
    + force(62)*y1**1*y2**3*y3**3& 
    + force(63)*y1**3*y2**2*y3**2& 
    + force(63)*y1**2*y2**3*y3**2& 
    + force(64)*y1**3*y2**3*y3**1& 
    + force(65)*y1**4*y2**0*y3**3& 
    + force(65)*y1**0*y2**4*y3**3& 
    + force(66)*y1**4*y2**1*y3**2& 
    + force(66)*y1**1*y2**4*y3**2& 
    + force(67)*y1**4*y2**2*y3**1& 
    + force(67)*y1**2*y2**4*y3**1& 
    + force(68)*y1**4*y2**3*y3**0& 
    + force(68)*y1**3*y2**4*y3**0& 
    + force(69)*y1**5*y2**0*y3**2& 
    + force(69)*y1**0*y2**5*y3**2& 
    + force(70)*y1**5*y2**1*y3**1& 
    + force(70)*y1**1*y2**5*y3**1& 
    + force(71)*y1**5*y2**2*y3**0& 
    + force(71)*y1**2*y2**5*y3**0& 
    + force(72)*y1**6*y2**0*y3**1& 
    + force(72)*y1**0*y2**6*y3**1& 
    + force(73)*y1**6*y2**1*y3**0& 
    + force(73)*y1**1*y2**6*y3**0& 
    + force(74)*y1**7*y2**0*y3**0& 
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then 
 v8 = force(75)*y1**0*y2**0*y3**8& 
    + force(76)*y1**1*y2**0*y3**7& 
    + force(76)*y1**0*y2**1*y3**7& 
    + force(77)*y1**1*y2**1*y3**6& 
    + force(78)*y1**2*y2**0*y3**6& 
    + force(78)*y1**0*y2**2*y3**6& 
    + force(79)*y1**2*y2**1*y3**5& 
    + force(79)*y1**1*y2**2*y3**5& 
    + force(80)*y1**2*y2**2*y3**4& 
    + force(81)*y1**3*y2**0*y3**5& 
    + force(81)*y1**0*y2**3*y3**5& 
    + force(82)*y1**3*y2**1*y3**4& 
    + force(82)*y1**1*y2**3*y3**4& 
    + force(83)*y1**3*y2**2*y3**3& 
    + force(83)*y1**2*y2**3*y3**3& 
    + force(84)*y1**3*y2**3*y3**2& 
    + force(85)*y1**4*y2**0*y3**4& 
    + force(85)*y1**0*y2**4*y3**4& 
    + force(86)*y1**4*y2**1*y3**3& 
    + force(86)*y1**1*y2**4*y3**3& 
    + force(87)*y1**4*y2**2*y3**2& 
    + force(87)*y1**2*y2**4*y3**2& 
    + force(88)*y1**4*y2**3*y3**1& 
    + force(88)*y1**3*y2**4*y3**1& 
    + force(89)*y1**4*y2**4*y3**0& 
    + force(90)*y1**5*y2**0*y3**3& 
    + force(90)*y1**0*y2**5*y3**3& 
    + force(91)*y1**5*y2**1*y3**2& 
    + force(91)*y1**1*y2**5*y3**2& 
    + force(92)*y1**5*y2**2*y3**1& 
    + force(92)*y1**2*y2**5*y3**1& 
    + force(93)*y1**5*y2**3*y3**0& 
    + force(93)*y1**3*y2**5*y3**0& 
    + force(94)*y1**6*y2**0*y3**2& 
    + force(94)*y1**0*y2**6*y3**2& 
    + force(95)*y1**6*y2**1*y3**1& 
    + force(95)*y1**1*y2**6*y3**1& 
    + force(96)*y1**6*y2**2*y3**0& 
    + force(96)*y1**2*y2**6*y3**0& 
    + force(97)*y1**7*y2**0*y3**1& 
    + force(97)*y1**0*y2**7*y3**1& 
    + force(98)*y1**7*y2**1*y3**0& 
    + force(98)*y1**1*y2**7*y3**0& 
    + force(99)*y1**8*y2**0*y3**0& 
    + force(99)*y1**0*y2**8*y3**0
endif


    f=(v0+v1+v2+v3+v4+v5+v6+v7+v8)*vdamp+vlong

return
end function potential_xy2_tyuterev_damp





  double precision function wifin_242_a(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r1,r2,th,reoh,thetae,b1,roh,alphaoh
      double precision            :: phh2,t0,ut,x0,ut2,x02,xs1,xs2,xst,rs,rm,rr1,rr2,xep1,xep2,xep3
      double precision            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2,v
      integer(4) :: i 
      !
      r1 = local(1)
      r2 = local(2)
      th = local(3)

      !
      
      !open(unit=32,status='old',form='formatted',file='fit.pot')
      !rewind 32
      !read(32,*) npropin
      !
      reoh  = xp(2) !  0.958649d0
      thetae= xp(3) ! 104.3475d0
      b1=2.0d0
      roh = xp(4) ! 0.951961d0
      alphaoh = xp(5) !2.587949757553683d0
      phh2= xp(6) !6.70164303995d0
      t0=0.01d0
      ut=20.d0
      x0=2.5d0
      ut2=20.d0
      x02=2.5d0
        
      thetae=thetae*.314159265358979312d01*.00555555555555555555d0

      xs1=(r1+r2)*0.5d0-reoh
      xs2=(r1-r2)*0.5d0
      xst=dcos(th)-dcos(thetae)

      rs=dsqrt(0.5d0)*(r1+r2)
      rm=dsqrt(0.5d0)*(r1-r2)

      rr1=r1-roh
      rr2=r2-roh

      xep1=dexp(-2.d0*alphaoh*rr1)-2.d0*dexp(-alphaoh*rr1)+1.d0
      xep2=dexp(-2.d0*alphaoh*rr2)-2.d0*dexp(-alphaoh*rr2)+1.d0
      xep3=dexp(-b1*(rr1**2+rr2**2))
      rhh=dsqrt(r1**2+r2**2-2.d0*r1*r2*dcos(th))
      vhh=xp(7)*dexp(-phh2*rhh) ! 0.900642240911285975d6
      vpb1=xp(8)*xep1  ! 0.518622556959170834d5
      vpb2=xp(8)*xep2 ! 0.518622556959170834d5


  v0= xp(  1)*xs1**0*xs2**0*xst**0
 vp1= xp(  9)*xs1**0*xs2**0*xst**1&
     +xp( 10)*xs1**1*xs2**0*xst**0&
     +xp( 11)*xs1**0*xs2**0*xst**2&
     +xp( 12)*xs1**0*xs2**2*xst**0&
     +xp( 13)*xs1**1*xs2**0*xst**1&
     +xp( 14)*xs1**2*xs2**0*xst**0&
     +xp( 15)*xs1**0*xs2**0*xst**3&
     +xp( 16)*xs1**0*xs2**2*xst**1&
     +xp( 17)*xs1**1*xs2**0*xst**2&
     +xp( 18)*xs1**1*xs2**2*xst**0&
     +xp( 19)*xs1**2*xs2**0*xst**1&
     +xp( 20)*xs1**3*xs2**0*xst**0&
     +xp( 21)*xs1**0*xs2**0*xst**4&
     +xp( 22)*xs1**0*xs2**2*xst**2&
     +xp( 23)*xs1**0*xs2**4*xst**0&
     +xp( 24)*xs1**1*xs2**0*xst**3&
     +xp( 25)*xs1**1*xs2**2*xst**1&
     +xp( 26)*xs1**2*xs2**0*xst**2&
     +xp( 27)*xs1**2*xs2**2*xst**0&
     +xp( 28)*xs1**3*xs2**0*xst**1&
     +xp( 29)*xs1**4*xs2**0*xst**0&
     +xp( 30)*xs1**0*xs2**0*xst**5&
     +xp( 31)*xs1**0*xs2**2*xst**3&
     +xp( 32)*xs1**0*xs2**4*xst**1&
     +xp( 33)*xs1**1*xs2**0*xst**4&
     +xp( 34)*xs1**1*xs2**2*xst**2&
     +xp( 35)*xs1**1*xs2**4*xst**0&
     +xp( 36)*xs1**2*xs2**0*xst**3&
     +xp( 37)*xs1**2*xs2**2*xst**1&
     +xp( 38)*xs1**3*xs2**0*xst**2&
     +xp( 39)*xs1**3*xs2**2*xst**0&
     +xp( 40)*xs1**4*xs2**0*xst**1&
     +xp( 41)*xs1**5*xs2**0*xst**0&
     +xp( 42)*xs1**0*xs2**0*xst**6&
     +xp( 43)*xs1**0*xs2**2*xst**4&
     +xp( 44)*xs1**0*xs2**4*xst**2&
     +xp( 45)*xs1**0*xs2**6*xst**0&
     +xp( 46)*xs1**1*xs2**0*xst**5&
     +xp( 47)*xs1**1*xs2**2*xst**3&
     +xp( 48)*xs1**1*xs2**4*xst**1&
     +xp( 49)*xs1**2*xs2**0*xst**4&
     +xp( 50)*xs1**2*xs2**2*xst**2&
     +xp( 51)*xs1**2*xs2**4*xst**0&
     +xp( 52)*xs1**3*xs2**0*xst**3&
     +xp( 53)*xs1**3*xs2**2*xst**1&
     +xp( 54)*xs1**4*xs2**0*xst**2&
     +xp( 55)*xs1**4*xs2**2*xst**0&
     +xp( 56)*xs1**5*xs2**0*xst**1&
     +xp( 57)*xs1**6*xs2**0*xst**0&
     +xp( 58)*xs1**0*xs2**0*xst**7&
     +xp( 59)*xs1**0*xs2**2*xst**5&
     +xp( 60)*xs1**0*xs2**4*xst**3&
     +xp( 61)*xs1**0*xs2**6*xst**1&
     +xp( 62)*xs1**1*xs2**0*xst**6&
     +xp( 63)*xs1**1*xs2**2*xst**4&
     +xp( 64)*xs1**1*xs2**4*xst**2&
     +xp( 65)*xs1**1*xs2**6*xst**0&
     +xp( 66)*xs1**2*xs2**0*xst**5&
     +xp( 67)*xs1**2*xs2**2*xst**3&
     +xp( 68)*xs1**2*xs2**4*xst**1&
     +xp( 69)*xs1**3*xs2**0*xst**4&
     +xp( 70)*xs1**3*xs2**2*xst**2&
     +xp( 71)*xs1**3*xs2**4*xst**0&
     +xp( 72)*xs1**4*xs2**0*xst**3&
     +xp( 73)*xs1**4*xs2**2*xst**1&
     +xp( 74)*xs1**5*xs2**0*xst**2&
     +xp( 75)*xs1**5*xs2**2*xst**0&
     +xp( 76)*xs1**6*xs2**0*xst**1&
     +xp( 77)*xs1**7*xs2**0*xst**0&
     +xp( 78)*xs1**0*xs2**0*xst**8&
     +xp( 79)*xs1**0*xs2**2*xst**6&
     +xp( 80)*xs1**0*xs2**4*xst**4&
     +xp( 81)*xs1**0*xs2**6*xst**2&
     +xp( 82)*xs1**0*xs2**8*xst**0&
     +xp( 83)*xs1**1*xs2**0*xst**7&
     +xp( 84)*xs1**1*xs2**2*xst**5&
     +xp( 85)*xs1**1*xs2**4*xst**3&
     +xp( 86)*xs1**1*xs2**6*xst**1&
     +xp( 87)*xs1**2*xs2**0*xst**6&
     +xp( 88)*xs1**2*xs2**2*xst**4&
     +xp( 89)*xs1**2*xs2**4*xst**2&
     +xp( 90)*xs1**2*xs2**6*xst**0&
     +xp( 91)*xs1**3*xs2**0*xst**5&
     +xp( 92)*xs1**3*xs2**2*xst**3&
     +xp( 93)*xs1**3*xs2**4*xst**1&
     +xp( 94)*xs1**4*xs2**0*xst**4&
     +xp( 95)*xs1**4*xs2**2*xst**2&
     +xp( 96)*xs1**4*xs2**4*xst**0&
     +xp( 97)*xs1**5*xs2**0*xst**3&
     +xp( 98)*xs1**5*xs2**2*xst**1&
     +xp( 99)*xs1**6*xs2**0*xst**2
 vp2= xp(100)*xs1**6*xs2**2*xst**0&
     +xp(101)*xs1**7*xs2**0*xst**1&
     +xp(102)*xs1**8*xs2**0*xst**0&
     +xp(103)*xs1**0*xs2**0*xst**9&
     +xp(104)*xs1**0*xs2**2*xst**7&
     +xp(105)*xs1**0*xs2**4*xst**5&
     +xp(106)*xs1**0*xs2**6*xst**3&
     +xp(107)*xs1**0*xs2**8*xst**1&
     +xp(108)*xs1**1*xs2**0*xst**8&
     +xp(109)*xs1**1*xs2**2*xst**6&
     +xp(110)*xs1**1*xs2**4*xst**4&
     +xp(111)*xs1**1*xs2**6*xst**2&
     +xp(112)*xs1**1*xs2**8*xst**0&
     +xp(113)*xs1**2*xs2**0*xst**7&
     +xp(114)*xs1**2*xs2**2*xst**5&
     +xp(115)*xs1**2*xs2**4*xst**3&
     +xp(116)*xs1**2*xs2**6*xst**1&
     +xp(117)*xs1**3*xs2**0*xst**6&
     +xp(118)*xs1**3*xs2**2*xst**4&
     +xp(119)*xs1**3*xs2**4*xst**2&
     +xp(120)*xs1**3*xs2**6*xst**0&
     +xp(121)*xs1**4*xs2**0*xst**5&
     +xp(122)*xs1**4*xs2**2*xst**3&
     +xp(123)*xs1**4*xs2**4*xst**1&
     +xp(124)*xs1**5*xs2**0*xst**4&
     +xp(125)*xs1**5*xs2**2*xst**2&
     +xp(126)*xs1**5*xs2**4*xst**0&
     +xp(127)*xs1**6*xs2**0*xst**3&
     +xp(128)*xs1**6*xs2**2*xst**1&
     +xp(129)*xs1**7*xs2**0*xst**2&
     +xp(130)*xs1**7*xs2**2*xst**0&
     +xp(131)*xs1**8*xs2**0*xst**1&
     +xp(132)*xs1**9*xs2**0*xst**0&
     +xp(133)*xs1**0*xs2**0*xst**10&
     +xp(134)*xs1**0*xs2**2*xst**8&
     +xp(135)*xs1**0*xs2**4*xst**6&
     +xp(136)*xs1**0*xs2**6*xst**4&
     +xp(137)*xs1**0*xs2**8*xst**2&
     +xp(138)*xs1**0*xs2**10*xst**0&
     +xp(139)*xs1**1*xs2**0*xst**9&
     +xp(140)*xs1**1*xs2**2*xst**7&
     +xp(141)*xs1**1*xs2**4*xst**5&
     +xp(142)*xs1**1*xs2**6*xst**3&
     +xp(143)*xs1**1*xs2**8*xst**1&
     +xp(144)*xs1**2*xs2**0*xst**8&
     +xp(145)*xs1**2*xs2**2*xst**6&
     +xp(146)*xs1**2*xs2**4*xst**4&
     +xp(147)*xs1**2*xs2**6*xst**2&
     +xp(148)*xs1**2*xs2**8*xst**0&
     +xp(149)*xs1**3*xs2**0*xst**7&
     +xp(150)*xs1**3*xs2**2*xst**5&
     +xp(151)*xs1**3*xs2**4*xst**3&
     +xp(152)*xs1**3*xs2**6*xst**1&
     +xp(153)*xs1**4*xs2**0*xst**6&
     +xp(154)*xs1**4*xs2**2*xst**4&
     +xp(155)*xs1**4*xs2**4*xst**2&
     +xp(156)*xs1**4*xs2**6*xst**0&
     +xp(157)*xs1**5*xs2**0*xst**5&
     +xp(158)*xs1**5*xs2**2*xst**3&
     +xp(159)*xs1**5*xs2**4*xst**1&
     +xp(160)*xs1**6*xs2**0*xst**4&
     +xp(161)*xs1**6*xs2**2*xst**2&
     +xp(162)*xs1**6*xs2**4*xst**0&
     +xp(163)*xs1**7*xs2**0*xst**3&
     +xp(164)*xs1**7*xs2**2*xst**1&
     +xp(165)*xs1**8*xs2**0*xst**2&
     +xp(166)*xs1**8*xs2**2*xst**0&
     +xp(167)*xs1**9*xs2**0*xst**1&
     +xp(168)*xs1**10*xs2**0*xst**0&
     +xp(169)*xs1**0*xs2**0*xst**11&
     +xp(170)*xs1**0*xs2**2*xst**9&
     +xp(171)*xs1**0*xs2**4*xst**7&
     +xp(172)*xs1**0*xs2**6*xst**5&
     +xp(173)*xs1**0*xs2**8*xst**3&
     +xp(174)*xs1**0*xs2**10*xst**1&
     +xp(175)*xs1**1*xs2**0*xst**10&
     +xp(176)*xs1**1*xs2**2*xst**8&
     +xp(177)*xs1**1*xs2**4*xst**6&
     +xp(178)*xs1**1*xs2**6*xst**4&
     +xp(179)*xs1**1*xs2**8*xst**2&
     +xp(180)*xs1**1*xs2**10*xst**0&
     +xp(181)*xs1**2*xs2**0*xst**9
 vp3= xp(182)*xs1**2*xs2**2*xst**7&
     +xp(183)*xs1**2*xs2**4*xst**5&
     +xp(184)*xs1**2*xs2**6*xst**3&
     +xp(185)*xs1**2*xs2**8*xst**1&
     +xp(186)*xs1**3*xs2**0*xst**8&
     +xp(187)*xs1**3*xs2**2*xst**6&
     +xp(188)*xs1**3*xs2**4*xst**4&
     +xp(189)*xs1**3*xs2**6*xst**2&
     +xp(190)*xs1**3*xs2**8*xst**0&
     +xp(191)*xs1**4*xs2**0*xst**7&
     +xp(192)*xs1**4*xs2**2*xst**5&
     +xp(193)*xs1**4*xs2**4*xst**3&
     +xp(194)*xs1**4*xs2**6*xst**1&
     +xp(195)*xs1**5*xs2**0*xst**6&
     +xp(196)*xs1**5*xs2**2*xst**4&
     +xp(197)*xs1**5*xs2**4*xst**2&
     +xp(198)*xs1**5*xs2**6*xst**0&
     +xp(199)*xs1**6*xs2**0*xst**5&
     +xp(200)*xs1**6*xs2**2*xst**3&
     +xp(201)*xs1**6*xs2**4*xst**1&
     +xp(202)*xs1**7*xs2**0*xst**4&
     +xp(203)*xs1**7*xs2**2*xst**2&
     +xp(204)*xs1**7*xs2**4*xst**0&
     +xp(205)*xs1**8*xs2**0*xst**3&
     +xp(206)*xs1**8*xs2**2*xst**1&
     +xp(207)*xs1**9*xs2**0*xst**2&
     +xp(208)*xs1**9*xs2**2*xst**0&
     +xp(209)*xs1**0*xs2**0*xst**12&
     +xp(210)*xs1**0*xs2**2*xst**10&
     +xp(211)*xs1**0*xs2**4*xst**8&
     +xp(212)*xs1**0*xs2**6*xst**6&
     +xp(213)*xs1**0*xs2**8*xst**4&
     +xp(214)*xs1**0*xs2**10*xst**2&
     +xp(215)*xs1**0*xs2**12*xst**0&
     +xp(216)*xs1**1*xs2**0*xst**11&    
     +xp(217)*xs1**1*xs2**2*xst**9&
     +xp(218)*xs1**1*xs2**4*xst**7&
     +xp(219)*xs1**1*xs2**6*xst**5&
     +xp(220)*xs1**1*xs2**8*xst**3&
     +xp(221)*xs1**1*xs2**10*xst**1&
     +xp(222)*xs1**2*xs2**0*xst**10&
     +xp(223)*xs1**2*xs2**2*xst**8&
     +xp(224)*xs1**2*xs2**4*xst**6&
     +xp(225)*xs1**2*xs2**6*xst**4&
     +xp(226)*xs1**2*xs2**8*xst**2&
     +xp(227)*xs1**2*xs2**10*xst**0&
     +xp(228)*xs1**3*xs2**0*xst**9&
     +xp(229)*xs1**3*xs2**2*xst**7&
     +xp(230)*xs1**3*xs2**4*xst**5&
     +xp(231)*xs1**3*xs2**6*xst**3&
     +xp(232)*xs1**3*xs2**8*xst**1&
     +xp(233)*xs1**4*xs2**0*xst**8&
     +xp(234)*xs1**4*xs2**2*xst**6&
     +xp(235)*xs1**4*xs2**4*xst**4&
     +xp(236)*xs1**4*xs2**6*xst**2&
     +xp(237)*xs1**4*xs2**8*xst**0&
     +xp(238)*xs1**5*xs2**0*xst**7&
     +xp(239)*xs1**5*xs2**2*xst**5&
     +xp(240)*xs1**5*xs2**4*xst**3&
     +xp(241)*xs1**5*xs2**6*xst**1&
     +xp(242)*xs1**6*xs2**0*xst**6&
     +xp(243)*xs1**6*xs2**2*xst**4&
     +xp(244)*xs1**6*xs2**4*xst**2&
     +xp(245)*xs1**6*xs2**6*xst**0&
     +xp(246)*xs1**7*xs2**0*xst**5&
     +xp(247)*xs1**7*xs2**2*xst**3



       vp=vp1+vp2+vp3
       !
       vps1=42395.535333d0*xep1
       vps2=42395.535333d0*xep2

       y1=1.d0/(1.d0+dexp(ut*(x0-r1)))
       y2=1.d0/(1.d0+dexp(ut*(x0-r2)))
       y12=1.d0/(1.d0+dexp(ut2*(x02-r1)))
       y22=1.d0/(1.d0+dexp(ut2*(x02-r2)))

       !vp=vp*xep3*(1-y12)*(1-y22)
       !voh1=vpb1*(1-y1)+y1*vps1
       !voh2=vpb2*(1-y2)+y2*vps2

       vp=vp*xep3
       voh1=vpb1 !*(1-y1)+y1*vps1
       voh2=vpb2 !*(1-y2)+y2*vps2

       wifin_242_a=v0+vp+voh1+voh2+vhh

        return
        end function wifin_242_a




  double precision function taylor_I(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,pi,re12,alphae,y1,y2,y3,v0,v2,v3,v4
      double precision            :: frr   ,frr1  ,faa   ,frrr  ,frrr1 ,frra  ,frra1 ,fraa  ,& 
                                     faaa  ,frrrr ,frrrr1,frrrr2,frrra ,frrra1,frraa ,frraa1,fraaa ,faaaa 
      integer(4) :: i 
      !
      r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)

      pi = 4.0d0 * atan2(1.0d0,1.0d0)

      re12      = xp(2)
      alphae    = xp(3)*pi/180.0d0
      !
      v0 = xp(1)
      !
      frr   = xp( 4)
      frr1  = xp( 5)
      faa   = xp( 6)
      frrr  = xp( 7)
      frrr1 = xp( 8)
      frra  = xp( 9)
      frra1 = xp(10)
      fraa  = xp(11)
      faaa  = xp(12)
      frrrr = xp(13)
      frrrr1= xp(14)
      frrrr2= xp(15)
      frrra = xp(16)
      frrra1= xp(17)
      frraa = xp(18)
      frraa1= xp(19)
      fraaa = xp(20)
      faaaa = xp(21)

      !
      y1 = (r12-re12)
      y2 = (r32-re12)
      y3 = (alpha-alphae)
      !
      v2 = frr*(y1**2+y2**2) + frr1*y1*y2 + faa*y3**2
      v3 = frrr*(y1**3+y2**3)+frrr1*(y1**2*y2+y1*y2**2)+ & 
           faaa*y3**4 + frra*(y1**2+y2**2)*y3+fraa*(y1+y2)*y3**2+frra1*(y1*y2*y3)
      v4 = frrrr*(y1**4+y2**4)+frrrr1*(y1**3*y2+y1*y2**3)+frrrr2*(y1**2*y2**2)+ & 
           frrra*(y1**3+y2**3)*y3+frrra1*(y1**2*y2+y1*y2**2)*y3+frraa*(y1**2+y2**2)*y3**2+&
           frraa1*y1*y2*y3**2+fraaa*(y1+y2)*y3**3+faaaa*y3**4


       taylor_I=v0+v2+v3+v4

        return
        end function taylor_I







  double precision function morse_tyuterev(local,ZPE,npropin,yp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  yp(npropin)
      double precision            :: r12,r32,alpha,aa1,pi,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3
      double precision            :: g1,g2,b1,b2,rhh,vhh,v1,v2,v3,v4,v5,v6,v7,v8,th1
      integer(4) :: i 
      !
      r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)

      pi = 4.0d0 * atan2(1.0d0,1.0d0)

      re12      = yp(2)
      alphae    = yp(3)*pi/180.0d0
      !
      aa1  = yp(4)
      b1   = yp(5)
      b2   = yp(6)
      g1   = yp(7)
      g2   = yp(8)
      !
      rhh=dsqrt(r12**2+r32**2-2.d0*r12*r32*dcos(alpha))
      vhh=b1*dexp(-g1*rhh)+b2*dexp(-g2*rhh**2)
      !
! calculate potential energy function values
!
      xst=cos(alphae)-cos(alpha)
!     coro=sin(rhoe)-sin(alpha)
      !
      xst=cotan(alpha*0.5d0)
      !
      y1=1.0d+00-exp(-aa1*(r12-re12))
      y2=1.0d+00-exp(-aa1*(r32-re12))
      !
      y3=-cos(alphae)+cos(alpha)
      !
      y3=cotan(alpha*0.5d0)
      !
 v0 = yp(  1)*y1**0*y2**0*y3**0
 v1 = yp( 10)*y1**0*y2**0*y3**1& 
    + yp( 11)*y1**1*y2**0*y3**0& 
    + yp( 11)*y1**0*y2**1*y3**0
 v2 = yp( 12)*y1**0*y2**0*y3**2& 
    + yp( 13)*y1**1*y2**0*y3**1& 
    + yp( 13)*y1**0*y2**1*y3**1& 
    + yp( 14)*y1**1*y2**1*y3**0& 
    + yp( 15)*y1**2*y2**0*y3**0& 
    + yp( 15)*y1**0*y2**2*y3**0
    !
    y3=-cos(alphae)+cos(alpha)
    !
 v3 = yp( 16)*y1**0*y2**0*y3**3& 
    + yp( 17)*y1**1*y2**0*y3**2& 
    + yp( 17)*y1**0*y2**1*y3**2& 
    + yp( 18)*y1**1*y2**1*y3**1& 
    + yp( 19)*y1**2*y2**0*y3**1& 
    + yp( 19)*y1**0*y2**2*y3**1& 
    + yp( 20)*y1**2*y2**1*y3**0& 
    + yp( 20)*y1**1*y2**2*y3**0& 
    + yp( 21)*y1**3*y2**0*y3**0& 
    + yp( 21)*y1**0*y2**3*y3**0
 v4 = yp( 22)*y1**0*y2**0*y3**4& 
    + yp( 23)*y1**1*y2**0*y3**3& 
    + yp( 23)*y1**0*y2**1*y3**3& 
    + yp( 24)*y1**1*y2**1*y3**2& 
    + yp( 25)*y1**2*y2**0*y3**2& 
    + yp( 25)*y1**0*y2**2*y3**2& 
    + yp( 26)*y1**2*y2**1*y3**1& 
    + yp( 26)*y1**1*y2**2*y3**1& 
    + yp( 27)*y1**2*y2**2*y3**0& 
    + yp( 28)*y1**3*y2**0*y3**1& 
    + yp( 28)*y1**0*y2**3*y3**1& 
    + yp( 29)*y1**3*y2**1*y3**0& 
    + yp( 29)*y1**1*y2**3*y3**0& 
    + yp( 30)*y1**4*y2**0*y3**0& 
    + yp( 30)*y1**0*y2**4*y3**0
 v5 = yp( 31)*y1**0*y2**0*y3**5& 
    + yp( 32)*y1**1*y2**0*y3**4& 
    + yp( 32)*y1**0*y2**1*y3**4& 
    + yp( 33)*y1**1*y2**1*y3**3& 
    + yp( 34)*y1**2*y2**0*y3**3& 
    + yp( 34)*y1**0*y2**2*y3**3& 
    + yp( 35)*y1**2*y2**1*y3**2& 
    + yp( 35)*y1**1*y2**2*y3**2& 
    + yp( 36)*y1**2*y2**2*y3**1& 
    + yp( 37)*y1**3*y2**0*y3**2& 
    + yp( 37)*y1**0*y2**3*y3**2& 
    + yp( 38)*y1**3*y2**1*y3**1& 
    + yp( 38)*y1**1*y2**3*y3**1& 
    + yp( 39)*y1**3*y2**2*y3**0& 
    + yp( 39)*y1**2*y2**3*y3**0& 
    + yp( 40)*y1**4*y2**0*y3**1& 
    + yp( 40)*y1**0*y2**4*y3**1& 
    + yp( 41)*y1**4*y2**1*y3**0& 
    + yp( 41)*y1**1*y2**4*y3**0& 
    + yp( 42)*y1**5*y2**0*y3**0& 
    + yp( 42)*y1**0*y2**5*y3**0
 v6 = yp( 43)*y1**0*y2**0*y3**6& 
    + yp( 44)*y1**1*y2**0*y3**5& 
    + yp( 44)*y1**0*y2**1*y3**5& 
    + yp( 45)*y1**1*y2**1*y3**4& 
    + yp( 46)*y1**2*y2**0*y3**4& 
    + yp( 46)*y1**0*y2**2*y3**4& 
    + yp( 47)*y1**2*y2**1*y3**3& 
    + yp( 47)*y1**1*y2**2*y3**3& 
    + yp( 48)*y1**2*y2**2*y3**2& 
    + yp( 49)*y1**3*y2**0*y3**3& 
    + yp( 49)*y1**0*y2**3*y3**3& 
    + yp( 50)*y1**3*y2**1*y3**2& 
    + yp( 50)*y1**1*y2**3*y3**2& 
    + yp( 51)*y1**3*y2**2*y3**1& 
    + yp( 51)*y1**2*y2**3*y3**1& 
    + yp( 52)*y1**3*y2**3*y3**0& 
    + yp( 53)*y1**4*y2**0*y3**2& 
    + yp( 53)*y1**0*y2**4*y3**2& 
    + yp( 54)*y1**4*y2**1*y3**1& 
    + yp( 54)*y1**1*y2**4*y3**1& 
    + yp( 55)*y1**4*y2**2*y3**0& 
    + yp( 55)*y1**2*y2**4*y3**0& 
    + yp( 56)*y1**5*y2**0*y3**1& 
    + yp( 56)*y1**0*y2**5*y3**1& 
    + yp( 57)*y1**5*y2**1*y3**0& 
    + yp( 57)*y1**1*y2**5*y3**0& 
    + yp( 58)*y1**6*y2**0*y3**0& 
    + yp( 58)*y1**0*y2**6*y3**0
 v7 = yp( 59)*y1**0*y2**0*y3**7& 
    + yp( 60)*y1**1*y2**0*y3**6& 
    + yp( 60)*y1**0*y2**1*y3**6& 
    + yp( 61)*y1**1*y2**1*y3**5& 
    + yp( 62)*y1**2*y2**0*y3**5& 
    + yp( 62)*y1**0*y2**2*y3**5& 
    + yp( 63)*y1**2*y2**1*y3**4& 
    + yp( 63)*y1**1*y2**2*y3**4& 
    + yp( 64)*y1**2*y2**2*y3**3& 
    + yp( 65)*y1**3*y2**0*y3**4& 
    + yp( 65)*y1**0*y2**3*y3**4& 
    + yp( 66)*y1**3*y2**1*y3**3& 
    + yp( 66)*y1**1*y2**3*y3**3& 
    + yp( 67)*y1**3*y2**2*y3**2& 
    + yp( 67)*y1**2*y2**3*y3**2& 
    + yp( 68)*y1**3*y2**3*y3**1& 
    + yp( 69)*y1**4*y2**0*y3**3& 
    + yp( 69)*y1**0*y2**4*y3**3& 
    + yp( 70)*y1**4*y2**1*y3**2& 
    + yp( 70)*y1**1*y2**4*y3**2& 
    + yp( 71)*y1**4*y2**2*y3**1& 
    + yp( 71)*y1**2*y2**4*y3**1& 
    + yp( 72)*y1**4*y2**3*y3**0& 
    + yp( 72)*y1**3*y2**4*y3**0& 
    + yp( 73)*y1**5*y2**0*y3**2& 
    + yp( 73)*y1**0*y2**5*y3**2& 
    + yp( 74)*y1**5*y2**1*y3**1& 
    + yp( 74)*y1**1*y2**5*y3**1& 
    + yp( 75)*y1**5*y2**2*y3**0& 
    + yp( 75)*y1**2*y2**5*y3**0& 
    + yp( 76)*y1**6*y2**0*y3**1& 
    + yp( 76)*y1**0*y2**6*y3**1& 
    + yp( 77)*y1**6*y2**1*y3**0& 
    + yp( 77)*y1**1*y2**6*y3**0& 
    + yp( 78)*y1**7*y2**0*y3**0& 
    + yp( 78)*y1**0*y2**7*y3**0
 v8 = yp( 79)*y1**0*y2**0*y3**8& 
    + yp( 80)*y1**1*y2**0*y3**7& 
    + yp( 80)*y1**0*y2**1*y3**7& 
    + yp( 81)*y1**1*y2**1*y3**6& 
    + yp( 82)*y1**2*y2**0*y3**6& 
    + yp( 82)*y1**0*y2**2*y3**6& 
    + yp( 83)*y1**2*y2**1*y3**5& 
    + yp( 83)*y1**1*y2**2*y3**5& 
    + yp( 84)*y1**2*y2**2*y3**4& 
    + yp( 85)*y1**3*y2**0*y3**5& 
    + yp( 85)*y1**0*y2**3*y3**5& 
    + yp( 86)*y1**3*y2**1*y3**4& 
    + yp( 86)*y1**1*y2**3*y3**4& 
    + yp( 87)*y1**3*y2**2*y3**3& 
    + yp( 87)*y1**2*y2**3*y3**3& 
    + yp( 88)*y1**3*y2**3*y3**2& 
    + yp( 89)*y1**4*y2**0*y3**4& 
    + yp( 89)*y1**0*y2**4*y3**4& 
    + yp( 90)*y1**4*y2**1*y3**3& 
    + yp( 90)*y1**1*y2**4*y3**3& 
    + yp( 91)*y1**4*y2**2*y3**2& 
    + yp( 91)*y1**2*y2**4*y3**2& 
    + yp( 92)*y1**4*y2**3*y3**1& 
    + yp( 92)*y1**3*y2**4*y3**1& 
    + yp( 93)*y1**4*y2**4*y3**0& 
    + yp( 94)*y1**5*y2**0*y3**3& 
    + yp( 94)*y1**0*y2**5*y3**3& 
    + yp( 95)*y1**5*y2**1*y3**2& 
    + yp( 95)*y1**1*y2**5*y3**2& 
    + yp( 96)*y1**5*y2**2*y3**1& 
    + yp( 96)*y1**2*y2**5*y3**1& 
    + yp( 97)*y1**5*y2**3*y3**0& 
    + yp( 97)*y1**3*y2**5*y3**0& 
    + yp( 98)*y1**6*y2**0*y3**2& 
    + yp( 98)*y1**0*y2**6*y3**2& 
    + yp( 99)*y1**6*y2**1*y3**1& 
    + yp( 99)*y1**1*y2**6*y3**1& 
    + yp(100)*y1**6*y2**2*y3**0& 
    + yp(100)*y1**2*y2**6*y3**0& 
    + yp(101)*y1**7*y2**0*y3**1& 
    + yp(101)*y1**0*y2**7*y3**1& 
    + yp(102)*y1**7*y2**1*y3**0& 
    + yp(102)*y1**1*y2**7*y3**0& 
    + yp(103)*y1**8*y2**0*y3**0& 
    + yp(103)*y1**0*y2**8*y3**0


    !th1 = 0.5d0*( 1.0d0-tanh( 0.0005d0*( v0+v1+v2-40000.0d0 ) ) )


    !f=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh
    !
    !morse_tyuterev=(v0+v1+v2)+(v3+v4+v5+v6+v7+v8)*th1+vhh
 

     morse_tyuterev=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh

        return
        end function morse_tyuterev




  double precision function morse_242(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,aa1,pi,re12,alphae,xst,y1,y2,xs1,xs2,v0,vp1,vp2,vp3
      integer(4) :: i 
      !
       r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
       aa1  = 1.4d0

       pi = 4.0d0 * atan2(1.0d0,1.0d0)

       re12      = xp(2)
       alphae      = xp(3)*pi/180.0d0


! calculate potential energy function values
!
      xst=cos(alphae)-cos(alpha)
!     coro=sin(rhoe)-sin(alpha)

      y1=1.0d+00-exp(-aa1*(r12-re12))
      y2=1.0d+00-exp(-aa1*(r32-re12))
      !
      xs1 = (y1+y2)*0.5d0
      xs2 = (y1-y2)*0.5d0
      !



  v0= xp(1)  *xs1**0*xs2**0*xst**0
 vp1= xp(4)  *xs1**0*xs2**0*xst**1&
     +xp(5)  *xs1**1*xs2**0*xst**0&
     +xp(6)  *xs1**0*xs2**0*xst**2&
     +xp(7)  *xs1**0*xs2**2*xst**0&
     +xp(8)  *xs1**1*xs2**0*xst**1&
     +xp(9)  *xs1**2*xs2**0*xst**0&
     +xp(10) *xs1**0*xs2**0*xst**3&
     +xp(11) *xs1**0*xs2**2*xst**1&
     +xp(12) *xs1**1*xs2**0*xst**2&
     +xp(13) *xs1**1*xs2**2*xst**0&
     +xp(14) *xs1**2*xs2**0*xst**1&
     +xp(15) *xs1**3*xs2**0*xst**0&
     +xp(16) *xs1**0*xs2**0*xst**4&
     +xp(17) *xs1**0*xs2**2*xst**2&
     +xp(18) *xs1**0*xs2**4*xst**0&
     +xp(19) *xs1**1*xs2**0*xst**3&
     +xp(20) *xs1**1*xs2**2*xst**1&
     +xp(21) *xs1**2*xs2**0*xst**2&
     +xp(22) *xs1**2*xs2**2*xst**0&
     +xp(23) *xs1**3*xs2**0*xst**1&
     +xp(24) *xs1**4*xs2**0*xst**0&
     +xp(25) *xs1**0*xs2**0*xst**5&
     +xp(26) *xs1**0*xs2**2*xst**3&
     +xp(27) *xs1**0*xs2**4*xst**1&
     +xp(28) *xs1**1*xs2**0*xst**4&
     +xp(29) *xs1**1*xs2**2*xst**2&
     +xp(30) *xs1**1*xs2**4*xst**0&
     +xp(31) *xs1**2*xs2**0*xst**3&
     +xp(32) *xs1**2*xs2**2*xst**1&
     +xp(33) *xs1**3*xs2**0*xst**2&
     +xp(34) *xs1**3*xs2**2*xst**0&
     +xp(35) *xs1**4*xs2**0*xst**1&
     +xp(36) *xs1**5*xs2**0*xst**0&
     +xp(37) *xs1**0*xs2**0*xst**6&
     +xp(38) *xs1**0*xs2**2*xst**4&
     +xp(39) *xs1**0*xs2**4*xst**2&
     +xp(40) *xs1**0*xs2**6*xst**0&
     +xp(41) *xs1**1*xs2**0*xst**5&
     +xp(42) *xs1**1*xs2**2*xst**3&
     +xp(43) *xs1**1*xs2**4*xst**1&
     +xp(44) *xs1**2*xs2**0*xst**4&
     +xp(45) *xs1**2*xs2**2*xst**2&
     +xp(46) *xs1**2*xs2**4*xst**0&
     +xp(47) *xs1**3*xs2**0*xst**3&
     +xp(48) *xs1**3*xs2**2*xst**1&
     +xp(49) *xs1**4*xs2**0*xst**2&
     +xp(50) *xs1**4*xs2**2*xst**0&
     +xp(51) *xs1**5*xs2**0*xst**1&
     +xp(52) *xs1**6*xs2**0*xst**0&
     +xp(53) *xs1**0*xs2**0*xst**7&
     +xp(54) *xs1**0*xs2**2*xst**5&
     +xp(55) *xs1**0*xs2**4*xst**3&
     +xp(56) *xs1**0*xs2**6*xst**1&
     +xp(57) *xs1**1*xs2**0*xst**6&
     +xp(58) *xs1**1*xs2**2*xst**4&
     +xp(59) *xs1**1*xs2**4*xst**2&
     +xp(60) *xs1**1*xs2**6*xst**0&
     +xp(61) *xs1**2*xs2**0*xst**5&
     +xp(62) *xs1**2*xs2**2*xst**3&
     +xp(63) *xs1**2*xs2**4*xst**1&
     +xp(64) *xs1**3*xs2**0*xst**4&
     +xp(65) *xs1**3*xs2**2*xst**2&
     +xp(66) *xs1**3*xs2**4*xst**0&
     +xp(67) *xs1**4*xs2**0*xst**3&
     +xp(68) *xs1**4*xs2**2*xst**1&
     +xp(69) *xs1**5*xs2**0*xst**2&
     +xp(70) *xs1**5*xs2**2*xst**0&
     +xp(71) *xs1**6*xs2**0*xst**1&
     +xp(72) *xs1**7*xs2**0*xst**0&
     +xp(73) *xs1**0*xs2**0*xst**8&
     +xp(74) *xs1**0*xs2**2*xst**6&
     +xp(75) *xs1**0*xs2**4*xst**4&
     +xp(76) *xs1**0*xs2**6*xst**2&
     +xp(77) *xs1**0*xs2**8*xst**0&
     +xp(78) *xs1**1*xs2**0*xst**7&
     +xp(79) *xs1**1*xs2**2*xst**5&
     +xp(80) *xs1**1*xs2**4*xst**3&
     +xp(81) *xs1**1*xs2**6*xst**1&
     +xp(82) *xs1**2*xs2**0*xst**6&
     +xp(83) *xs1**2*xs2**2*xst**4&
     +xp(84) *xs1**2*xs2**4*xst**2&
     +xp(85) *xs1**2*xs2**6*xst**0&
     +xp(86) *xs1**3*xs2**0*xst**5&
     +xp(87) *xs1**3*xs2**2*xst**3&
     +xp(88) *xs1**3*xs2**4*xst**1&
     +xp(89) *xs1**4*xs2**0*xst**4&
     +xp(90) *xs1**4*xs2**2*xst**2&
     +xp(91) *xs1**4*xs2**4*xst**0&
     +xp(92) *xs1**5*xs2**0*xst**3&
     +xp(93) *xs1**5*xs2**2*xst**1&
     +xp(94) *xs1**6*xs2**0*xst**2
 vp2= xp(95) *xs1**6*xs2**2*xst**0&
     +xp(96) *xs1**7*xs2**0*xst**1&
     +xp(97) *xs1**8*xs2**0*xst**0&
     +xp(98) *xs1**0*xs2**0*xst**9&
     +xp(99) *xs1**0*xs2**2*xst**7&
     +xp(100)*xs1**0*xs2**4*xst**5&
     +xp(101)*xs1**0*xs2**6*xst**3&
     +xp(102)*xs1**0*xs2**8*xst**1&
     +xp(103)*xs1**1*xs2**0*xst**8&
     +xp(104)*xs1**1*xs2**2*xst**6&
     +xp(105)*xs1**1*xs2**4*xst**4&
     +xp(106)*xs1**1*xs2**6*xst**2&
     +xp(107)*xs1**1*xs2**8*xst**0&
     +xp(108)*xs1**2*xs2**0*xst**7&
     +xp(109)*xs1**2*xs2**2*xst**5&
     +xp(110)*xs1**2*xs2**4*xst**3&
     +xp(111)*xs1**2*xs2**6*xst**1&
     +xp(112)*xs1**3*xs2**0*xst**6&
     +xp(113)*xs1**3*xs2**2*xst**4&
     +xp(114)*xs1**3*xs2**4*xst**2&
     +xp(115)*xs1**3*xs2**6*xst**0&
     +xp(116)*xs1**4*xs2**0*xst**5&
     +xp(117)*xs1**4*xs2**2*xst**3&
     +xp(118)*xs1**4*xs2**4*xst**1&
     +xp(119)*xs1**5*xs2**0*xst**4&
     +xp(120)*xs1**5*xs2**2*xst**2&
     +xp(121)*xs1**5*xs2**4*xst**0&
     +xp(122)*xs1**6*xs2**0*xst**3&
     +xp(123)*xs1**6*xs2**2*xst**1&
     +xp(124)*xs1**7*xs2**0*xst**2&
     +xp(125)*xs1**7*xs2**2*xst**0&
     +xp(126)*xs1**8*xs2**0*xst**1&
     +xp(127)*xs1**9*xs2**0*xst**0&
     +xp(128)*xs1**0*xs2**0*xst**10&
     +xp(129)*xs1**0*xs2**2*xst**8&
     +xp(130)*xs1**0*xs2**4*xst**6&
     +xp(131)*xs1**0*xs2**6*xst**4&
     +xp(132)*xs1**0*xs2**8*xst**2&
     +xp(133)*xs1**0*xs2**10*xst**0&
     +xp(134)*xs1**1*xs2**0*xst**9&
     +xp(135)*xs1**1*xs2**2*xst**7&
     +xp(136)*xs1**1*xs2**4*xst**5&
     +xp(137)*xs1**1*xs2**6*xst**3&
     +xp(138)*xs1**1*xs2**8*xst**1&
     +xp(139)*xs1**2*xs2**0*xst**8&
     +xp(140)*xs1**2*xs2**2*xst**6&
     +xp(141)*xs1**2*xs2**4*xst**4&
     +xp(142)*xs1**2*xs2**6*xst**2&
     +xp(143)*xs1**2*xs2**8*xst**0&
     +xp(144)*xs1**3*xs2**0*xst**7&
     +xp(145)*xs1**3*xs2**2*xst**5&
     +xp(146)*xs1**3*xs2**4*xst**3&
     +xp(147)*xs1**3*xs2**6*xst**1&
     +xp(148)*xs1**4*xs2**0*xst**6&
     +xp(149)*xs1**4*xs2**2*xst**4&
     +xp(150)*xs1**4*xs2**4*xst**2&
     +xp(151)*xs1**4*xs2**6*xst**0&
     +xp(152)*xs1**5*xs2**0*xst**5&
     +xp(153)*xs1**5*xs2**2*xst**3&
     +xp(154)*xs1**5*xs2**4*xst**1&
     +xp(155)*xs1**6*xs2**0*xst**4&
     +xp(156)*xs1**6*xs2**2*xst**2&
     +xp(157)*xs1**6*xs2**4*xst**0&
     +xp(158)*xs1**7*xs2**0*xst**3&
     +xp(159)*xs1**7*xs2**2*xst**1&
     +xp(160)*xs1**8*xs2**0*xst**2&
     +xp(161)*xs1**8*xs2**2*xst**0&
     +xp(162)*xs1**9*xs2**0*xst**1&
     +xp(163)*xs1**10*xs2**0*xst**0&
     +xp(164)*xs1**0*xs2**0*xst**11&
     +xp(165)*xs1**0*xs2**2*xst**9&
     +xp(166)*xs1**0*xs2**4*xst**7&
     +xp(167)*xs1**0*xs2**6*xst**5&
     +xp(168)*xs1**0*xs2**8*xst**3&
     +xp(169)*xs1**0*xs2**10*xst**1&
     +xp(170)*xs1**1*xs2**0*xst**10&
     +xp(171)*xs1**1*xs2**2*xst**8&
     +xp(172)*xs1**1*xs2**4*xst**6&
     +xp(173)*xs1**1*xs2**6*xst**4&
     +xp(174)*xs1**1*xs2**8*xst**2&
     +xp(175)*xs1**1*xs2**10*xst**0&
     +xp(176)*xs1**2*xs2**0*xst**9
 vp3= xp(177)*xs1**2*xs2**2*xst**7&
     +xp(178)*xs1**2*xs2**4*xst**5&
     +xp(179)*xs1**2*xs2**6*xst**3&
     +xp(180)*xs1**2*xs2**8*xst**1&
     +xp(181)*xs1**3*xs2**0*xst**8&
     +xp(182)*xs1**3*xs2**2*xst**6&
     +xp(183)*xs1**3*xs2**4*xst**4&
     +xp(184)*xs1**3*xs2**6*xst**2&
     +xp(185)*xs1**3*xs2**8*xst**0&
     +xp(186)*xs1**4*xs2**0*xst**7&
     +xp(187)*xs1**4*xs2**2*xst**5&
     +xp(188)*xs1**4*xs2**4*xst**3&
     +xp(189)*xs1**4*xs2**6*xst**1&
     +xp(190)*xs1**5*xs2**0*xst**6&
     +xp(191)*xs1**5*xs2**2*xst**4&
     +xp(192)*xs1**5*xs2**4*xst**2&
     +xp(193)*xs1**5*xs2**6*xst**0&
     +xp(194)*xs1**6*xs2**0*xst**5&
     +xp(195)*xs1**6*xs2**2*xst**3&
     +xp(196)*xs1**6*xs2**4*xst**1&
     +xp(197)*xs1**7*xs2**0*xst**4&
     +xp(198)*xs1**7*xs2**2*xst**2&
     +xp(199)*xs1**7*xs2**4*xst**0&
     +xp(200)*xs1**8*xs2**0*xst**3&
     +xp(201)*xs1**8*xs2**2*xst**1&
     +xp(202)*xs1**9*xs2**0*xst**2&
     +xp(203)*xs1**9*xs2**2*xst**0&
     +xp(204)*xs1**0*xs2**0*xst**12&
     +xp(205)*xs1**0*xs2**2*xst**10&
     +xp(206)*xs1**0*xs2**4*xst**8&
     +xp(207)*xs1**0*xs2**6*xst**6&
     +xp(208)*xs1**0*xs2**8*xst**4&
     +xp(209)*xs1**0*xs2**10*xst**2&
     +xp(210)*xs1**0*xs2**12*xst**0&
     +xp(211)*xs1**1*xs2**0*xst**11&
     +xp(212)*xs1**1*xs2**2*xst**9&
     +xp(213)*xs1**1*xs2**4*xst**7&
     +xp(214)*xs1**1*xs2**6*xst**5&
     +xp(215)*xs1**1*xs2**8*xst**3&
     +xp(216)*xs1**1*xs2**10*xst**1&
     +xp(217)*xs1**2*xs2**0*xst**10&
     +xp(218)*xs1**2*xs2**2*xst**8&
     +xp(219)*xs1**2*xs2**4*xst**6&
     +xp(220)*xs1**2*xs2**6*xst**4&
     +xp(221)*xs1**2*xs2**8*xst**2&
     +xp(222)*xs1**2*xs2**10*xst**0&
     +xp(223)*xs1**3*xs2**0*xst**9&
     +xp(224)*xs1**3*xs2**2*xst**7&
     +xp(225)*xs1**3*xs2**4*xst**5&
     +xp(226)*xs1**3*xs2**6*xst**3&
     +xp(227)*xs1**3*xs2**8*xst**1&
     +xp(228)*xs1**4*xs2**0*xst**8&
     +xp(229)*xs1**4*xs2**2*xst**6&
     +xp(230)*xs1**4*xs2**4*xst**4&
     +xp(231)*xs1**4*xs2**6*xst**2&
     +xp(232)*xs1**4*xs2**8*xst**0&
     +xp(233)*xs1**5*xs2**0*xst**7&
     +xp(234)*xs1**5*xs2**2*xst**5&
     +xp(235)*xs1**5*xs2**4*xst**3&
     +xp(236)*xs1**5*xs2**6*xst**1&
     +xp(237)*xs1**6*xs2**0*xst**6&
     +xp(238)*xs1**6*xs2**2*xst**4&
     +xp(239)*xs1**6*xs2**4*xst**2&
     +xp(240)*xs1**6*xs2**6*xst**0&
     +xp(241)*xs1**7*xs2**0*xst**5&
     +xp(242)*xs1**7*xs2**2*xst**3




       morse_242=v0+vp1+vp2+vp3

        return
        end function morse_242



double precision function R_r12_morbid_switch(local,ZPE,npropin,xp)  result (v)

   double precision,intent(in) ::  local(3)
   double precision,intent(in) ::  ZPE
   integer,intent(in)          ::  npropin
   double precision,intent(in) ::  xp(npropin)
   double precision            ::  r12,r32,alpha,pi,r13,R,alpha1s,alpha2s,gamma1,gamma2,v0,v1,v2,v3
   !
   pi = 4.0d0 * atan2(1.0d0,1.0d0)
   !
   r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
   !
   r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
   !
   R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
   !
   v0       = xp(1)
   alpha1s  = xp(2)*pi/180.0d0
   alpha2s  = xp(3)*pi/180.0d0
   gamma1   = xp(4)
   gamma2   = xp(5)
   !
   v1 =  R_r12_dr_morse_pol(local,ZPE,253,xp(6:258))
   v2 =  morbid_mep2(local,ZPE,76,xp(259:334))
   v3 =  morbid_mep2(local,ZPE,npropin,xp(335:410))
   !
   v = v0
   v = v  +  v1*0.5d0*(1.0d0+tanh(-gamma1*(alpha-alpha1s)))
   v = v  +  v2*0.5d0*(1.0d0+tanh( gamma1*(alpha-alpha1s)))*0.5d0*(1.0d0+tanh(-gamma2*(alpha-alpha2s)))
   v = v  +  v3*0.5d0*(1.0d0+tanh( gamma2*(alpha-alpha2s)))
   !
end function R_r12_morbid_switch



  double precision function R_r12_dr_morse_pol(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,pi,r13e,Re,b1,b2,y,r13,R,a0,a1,a2,a3,a4,a5,a6,alphae,xs3,xs1,xs2,v0,vp1,vp2,vp3
      double precision            :: t1,t2,t3,t4,t5,t6,f1,f2,f3,f4,f5,f6,r13_inf,a13,r12e,De1,De2,y1,y2,vhh,bR1,bR2,bR3,coro,bR0,f0,frr
      integer(4) :: i 
       !
       !
       pi = 4.0d0 * atan2(1.0d0,1.0d0)
       !
       r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
       !
       r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
       !
	   R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
       !
       r12e   = xp(2)
       alphae = xp(3)*pi/180.0d0
       b1  = xp(4)

       bR0 = xp(5)
       bR1 = xp(6)
       bR2 = xp(7)
       bR3 = xp(8)
       !
       coro=cos(alpha)-cos(alphae)
       !
       b2 = bR0 + bR1*coro+bR2*coro**2+bR3*coro**3
       !
	   Re = 0.5d0*sqrt(r12e**2+r12e**2+2.0d0*r12e*r12e*cos(alphae))
       !
       a0  = xp( 9)
       a1  = xp(10)
       a2  = xp(11)
       a3  = xp(12)

       !f0  = xp(13)
       !f1  = xp(14)
       !f2  = xp(15)
       !f3  = xp(16)
       !
       !frr = f0 + f1*coro+f2*coro**2+f3*coro**3
       !
       !y = 1.0d0/(alpha+1e-14)
       !
       coro=cos(alpha)-cos(alphae)
       !
       f1 = coro**1
       f2 = coro**2
       f3 = coro**3
       f4 = coro**4
       f5 = coro**5
       f6 = coro**6
       !
       !f1 = exp(-a13*t1)
       !f2 = exp(-a13*t2)
       !f3 = exp(-a13*t3)
       !f4 = exp(-a13*t4)
       !f5 = exp(-a13*t5)
       !
       !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
       !
       r13e = a0+a1*f1+a2*f2+a3*f3 ! +a4*f4+a5*f5+a6*f6
       !
       !r12e = sqrt(0.25d0*r13e**2+R**2)
       !
       !r13e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
       !
       y1 = 1.0d+00-exp(-b1*(R-Re))
       y2 = 1.0d+00-exp(-b2*(r13-r13e))
       !
       !y1 = 1.0d+00-exp(-b1*((alpha)-alphae))
       !
       vhh=0   !De1*y1**2
       !
       xs1 = y1 ! (R-Re)
       xs3 = y2 ! (r13-r13e)
       xs2 = (r12-r32)*0.5d0
       !

  v0= xp( 1) *xs1**0*xs2**0*xs3**0  ! +frr*xs3**2
 vp1= xp( 15)*xs1**0*xs2**0*xs3**1&
     +xp( 16)*xs1**1*xs2**0*xs3**0&
     +xp( 17)*xs1**0*xs2**0*xs3**2&
     +xp( 18)*xs1**0*xs2**2*xs3**0&
     +xp( 19)*xs1**1*xs2**0*xs3**1&
     +xp( 20)*xs1**2*xs2**0*xs3**0&
     +xp( 21)*xs1**0*xs2**0*xs3**3&
     +xp( 22)*xs1**0*xs2**2*xs3**1&
     +xp( 23)*xs1**1*xs2**0*xs3**2&
     +xp( 24)*xs1**1*xs2**2*xs3**0&
     +xp( 25)*xs1**2*xs2**0*xs3**1&
     +xp( 26)*xs1**3*xs2**0*xs3**0&
     +xp( 27)*xs1**0*xs2**0*xs3**4&
     +xp( 28)*xs1**0*xs2**2*xs3**2&
     +xp( 29)*xs1**0*xs2**4*xs3**0&
     +xp( 30)*xs1**1*xs2**0*xs3**3&
     +xp( 31)*xs1**1*xs2**2*xs3**1&
     +xp( 32)*xs1**2*xs2**0*xs3**2&
     +xp( 33)*xs1**2*xs2**2*xs3**0&
     +xp( 34)*xs1**3*xs2**0*xs3**1&
     +xp( 35)*xs1**4*xs2**0*xs3**0&
     +xp( 36)*xs1**0*xs2**0*xs3**5&
     +xp( 37)*xs1**0*xs2**2*xs3**3&
     +xp( 38)*xs1**0*xs2**4*xs3**1&
     +xp( 39)*xs1**1*xs2**0*xs3**4&
     +xp( 40)*xs1**1*xs2**2*xs3**2&
     +xp( 41)*xs1**1*xs2**4*xs3**0&
     +xp( 42)*xs1**2*xs2**0*xs3**3&
     +xp( 43)*xs1**2*xs2**2*xs3**1&
     +xp( 44)*xs1**3*xs2**0*xs3**2&
     +xp( 45)*xs1**3*xs2**2*xs3**0&
     +xp( 46)*xs1**4*xs2**0*xs3**1&
     +xp( 47)*xs1**5*xs2**0*xs3**0&
     +xp( 48)*xs1**0*xs2**0*xs3**6&
     +xp( 49)*xs1**0*xs2**2*xs3**4&
     +xp( 50)*xs1**0*xs2**4*xs3**2&
     +xp( 51)*xs1**0*xs2**6*xs3**0&
     +xp( 52)*xs1**1*xs2**0*xs3**5&
     +xp( 53)*xs1**1*xs2**2*xs3**3&
     +xp( 54)*xs1**1*xs2**4*xs3**1&
     +xp( 55)*xs1**2*xs2**0*xs3**4&
     +xp( 56)*xs1**2*xs2**2*xs3**2&
     +xp( 57)*xs1**2*xs2**4*xs3**0&
     +xp( 58)*xs1**3*xs2**0*xs3**3&
     +xp( 59)*xs1**3*xs2**2*xs3**1&
     +xp( 60)*xs1**4*xs2**0*xs3**2&
     +xp( 61)*xs1**4*xs2**2*xs3**0&
     +xp( 62)*xs1**5*xs2**0*xs3**1&
     +xp( 63)*xs1**6*xs2**0*xs3**0&
     +xp( 64)*xs1**0*xs2**0*xs3**7&
     +xp( 65)*xs1**0*xs2**2*xs3**5&
     +xp( 66)*xs1**0*xs2**4*xs3**3&
     +xp( 67)*xs1**0*xs2**6*xs3**1&
     +xp( 68)*xs1**1*xs2**0*xs3**6&
     +xp( 69)*xs1**1*xs2**2*xs3**4&
     +xp( 70)*xs1**1*xs2**4*xs3**2&
     +xp( 71)*xs1**1*xs2**6*xs3**0&
     +xp( 72)*xs1**2*xs2**0*xs3**5&
     +xp( 73)*xs1**2*xs2**2*xs3**3&
     +xp( 74)*xs1**2*xs2**4*xs3**1&
     +xp( 75)*xs1**3*xs2**0*xs3**4&
     +xp( 76)*xs1**3*xs2**2*xs3**2&
     +xp( 77)*xs1**3*xs2**4*xs3**0&
     +xp( 78)*xs1**4*xs2**0*xs3**3&
     +xp( 79)*xs1**4*xs2**2*xs3**1&
     +xp( 80)*xs1**5*xs2**0*xs3**2&
     +xp( 81)*xs1**5*xs2**2*xs3**0&
     +xp( 82)*xs1**6*xs2**0*xs3**1&
     +xp( 83)*xs1**7*xs2**0*xs3**0&
     +xp( 84)*xs1**0*xs2**0*xs3**8&
     +xp( 85)*xs1**0*xs2**2*xs3**6&
     +xp( 86)*xs1**0*xs2**4*xs3**4&
     +xp( 87)*xs1**0*xs2**6*xs3**2&
     +xp( 88)*xs1**0*xs2**8*xs3**0&
     +xp( 89)*xs1**1*xs2**0*xs3**7&
     +xp( 90)*xs1**1*xs2**2*xs3**5&
     +xp( 91)*xs1**1*xs2**4*xs3**3&
     +xp( 92)*xs1**1*xs2**6*xs3**1&
     +xp( 93)*xs1**2*xs2**0*xs3**6&
     +xp( 94)*xs1**2*xs2**2*xs3**4&
     +xp( 95)*xs1**2*xs2**4*xs3**2&
     +xp( 96)*xs1**2*xs2**6*xs3**0&
     +xp( 97)*xs1**3*xs2**0*xs3**5&
     +xp( 98)*xs1**3*xs2**2*xs3**3&
     +xp( 99)*xs1**3*xs2**4*xs3**1&
     +xp(100)*xs1**4*xs2**0*xs3**4&
     +xp(101)*xs1**4*xs2**2*xs3**2&
     +xp(102)*xs1**4*xs2**4*xs3**0&
     +xp(103)*xs1**5*xs2**0*xs3**3&
     +xp(104)*xs1**5*xs2**2*xs3**1&
     +xp(105)*xs1**6*xs2**0*xs3**2
 vp2= xp(106)*xs1**6*xs2**2*xs3**0&
     +xp(107)*xs1**7*xs2**0*xs3**1&
     +xp(108)*xs1**8*xs2**0*xs3**0&
     +xp(109)*xs1**0*xs2**0*xs3**9&
     +xp(110)*xs1**0*xs2**2*xs3**7&
     +xp(111)*xs1**0*xs2**4*xs3**5&
     +xp(112)*xs1**0*xs2**6*xs3**3&
     +xp(113)*xs1**0*xs2**8*xs3**1&
     +xp(114)*xs1**1*xs2**0*xs3**8&
     +xp(115)*xs1**1*xs2**2*xs3**6&
     +xp(116)*xs1**1*xs2**4*xs3**4&
     +xp(117)*xs1**1*xs2**6*xs3**2&
     +xp(118)*xs1**1*xs2**8*xs3**0&
     +xp(119)*xs1**2*xs2**0*xs3**7&
     +xp(120)*xs1**2*xs2**2*xs3**5&
     +xp(121)*xs1**2*xs2**4*xs3**3&
     +xp(122)*xs1**2*xs2**6*xs3**1&
     +xp(123)*xs1**3*xs2**0*xs3**6&
     +xp(124)*xs1**3*xs2**2*xs3**4&
     +xp(125)*xs1**3*xs2**4*xs3**2&
     +xp(126)*xs1**3*xs2**6*xs3**0&
     +xp(127)*xs1**4*xs2**0*xs3**5&
     +xp(128)*xs1**4*xs2**2*xs3**3&
     +xp(129)*xs1**4*xs2**4*xs3**1&
     +xp(130)*xs1**5*xs2**0*xs3**4&
     +xp(131)*xs1**5*xs2**2*xs3**2&
     +xp(132)*xs1**5*xs2**4*xs3**0&
     +xp(133)*xs1**6*xs2**0*xs3**3&
     +xp(134)*xs1**6*xs2**2*xs3**1&
     +xp(135)*xs1**7*xs2**0*xs3**2&
     +xp(136)*xs1**7*xs2**2*xs3**0&
     +xp(137)*xs1**8*xs2**0*xs3**1&
     +xp(138)*xs1**9*xs2**0*xs3**0&
     +xp(139)*xs1**0*xs2**0*xs3**10&
     +xp(140)*xs1**0*xs2**2*xs3**8&
     +xp(141)*xs1**0*xs2**4*xs3**6&
     +xp(142)*xs1**0*xs2**6*xs3**4&
     +xp(143)*xs1**0*xs2**8*xs3**2&
     +xp(144)*xs1**0*xs2**10*xs3**0&
     +xp(145)*xs1**1*xs2**0*xs3**9&
     +xp(146)*xs1**1*xs2**2*xs3**7&
     +xp(147)*xs1**1*xs2**4*xs3**5&
     +xp(148)*xs1**1*xs2**6*xs3**3&
     +xp(149)*xs1**1*xs2**8*xs3**1&
     +xp(150)*xs1**2*xs2**0*xs3**8&
     +xp(151)*xs1**2*xs2**2*xs3**6&
     +xp(152)*xs1**2*xs2**4*xs3**4&
     +xp(153)*xs1**2*xs2**6*xs3**2&
     +xp(154)*xs1**2*xs2**8*xs3**0&
     +xp(155)*xs1**3*xs2**0*xs3**7&
     +xp(156)*xs1**3*xs2**2*xs3**5&
     +xp(157)*xs1**3*xs2**4*xs3**3&
     +xp(158)*xs1**3*xs2**6*xs3**1&
     +xp(159)*xs1**4*xs2**0*xs3**6&
     +xp(160)*xs1**4*xs2**2*xs3**4&
     +xp(161)*xs1**4*xs2**4*xs3**2&
     +xp(162)*xs1**4*xs2**6*xs3**0&
     +xp(163)*xs1**5*xs2**0*xs3**5&
     +xp(164)*xs1**5*xs2**2*xs3**3&
     +xp(165)*xs1**5*xs2**4*xs3**1&
     +xp(166)*xs1**6*xs2**0*xs3**4&
     +xp(167)*xs1**6*xs2**2*xs3**2&
     +xp(168)*xs1**6*xs2**4*xs3**0&
     +xp(169)*xs1**7*xs2**0*xs3**3&
     +xp(170)*xs1**7*xs2**2*xs3**1&
     +xp(171)*xs1**8*xs2**0*xs3**2&
     +xp(172)*xs1**8*xs2**2*xs3**0&
     +xp(173)*xs1**9*xs2**0*xs3**1&
     +xp(174)*xs1**10*xs2**0*xs3**0&
     +xp(175)*xs1**0*xs2**0*xs3**11&
     +xp(176)*xs1**0*xs2**2*xs3**9&
     +xp(177)*xs1**0*xs2**4*xs3**7&
     +xp(178)*xs1**0*xs2**6*xs3**5&
     +xp(179)*xs1**0*xs2**8*xs3**3&
     +xp(180)*xs1**0*xs2**10*xs3**1&
     +xp(181)*xs1**1*xs2**0*xs3**10&
     +xp(182)*xs1**1*xs2**2*xs3**8&
     +xp(183)*xs1**1*xs2**4*xs3**6&
     +xp(184)*xs1**1*xs2**6*xs3**4&
     +xp(185)*xs1**1*xs2**8*xs3**2&
     +xp(186)*xs1**1*xs2**10*xs3**0&
     +xp(187)*xs1**2*xs2**0*xs3**9
 vp3= xp(188)*xs1**2*xs2**2*xs3**7&
     +xp(189)*xs1**2*xs2**4*xs3**5&
     +xp(190)*xs1**2*xs2**6*xs3**3&
     +xp(191)*xs1**2*xs2**8*xs3**1&
     +xp(192)*xs1**3*xs2**0*xs3**8&
     +xp(193)*xs1**3*xs2**2*xs3**6&
     +xp(194)*xs1**3*xs2**4*xs3**4&
     +xp(195)*xs1**3*xs2**6*xs3**2&
     +xp(196)*xs1**3*xs2**8*xs3**0&
     +xp(197)*xs1**4*xs2**0*xs3**7&
     +xp(198)*xs1**4*xs2**2*xs3**5&
     +xp(199)*xs1**4*xs2**4*xs3**3&
     +xp(200)*xs1**4*xs2**6*xs3**1&
     +xp(201)*xs1**5*xs2**0*xs3**6&
     +xp(202)*xs1**5*xs2**2*xs3**4&
     +xp(203)*xs1**5*xs2**4*xs3**2&
     +xp(204)*xs1**5*xs2**6*xs3**0&
     +xp(205)*xs1**6*xs2**0*xs3**5&
     +xp(206)*xs1**6*xs2**2*xs3**3&
     +xp(207)*xs1**6*xs2**4*xs3**1&
     +xp(208)*xs1**7*xs2**0*xs3**4&
     +xp(209)*xs1**7*xs2**2*xs3**2&
     +xp(210)*xs1**7*xs2**4*xs3**0&
     +xp(211)*xs1**8*xs2**0*xs3**3&
     +xp(212)*xs1**8*xs2**2*xs3**1&
     +xp(213)*xs1**9*xs2**0*xs3**2&
     +xp(214)*xs1**9*xs2**2*xs3**0&
     +xp(215)*xs1**0*xs2**0*xs3**12&
     +xp(216)*xs1**0*xs2**2*xs3**10&
     +xp(217)*xs1**0*xs2**4*xs3**8&
     +xp(218)*xs1**0*xs2**6*xs3**6&
     +xp(219)*xs1**0*xs2**8*xs3**4&
     +xp(220)*xs1**0*xs2**10*xs3**2&
     +xp(221)*xs1**0*xs2**12*xs3**0&
     +xp(222)*xs1**1*xs2**0*xs3**11&
     +xp(223)*xs1**1*xs2**2*xs3**9&
     +xp(224)*xs1**1*xs2**4*xs3**7&
     +xp(225)*xs1**1*xs2**6*xs3**5&
     +xp(226)*xs1**1*xs2**8*xs3**3&
     +xp(227)*xs1**1*xs2**10*xs3**1&
     +xp(228)*xs1**2*xs2**0*xs3**10&
     +xp(229)*xs1**2*xs2**2*xs3**8&
     +xp(230)*xs1**2*xs2**4*xs3**6&
     +xp(231)*xs1**2*xs2**6*xs3**4&
     +xp(232)*xs1**2*xs2**8*xs3**2&
     +xp(233)*xs1**2*xs2**10*xs3**0&
     +xp(234)*xs1**3*xs2**0*xs3**9&
     +xp(235)*xs1**3*xs2**2*xs3**7&
     +xp(236)*xs1**3*xs2**4*xs3**5&
     +xp(237)*xs1**3*xs2**6*xs3**3&
     +xp(238)*xs1**3*xs2**8*xs3**1&
     +xp(239)*xs1**4*xs2**0*xs3**8&
     +xp(240)*xs1**4*xs2**2*xs3**6&
     +xp(241)*xs1**4*xs2**4*xs3**4&
     +xp(242)*xs1**4*xs2**6*xs3**2&
     +xp(243)*xs1**4*xs2**8*xs3**0&
     +xp(244)*xs1**5*xs2**0*xs3**7&
     +xp(245)*xs1**5*xs2**2*xs3**5&
     +xp(246)*xs1**5*xs2**4*xs3**3&
     +xp(247)*xs1**5*xs2**6*xs3**1&
     +xp(248)*xs1**6*xs2**0*xs3**6&
     +xp(249)*xs1**6*xs2**2*xs3**4&
     +xp(250)*xs1**6*xs2**4*xs3**2&
     +xp(251)*xs1**6*xs2**6*xs3**0&
     +xp(252)*xs1**7*xs2**0*xs3**5&
     +xp(253)*xs1**7*xs2**2*xs3**3




       R_r12_dr_morse_pol=(v0+vp1+vp2+vp3)+vhh

        return
        end function R_r12_dr_morse_pol




  double precision function R_r12_dr_morse_pol_x(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,pi,r13e,Re,b1,b2,y,r13,R,a0,a1,a2,a3,a4,a5,a6,alphae,xs3,xs1,xs2,v0,vp1,vp2,vp3
      double precision            :: t1,t2,t3,t4,t5,t6,f1,f2,f3,f4,f5,f6,r13_inf,a13,r12e,De1,De2,y1,y2,vhh,bR1,bR2,bR3,coro
      integer(4) :: i 
       !
       !
       pi = 4.0d0 * atan2(1.0d0,1.0d0)
       !
       r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
       !
       r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
       !
	   R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
       !
       r12e   = xp(2)
       alphae = xp(3)*pi/180.0d0
       !
       coro=cos(alpha)-cos(alphae)
       !
       b1  = xp(4)
       bR1 = xp(5)
       bR2 = xp(6)
       bR3 = xp(7)
       !
       b2 = bR1*coro**0+bR2*coro**1+bR3*coro**2
       !
	   Re = 0.5d0*sqrt(r12e**2+r12e**2+2.0d0*r12e*r12e*cos(alphae))
       !
       a0  = xp( 8)
       a1  = xp( 9)
       a2  = xp(10)
       a3  = xp(11)
       a4  = xp(12)
       a5  = xp(13)
       a6  = xp(14)
       !
       !y = 1.0d0/(alpha+1e-14)
       !
       f1 = coro**1
       f2 = coro**2
       f3 = coro**3
       f4 = coro**4
       f5 = coro**5
       f6 = coro**6
       !
       !f1 = exp(-a13*t1)
       !f2 = exp(-a13*t2)
       !f3 = exp(-a13*t3)
       !f4 = exp(-a13*t4)
       !f5 = exp(-a13*t5)
       !
       !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
       !
       r13e = a0+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5+a6*f6
       !
       !r12e = sqrt(0.25d0*r13e**2+R**2)
       !
       !r13e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
       !
       R = 0.5d0*r13e/tan(0.5d0*alpha)
       !
       y1 = 1.0d+00-exp(-b1*(R-Re))
       y2 = 1.0d+00-exp(-b2*(r13-r13e))
       !
       !y1 = 1.0d+00-exp(-b1*((alpha)-alphae))
       !
       vhh=0   !De1*y1**2
       !
       xs1 = y1 ! (R-Re)
       xs3 = y2 ! (r13-r13e)
       xs2 = (r12-r32)*0.5d0
       !
       xs1 = coro


  v0= xp( 1) *xs1**0*xs2**0*xs3**0
 vp1= xp( 15)*xs1**0*xs2**0*xs3**1&
     +xp( 16)*xs1**1*xs2**0*xs3**0&
     +xp( 17)*xs1**0*xs2**0*xs3**2&
     +xp( 18)*xs1**0*xs2**2*xs3**0&
     +xp( 19)*xs1**1*xs2**0*xs3**1&
     !+xp( 20)*xs1**2*xs2**0*xs3**0&
     +xp( 20)*y1**2*xs2**0*xs3**0&
     +xp( 21)*xs1**0*xs2**0*xs3**3&
     +xp( 22)*xs1**0*xs2**2*xs3**1&
     +xp( 23)*xs1**1*xs2**0*xs3**2&
     +xp( 24)*xs1**1*xs2**2*xs3**0&
     +xp( 25)*xs1**2*xs2**0*xs3**1&
     +xp( 26)*xs1**3*xs2**0*xs3**0&
     +xp( 27)*xs1**0*xs2**0*xs3**4&
     +xp( 28)*xs1**0*xs2**2*xs3**2&
     +xp( 29)*xs1**0*xs2**4*xs3**0&
     +xp( 30)*xs1**1*xs2**0*xs3**3&
     +xp( 31)*xs1**1*xs2**2*xs3**1&
     +xp( 32)*xs1**2*xs2**0*xs3**2&
     +xp( 33)*xs1**2*xs2**2*xs3**0&
     +xp( 34)*xs1**3*xs2**0*xs3**1&
     +xp( 35)*xs1**4*xs2**0*xs3**0&
     +xp( 36)*xs1**0*xs2**0*xs3**5&
     +xp( 37)*xs1**0*xs2**2*xs3**3&
     +xp( 38)*xs1**0*xs2**4*xs3**1&
     +xp( 39)*xs1**1*xs2**0*xs3**4&
     +xp( 40)*xs1**1*xs2**2*xs3**2&
     +xp( 41)*xs1**1*xs2**4*xs3**0&
     +xp( 42)*xs1**2*xs2**0*xs3**3&
     +xp( 43)*xs1**2*xs2**2*xs3**1&
     +xp( 44)*xs1**3*xs2**0*xs3**2&
     +xp( 45)*xs1**3*xs2**2*xs3**0&
     +xp( 46)*xs1**4*xs2**0*xs3**1&
     +xp( 47)*xs1**5*xs2**0*xs3**0&
     +xp( 48)*xs1**0*xs2**0*xs3**6&
     +xp( 49)*xs1**0*xs2**2*xs3**4&
     +xp( 50)*xs1**0*xs2**4*xs3**2&
     +xp( 51)*xs1**0*xs2**6*xs3**0&
     +xp( 52)*xs1**1*xs2**0*xs3**5&
     +xp( 53)*xs1**1*xs2**2*xs3**3&
     +xp( 54)*xs1**1*xs2**4*xs3**1&
     +xp( 55)*xs1**2*xs2**0*xs3**4&
     +xp( 56)*xs1**2*xs2**2*xs3**2&
     +xp( 57)*xs1**2*xs2**4*xs3**0&
     +xp( 58)*xs1**3*xs2**0*xs3**3&
     +xp( 59)*xs1**3*xs2**2*xs3**1&
     +xp( 60)*xs1**4*xs2**0*xs3**2&
     +xp( 61)*xs1**4*xs2**2*xs3**0&
     +xp( 62)*xs1**5*xs2**0*xs3**1&
     +xp( 63)*xs1**6*xs2**0*xs3**0&
     +xp( 64)*xs1**0*xs2**0*xs3**7&
     +xp( 65)*xs1**0*xs2**2*xs3**5&
     +xp( 66)*xs1**0*xs2**4*xs3**3&
     +xp( 67)*xs1**0*xs2**6*xs3**1&
     +xp( 68)*xs1**1*xs2**0*xs3**6&
     +xp( 69)*xs1**1*xs2**2*xs3**4&
     +xp( 70)*xs1**1*xs2**4*xs3**2&
     +xp( 71)*xs1**1*xs2**6*xs3**0&
     +xp( 72)*xs1**2*xs2**0*xs3**5&
     +xp( 73)*xs1**2*xs2**2*xs3**3&
     +xp( 74)*xs1**2*xs2**4*xs3**1&
     +xp( 75)*xs1**3*xs2**0*xs3**4&
     +xp( 76)*xs1**3*xs2**2*xs3**2&
     +xp( 77)*xs1**3*xs2**4*xs3**0&
     +xp( 78)*xs1**4*xs2**0*xs3**3&
     +xp( 79)*xs1**4*xs2**2*xs3**1&
     +xp( 80)*xs1**5*xs2**0*xs3**2&
     +xp( 81)*xs1**5*xs2**2*xs3**0&
     +xp( 82)*xs1**6*xs2**0*xs3**1&
     +xp( 83)*xs1**7*xs2**0*xs3**0&
     +xp( 84)*xs1**0*xs2**0*xs3**8&
     +xp( 85)*xs1**0*xs2**2*xs3**6&
     +xp( 86)*xs1**0*xs2**4*xs3**4&
     +xp( 87)*xs1**0*xs2**6*xs3**2&
     +xp( 88)*xs1**0*xs2**8*xs3**0&
     +xp( 89)*xs1**1*xs2**0*xs3**7&
     +xp( 90)*xs1**1*xs2**2*xs3**5&
     +xp( 91)*xs1**1*xs2**4*xs3**3&
     +xp( 92)*xs1**1*xs2**6*xs3**1&
     +xp( 93)*xs1**2*xs2**0*xs3**6&
     +xp( 94)*xs1**2*xs2**2*xs3**4&
     +xp( 95)*xs1**2*xs2**4*xs3**2&
     +xp( 96)*xs1**2*xs2**6*xs3**0&
     +xp( 97)*xs1**3*xs2**0*xs3**5&
     +xp( 98)*xs1**3*xs2**2*xs3**3&
     +xp( 99)*xs1**3*xs2**4*xs3**1&
     +xp(100)*xs1**4*xs2**0*xs3**4&
     +xp(101)*xs1**4*xs2**2*xs3**2&
     +xp(102)*xs1**4*xs2**4*xs3**0&
     +xp(103)*xs1**5*xs2**0*xs3**3&
     +xp(104)*xs1**5*xs2**2*xs3**1&
     +xp(105)*xs1**6*xs2**0*xs3**2
 vp2= xp(106)*xs1**6*xs2**2*xs3**0&
     +xp(107)*xs1**7*xs2**0*xs3**1&
     +xp(108)*xs1**8*xs2**0*xs3**0&
     +xp(109)*xs1**0*xs2**0*xs3**9&
     +xp(110)*xs1**0*xs2**2*xs3**7&
     +xp(111)*xs1**0*xs2**4*xs3**5&
     +xp(112)*xs1**0*xs2**6*xs3**3&
     +xp(113)*xs1**0*xs2**8*xs3**1&
     +xp(114)*xs1**1*xs2**0*xs3**8&
     +xp(115)*xs1**1*xs2**2*xs3**6&
     +xp(116)*xs1**1*xs2**4*xs3**4&
     +xp(117)*xs1**1*xs2**6*xs3**2&
     +xp(118)*xs1**1*xs2**8*xs3**0&
     +xp(119)*xs1**2*xs2**0*xs3**7&
     +xp(120)*xs1**2*xs2**2*xs3**5&
     +xp(121)*xs1**2*xs2**4*xs3**3&
     +xp(122)*xs1**2*xs2**6*xs3**1&
     +xp(123)*xs1**3*xs2**0*xs3**6&
     +xp(124)*xs1**3*xs2**2*xs3**4&
     +xp(125)*xs1**3*xs2**4*xs3**2&
     +xp(126)*xs1**3*xs2**6*xs3**0&
     +xp(127)*xs1**4*xs2**0*xs3**5&
     +xp(128)*xs1**4*xs2**2*xs3**3&
     +xp(129)*xs1**4*xs2**4*xs3**1&
     +xp(130)*xs1**5*xs2**0*xs3**4&
     +xp(131)*xs1**5*xs2**2*xs3**2&
     +xp(132)*xs1**5*xs2**4*xs3**0&
     +xp(133)*xs1**6*xs2**0*xs3**3&
     +xp(134)*xs1**6*xs2**2*xs3**1&
     +xp(135)*xs1**7*xs2**0*xs3**2&
     +xp(136)*xs1**7*xs2**2*xs3**0&
     +xp(137)*xs1**8*xs2**0*xs3**1&
     +xp(138)*xs1**9*xs2**0*xs3**0&
     +xp(139)*xs1**0*xs2**0*xs3**10&
     +xp(140)*xs1**0*xs2**2*xs3**8&
     +xp(141)*xs1**0*xs2**4*xs3**6&
     +xp(142)*xs1**0*xs2**6*xs3**4&
     +xp(143)*xs1**0*xs2**8*xs3**2&
     +xp(144)*xs1**0*xs2**10*xs3**0&
     +xp(145)*xs1**1*xs2**0*xs3**9&
     +xp(146)*xs1**1*xs2**2*xs3**7&
     +xp(147)*xs1**1*xs2**4*xs3**5&
     +xp(148)*xs1**1*xs2**6*xs3**3&
     +xp(149)*xs1**1*xs2**8*xs3**1&
     +xp(150)*xs1**2*xs2**0*xs3**8&
     +xp(151)*xs1**2*xs2**2*xs3**6&
     +xp(152)*xs1**2*xs2**4*xs3**4&
     +xp(153)*xs1**2*xs2**6*xs3**2&
     +xp(154)*xs1**2*xs2**8*xs3**0&
     +xp(155)*xs1**3*xs2**0*xs3**7&
     +xp(156)*xs1**3*xs2**2*xs3**5&
     +xp(157)*xs1**3*xs2**4*xs3**3&
     +xp(158)*xs1**3*xs2**6*xs3**1&
     +xp(159)*xs1**4*xs2**0*xs3**6&
     +xp(160)*xs1**4*xs2**2*xs3**4&
     +xp(161)*xs1**4*xs2**4*xs3**2&
     +xp(162)*xs1**4*xs2**6*xs3**0&
     +xp(163)*xs1**5*xs2**0*xs3**5&
     +xp(164)*xs1**5*xs2**2*xs3**3&
     +xp(165)*xs1**5*xs2**4*xs3**1&
     +xp(166)*xs1**6*xs2**0*xs3**4&
     +xp(167)*xs1**6*xs2**2*xs3**2&
     +xp(168)*xs1**6*xs2**4*xs3**0&
     +xp(169)*xs1**7*xs2**0*xs3**3&
     +xp(170)*xs1**7*xs2**2*xs3**1&
     +xp(171)*xs1**8*xs2**0*xs3**2&
     +xp(172)*xs1**8*xs2**2*xs3**0&
     +xp(173)*xs1**9*xs2**0*xs3**1&
     +xp(174)*xs1**10*xs2**0*xs3**0&
     +xp(175)*xs1**0*xs2**0*xs3**11&
     +xp(176)*xs1**0*xs2**2*xs3**9&
     +xp(177)*xs1**0*xs2**4*xs3**7&
     +xp(178)*xs1**0*xs2**6*xs3**5&
     +xp(179)*xs1**0*xs2**8*xs3**3&
     +xp(180)*xs1**0*xs2**10*xs3**1&
     +xp(181)*xs1**1*xs2**0*xs3**10&
     +xp(182)*xs1**1*xs2**2*xs3**8&
     +xp(183)*xs1**1*xs2**4*xs3**6&
     +xp(184)*xs1**1*xs2**6*xs3**4&
     +xp(185)*xs1**1*xs2**8*xs3**2&
     +xp(186)*xs1**1*xs2**10*xs3**0&
     +xp(187)*xs1**2*xs2**0*xs3**9
 vp3= xp(188)*xs1**2*xs2**2*xs3**7&
     +xp(189)*xs1**2*xs2**4*xs3**5&
     +xp(190)*xs1**2*xs2**6*xs3**3&
     +xp(191)*xs1**2*xs2**8*xs3**1&
     +xp(192)*xs1**3*xs2**0*xs3**8&
     +xp(193)*xs1**3*xs2**2*xs3**6&
     +xp(194)*xs1**3*xs2**4*xs3**4&
     +xp(195)*xs1**3*xs2**6*xs3**2&
     +xp(196)*xs1**3*xs2**8*xs3**0&
     +xp(197)*xs1**4*xs2**0*xs3**7&
     +xp(198)*xs1**4*xs2**2*xs3**5&
     +xp(199)*xs1**4*xs2**4*xs3**3&
     +xp(200)*xs1**4*xs2**6*xs3**1&
     +xp(201)*xs1**5*xs2**0*xs3**6&
     +xp(202)*xs1**5*xs2**2*xs3**4&
     +xp(203)*xs1**5*xs2**4*xs3**2&
     +xp(204)*xs1**5*xs2**6*xs3**0&
     +xp(205)*xs1**6*xs2**0*xs3**5&
     +xp(206)*xs1**6*xs2**2*xs3**3&
     +xp(207)*xs1**6*xs2**4*xs3**1&
     +xp(208)*xs1**7*xs2**0*xs3**4&
     +xp(209)*xs1**7*xs2**2*xs3**2&
     +xp(210)*xs1**7*xs2**4*xs3**0&
     +xp(211)*xs1**8*xs2**0*xs3**3&
     +xp(212)*xs1**8*xs2**2*xs3**1&
     +xp(213)*xs1**9*xs2**0*xs3**2&
     +xp(214)*xs1**9*xs2**2*xs3**0&
     +xp(215)*xs1**0*xs2**0*xs3**12&
     +xp(216)*xs1**0*xs2**2*xs3**10&
     +xp(217)*xs1**0*xs2**4*xs3**8&
     +xp(218)*xs1**0*xs2**6*xs3**6&
     +xp(219)*xs1**0*xs2**8*xs3**4&
     +xp(220)*xs1**0*xs2**10*xs3**2&
     +xp(221)*xs1**0*xs2**12*xs3**0&
     +xp(222)*xs1**1*xs2**0*xs3**11&
     +xp(223)*xs1**1*xs2**2*xs3**9&
     +xp(224)*xs1**1*xs2**4*xs3**7&
     +xp(225)*xs1**1*xs2**6*xs3**5&
     +xp(226)*xs1**1*xs2**8*xs3**3&
     +xp(227)*xs1**1*xs2**10*xs3**1&
     +xp(228)*xs1**2*xs2**0*xs3**10&
     +xp(229)*xs1**2*xs2**2*xs3**8&
     +xp(230)*xs1**2*xs2**4*xs3**6&
     +xp(231)*xs1**2*xs2**6*xs3**4&
     +xp(232)*xs1**2*xs2**8*xs3**2&
     +xp(233)*xs1**2*xs2**10*xs3**0&
     +xp(234)*xs1**3*xs2**0*xs3**9&
     +xp(235)*xs1**3*xs2**2*xs3**7&
     +xp(236)*xs1**3*xs2**4*xs3**5&
     +xp(237)*xs1**3*xs2**6*xs3**3&
     +xp(238)*xs1**3*xs2**8*xs3**1&
     +xp(239)*xs1**4*xs2**0*xs3**8&
     +xp(240)*xs1**4*xs2**2*xs3**6&
     +xp(241)*xs1**4*xs2**4*xs3**4&
     +xp(242)*xs1**4*xs2**6*xs3**2&
     +xp(243)*xs1**4*xs2**8*xs3**0&
     +xp(244)*xs1**5*xs2**0*xs3**7&
     +xp(245)*xs1**5*xs2**2*xs3**5&
     +xp(246)*xs1**5*xs2**4*xs3**3&
     +xp(247)*xs1**5*xs2**6*xs3**1&
     +xp(248)*xs1**6*xs2**0*xs3**6&
     +xp(249)*xs1**6*xs2**2*xs3**4&
     +xp(250)*xs1**6*xs2**4*xs3**2&
     +xp(251)*xs1**6*xs2**6*xs3**0&
     +xp(252)*xs1**7*xs2**0*xs3**5&
     +xp(253)*xs1**7*xs2**2*xs3**3




       R_r12_dr_morse_pol_x=(v0+vp1+vp2+vp3)+vhh

        return
        end function R_r12_dr_morse_pol_x



  double precision function V_rep_disp(local,ZPE,npropin,xp)

      implicit real*8 (a-h,o-z)

      double precision,intent(in) ::  local(3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          ::  npropin
      double precision,intent(in) ::  xp(npropin)
      double precision            :: r12,r32,alpha,pi,r13e,Re,b1,b2,y,r13,R,a0,a1,a2,a3,a4,a5,a6,alphae,xs3,xs1,xs2,v0,vp1,vp2,vp3
      double precision            :: t1,t2,t3,t4,t5,t6,f1,f2,f3,f4,f5,f6,r13_inf,a13,r12e,De1,De2,y1,y2,vhh,bR1,bR2,bR3,coro,bR0,f0,frr
      double precision            :: A,C6,C7,C8,beta,Vrep,Vdisp,DHH,v2
      integer(4) :: i 
       !
       !
       pi = 4.0d0 * atan2(1.0d0,1.0d0)
       !
       r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
       !
       r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
       !
	   R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
       !
       r12e   = xp(2)
       alphae = xp(3)*pi/180.0d0
       b1  = xp(4)

       bR0 = xp(5)
       bR1 = xp(6)
       bR2 = xp(7)
       bR3 = xp(8)
       !
       coro=cos(alpha)-cos(alphae)
       !
       b2 = bR0 + bR1*coro+bR2*coro**2+bR3*coro**3
       !
	   Re = 0.5d0*sqrt(r12e**2+r12e**2+2.0d0*r12e*r12e*cos(alphae))
       !
       a0  = xp( 9)
       a1  = xp(10)
       a2  = xp(11)
       a3  = xp(12)

       !f0  = xp(13)
       !f1  = xp(14)
       !f2  = xp(15)
       !f3  = xp(16)
       !
       !frr = f0 + f1*coro+f2*coro**2+f3*coro**3
       !
       !y = 1.0d0/(alpha+1e-14)
       !
       coro=cos(alpha)-cos(alphae)
       !
       f1 = coro**1
       f2 = coro**2
       f3 = coro**3
       f4 = coro**4
       f5 = coro**5
       f6 = coro**6
       !
       !f1 = exp(-a13*t1)
       !f2 = exp(-a13*t2)
       !f3 = exp(-a13*t3)
       !f4 = exp(-a13*t4)
       !f5 = exp(-a13*t5)
       !
       !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
       !
       r13e = a0+a1*f1+a2*f2+a3*f3 ! +a4*f4+a5*f5+a6*f6
       !
       !r12e = sqrt(0.25d0*r13e**2+R**2)
       !
       !r13e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
       !
       y1 = 1.0d+00-exp(-b1*(R-Re))
       y2 = 1.0d+00-exp(-b2*(r13-r13e))
       !
       !y1 = 1.0d+00-exp(-b1*((alpha)-alphae))
       !
       vhh=0   !De1*y1**2
       !
       xs1 = y1 ! (R-Re)
       xs3 = y2 ! (r13-r13e)
       xs2 = (r12-r32)*0.5d0
       !
       A = xp(13)
       !
       beta = xp(14)
       C6  = xp(15)
       C7  = xp(16)
       C8  = xp(17)
       !
       Vrep =A*exp(-beta*(R-Re))
       !
       Vdisp = -( C6*D_toennies(R,beta,6)/R**6+C7*D_toennies(R,beta,7)/R**7+C8*D_toennies(R,beta,8)/R**8 )
       !
       v0 = xp(1)
       !
       DHH = xp(2)
       !
       v2 = DHH*y2**2
       !
       V_rep_disp=v0+v2+Vdisp+Vrep

        return
        end function V_rep_disp


 double precision function D_toennies(R,beta,n)
   !
   double precision, intent(in) :: R,beta
   integer, intent(in) :: n
   double precision :: n_,temp
   integer :: k
   !
   temp = 0
   !
   n_ = 1.d0
   !
   do k = 0,n
     !
     n_ = n_*real(max(k,1),8)
     !
     temp = temp + (beta*R)**k/n_
     !
   enddo
   !
   D_toennies = 1.0d0-exp(-beta*R)*temp
   !
 end function D_toennies


   function fakt(a) result (f)

      real(8),intent(in) :: a
      real(8)            :: ax,f
      integer         :: i,ic
!
      ax=a
      f=1.0d0
      if(abs(ax)<1.d-24) return
      f=.1d0
      if(ax.lt.0.d0) then 
         write (*,"(1h0,' fkt.err  negative argument for functi on fakt. argument = ',e12.5)") ax
         stop 'fkt.err  negative argument'
      endif 
      !
      ic=idnint(ax)
      ax=ax/10.0d0
      f=ax
      do  i=1,ic-1
        f=f*(ax-real(i,8)*0.1d0)
      enddo

    end function fakt


 double precision function morse(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 
        integer            :: i,i1,i2,i3
        double precision   :: ve    ,re12  ,aa1, alpha

        double precision   :: pi,y1,y2,y3,v0,coro,alphae,r12,r32

        double precision   :: f(0:12,0:10,0:10)
        !
        pi = 4.0d0 * atan2(1.0d0,1.0d0)

        VE        = param( 1)
        RE12      = param( 2)
        alphae      = param( 3)*pi/180.0d0

        AA1       = param(4)
        !
        F = 0 

        i = 4

        F( 1:12,0,0)  = param( i+1:i+12) ; i = i + 12 

        F( 1:11,0,1)  = param( i+1:i+11) 
        F( 1:11,1,0)  = param( i+1:i+11) ; i = i + 11

        F( 0:10,0,2)  = param( i+1:i+11) 
        F( 0:10,2,0)  = param( i+1:i+11) ; i = i + 11
        F( 0:10,1,1)  = param( i+1:i+11) ; i = i + 11

        F( 0: 9,0,3)  = param( i+1:i+10) 
        F( 0: 9,3,0)  = param( i+1:i+10) ; i = i + 10
        F( 0: 9,1,2)  = param( i+1:i+10) 
        F( 0: 9,2,1)  = param( i+1:i+10) ; i = i + 10

        F( 0: 8,0,4)  = param( i+1:i+ 9) 
        F( 0: 8,4,0)  = param( i+1:i+ 9) ; i = i + 9 
        F( 0: 8,1,3)  = param( i+1:i+ 9) 
        F( 0: 8,3,1)  = param( i+1:i+ 9) ; i = i + 9 
        F( 0: 8,2,2)  = param( i+1:i+ 9) ; i = i + 9 

        F( 0: 7,0,5)  = param( i+1:i+ 8) 
        F( 0: 7,5,0)  = param( i+1:i+ 8) ; i = i + 8 
        F( 0: 7,1,4)  = param( i+1:i+ 8) 
        F( 0: 7,4,1)  = param( i+1:i+ 8) ; i = i + 8 
        F( 0: 7,2,3)  = param( i+1:i+ 8) 
        F( 0: 7,3,2)  = param( i+1:i+ 8) ; i = i + 8 

        F( 0: 6,0,6)  = param( i+1:i+ 7) 
        F( 0: 6,6,0)  = param( i+1:i+ 7) ; i = i + 7 
        F( 0: 6,1,5)  = param( i+1:i+ 7) 
        F( 0: 6,5,1)  = param( i+1:i+ 7) ; i = i + 7 
        F( 0: 6,2,4)  = param( i+1:i+ 7) 
        F( 0: 6,4,2)  = param( i+1:i+ 7) ; i = i + 7 
        F( 0: 6,3,3)  = param( i+1:i+ 7) ; i = i + 7 

        F( 0: 5,0,7)  = param( i+1:i+ 6) 
        F( 0: 5,7,0)  = param( i+1:i+ 6) ; i = i + 6 
        F( 0: 5,1,6)  = param( i+1:i+ 6) 
        F( 0: 5,6,1)  = param( i+1:i+ 6) ; i = i + 6 
        F( 0: 5,2,5)  = param( i+1:i+ 6) 
        F( 0: 5,5,2)  = param( i+1:i+ 6) ; i = i + 6 
        F( 0: 5,3,4)  = param( i+1:i+ 6) 
        F( 0: 5,4,3)  = param( i+1:i+ 6) ; i = i + 6 

        F( 0: 4,0,8)  = param( i+1:i+ 5) 
        F( 0: 4,8,0)  = param( i+1:i+ 5) ; i = i + 5 
        F( 0: 4,1,7)  = param( i+1:i+ 5) 
        F( 0: 4,7,1)  = param( i+1:i+ 5) ; i = i + 5 
        F( 0: 4,2,6)  = param( i+1:i+ 5) 
        F( 0: 4,6,2)  = param( i+1:i+ 5) ; i = i + 5 
        F( 0: 4,3,5)  = param( i+1:i+ 5) 
        F( 0: 4,5,3)  = param( i+1:i+ 5) ; i = i + 5 
        F( 0: 4,4,4)  = param( i+1:i+ 5) ; i = i + 5 

        F( 0: 3,0,9)  = param( i+1:i+ 4) 
        F( 0: 3,9,0)  = param( i+1:i+ 4) ; i = i + 4 
        F( 0: 3,1,8)  = param( i+1:i+ 4) 
        F( 0: 3,8,1)  = param( i+1:i+ 4) ; i = i + 4 
        F( 0: 3,2,7)  = param( i+1:i+ 4) 
        F( 0: 3,7,2)  = param( i+1:i+ 4) ; i = i + 4
        F( 0: 3,3,6)  = param( i+1:i+ 4) 
        F( 0: 3,6,3)  = param( i+1:i+ 4) ; i = i + 4 
        F( 0: 3,4,5)  = param( i+1:i+ 4) 
        F( 0: 3,5,4)  = param( i+1:i+ 4) ; i = i + 4

        F( 0:2, 0,10) = param( i+1:i+ 3) 
        F( 0:2,10, 0) = param( i+1:i+ 3) ; i = i + 3
        F( 0:2, 1, 9) = param( i+1:i+ 3) 
        F( 0:2, 9, 1) = param( i+1:i+ 3) ; i = i + 3 
        F( 0:2, 2, 8) = param( i+1:i+ 3) 
        F( 0:2, 8, 2) = param( i+1:i+ 3) ; i = i + 3 
        F( 0:2, 3 ,7) = param( i+1:i+ 3) 
        F( 0:2, 7, 3) = param( i+1:i+ 3) ; i = i + 3 
        F( 0:2, 4, 6) = param( i+1:i+ 3) 
        F( 0:2, 6, 4) = param( i+1:i+ 3) ; i = i + 3 
        F( 0:2, 5, 5) = param( i+1:i+ 3) ; i = i + 3 


!
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)


! calculate potential energy function values
!
      y1=cos(alpha)-cos(alphae)
!     coro=sin(rhoe)-sin(alpha)

      y2=1.0d+00-exp(-aa1*(r12-re12))
      y3=1.0d+00-exp(-aa1*(r32-re12))

!
! calculate potential energy function values
!     
      v0 = 0 
      do i1 = 0,12
        do i2 = 0,10
          do i3 = 0,i2
             v0 = v0 + 0.5d0*f(i1,i2,i3)*y1**i1*(y2**i2*y3**i3+y2**i3*y3**i2)
             !
             continue
             !
         enddo
        enddo
      enddo




       morse = ve+v0


 end function morse




 double precision function morbid(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 

      double precision   ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,  &
            f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,    & 
            fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
            fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
            f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
            f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133,r12,r32,alpha


      double precision  :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,   &    
            fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,    &
            f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,   &
            f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,       & 
            re32, aa3  ,f3,f33  ,f333 ,f133 ,f3333,f1333      
      double precision ::  &
      fe11111,f11111,fe11113,f11113,fe11133,f11133,fe11333,f11333,&
      fe13333,f13333,fe33333,f33333,fe111111,f111111,fe111113,f111113,&
      fe111133,f111133,fe111333,f111333,fe113333,f113333,fe133333,f133333,&
      fe333333,f333333
      double precision ::  &
      f2a1111,f2a3333,f2a1113,f2a1333,f2a1133,f1a11111,f2a11111,&
      f1a11113,f2a11113,f1a11133,f2a11133,&
      f1a11333,f2a11333,f1a13333,f2a13333,&
      f1a33333,f2a33333,f1a111111,f1a111113,f1a111133,f1a111333,f1a113333,&
      f1a133333,f1a333333      

      double precision ::  &
      f1111111,f1111113,f1111133,f11111111,f11111113,&
      f11111133,f11111333,f11113333,f1111333,f3333333,f1333333,f1133333,f1113333,&
      f33333333,f13333333,f11333333,f11133333

      double precision ::  &
      fe1111111,fe1111113,fe1111133,fe1111333,fe1113333,fe1133333,fe1333333,fe3333333,&
      fe11111111,fe11111113,fe11111133,fe11111333,fe11113333,fe11133333,fe13333333,fe33333333,fe11333333


        double precision   :: pi,y1,y3,v,coro,v0,RHOE,v5,v6,v7,v8



        pi = 4.0d0 * atan2(1.0d0,1.0d0)

        VE        = param( 1)
        RHOE      = param( 2)*pi/180.0d0
        FA1       = param( 3)
        FA2       = param( 4)
        FA3       = param( 5)
        FA4       = param( 6)
        FA5       = param( 7)
        FA6       = param( 8)
        FA7       = param( 9)
        FA8       = param(10)
        RE12      = param(11)
        AA1       = param(12)
        F1A1      = param(13)
        F2A1      = param(14)
        F3A1      = param(15)
        F4A1      = param(16)
        F11       = param(17)
        F1A11     = param(18)
        F2A11     = param(19)
        F3A11     = param(20)
        F13       = param(21)
        F1A13     = param(22)
        F2A13     = param(23)
        F3A13     = param(24)
        F111      = param(25)
        F1A111    = param(26)
        F2A111    = param(27)
        F113      = param(28)
        F1A113    = param(29)
        F2A113    = param(30)
        F1111     = param(31)
        F1A1111   = param(32)
        F1113     = param(33)
        F1A1113   = param(34)
        F1133     = param(35)
        F1A1133   = param(36)


      f11111   = param(37)
      f1a11111 = param(38)
      f2a11111 = param(39)

      f11113   = param(40)
      f1a11113 = param(41)
      f2a11113 = param(42)
      f11133   = param(43)
      f1a11133 = param(44)
      f2a11133 = param(45)


      f111111  = param(46)
      f1a111111= param(47)

      f111113  = param(48)
      f1a111113= param(49)
      f111133  = param(50)
      f1a111133= param(51)
      f111333  = param(52)
      f1a111333= param(53)



        f1 = 0.0d0
        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f1a1333 = f1a1113

        f33333= f11111
        f13333= f11113
        f11333= f11133
        f1a33333= f1a11111
        f1a13333= f1a11113
        f1a11333= f1a11133
        f2a33333= f2a11111
        f2a13333= f2a11113
        f2a11333= f2a11133

        f333333= f111111
        f133333= f111113
        f113333= f111133
        f1a333333= f1a111111
        f1a133333= f1a111113
        f1a113333= f1a111133



!
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)


! calculate potential energy function values
!
       coro=cos(rhoe)+cos(alpha)
       !
       !coro=alpha-(pi-rhoe)
       !
       !coro=sin(rhoe)-sin(alpha)

      y1=1.0d+00-exp(-aa1*(r12-re12))
      y3=1.0d+00-exp(-aa1*(r32-re12))

!
! calculate potential energy function values
!

      v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 +fa6*coro**6+fa7*coro**7+fa8*coro**8


      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro !!!!!+f2a1111*coro**2
      fe3333= f3333+f1a3333*coro !!!!!+f2a3333*coro**2
      fe1113= f1113+f1a1113*coro !!!!!+f2a1113*coro**2
      fe1333= f1333+f1a1333*coro !!!!!+f2a1333*coro**2
      fe1133= f1133+f1a1133*coro !!!!!+f2a1133*coro**2

      fe11111= f11111+f1a11111*coro+f2a11111*coro**2
      fe11113= f11113+f1a11113*coro+f2a11113*coro**2
      fe11133= f11133+f1a11133*coro+f2a11133*coro**2
      fe11333= f11333+f1a11333*coro+f2a11333*coro**2
      fe13333= f13333+f1a13333*coro+f2a13333*coro**2
      fe33333= f33333+f1a33333*coro+f2a33333*coro**2

      fe111111= f111111+f1a111111*coro
      fe111113= f111113+f1a111113*coro
      fe111133= f111133+f1a111133*coro
      fe111333= f111333+f1a111333*coro
      fe113333= f113333+f1a113333*coro
      fe133333= f133333+f1a133333*coro
      fe333333= f333333+f1a333333*coro

      !fe1111111= f1111111 !+f1a1111111*coro
      !fe1111113= f1111113 !+f1a1111113*coro
      !fe1111133= f1111133 !+f1a1111133*coro
      !fe1111333= f1111333 !+f1a1111333*coro
      !fe1113333= f1113333 !+f1a1113333*coro
      !fe1133333= f1133333 !+f1a1133333*coro
      !fe1333333= f1333333 !+f1a1333333*coro
      !fe3333333= f3333333 !+f1a3333333*coro

      !fe11111111= f11111111 !+f1a11111111*coro
      !fe11111113= f11111113 !+f1a11111113*coro
      !fe11111133= f11111133 !+f1a11111133*coro
      !fe11111333= f11111333 !+f1a11111333*coro
      !fe11113333= f11113333 !+f1a11113333*coro
      !fe11133333= f11133333 !!+f1a11333333*coro
      !fe11333333= f11333333 !+f1a13333333*coro
      !fe13333333= f13333333 !+f1a13333333*coro
      !fe33333333= f33333333 !+f1a33333333*coro



      v     =  v0+fe1*y1+fe3*y3                          &  
              +fe11*y1**2+fe33*y3**2+fe13*y1*y3          & 
              +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3    &
              +fe133*y1*y3**2                            & 
              +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3 & 
              +fe1333*y1*y3**3+fe1133*y1**2*y3**2

      v5 = fe11111*y1**5+fe11113*y1**4*y3**1+fe11133*y1**3*y3**2+fe11333*y1**2*y3**3+fe13333*y1**1*y3**4+fe33333*y3**5
      v6 = fe111111*y1**6+fe111113*y1**5*y3**1+fe111133*y1**4*y3**2+fe111333*y1**3*y3**3+fe113333*y1**2*y3**4+fe133333*y1**1*y3**5+fe333333*y3**6
      !v7 = fe1111111*y1**7+fe1111113*y1**6*y3**1+fe1111133*y1**5*y3**2+fe1111333*y1**4*y3**3+&
      !     fe1113333*y1**3*y3**4+fe1133333*y1**2*y3**5+fe1333333*y1**1*y3**6+fe3333333*y1**0*y3**7
      !v8 = fe11111111*y1**8+fe11111113*y1**7*y3**1+fe11111133*y1**6*y3**2+fe11111333*y1**5*y3**3+&
      !     fe11113333*y1**4*y3**4+fe11133333*y1**3*y3**5+fe11333333*y1**2*y3**6+fe13333333*y1**1*y3**7+fe33333333*y1**0*y3**8

       morbid = v+v5+v6


 end function morbid





 double precision function morbid_MEP(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 

      double precision   ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,  &
            f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,    & 
            fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
            fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
            f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
            f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133,r12,r32,alpha


      double precision  :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,   &    
            fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,    &
            f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,   &
            f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,       & 
            re32, aa3  ,f3,f33  ,f333 ,f133 ,f3333,f1333      
      double precision ::  &
      fe11111,f11111,fe11113,f11113,fe11133,f11133,fe11333,f11333,&
      fe13333,f13333,fe33333,f33333,fe111111,f111111,fe111113,f111113,&
      fe111133,f111133,fe111333,f111333,fe113333,f113333,fe133333,f133333,&
      fe333333,f333333
      double precision ::  &
      f2a1111,f2a3333,f2a1113,f2a1333,f2a1133,f1a11111,f2a11111,&
      f1a11113,f2a11113,f1a11133,f2a11133,&
      f1a11333,f2a11333,f1a13333,f2a13333,&
      f1a33333,f2a33333,f1a111111,f1a111113,f1a111133,f1a111333,f1a113333,&
      f1a133333,f1a333333      

      double precision ::  &
      f1111111,f1111113,f1111133,f11111111,f11111113,&
      f11111133,f11111333,f11113333,f1111333,f3333333,f1333333,f1133333,f1113333,&
      f33333333,f13333333,f11333333,f11133333

      double precision ::  &
      fe1111111,fe1111113,fe1111133,fe1111333,fe1113333,fe1133333,fe1333333,fe3333333,&
      fe11111111,fe11111113,fe11111133,fe11111333,fe11113333,fe11133333,fe13333333,fe33333333,fe11333333,&
      F4A11,f4a33


        double precision   :: pi,y1,y3,v,coro,v0,RHOE,v5,v6,v7,v8,alphae,r13e,r12e,r13_0,r13_inf,a13,r,r13,re
        double precision   :: a0,a1,a2,a3,a4,a5,a6,y,vhh,fhh,ahh,g1,g2,g3,g4,g5,g6



        !
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
        !
        r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
        !
	    R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
        !
        pi = 4.0d0 * atan2(1.0d0,1.0d0)
        !
        VE        = param( 1)
        !
        a0      = param(2)
        a1      = param(3)
        a2      = param(4)
        a3      = param(5)
        a4      = param(6)
        a5      = param(7)
        a6      = param(8)
        !
        !a0      = param(2)
        !a13     = param(3)
        !a1      = param(4)
        !a2      = param(5)
        !a3      = param(6)
        !a4      = param(7)
        !a5      = param(8)
        !
        !r13e = a0+a1*exp(-a13*R**2)+a2*exp(-a13*R**4)+a3*exp(-a13*R**6)+a4*exp(-a13*R**8)+a5*exp(-a13*R**10)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        g1 = alpha**1
        g2 = alpha**2
        g3 = alpha**3
        g4 = alpha**4
        g5 = alpha**5
        g6 = alpha**6
        !
        !f1 = exp(-a13*t1)
        !f2 = exp(-a13*t2)
        !f3 = exp(-a13*t3)
        !f4 = exp(-a13*t4)
        !f5 = exp(-a13*t5)
        !
        !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
        !
        r13e = a0+a1*g1+a2*g2+a3*g3+a4*g4+a5*g5+a6*g6
        !
        !r13e = r13_inf+(r13_0-r13_inf)*exp(-a13*sqrt(R)**9)
        !
        r12e = 0.5d0*r13e/sin(alpha*0.5d0)  !  sqrt(0.25d0*r13e**2+R**2)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        !r12e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
        !r13e = 2.0d0*r12e*sin(0.5d0*alpha)
        !
        !alphae = atan(0.5d0*r13e/R)*2.d0
        !
        !
        !RE12      = param(14)
        !RHOE      = param( 5)*pi/180.0d0
        !
        !Re = 0.5d0*sqrt(2.0d0*r12e**2-2.0d0*re12*re12*cos(rhoe))
        !
        RHOE      = param( 9)*pi/180.0d0
        RE12      = param(18)
        A0        = param(19)
        a1        = param(20)
        a2        = param(21)
        a3        = param(22)
        !
        aa1 = a1*alpha+a2*alpha**2+a3*alpha**3
        aa1 = a0
        !
        coro=cos(alpha)+cos(rhoe)
        !
        !coro=sin(2.0d0*alpha) ! +cos(2.0d0*rhoe)
        !
        !coro=1.0d0-exp(-a0*(alpha-(pi-rhoe)))
        !
        !coro=alpha-(pi-rhoe)
        !
        !coro=sin(alpha)+sin(rhoe)
        !
        !coro = 1.0d+00-exp(-1.0d0*(pi-rhoe-alpha))
        !
        !a13 = a1
        !
        !coro = 1.0d+00-exp(-a13*(r13-r13e))
        !
        !
        !aa1 = a0
        !
        y1=1.0d+00-exp(-aa1*(r12-re12))
        y3=1.0d+00-exp(-aa1*(r32-re12))
        !
        ahh = param(67)
        fhh = param(68)
        !
        vhh = fhh*exp(-ahh*(pi-alpha)) !  0.9e5*exp(-10d0*r13) !+0.518622e5*y1**2+0.518622e5*y3**2
        !
        !y1=1.0d+00-exp(-aa1*(r12-r12e))
        !y3=1.0d+00-exp(-aa1*(r32-r12e))
        FA1       = param(10)
        FA2       = param(11)
        FA3       = param(12)
        FA4       = param(13)
        FA5       = param(14)
        FA6       = param(15)
        FA7       = param(16)
        FA8       = param(17)
        F1        = param(23)
        F1A1      = param(24)
        F2A1      = param(25)
        F3A1      = param(26)
        F4A1      = param(27)
        F11       = param(28)
        F1A11     = param(29)
        F2A11     = param(30)
        F3A11     = param(31)
        F4A11     = param(32)
        F13       = param(33)
        F1A13     = param(34)
        F2A13     = param(35)
        F3A13     = param(36)
        F111      = param(37)
        F1A111    = param(38)
        F2A111    = param(39)
        F113      = param(40)
        F1A113    = param(41)
        F2A113    = param(42)
        F1111     = param(43)
        F1A1111   = param(44)
        F2A1111   = param(45)
        F1113     = param(46)
        F1A1113   = param(47)
        F1133     = param(48)
        F1A1133   = param(49)
        f11111    = param(50)
        f1a11111  = param(51)
        f2a11111  = param(52)
        f11113    = param(53)
        f1a11113  = param(54)
        f2a11113  = param(55)
        f11133    = param(56)
        f1a11133  = param(57)
        f2a11133  = param(58)
        f111111   = param(59)
        f1a111111 = param(60)
        f111113   = param(61)
        f1a111113 = param(62)
        f111133   = param(63)
        f1a111133 = param(64)
        f111333   = param(65)
        f1a111333 = param(66)


        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f4a33   = f4a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f2a3333 = f2a1111
        f1a1333 = f1a1113

        f33333= f11111
        f13333= f11113
        f11333= f11133
        f1a33333= f1a11111
        f1a13333= f1a11113
        f1a11333= f1a11133
        f2a33333= f2a11111
        f2a13333= f2a11113
        f2a11333= f2a11133

        f333333= f111111
        f133333= f111113
        f113333= f111133
        f1a333333= f1a111111
        f1a133333= f1a111113
        f1a113333= f1a111133

!
! calculate potential energy function values
!

      v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 +fa6*coro**6+fa7*coro**7+fa8*coro**8
      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3+f4a33*coro**4
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro+f2a1111*coro**2
      fe3333= f3333+f1a3333*coro+f2a3333*coro**2
      fe1113= f1113+f1a1113*coro !!!!!+f2a1113*coro**2
      fe1333= f1333+f1a1333*coro !!!!!+f2a1333*coro**2
      fe1133= f1133+f1a1133*coro !!!!!+f2a1133*coro**2

      fe11111= f11111+f1a11111*coro+f2a11111*coro**2
      fe11113= f11113+f1a11113*coro+f2a11113*coro**2
      fe11133= f11133+f1a11133*coro+f2a11133*coro**2
      fe11333= f11333+f1a11333*coro+f2a11333*coro**2
      fe13333= f13333+f1a13333*coro+f2a13333*coro**2
      fe33333= f33333+f1a33333*coro+f2a33333*coro**2

      fe111111= f111111+f1a111111*coro
      fe111113= f111113+f1a111113*coro
      fe111133= f111133+f1a111133*coro
      fe111333= f111333+f1a111333*coro
      fe113333= f113333+f1a113333*coro
      fe133333= f133333+f1a133333*coro
      fe333333= f333333+f1a333333*coro

      !fe1111111= f1111111 !+f1a1111111*coro
      !fe1111113= f1111113 !+f1a1111113*coro
      !fe1111133= f1111133 !+f1a1111133*coro
      !fe1111333= f1111333 !+f1a1111333*coro
      !fe1113333= f1113333 !+f1a1113333*coro
      !fe1133333= f1133333 !+f1a1133333*coro
      !fe1333333= f1333333 !+f1a1333333*coro
      !fe3333333= f3333333 !+f1a3333333*coro

      !fe11111111= f11111111 !+f1a11111111*coro
      !fe11111113= f11111113 !+f1a11111113*coro
      !fe11111133= f11111133 !+f1a11111133*coro
      !fe11111333= f11111333 !+f1a11111333*coro
      !fe11113333= f11113333 !+f1a11113333*coro
      !fe11133333= f11133333 !!+f1a11333333*coro
      !fe11333333= f11333333 !+f1a13333333*coro
      !fe13333333= f13333333 !+f1a13333333*coro
      !fe33333333= f33333333 !+f1a33333333*coro



      v     =  v0+fe1*y1+fe3*y3                          &  
              +fe11*y1**2+fe33*y3**2+fe13*y1*y3          & 
              +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3    &
              +fe133*y1*y3**2                            & 
              +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3 & 
              +fe1333*y1*y3**3+fe1133*y1**2*y3**2

      v5 = fe11111*y1**5+fe11113*y1**4*y3**1+fe11133*y1**3*y3**2+fe11333*y1**2*y3**3+fe13333*y1**1*y3**4+fe33333*y3**5
      v6 = fe111111*y1**6+fe111113*y1**5*y3**1+fe111133*y1**4*y3**2+fe111333*y1**3*y3**3+fe113333*y1**2*y3**4+fe133333*y1**1*y3**5+fe333333*y3**6
      !v7 = fe1111111*y1**7+fe1111113*y1**6*y3**1+fe1111133*y1**5*y3**2+fe1111333*y1**4*y3**3+&
      !     fe1113333*y1**3*y3**4+fe1133333*y1**2*y3**5+fe1333333*y1**1*y3**6+fe3333333*y1**0*y3**7
      !v8 = fe11111111*y1**8+fe11111113*y1**7*y3**1+fe11111133*y1**6*y3**2+fe11111333*y1**5*y3**3+&
      !     fe11113333*y1**4*y3**4+fe11133333*y1**3*y3**5+fe11333333*y1**2*y3**6+fe13333333*y1**1*y3**7+fe33333333*y1**0*y3**8

       morbid_mep = v+v5+v6+vhh


 end function morbid_MEP





 double precision function morbid_MEP2(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 

      double precision   ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,  &
            f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,    & 
            fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
            fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
            f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
            f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133,r12,r32,alpha


      double precision  :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,   &    
            fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,    &
            f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,   &
            f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,       & 
            re32, aa3  ,f3,f33  ,f333 ,f133 ,f3333,f1333      
      double precision ::  &
      fe11111,f11111,fe11113,f11113,fe11133,f11133,fe11333,f11333,&
      fe13333,f13333,fe33333,f33333,fe111111,f111111,fe111113,f111113,&
      fe111133,f111133,fe111333,f111333,fe113333,f113333,fe133333,f133333,&
      fe333333,f333333
      double precision ::  &
      f2a1111,f2a3333,f2a1113,f2a1333,f2a1133,f1a11111,f2a11111,&
      f1a11113,f2a11113,f1a11133,f2a11133,&
      f1a11333,f2a11333,f1a13333,f2a13333,&
      f1a33333,f2a33333,f1a111111,f1a111113,f1a111133,f1a111333,f1a113333,&
      f1a133333,f1a333333      

      double precision ::  &
      f1111111,f1111113,f1111133,f11111111,f11111113,&
      f11111133,f11111333,f11113333,f1111333,f3333333,f1333333,f1133333,f1113333,&
      f33333333,f13333333,f11333333,f11133333

      double precision ::  &
      fe1111111,fe1111113,fe1111133,fe1111333,fe1113333,fe1133333,fe1333333,fe3333333,&
      fe11111111,fe11111113,fe11111133,fe11111333,fe11113333,fe11133333,fe13333333,fe33333333,fe11333333,&
      F4A11,f4a33,F5A11,F6A11,F5A33,F6A33,F4A13,F5A13,F6A13


        double precision   :: pi,y1,y3,v,coro,v0,RHOE,v5,v6,v7,v8,alphae,r13e,r12e,r13_0,r13_inf,a13,r,r13,re
        double precision   :: a0,a1,a2,a3,a4,a5,a6,y,vhh,fhh,ahh,g1,g2,g3,g4,g5,g6



        !
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
        !
        r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
        !
	    R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
        !
        pi = 4.0d0 * atan2(1.0d0,1.0d0)
        !
        VE        = param( 1)
        RHOE      = param( 9)*pi/180.0d0
        !
        a0      = param(2)
        a1      = param(3)
        a2      = param(4)
        a3      = param(5)
        a4      = param(6)
        a5      = param(7)
        a6      = param(8)
        !
        !a0      = param(2)
        !a13     = param(3)
        !a1      = param(4)
        !a2      = param(5)
        !a3      = param(6)
        !a4      = param(7)
        !a5      = param(8)
        !
        !r13e = a0+a1*exp(-a13*R**2)+a2*exp(-a13*R**4)+a3*exp(-a13*R**6)+a4*exp(-a13*R**8)+a5*exp(-a13*R**10)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        coro=cos(alpha)+cos(rhoe)
        !
        g1 = coro**1
        g2 = coro**2
        g3 = coro**3
        g4 = coro**4
        g5 = coro**5
        g6 = coro**6
        !
        !f1 = exp(-a13*t1)
        !f2 = exp(-a13*t2)
        !f3 = exp(-a13*t3)
        !f4 = exp(-a13*t4)
        !f5 = exp(-a13*t5)
        !
        !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
        !
        r13e = a0+a1*g1+a2*g2+a3*g3+a4*g4+a5*g5+a6*g6
        !
        !r13e = r13_inf+(r13_0-r13_inf)*exp(-a13*sqrt(R)**9)
        !
        r12e = 0.5d0*r13e/sin(alpha*0.5d0)  !  sqrt(0.25d0*r13e**2+R**2)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        !r12e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
        !r13e = 2.0d0*r12e*sin(0.5d0*alpha)
        !
        !alphae = atan(0.5d0*r13e/R)*2.d0
        !
        !
        !RE12      = param(14)
        !RHOE      = param( 5)*pi/180.0d0
        !
        !Re = 0.5d0*sqrt(2.0d0*r12e**2-2.0d0*re12*re12*cos(rhoe))
        !
        RE12      = param(18)
        A0        = param(19)
        a1        = param(20)
        a2        = param(21)
        a3        = param(22)
        a4        = param(23)
        a5        = param(24)
        a6        = param(25)
        !
        aa1 = a0+a1*coro+a2*coro**2+a3*coro**3+a4*coro**4+a5*coro**5+a6*coro**6
        !aa1 = a0
        !
        coro=cos(alpha)+cos(rhoe)
        !
        !coro=sin(2.0d0*alpha) ! +cos(2.0d0*rhoe)
        !
        !coro=1.0d0-exp(-a0*(alpha-(pi-rhoe)))
        !
        !coro=alpha-(pi-rhoe)
        !
        !coro=sin(alpha)+sin(rhoe)
        !
        !coro = 1.0d+00-exp(-1.0d0*(pi-rhoe-alpha))
        !
        !a13 = a1
        !
        !coro = 1.0d+00-exp(-a13*(r13-r13e))
        !
        !
        !aa1 = a0
        !
        y1=1.0d+00-exp(-aa1*(r12-r12e))
        y3=1.0d+00-exp(-aa1*(r32-r12e))
        !
        ahh = param(75)
        fhh = param(76)
        !
        !y = 0.5d0*(r12-r32) ! *exp(-ahh*(r12-r32)**2)
        y = 0.5d0*(y1-y3) ! *exp(-ahh*(r12-r32)**2)
        !
        vhh = fhh*(1.d0-exp(-ahh*(r13-r13e) ) )**2 !  0.9e5*exp(-10d0*r13) !+0.518622e5*y1**2+0.518622e5*y3**2
        !
        !y1=1.0d+00-exp(-aa1*(r12-r12e))
        !y3=1.0d+00-exp(-aa1*(r32-r12e))
        FA1       = param(10)
        FA2       = param(11)
        FA3       = param(12)
        FA4       = param(13)
        FA5       = param(14)
        FA6       = param(15)
        FA7       = param(16)
        FA8       = param(17)
        F1        = param(26)
        F1A1      = param(27)
        F2A1      = param(28)
        F3A1      = param(29)
        F4A1      = param(30)
        F11       = param(31)
        F1A11     = param(32)
        F2A11     = param(33)
        F3A11     = param(34)
        F4A11     = param(35)
        F5A11     = param(36)
        F6A11     = param(37)
        F13       = param(38)
        F1A13     = param(39)
        F2A13     = param(40)
        F3A13     = param(41)
        F4A13     = param(42)
        F5A13     = param(43)
        F6A13     = param(44)
        F111      = param(45)
        F1A111    = param(46)
        F2A111    = param(47)
        F113      = param(48)
        F1A113    = param(49)
        F2A113    = param(50)
        F1111     = param(51)
        F1A1111   = param(52)
        F2A1111   = param(53)
        F1113     = param(54)
        F1A1113   = param(55)
        F1133     = param(56)
        F1A1133   = param(57)
        f11111    = param(58)
        f1a11111  = param(59)
        f2a11111  = param(60)
        f11113    = param(61)
        f1a11113  = param(62)
        f2a11113  = param(63)
        f11133    = param(64)
        f1a11133  = param(65)
        f2a11133  = param(66)
        f111111   = param(67)
        f1a111111 = param(68)
        f111113   = param(69)
        f1a111113 = param(70)
        f111133   = param(71)
        f1a111133 = param(72)
        f111333   = param(73)
        f1a111333 = param(74)


        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f4a33   = f4a11
        f5a33   = f5a11
        f6a33   = f6a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f2a3333 = f2a1111
        f1a1333 = f1a1113

        f33333= f11111
        f13333= f11113
        f11333= f11133
        f1a33333= f1a11111
        f1a13333= f1a11113
        f1a11333= f1a11133
        f2a33333= f2a11111
        f2a13333= f2a11113
        f2a11333= f2a11133

        f333333= f111111
        f133333= f111113
        f113333= f111133
        f1a333333= f1a111111
        f1a133333= f1a111113
        f1a113333= f1a111133

!
! calculate potential energy function values
!

      v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 +fa6*coro**6+fa7*coro**7+fa8*coro**8
      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4+f5a11*coro**5+f6a11*coro**6
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3+f4a33*coro**4+f5a33*coro**5+f6a33*coro**6
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3+f4a13*coro**4+f5a13*coro**5+f6a13*coro**6
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro+f2a1111*coro**2
      fe3333= f3333+f1a3333*coro+f2a3333*coro**2
      fe1113= f1113+f1a1113*coro !!!!!+f2a1113*coro**2
      fe1333= f1333+f1a1333*coro !!!!!+f2a1333*coro**2
      fe1133= f1133+f1a1133*coro !!!!!+f2a1133*coro**2

      fe11111= f11111+f1a11111*coro+f2a11111*coro**2
      fe11113= f11113+f1a11113*coro+f2a11113*coro**2
      fe11133= f11133+f1a11133*coro+f2a11133*coro**2
      fe11333= f11333+f1a11333*coro+f2a11333*coro**2
      fe13333= f13333+f1a13333*coro+f2a13333*coro**2
      fe33333= f33333+f1a33333*coro+f2a33333*coro**2

      fe111111= f111111+f1a111111*coro
      fe111113= f111113+f1a111113*coro
      fe111133= f111133+f1a111133*coro
      fe111333= f111333+f1a111333*coro
      fe113333= f113333+f1a113333*coro
      fe133333= f133333+f1a133333*coro
      fe333333= f333333+f1a333333*coro

      !fe1111111= f1111111 !+f1a1111111*coro
      !fe1111113= f1111113 !+f1a1111113*coro
      !fe1111133= f1111133 !+f1a1111133*coro
      !fe1111333= f1111333 !+f1a1111333*coro
      !fe1113333= f1113333 !+f1a1113333*coro
      !fe1133333= f1133333 !+f1a1133333*coro
      !fe1333333= f1333333 !+f1a1333333*coro
      !fe3333333= f3333333 !+f1a3333333*coro

      !fe11111111= f11111111 !+f1a11111111*coro
      !fe11111113= f11111113 !+f1a11111113*coro
      !fe11111133= f11111133 !+f1a11111133*coro
      !fe11111333= f11111333 !+f1a11111333*coro
      !fe11113333= f11113333 !+f1a11113333*coro
      !fe11133333= f11133333 !!+f1a11333333*coro
      !fe11333333= f11333333 !+f1a13333333*coro
      !fe13333333= f13333333 !+f1a13333333*coro
      !fe33333333= f33333333 !+f1a33333333*coro



      v     =  v0+fe1*y1+fe3*y3                          &  
!              +fe11*y1**2+fe33*y3**2+fe13*y1*y3          & 
               +fe11*y1**2+fe33*y3**2+fe13*y**2          & 
!               +fe11*y1**2+fe33*y3**2 &
!               +fe13*y**2          & 
!               +fe13*(r12-r12e)*(r32-r12e) &
              +fe111*(y1**3+y3**3)+fe113*(y1+y3)*y**2 & 
              +fe1111*(y1**4+y3**4)+fe1133*(y1**2+y3**2)*y**2+fe1113*y**4

      v5 = fe11111*y1**5+fe11113*y1**4*y3**1+fe11133*y1**3*y3**2+fe11333*y1**2*y3**3+fe13333*y1**1*y3**4+fe33333*y3**5
      v6 = fe111111*y1**6+fe111113*y1**5*y3**1+fe111133*y1**4*y3**2+fe111333*y1**3*y3**3+fe113333*y1**2*y3**4+fe133333*y1**1*y3**5+fe333333*y3**6
      !v7 = fe1111111*y1**7+fe1111113*y1**6*y3**1+fe1111133*y1**5*y3**2+fe1111333*y1**4*y3**3+&
      !     fe1113333*y1**3*y3**4+fe1133333*y1**2*y3**5+fe1333333*y1**1*y3**6+fe3333333*y1**0*y3**7
      !v8 = fe11111111*y1**8+fe11111113*y1**7*y3**1+fe11111133*y1**6*y3**2+fe11111333*y1**5*y3**3+&
      !     fe11113333*y1**4*y3**4+fe11133333*y1**3*y3**5+fe11333333*y1**2*y3**6+fe13333333*y1**1*y3**7+fe33333333*y1**0*y3**8

       morbid_MEP2 = v+v5+v6+vhh


 end function morbid_MEP2




 double precision function morbid_MEP3(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 

      double precision   ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,  &
            f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,    & 
            fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
            fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
            f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
            f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133,r12,r32,alpha


      double precision  :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,   &    
            fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,    &
            f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,   &
            f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,       & 
            re32, aa3  ,f3,f33  ,f333 ,f133 ,f3333,f1333      
      double precision ::  &
      fe11111,f11111,fe11113,f11113,fe11133,f11133,fe11333,f11333,&
      fe13333,f13333,fe33333,f33333,fe111111,f111111,fe111113,f111113,&
      fe111133,f111133,fe111333,f111333,fe113333,f113333,fe133333,f133333,&
      fe333333,f333333
      double precision ::  &
      f2a1111,f2a3333,f2a1113,f2a1333,f2a1133,f1a11111,f2a11111,&
      f1a11113,f2a11113,f1a11133,f2a11133,&
      f1a11333,f2a11333,f1a13333,f2a13333,&
      f1a33333,f2a33333,f1a111111,f1a111113,f1a111133,f1a111333,f1a113333,&
      f1a133333,f1a333333      

      double precision ::  &
      f1111111,f1111113,f1111133,f11111111,f11111113,&
      f11111133,f11111333,f11113333,f1111333,f3333333,f1333333,f1133333,f1113333,&
      f33333333,f13333333,f11333333,f11133333

      double precision ::  &
      fe1111111,fe1111113,fe1111133,fe1111333,fe1113333,fe1133333,fe1333333,fe3333333,&
      fe11111111,fe11111113,fe11111133,fe11111333,fe11113333,fe11133333,fe13333333,fe33333333,fe11333333,&
      F4A11,f4a33,F5A11,F6A11,F5A33,F6A33,F4A13,F5A13,F6A13


        double precision   :: pi,y1,y3,v,coro,v0,RHOE,v5,v6,v7,v8,alphae,r13e,r12e,r13_0,r13_inf,a13,r,r13,re
        double precision   :: a0,a1,a2,a3,a4,a5,a6,y,vhh,fhh,ahh,g1,g2,g3,g4,g5,g6



        !
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
        !
        r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
        !
	    R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
        !
        pi = 4.0d0 * atan2(1.0d0,1.0d0)
        !
        VE        = param( 1)
        RHOE      = param( 9)*pi/180.0d0
        !
        a0      = param(2)
        a1      = param(3)
        a2      = param(4)
        a3      = param(5)
        a4      = param(6)
        a5      = param(7)
        a6      = param(8)
        !
        !a0      = param(2)
        !a13     = param(3)
        !a1      = param(4)
        !a2      = param(5)
        !a3      = param(6)
        !a4      = param(7)
        !a5      = param(8)
        !
        !r13e = a0+a1*exp(-a13*R**2)+a2*exp(-a13*R**4)+a3*exp(-a13*R**6)+a4*exp(-a13*R**8)+a5*exp(-a13*R**10)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        coro=cos(alpha)+cos(rhoe)
        !
        g1 = coro**1
        g2 = coro**2
        g3 = coro**3
        g4 = coro**4
        g5 = coro**5
        g6 = coro**6
        !
        !f1 = exp(-a13*t1)
        !f2 = exp(-a13*t2)
        !f3 = exp(-a13*t3)
        !f4 = exp(-a13*t4)
        !f5 = exp(-a13*t5)
        !
        !r13e = r13_inf+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5
        !
        r13e = a0+a1*g1+a2*g2+a3*g3+a4*g4+a5*g5+a6*g6
        !
        !r13e = r13_inf+(r13_0-r13_inf)*exp(-a13*sqrt(R)**9)
        !
        r12e = 0.5d0*r13e/sin(alpha*0.5d0)  !  sqrt(0.25d0*r13e**2+R**2)
        !
        !y = 1.0d0/(alpha+1e-14)
        !
        !r12e = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
        !r13e = 2.0d0*r12e*sin(0.5d0*alpha)
        !
        !alphae = atan(0.5d0*r13e/R)*2.d0
        !
        !
        !RE12      = param(14)
        !RHOE      = param( 5)*pi/180.0d0
        !
        !Re = 0.5d0*sqrt(2.0d0*r12e**2-2.0d0*re12*re12*cos(rhoe))
        !
        RE12      = param(18)
        A0        = param(19)
        a1        = param(20)
        a2        = param(21)
        a3        = param(22)
        a4        = param(23)
        a5        = param(24)
        a6        = param(25)
        !
        aa1 = a0+a1*coro+a2*coro**2+a3*coro**3+a4*coro**4+a5*coro**5+a6*coro**6
        !aa1 = a0
        !
        coro=cos(alpha)+cos(rhoe)
        !
        !coro=sin(2.0d0*alpha) ! +cos(2.0d0*rhoe)
        !
        !coro=1.0d0-exp(-a0*(alpha-(pi-rhoe)))
        !
        !coro=alpha-(pi-rhoe)
        !
        !coro=sin(alpha)+sin(rhoe)
        !
        !coro = 1.0d+00-exp(-1.0d0*(pi-rhoe-alpha))
        !
        !a13 = a1
        !
        !coro = 1.0d+00-exp(-a13*(r13-r13e))
        !
        !
        !aa1 = a0
        !
        y1=1.0d+00-exp(-aa1*(r12-r12e))
        y3=1.0d+00-exp(-aa1*(r32-r12e))
        !
        ahh = param(75)
        fhh = param(76)
        !
        !y = 0.5d0*(r12-r32) ! *exp(-ahh*(r12-r32)**2)
        y = 0.5d0*(y1-y3) ! *exp(-ahh*(r12-r32)**2)
        !
        vhh = fhh*(1.d0-exp(-ahh*(r13-r13e) ) )**2 !  0.9e5*exp(-10d0*r13) !+0.518622e5*y1**2+0.518622e5*y3**2
        !
        !y1=1.0d+00-exp(-aa1*(r12-r12e))
        !y3=1.0d+00-exp(-aa1*(r32-r12e))
        FA1       = param(10)
        FA2       = param(11)
        FA3       = param(12)
        FA4       = param(13)
        FA5       = param(14)
        FA6       = param(15)
        FA7       = param(16)
        FA8       = param(17)
        F1        = param(26)
        F1A1      = param(27)
        F2A1      = param(28)
        F3A1      = param(29)
        F4A1      = param(30)
        F11       = param(31)
        F1A11     = param(32)
        F2A11     = param(33)
        F3A11     = param(34)
        F4A11     = param(35)
        F5A11     = param(36)
        F6A11     = param(37)
        F13       = param(38)
        F1A13     = param(39)
        F2A13     = param(40)
        F3A13     = param(41)
        F4A13     = param(42)
        F5A13     = param(43)
        F6A13     = param(44)
        F111      = param(45)
        F1A111    = param(46)
        F2A111    = param(47)
        F113      = param(48)
        F1A113    = param(49)
        F2A113    = param(50)
        F1111     = param(51)
        F1A1111   = param(52)
        F2A1111   = param(53)
        F1113     = param(54)
        F1A1113   = param(55)
        F1133     = param(56)
        F1A1133   = param(57)
        f11111    = param(58)
        f1a11111  = param(59)
        f2a11111  = param(60)
        f11113    = param(61)
        f1a11113  = param(62)
        f2a11113  = param(63)
        f11133    = param(64)
        f1a11133  = param(65)
        f2a11133  = param(66)
        f111111   = param(67)
        f1a111111 = param(68)
        f111113   = param(69)
        f1a111113 = param(70)
        f111133   = param(71)
        f1a111133 = param(72)
        f111333   = param(73)
        f1a111333 = param(74)


        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f4a33   = f4a11
        f5a33   = f5a11
        f6a33   = f6a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f2a3333 = f2a1111
        f1a1333 = f1a1113

        f33333= f11111
        f13333= f11113
        f11333= f11133
        f1a33333= f1a11111
        f1a13333= f1a11113
        f1a11333= f1a11133
        f2a33333= f2a11111
        f2a13333= f2a11113
        f2a11333= f2a11133

        f333333= f111111
        f133333= f111113
        f113333= f111133
        f1a333333= f1a111111
        f1a133333= f1a111113
        f1a113333= f1a111133

!
! calculate potential energy function values
!

      v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 +fa6*coro**6+fa7*coro**7+fa8*coro**8
      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4+f5a11*coro**5+f6a11*coro**6
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3+f4a33*coro**4+f5a33*coro**5+f6a33*coro**6
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3+f4a13*coro**4+f5a13*coro**5+f6a13*coro**6
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro+f2a1111*coro**2
      fe3333= f3333+f1a3333*coro+f2a3333*coro**2
      fe1113= f1113+f1a1113*coro !!!!!+f2a1113*coro**2
      fe1333= f1333+f1a1333*coro !!!!!+f2a1333*coro**2
      fe1133= f1133+f1a1133*coro !!!!!+f2a1133*coro**2

      fe11111= f11111+f1a11111*coro+f2a11111*coro**2
      fe11113= f11113+f1a11113*coro+f2a11113*coro**2
      fe11133= f11133+f1a11133*coro+f2a11133*coro**2
      fe11333= f11333+f1a11333*coro+f2a11333*coro**2
      fe13333= f13333+f1a13333*coro+f2a13333*coro**2
      fe33333= f33333+f1a33333*coro+f2a33333*coro**2

      fe111111= f111111+f1a111111*coro
      fe111113= f111113+f1a111113*coro
      fe111133= f111133+f1a111133*coro
      fe111333= f111333+f1a111333*coro
      fe113333= f113333+f1a113333*coro
      fe133333= f133333+f1a133333*coro
      fe333333= f333333+f1a333333*coro

      !fe1111111= f1111111 !+f1a1111111*coro
      !fe1111113= f1111113 !+f1a1111113*coro
      !fe1111133= f1111133 !+f1a1111133*coro
      !fe1111333= f1111333 !+f1a1111333*coro
      !fe1113333= f1113333 !+f1a1113333*coro
      !fe1133333= f1133333 !+f1a1133333*coro
      !fe1333333= f1333333 !+f1a1333333*coro
      !fe3333333= f3333333 !+f1a3333333*coro

      !fe11111111= f11111111 !+f1a11111111*coro
      !fe11111113= f11111113 !+f1a11111113*coro
      !fe11111133= f11111133 !+f1a11111133*coro
      !fe11111333= f11111333 !+f1a11111333*coro
      !fe11113333= f11113333 !+f1a11113333*coro
      !fe11133333= f11133333 !!+f1a11333333*coro
      !fe11333333= f11333333 !+f1a13333333*coro
      !fe13333333= f13333333 !+f1a13333333*coro
      !fe33333333= f33333333 !+f1a33333333*coro



      v     =  v0+fe1*y1+fe3*y3                          &  
!               +fe11*y1**2+fe33*y3**2+fe13*y1*y3          & 
               +fe11*y1**2+fe33*y3**2+fe13*y**2          & 
!               +fe11*y1**2+fe33*y3**2 &
!               +fe13*y**2          & 
!               +fe13*(r12-r12e)*(r32-r12e) &
              +fe111*(y1**3+y3**3)+fe113*(y1+y3)*y**2 & 
              +fe1111*(y1**4+y3**4)+fe1133*(y1**2+y3**2)*y**2+fe1113*y**4

      v5 = fe11111*y1**5+fe11113*y1**4*y3**1+fe11133*y1**3*y3**2+fe11333*y1**2*y3**3+fe13333*y1**1*y3**4+fe33333*y3**5
      v6 = fe111111*y1**6+fe111113*y1**5*y3**1+fe111133*y1**4*y3**2+fe111333*y1**3*y3**3+fe113333*y1**2*y3**4+fe133333*y1**1*y3**5+fe333333*y3**6
      !v7 = fe1111111*y1**7+fe1111113*y1**6*y3**1+fe1111133*y1**5*y3**2+fe1111333*y1**4*y3**3+&
      !     fe1113333*y1**3*y3**4+fe1133333*y1**2*y3**5+fe1333333*y1**1*y3**6+fe3333333*y1**0*y3**7
      !v8 = fe11111111*y1**8+fe11111113*y1**7*y3**1+fe11111133*y1**6*y3**2+fe11111333*y1**5*y3**3+&
      !     fe11113333*y1**4*y3**4+fe11133333*y1**3*y3**5+fe11333333*y1**2*y3**6+fe13333333*y1**1*y3**7+fe33333333*y1**0*y3**8

       morbid_MEP3 = v+v5+v6+vhh


 end function morbid_MEP3






 double precision function morbid_MEP4(local,ZPE,parmax,param)

 double precision,intent(in) ::  local(3)
 double precision,intent(in) ::  ZPE
 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)

                 

      double precision   ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,  &
            f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,    & 
            fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
            fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
            f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
            f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133,r12,r32,alpha


      double precision  :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,   &    
            fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,    &
            f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,   &
            f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,       & 
            re32, aa3  ,f3,f33  ,f333 ,f133 ,f3333,f1333      
      double precision ::  &
      fe11111,f11111,fe11113,f11113,fe11133,f11133,fe11333,f11333,&
      fe13333,f13333,fe33333,f33333,fe111111,f111111,fe111113,f111113,&
      fe111133,f111133,fe111333,f111333,fe113333,f113333,fe133333,f133333,&
      fe333333,f333333
      double precision ::  &
      f2a1111,f2a3333,f2a1113,f2a1333,f2a1133,f1a11111,f2a11111,&
      f1a11113,f2a11113,f1a11133,f2a11133,&
      f1a11333,f2a11333,f1a13333,f2a13333,&
      f1a33333,f2a33333,f1a111111,f1a111113,f1a111133,f1a111333,f1a113333,&
      f1a133333,f1a333333      

      double precision ::  &
      f1111111,f1111113,f1111133,f11111111,f11111113,&
      f11111133,f11111333,f11113333,f1111333,f3333333,f1333333,f1133333,f1113333,&
      f33333333,f13333333,f11333333,f11133333

      double precision ::  &
      fe1111111,fe1111113,fe1111133,fe1111333,fe1113333,fe1133333,fe1333333,fe3333333,&
      fe11111111,fe11111113,fe11111133,fe11111333,fe11113333,fe11133333,fe13333333,fe33333333,fe11333333,&
      F4A11,f4a33,F5A11,F6A11,F5A33,F6A33,F4A13,F5A13,F6A13


        double precision   :: pi,y1,y3,v,coro,v0,RHOE,v5,v6,v7,v8,alphae,r13long,r12e,r13_0,r13_inf,a13,r,r13,re
        double precision   :: a0,a1,a2,a3,a4,a5,a6,y,vhh,fhh,ahh,g1,g2,g3,g4,g5,g6,a12,b13,re13,De,De13,damp2,damp4,vlong,vdamp,x1,x3,dampR2



        !
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
        !
        r13 = sqrt(r12**2+r32**2-2.0d0*r12*r32*cos(alpha))
        !
	    R = 0.5d0*sqrt(r12**2+r32**2+2.0d0*r12*r32*cos(alpha))
        !
        pi = 4.0d0 * atan2(1.0d0,1.0d0)
        !
        VE        = param( 1)
        !
        re12      = param( 2)
        alphae    = param( 3)*pi/180.0d0
        !
        re13 = re12*2.0d0*sin(alphae*0.5d0)
        !
        Re = sqrt(re12**2-0.25d0*re13**2)
        !
        r13long   = param(4)
        !
        a12       = param(5)
        a13       = param(6)
        !
        b13       = param(7)
        !
        De        = param(8)
        De13      = param(9)
        !
        damp2     = param(10)
        damp4     = param(11)
        dampR2    = param(12)
        !
        !coro=cos(alpha)-cos(alphae)
        !
        y1=1.0d+00-exp(-a12*(r12-re12))
        y3=1.0d+00-exp(-a12*(r32-re12))
        !
        x1 = r12-re12
        x3 = r32-re12
        !
        coro = 1.0d+00-exp(-a13*(r13-re13))
        !
        !vlong = De13*(1.0d0-exp(-b13*(r13-r13long) ) )**2 + De
        !
        vdamp = exp( -damp2*( x1**2+x3**2)-damp4*( x1**4+x3**4 )-dampR2*(R-Re)**2)
        !
        vlong = De*(1.0d0-exp(-b13*(R-Re) ) )**2+De13*(1.0d0-exp(-b13*(r13-r13long) ) )**2*(1.0d0-vdamp)
        !
        FA1       = param(13)
        FA2       = param(14)
        FA3       = param(15)
        FA4       = param(16)
        FA5       = param(17)
        FA6       = param(18)
        FA7       = param(19)
        FA8       = param(20)
        F1        = param(21)
        F1A1      = param(22)
        F2A1      = param(23)
        F3A1      = param(24)
        F4A1      = param(25)
        F11       = param(26)
        F1A11     = param(27)
        F2A11     = param(28)
        F3A11     = param(29)
        F4A11     = param(30)
        F5A11     = param(31)
        F6A11     = param(32)
        F13       = param(33)
        F1A13     = param(34)
        F2A13     = param(35)
        F3A13     = param(36)
        F4A13     = param(37)
        F5A13     = param(38)
        F6A13     = param(39)
        F111      = param(40)
        F1A111    = param(41)
        F2A111    = param(42)
        F113      = param(43)
        F1A113    = param(44)
        F2A113    = param(45)
        F1111     = param(46)
        F1A1111   = param(47)
        F2A1111   = param(48)
        F1113     = param(49)
        F1A1113   = param(50)
        F1133     = param(51)
        F1A1133   = param(52)
        f11111    = param(53)
        f1a11111  = param(54)
        f2a11111  = param(55)
        f11113    = param(56)
        f1a11113  = param(57)
        f2a11113  = param(58)
        f11133    = param(59)
        f1a11133  = param(60)
        f2a11133  = param(61)
        f111111   = param(62)
        f1a111111 = param(63)
        f111113   = param(64)
        f1a111113 = param(65)
        f111133   = param(66)
        f1a111133 = param(67)
        f111333   = param(68)
        f1a111333 = param(69)


        f3  = f1
        f33  = f11
        f333  = f111
        f133  = f113
        f3333  = f1111
        f1333  = f1113

        f1a3    = f1a1
        f2a3    = f2a1
        f3a3    = f3a1
        f4a3    = f4a1
        f1a33   = f1a11 
        f2a33   = f2a11
        f3a33   = f3a11
        f4a33   = f4a11
        f5a33   = f5a11
        f6a33   = f6a11
        f1a333  = f1a111
        f2a333  = f2a111
        f1a133  = f1a113
        f2a133  = f2a113
        f1a3333 = f1a1111
        f2a3333 = f2a1111
        f1a1333 = f1a1113

        f33333= f11111
        f13333= f11113
        f11333= f11133
        f1a33333= f1a11111
        f1a13333= f1a11113
        f1a11333= f1a11133
        f2a33333= f2a11111
        f2a13333= f2a11113
        f2a11333= f2a11133

        f333333= f111111
        f133333= f111113
        f113333= f111133
        f1a333333= f1a111111
        f1a133333= f1a111113
        f1a113333= f1a111133

!
! calculate potential energy function values
!

      v0= fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5 +fa6*coro**6+fa7*coro**7+fa8*coro**8
      fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
      fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
      fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4+f5a11*coro**5+f6a11*coro**6
      fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3+f4a33*coro**4+f5a33*coro**5+f6a33*coro**6
      fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3+f4a13*coro**4+f5a13*coro**5+f6a13*coro**6
      fe111= f111+f1a111*coro+f2a111*coro**2
      fe333= f333+f1a333*coro+f2a333*coro**2
      fe113= f113+f1a113*coro+f2a113*coro**2
      fe133= f133+f1a133*coro+f2a133*coro**2
      fe1111= f1111+f1a1111*coro+f2a1111*coro**2
      fe3333= f3333+f1a3333*coro+f2a3333*coro**2
      fe1113= f1113+f1a1113*coro !!!!!+f2a1113*coro**2
      fe1333= f1333+f1a1333*coro !!!!!+f2a1333*coro**2
      fe1133= f1133+f1a1133*coro !!!!!+f2a1133*coro**2

      fe11111= f11111+f1a11111*coro+f2a11111*coro**2
      fe11113= f11113+f1a11113*coro+f2a11113*coro**2
      fe11133= f11133+f1a11133*coro+f2a11133*coro**2
      fe11333= f11333+f1a11333*coro+f2a11333*coro**2
      fe13333= f13333+f1a13333*coro+f2a13333*coro**2
      fe33333= f33333+f1a33333*coro+f2a33333*coro**2

      fe111111= f111111+f1a111111*coro
      fe111113= f111113+f1a111113*coro
      fe111133= f111133+f1a111133*coro
      fe111333= f111333+f1a111333*coro
      fe113333= f113333+f1a113333*coro
      fe133333= f133333+f1a133333*coro
      fe333333= f333333+f1a333333*coro

      !fe1111111= f1111111 !+f1a1111111*coro
      !fe1111113= f1111113 !+f1a1111113*coro
      !fe1111133= f1111133 !+f1a1111133*coro
      !fe1111333= f1111333 !+f1a1111333*coro
      !fe1113333= f1113333 !+f1a1113333*coro
      !fe1133333= f1133333 !+f1a1133333*coro
      !fe1333333= f1333333 !+f1a1333333*coro
      !fe3333333= f3333333 !+f1a3333333*coro

      !fe11111111= f11111111 !+f1a11111111*coro
      !fe11111113= f11111113 !+f1a11111113*coro
      !fe11111133= f11111133 !+f1a11111133*coro
      !fe11111333= f11111333 !+f1a11111333*coro
      !fe11113333= f11113333 !+f1a11113333*coro
      !fe11133333= f11133333 !!+f1a11333333*coro
      !fe11333333= f11333333 !+f1a13333333*coro
      !fe13333333= f13333333 !+f1a13333333*coro
      !fe33333333= f33333333 !+f1a33333333*coro



      v     =  v0+fe1*y1+fe3*y3                          &  
              +fe11*y1**2+fe33*y3**2+fe13*y1*y3          & 
              +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3    &
              +fe133*y1*y3**2                            & 
              +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3 & 
              +fe1333*y1*y3**3+fe1133*y1**2*y3**2

      v5 = fe11111*y1**5+fe11113*y1**4*y3**1+fe11133*y1**3*y3**2+fe11333*y1**2*y3**3+fe13333*y1**1*y3**4+fe33333*y3**5
      v6 = fe111111*y1**6+fe111113*y1**5*y3**1+fe111133*y1**4*y3**2+fe111333*y1**3*y3**3+fe113333*y1**2*y3**4+fe133333*y1**1*y3**5+fe333333*y3**6
      !
      !v7 = fe1111111*y1**7+fe1111113*y1**6*y3**1+fe1111133*y1**5*y3**2+fe1111333*y1**4*y3**3+&
      !     fe1113333*y1**3*y3**4+fe1133333*y1**2*y3**5+fe1333333*y1**1*y3**6+fe3333333*y1**0*y3**7
      !v8 = fe11111111*y1**8+fe11111113*y1**7*y3**1+fe11111133*y1**6*y3**2+fe11111333*y1**5*y3**3+&
      !     fe11113333*y1**4*y3**4+fe11133333*y1**3*y3**5+fe11333333*y1**2*y3**6+fe13333333*y1**1*y3**7+fe33333333*y1**0*y3**8

       morbid_MEP4 = ve+ (v+v5+v6)*vdamp+vlong


 end function morbid_MEP4


  double precision function dipol_xy2_p_dump(local,ZPE,N ,force) result(f)
      !
      double precision,intent(in) :: local(1:3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          :: N
      double precision,intent(in) :: force(1:N)
      !
      double precision            :: y1,y2,y3,r1,r2,alpha
      double precision            :: ae,re,v0,v1,v2,v3,v4,v5,v6,v7,v8,b0
      !
      double precision,parameter  :: pi=3.1415926535897932385d0
      !
      re     = force(1)
      ae     = force(2)/180.0d0*pi
      !
      r1 = local(1)
      r2 = local(2)
      alpha = local(3)
      !
      b0 = 1.0d0
      !
      y1 = (r1 - re)*exp(-b0*(r1-re)**2)
      y2 = (r2 - re)*exp(-b0*(r2-re)**2)
      y3 = cos(alpha) - cos(ae)
      !
 v1 = force( 3)*y1**1*y2**0*y3**0& 
    - force( 3)*y1**0*y2**1*y3**0
 v2 = force( 4)*y1**1*y2**0*y3**1& 
    - force( 4)*y1**0*y2**1*y3**1&  
    + force( 5)*y1**2*y2**0*y3**0& 
    - force( 5)*y1**0*y2**2*y3**0

 v3 = force( 6)*y1**1*y2**0*y3**2& 
    - force( 6)*y1**0*y2**1*y3**2& 
    + force( 7)*y1**2*y2**0*y3**1& 
    - force( 7)*y1**0*y2**2*y3**1& 
    + force( 8)*y1**2*y2**1*y3**0& 
    - force( 8)*y1**1*y2**2*y3**0& 
    + force( 9)*y1**3*y2**0*y3**0& 
    - force( 9)*y1**0*y2**3*y3**0

 if (N> 9) then 
  v4 =force(10)*y1**1*y2**0*y3**3& 
    - force(10)*y1**0*y2**1*y3**3& 
    + force(11)*y1**2*y2**0*y3**2& 
    - force(11)*y1**0*y2**2*y3**2& 
    + force(12)*y1**2*y2**1*y3**1& 
    - force(12)*y1**1*y2**2*y3**1& 
    + force(13)*y1**3*y2**0*y3**1& 
    - force(13)*y1**0*y2**3*y3**1& 
    + force(14)*y1**3*y2**1*y3**0& 
    - force(14)*y1**1*y2**3*y3**0& 
    + force(15)*y1**4*y2**0*y3**0& 
    - force(15)*y1**0*y2**4*y3**0
endif

 if (N>15) then 
  v5 =force(16)*y1**1*y2**0*y3**4& 
    - force(16)*y1**0*y2**1*y3**4& 
    + force(17)*y1**2*y2**0*y3**3& 
    - force(17)*y1**0*y2**2*y3**3& 
    + force(18)*y1**2*y2**1*y3**2& 
    - force(18)*y1**1*y2**2*y3**2& 
    + force(19)*y1**3*y2**0*y3**2& 
    - force(19)*y1**0*y2**3*y3**2& 
    + force(20)*y1**3*y2**1*y3**1& 
    - force(20)*y1**1*y2**3*y3**1& 
    + force(21)*y1**3*y2**2*y3**0& 
    - force(21)*y1**2*y2**3*y3**0& 
    + force(22)*y1**4*y2**0*y3**1& 
    - force(22)*y1**0*y2**4*y3**1& 
    + force(23)*y1**4*y2**1*y3**0& 
    - force(23)*y1**1*y2**4*y3**0& 
    + force(24)*y1**5*y2**0*y3**0& 
    - force(24)*y1**0*y2**5*y3**0
endif

 if (N>24) then 
  v6 =force(25)*y1**1*y2**0*y3**5& 
    - force(25)*y1**0*y2**1*y3**5& 
    + force(26)*y1**2*y2**0*y3**4& 
    - force(26)*y1**0*y2**2*y3**4& 
    + force(27)*y1**2*y2**1*y3**3& 
    - force(27)*y1**1*y2**2*y3**3& 
    + force(28)*y1**3*y2**0*y3**3& 
    - force(28)*y1**0*y2**3*y3**3& 
    + force(29)*y1**3*y2**1*y3**2& 
    - force(29)*y1**1*y2**3*y3**2& 
    + force(30)*y1**3*y2**2*y3**1& 
    - force(30)*y1**2*y2**3*y3**1& 
    + force(31)*y1**4*y2**0*y3**2& 
    - force(31)*y1**0*y2**4*y3**2& 
    + force(32)*y1**4*y2**1*y3**1& 
    - force(32)*y1**1*y2**4*y3**1& 
    + force(33)*y1**4*y2**2*y3**0& 
    - force(33)*y1**2*y2**4*y3**0& 
    + force(34)*y1**5*y2**0*y3**1& 
    - force(34)*y1**0*y2**5*y3**1& 
    + force(35)*y1**5*y2**1*y3**0& 
    - force(35)*y1**1*y2**5*y3**0& 
    + force(36)*y1**6*y2**0*y3**0& 
    - force(36)*y1**0*y2**6*y3**0
 endif

 if (N>36) then 
 v7 = force(37)*y1**1*y2**0*y3**6& 
    - force(37)*y1**0*y2**1*y3**6& 
    + force(38)*y1**2*y2**0*y3**5& 
    - force(38)*y1**0*y2**2*y3**5& 
    + force(39)*y1**2*y2**1*y3**4& 
    - force(39)*y1**1*y2**2*y3**4& 
    + force(40)*y1**3*y2**0*y3**4& 
    - force(40)*y1**0*y2**3*y3**4& 
    + force(41)*y1**3*y2**1*y3**3& 
    - force(41)*y1**1*y2**3*y3**3& 
    + force(42)*y1**3*y2**2*y3**2& 
    - force(42)*y1**2*y2**3*y3**2& 
    + force(43)*y1**4*y2**0*y3**3& 
    - force(43)*y1**0*y2**4*y3**3& 
    + force(44)*y1**4*y2**1*y3**2& 
    - force(44)*y1**1*y2**4*y3**2& 
    + force(45)*y1**4*y2**2*y3**1& 
    - force(45)*y1**2*y2**4*y3**1& 
    + force(46)*y1**4*y2**3*y3**0& 
    - force(46)*y1**3*y2**4*y3**0& 
    + force(47)*y1**5*y2**0*y3**2& 
    - force(47)*y1**0*y2**5*y3**2& 
    + force(48)*y1**5*y2**1*y3**1& 
    - force(48)*y1**1*y2**5*y3**1& 
    + force(49)*y1**5*y2**2*y3**0& 
    - force(49)*y1**2*y2**5*y3**0& 
    + force(50)*y1**6*y2**0*y3**1& 
    - force(50)*y1**0*y2**6*y3**1& 
    + force(51)*y1**6*y2**1*y3**0& 
    - force(51)*y1**1*y2**6*y3**0& 
    + force(52)*y1**7*y2**0*y3**0& 
    - force(52)*y1**0*y2**7*y3**0
 endif

 if (N>52) then 
 v8 = force(53)*y1**1*y2**0*y3**7& 
    - force(53)*y1**0*y2**1*y3**7& 
    + force(54)*y1**2*y2**0*y3**6& 
    - force(54)*y1**0*y2**2*y3**6& 
    + force(55)*y1**2*y2**1*y3**5& 
    - force(55)*y1**1*y2**2*y3**5& 
    + force(56)*y1**3*y2**0*y3**5& 
    - force(56)*y1**0*y2**3*y3**5& 
    + force(57)*y1**3*y2**1*y3**4& 
    - force(57)*y1**1*y2**3*y3**4& 
    + force(58)*y1**3*y2**2*y3**3& 
    - force(58)*y1**2*y2**3*y3**3& 
    + force(59)*y1**4*y2**0*y3**4& 
    - force(59)*y1**0*y2**4*y3**4& 
    + force(60)*y1**4*y2**1*y3**3& 
    - force(60)*y1**1*y2**4*y3**3& 
    + force(61)*y1**4*y2**2*y3**2& 
    - force(61)*y1**2*y2**4*y3**2& 
    + force(62)*y1**4*y2**3*y3**1& 
    - force(62)*y1**3*y2**4*y3**1& 
    + force(63)*y1**5*y2**0*y3**3& 
    - force(63)*y1**0*y2**5*y3**3& 
    + force(64)*y1**5*y2**1*y3**2& 
    - force(64)*y1**1*y2**5*y3**2& 
    + force(65)*y1**5*y2**2*y3**1& 
    - force(65)*y1**2*y2**5*y3**1& 
    + force(66)*y1**5*y2**3*y3**0& 
    - force(66)*y1**3*y2**5*y3**0& 
    + force(67)*y1**6*y2**0*y3**2& 
    - force(67)*y1**0*y2**6*y3**2& 
    + force(68)*y1**6*y2**1*y3**1& 
    - force(68)*y1**1*y2**6*y3**1& 
    + force(69)*y1**6*y2**2*y3**0& 
    - force(69)*y1**2*y2**6*y3**0& 
    + force(70)*y1**7*y2**0*y3**1& 
    - force(70)*y1**0*y2**7*y3**1& 
    + force(71)*y1**7*y2**1*y3**0& 
    - force(71)*y1**1*y2**7*y3**0& 
    + force(72)*y1**8*y2**0*y3**0& 
    - force(72)*y1**0*y2**8*y3**0
endif


    f=(v1+v2+v3+v4+v5+v6+v7+v8)
      !
  end function dipol_xy2_p_dump





  double precision function dipol_xy2_p(local,ZPE,N ,force) result(f)
      !
      double precision,intent(in) :: local(1:3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          :: N
      double precision,intent(in) :: force(1:N)
      !
      double precision            :: y1,y2,y3,r1,r2,alpha
      double precision            :: ae,re,v0,v1,v2,v3,v4,v5,v6,v7,v8
      !
      double precision,parameter  :: pi=3.1415926535897932385d0
      !
      re     = force(1)
      ae     = force(2)/180.0d0*pi
      !
      r1 = local(1)
      r2 = local(2)
      alpha = local(3)
      !
      y1 = r1 - re
      y2 = r2 - re
      y3 = cos(alpha) - cos(ae)
      !
 v1 = force( 3)*y1**1*y2**0*y3**0& 
    - force( 3)*y1**0*y2**1*y3**0
 v2 = force( 4)*y1**1*y2**0*y3**1& 
    - force( 4)*y1**0*y2**1*y3**1&  
    + force( 5)*y1**2*y2**0*y3**0& 
    - force( 5)*y1**0*y2**2*y3**0

 v3 = force( 6)*y1**1*y2**0*y3**2& 
    - force( 6)*y1**0*y2**1*y3**2& 
    + force( 7)*y1**2*y2**0*y3**1& 
    - force( 7)*y1**0*y2**2*y3**1& 
    + force( 8)*y1**2*y2**1*y3**0& 
    - force( 8)*y1**1*y2**2*y3**0& 
    + force( 9)*y1**3*y2**0*y3**0& 
    - force( 9)*y1**0*y2**3*y3**0

 if (N> 9) then 
  v4 =force(10)*y1**1*y2**0*y3**3& 
    - force(10)*y1**0*y2**1*y3**3& 
    + force(11)*y1**2*y2**0*y3**2& 
    - force(11)*y1**0*y2**2*y3**2& 
    + force(12)*y1**2*y2**1*y3**1& 
    - force(12)*y1**1*y2**2*y3**1& 
    + force(13)*y1**3*y2**0*y3**1& 
    - force(13)*y1**0*y2**3*y3**1& 
    + force(14)*y1**3*y2**1*y3**0& 
    - force(14)*y1**1*y2**3*y3**0& 
    + force(15)*y1**4*y2**0*y3**0& 
    - force(15)*y1**0*y2**4*y3**0
endif

 if (N>15) then 
  v5 =force(16)*y1**1*y2**0*y3**4& 
    - force(16)*y1**0*y2**1*y3**4& 
    + force(17)*y1**2*y2**0*y3**3& 
    - force(17)*y1**0*y2**2*y3**3& 
    + force(18)*y1**2*y2**1*y3**2& 
    - force(18)*y1**1*y2**2*y3**2& 
    + force(19)*y1**3*y2**0*y3**2& 
    - force(19)*y1**0*y2**3*y3**2& 
    + force(20)*y1**3*y2**1*y3**1& 
    - force(20)*y1**1*y2**3*y3**1& 
    + force(21)*y1**3*y2**2*y3**0& 
    - force(21)*y1**2*y2**3*y3**0& 
    + force(22)*y1**4*y2**0*y3**1& 
    - force(22)*y1**0*y2**4*y3**1& 
    + force(23)*y1**4*y2**1*y3**0& 
    - force(23)*y1**1*y2**4*y3**0& 
    + force(24)*y1**5*y2**0*y3**0& 
    - force(24)*y1**0*y2**5*y3**0
endif

 if (N>24) then 
  v6 =force(25)*y1**1*y2**0*y3**5& 
    - force(25)*y1**0*y2**1*y3**5& 
    + force(26)*y1**2*y2**0*y3**4& 
    - force(26)*y1**0*y2**2*y3**4& 
    + force(27)*y1**2*y2**1*y3**3& 
    - force(27)*y1**1*y2**2*y3**3& 
    + force(28)*y1**3*y2**0*y3**3& 
    - force(28)*y1**0*y2**3*y3**3& 
    + force(29)*y1**3*y2**1*y3**2& 
    - force(29)*y1**1*y2**3*y3**2& 
    + force(30)*y1**3*y2**2*y3**1& 
    - force(30)*y1**2*y2**3*y3**1& 
    + force(31)*y1**4*y2**0*y3**2& 
    - force(31)*y1**0*y2**4*y3**2& 
    + force(32)*y1**4*y2**1*y3**1& 
    - force(32)*y1**1*y2**4*y3**1& 
    + force(33)*y1**4*y2**2*y3**0& 
    - force(33)*y1**2*y2**4*y3**0& 
    + force(34)*y1**5*y2**0*y3**1& 
    - force(34)*y1**0*y2**5*y3**1& 
    + force(35)*y1**5*y2**1*y3**0& 
    - force(35)*y1**1*y2**5*y3**0& 
    + force(36)*y1**6*y2**0*y3**0& 
    - force(36)*y1**0*y2**6*y3**0
 endif

 if (N>36) then 
 v7 = force(37)*y1**1*y2**0*y3**6& 
    - force(37)*y1**0*y2**1*y3**6& 
    + force(38)*y1**2*y2**0*y3**5& 
    - force(38)*y1**0*y2**2*y3**5& 
    + force(39)*y1**2*y2**1*y3**4& 
    - force(39)*y1**1*y2**2*y3**4& 
    + force(40)*y1**3*y2**0*y3**4& 
    - force(40)*y1**0*y2**3*y3**4& 
    + force(41)*y1**3*y2**1*y3**3& 
    - force(41)*y1**1*y2**3*y3**3& 
    + force(42)*y1**3*y2**2*y3**2& 
    - force(42)*y1**2*y2**3*y3**2& 
    + force(43)*y1**4*y2**0*y3**3& 
    - force(43)*y1**0*y2**4*y3**3& 
    + force(44)*y1**4*y2**1*y3**2& 
    - force(44)*y1**1*y2**4*y3**2& 
    + force(45)*y1**4*y2**2*y3**1& 
    - force(45)*y1**2*y2**4*y3**1& 
    + force(46)*y1**4*y2**3*y3**0& 
    - force(46)*y1**3*y2**4*y3**0& 
    + force(47)*y1**5*y2**0*y3**2& 
    - force(47)*y1**0*y2**5*y3**2& 
    + force(48)*y1**5*y2**1*y3**1& 
    - force(48)*y1**1*y2**5*y3**1& 
    + force(49)*y1**5*y2**2*y3**0& 
    - force(49)*y1**2*y2**5*y3**0& 
    + force(50)*y1**6*y2**0*y3**1& 
    - force(50)*y1**0*y2**6*y3**1& 
    + force(51)*y1**6*y2**1*y3**0& 
    - force(51)*y1**1*y2**6*y3**0& 
    + force(52)*y1**7*y2**0*y3**0& 
    - force(52)*y1**0*y2**7*y3**0
 endif

 if (N>52) then 
 v8 = force(53)*y1**1*y2**0*y3**7& 
    - force(53)*y1**0*y2**1*y3**7& 
    + force(54)*y1**2*y2**0*y3**6& 
    - force(54)*y1**0*y2**2*y3**6& 
    + force(55)*y1**2*y2**1*y3**5& 
    - force(55)*y1**1*y2**2*y3**5& 
    + force(56)*y1**3*y2**0*y3**5& 
    - force(56)*y1**0*y2**3*y3**5& 
    + force(57)*y1**3*y2**1*y3**4& 
    - force(57)*y1**1*y2**3*y3**4& 
    + force(58)*y1**3*y2**2*y3**3& 
    - force(58)*y1**2*y2**3*y3**3& 
    + force(59)*y1**4*y2**0*y3**4& 
    - force(59)*y1**0*y2**4*y3**4& 
    + force(60)*y1**4*y2**1*y3**3& 
    - force(60)*y1**1*y2**4*y3**3& 
    + force(61)*y1**4*y2**2*y3**2& 
    - force(61)*y1**2*y2**4*y3**2& 
    + force(62)*y1**4*y2**3*y3**1& 
    - force(62)*y1**3*y2**4*y3**1& 
    + force(63)*y1**5*y2**0*y3**3& 
    - force(63)*y1**0*y2**5*y3**3& 
    + force(64)*y1**5*y2**1*y3**2& 
    - force(64)*y1**1*y2**5*y3**2& 
    + force(65)*y1**5*y2**2*y3**1& 
    - force(65)*y1**2*y2**5*y3**1& 
    + force(66)*y1**5*y2**3*y3**0& 
    - force(66)*y1**3*y2**5*y3**0& 
    + force(67)*y1**6*y2**0*y3**2& 
    - force(67)*y1**0*y2**6*y3**2& 
    + force(68)*y1**6*y2**1*y3**1& 
    - force(68)*y1**1*y2**6*y3**1& 
    + force(69)*y1**6*y2**2*y3**0& 
    - force(69)*y1**2*y2**6*y3**0& 
    + force(70)*y1**7*y2**0*y3**1& 
    - force(70)*y1**0*y2**7*y3**1& 
    + force(71)*y1**7*y2**1*y3**0& 
    - force(71)*y1**1*y2**7*y3**0& 
    + force(72)*y1**8*y2**0*y3**0& 
    - force(72)*y1**0*y2**8*y3**0
endif


    f=(v1+v2+v3+v4+v5+v6+v7+v8)
      !
  end function dipol_xy2_p


  double precision function dipol_xy2_q(local,ZPE,N ,force) result(f)
      !
      double precision,intent(in) :: local(1:3)
      double precision,intent(in) ::  ZPE
      integer,intent(in)          :: N
      double precision,intent(in) :: force(1:N)
      !
      double precision            :: y1,y2,y3,r1,r2,alpha
      double precision            :: ae,re,v0,v1,v2,v3,v4,v5,v6,v7,v8,b0
      !
      double precision,parameter  :: pi=3.1415926535897932385d0
      !
      re     = force(2)
      ae     = force(3)/180.0d0*pi
      !
      b0 = force(4)
      !
      r1 = local(1)
      r2 = local(2)
      alpha = local(3)
      !
      y1 = (r1 - re)*exp(-b0*(r1-re)**2)
      y2 = (r2 - re)*exp(-b0*(r2-re)**2)
      y3 = cos(alpha) - cos(ae)
      !
 v0 = force(1)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1& 
    + force(7)*y1**1*y2**0*y3**0& 
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2& 
    + force(9)*y1**1*y2**0*y3**1& 
    + force(9)*y1**0*y2**1*y3**1& 
    + force(10)*y1**1*y2**1*y3**0& 
    + force(11)*y1**2*y2**0*y3**0& 
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3& 
    + force(13)*y1**1*y2**0*y3**2& 
    + force(13)*y1**0*y2**1*y3**2& 
    + force(14)*y1**1*y2**1*y3**1& 
    + force(15)*y1**2*y2**0*y3**1& 
    + force(15)*y1**0*y2**2*y3**1& 
    + force(16)*y1**2*y2**1*y3**0& 
    + force(16)*y1**1*y2**2*y3**0& 
    + force(17)*y1**3*y2**0*y3**0& 
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then 
  v4 = force(18)*y1**0*y2**0*y3**4& 
    + force(19)*y1**1*y2**0*y3**3& 
    + force(19)*y1**0*y2**1*y3**3& 
    + force(20)*y1**1*y2**1*y3**2& 
    + force(21)*y1**2*y2**0*y3**2& 
    + force(21)*y1**0*y2**2*y3**2& 
    + force(22)*y1**2*y2**1*y3**1& 
    + force(22)*y1**1*y2**2*y3**1& 
    + force(23)*y1**2*y2**2*y3**0& 
    + force(24)*y1**3*y2**0*y3**1& 
    + force(24)*y1**0*y2**3*y3**1& 
    + force(25)*y1**3*y2**1*y3**0& 
    + force(25)*y1**1*y2**3*y3**0& 
    + force(26)*y1**4*y2**0*y3**0& 
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then 
  v5 = force(27)*y1**0*y2**0*y3**5& 
    + force(28)*y1**1*y2**0*y3**4& 
    + force(28)*y1**0*y2**1*y3**4& 
    + force(29)*y1**1*y2**1*y3**3& 
    + force(30)*y1**2*y2**0*y3**3& 
    + force(30)*y1**0*y2**2*y3**3& 
    + force(31)*y1**2*y2**1*y3**2& 
    + force(31)*y1**1*y2**2*y3**2& 
    + force(32)*y1**2*y2**2*y3**1& 
    + force(33)*y1**3*y2**0*y3**2& 
    + force(33)*y1**0*y2**3*y3**2& 
    + force(34)*y1**3*y2**1*y3**1& 
    + force(34)*y1**1*y2**3*y3**1& 
    + force(35)*y1**3*y2**2*y3**0& 
    + force(35)*y1**2*y2**3*y3**0& 
    + force(36)*y1**4*y2**0*y3**1& 
    + force(36)*y1**0*y2**4*y3**1& 
    + force(37)*y1**4*y2**1*y3**0& 
    + force(37)*y1**1*y2**4*y3**0& 
    + force(38)*y1**5*y2**0*y3**0& 
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then 
  v6 = force(39)*y1**0*y2**0*y3**6& 
    + force(40)*y1**1*y2**0*y3**5& 
    + force(40)*y1**0*y2**1*y3**5& 
    + force(41)*y1**1*y2**1*y3**4& 
    + force(42)*y1**2*y2**0*y3**4& 
    + force(42)*y1**0*y2**2*y3**4& 
    + force(43)*y1**2*y2**1*y3**3& 
    + force(43)*y1**1*y2**2*y3**3& 
    + force(44)*y1**2*y2**2*y3**2& 
    + force(45)*y1**3*y2**0*y3**3& 
    + force(45)*y1**0*y2**3*y3**3& 
    + force(46)*y1**3*y2**1*y3**2& 
    + force(46)*y1**1*y2**3*y3**2& 
    + force(47)*y1**3*y2**2*y3**1& 
    + force(47)*y1**2*y2**3*y3**1& 
    + force(48)*y1**3*y2**3*y3**0& 
    + force(49)*y1**4*y2**0*y3**2& 
    + force(49)*y1**0*y2**4*y3**2& 
    + force(50)*y1**4*y2**1*y3**1& 
    + force(50)*y1**1*y2**4*y3**1& 
    + force(51)*y1**4*y2**2*y3**0& 
    + force(51)*y1**2*y2**4*y3**0& 
    + force(52)*y1**5*y2**0*y3**1& 
    + force(52)*y1**0*y2**5*y3**1& 
    + force(53)*y1**5*y2**1*y3**0& 
    + force(53)*y1**1*y2**5*y3**0& 
    + force(54)*y1**6*y2**0*y3**0& 
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then 
 v7 = force(55)*y1**0*y2**0*y3**7& 
    + force(56)*y1**1*y2**0*y3**6& 
    + force(56)*y1**0*y2**1*y3**6& 
    + force(57)*y1**1*y2**1*y3**5& 
    + force(58)*y1**2*y2**0*y3**5& 
    + force(58)*y1**0*y2**2*y3**5& 
    + force(59)*y1**2*y2**1*y3**4& 
    + force(59)*y1**1*y2**2*y3**4& 
    + force(60)*y1**2*y2**2*y3**3& 
    + force(61)*y1**3*y2**0*y3**4& 
    + force(61)*y1**0*y2**3*y3**4& 
    + force(62)*y1**3*y2**1*y3**3& 
    + force(62)*y1**1*y2**3*y3**3& 
    + force(63)*y1**3*y2**2*y3**2& 
    + force(63)*y1**2*y2**3*y3**2& 
    + force(64)*y1**3*y2**3*y3**1& 
    + force(65)*y1**4*y2**0*y3**3& 
    + force(65)*y1**0*y2**4*y3**3& 
    + force(66)*y1**4*y2**1*y3**2& 
    + force(66)*y1**1*y2**4*y3**2& 
    + force(67)*y1**4*y2**2*y3**1& 
    + force(67)*y1**2*y2**4*y3**1& 
    + force(68)*y1**4*y2**3*y3**0& 
    + force(68)*y1**3*y2**4*y3**0& 
    + force(69)*y1**5*y2**0*y3**2& 
    + force(69)*y1**0*y2**5*y3**2& 
    + force(70)*y1**5*y2**1*y3**1& 
    + force(70)*y1**1*y2**5*y3**1& 
    + force(71)*y1**5*y2**2*y3**0& 
    + force(71)*y1**2*y2**5*y3**0& 
    + force(72)*y1**6*y2**0*y3**1& 
    + force(72)*y1**0*y2**6*y3**1& 
    + force(73)*y1**6*y2**1*y3**0& 
    + force(73)*y1**1*y2**6*y3**0& 
    + force(74)*y1**7*y2**0*y3**0& 
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then 
 v8 = force(75)*y1**0*y2**0*y3**8& 
    + force(76)*y1**1*y2**0*y3**7& 
    + force(76)*y1**0*y2**1*y3**7& 
    + force(77)*y1**1*y2**1*y3**6& 
    + force(78)*y1**2*y2**0*y3**6& 
    + force(78)*y1**0*y2**2*y3**6& 
    + force(79)*y1**2*y2**1*y3**5& 
    + force(79)*y1**1*y2**2*y3**5& 
    + force(80)*y1**2*y2**2*y3**4& 
    + force(81)*y1**3*y2**0*y3**5& 
    + force(81)*y1**0*y2**3*y3**5& 
    + force(82)*y1**3*y2**1*y3**4& 
    + force(82)*y1**1*y2**3*y3**4& 
    + force(83)*y1**3*y2**2*y3**3& 
    + force(83)*y1**2*y2**3*y3**3& 
    + force(84)*y1**3*y2**3*y3**2& 
    + force(85)*y1**4*y2**0*y3**4& 
    + force(85)*y1**0*y2**4*y3**4& 
    + force(86)*y1**4*y2**1*y3**3& 
    + force(86)*y1**1*y2**4*y3**3& 
    + force(87)*y1**4*y2**2*y3**2& 
    + force(87)*y1**2*y2**4*y3**2& 
    + force(88)*y1**4*y2**3*y3**1& 
    + force(88)*y1**3*y2**4*y3**1& 
    + force(89)*y1**4*y2**4*y3**0& 
    + force(90)*y1**5*y2**0*y3**3& 
    + force(90)*y1**0*y2**5*y3**3& 
    + force(91)*y1**5*y2**1*y3**2& 
    + force(91)*y1**1*y2**5*y3**2& 
    + force(92)*y1**5*y2**2*y3**1& 
    + force(92)*y1**2*y2**5*y3**1& 
    + force(93)*y1**5*y2**3*y3**0& 
    + force(93)*y1**3*y2**5*y3**0& 
    + force(94)*y1**6*y2**0*y3**2& 
    + force(94)*y1**0*y2**6*y3**2& 
    + force(95)*y1**6*y2**1*y3**1& 
    + force(95)*y1**1*y2**6*y3**1& 
    + force(96)*y1**6*y2**2*y3**0& 
    + force(96)*y1**2*y2**6*y3**0& 
    + force(97)*y1**7*y2**0*y3**1& 
    + force(97)*y1**0*y2**7*y3**1& 
    + force(98)*y1**7*y2**1*y3**0& 
    + force(98)*y1**1*y2**7*y3**0& 
    + force(99)*y1**8*y2**0*y3**0& 
    + force(99)*y1**0*y2**8*y3**0
endif


    f=(v0+v1+v2+v3+v4+v5+v6+v7+v8)*sin(alpha)
      !
  end function dipol_xy2_q




!
!   Skip n lines  in the input file 
!
  subroutine skiplines( inpunit,n )
  integer,intent(in) :: n,inpunit
  character(len=80)  :: label
  integer            :: i0

    do i0=1,n
       read  (inpunit,"(a80)") label 
    enddo

  end subroutine skiplines


  subroutine linur(dimen,npar,coeff,constant,solution,error)

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur     

!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  integer,intent(in)           :: npar,dimen
  double precision,intent(in)  :: al(npar,npar)
  double precision,intent(out) :: ai(npar,npar)
  double precision             :: h(npar),p,q
  integer                      :: i1,i2,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine invmat
!------------------------------------------!


! subroutine allocate_array(name,size,string)
!  integer,intent(in) :: size
!  character(len=*),intent(in) :: string
!  double precision,intent(in,out) :: name(:) 
!  integer :: alloc
!
!   if (.not.allocated(name)) then
!   allocate (name(size),stat=alloc)
!   if (alloc/=0) then 
!    write(6,"('out of memory - ',A,' size = ',I)") trim(string),size
! stop
!   endif 
!   else 
!     write(6,"('It has been allocated :',A)") trim(string)
!     stop
!   endif 
!  end subroutine allocate_array


  end module xy2_globfit

  program driver
    use xy2_globfit
    character(len=20)  :: ch_t
    integer :: i

    do i=1,9
      read  (5,"(a20)") ch_t
    enddo
    !
    read  (5,"(a20)") ch_t

   select case (trim(ch_t))
     case default
       write (6,"('uknown potential function ',a20)") trim(ch_t)
       stop 'uknown potential function '
     case ('morse_tyuterev','MORSE_TYUETEREV')
       call fitting(morse_tyuterev)
     case ('morse_tyuterev_II','MORSE_TYUETEREV_II')
       call fitting(morse_tyuterev_II)
     case ('morse_tyuterev_damp')
       call fitting(potential_xy2_tyuterev_damp)
     case ('morse_242','MORSE_242')
       call fitting(morse_242)
     case ('wifin_242','WIFIN_242')
       call fitting(wifin_242)
     case ('MORBID','morbid')
       call fitting(morbid)
     case ('MORBID_MEP','morbid_mep')
       call fitting(morbid_mep)
     case ('MORSE','morse')
       call fitting(morse)
     case ('TAYLOR_I','taylor_I')
       call fitting(taylor_i)
     case ('r_r12_dr','R_R12_DR','R_r12_dr')
       call fitting(R_r12_dr_morse_pol)
     case ('r_r12_dr_x')
       call fitting(R_r12_dr_morse_pol_x)
     case ('MORBID_MEP2','morbid_mep2')
       call fitting(morbid_mep2)
     case ('MORBID_MEP3','morbid_mep3')
       call fitting(morbid_mep3)
     case ('MORBID_MEP4','morbid_mep4')
       call fitting(morbid_mep4)
     case ('V_REP_DISP','V_rep_disp')
       call fitting(V_rep_disp)
     case ('R12_MORBID_SWITCH','r12_morbid_switch')
       call fitting(R_r12_morbid_switch)
     case ('dipol_xy2_p')
       call fitting(dipol_xy2_p)
     case ('dipol_xy2_q')
       call fitting(dipol_xy2_q)
     case ('dipol_xy2_p_dump')
       call fitting(dipol_xy2_p_dump)
   end select
    !
    
  end program driver
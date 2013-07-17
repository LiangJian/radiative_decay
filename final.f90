! ftn -fastsse -Mbyteswapio -o final res_func.f90 res_chi.f90 final.f90

      	program three_point_funcion

	implicit none

        integer, parameter :: nt=96, mt=35, nt1=nt/2+1
        integer, parameter :: lbin = 50, nbin = 100
	integer, parameter :: nfac=10,nin=2,nq=1,nout=2,nom=3,nmax=100
        integer, parameter :: nf=3,np=3,nm=8

        integer iq,i,ii,it,ic,itmp,mfac,itq,im,i_i,i_q,i_f

        real*8  factor(7,3),dmass(27,2),mass1mm,atm1
	complex*16 Goper(nt,27,3,7,lbin)
        complex*16 corre(nt,12,27,lbin)
        complex*16 tco(nt1,27,3,4,7,nbin)
        complex*16 tcot(nt1,27,3,4,7,nbin)
        complex*16 tcoaver(nt1,27,3,4,7)
	complex*16 psi(nt1,27,3)
        complex*16 psit(nt1,27,3)
        complex*16 psiaver(nt1,27,3)
        real*8 C2new(nt,3,3,7,27)
        real*8 C2g(nt1,7,27),C2gaver(nt1,7,27),C2gt(nt1,7,27)
	real*8 sigma(2,nt1,27,3,4,7)
	real*8 parameter(27,3,4,7,10,5)

        complex*16 ressort(nfac,(2*nin+1)*(2*nq+2)*(2*nout+1)*nom*nom*nom)
        integer resind(9,(2*nin+1)*(2*nq+2)*(2*nout+1)*nom*nom*nom)
	integer in,out,dp,flag,ns,indall
	integer ind(10),imom(3,10)
        real*8  m1,m2,p(4,2),pq(4),qsq(10)
	complex*16 tconew(nbin,nt1,nmax,10)
	complex*16 tcoavernew(nt1,nmax,10)
	real*8 sigmanew(2,nt1,nmax,10)
	complex*16 coef(nfac,nmax,10)
	complex*16 resfac(nfac,nbin,nt1,10)
	complex*16 facaver(nfac,nt1,10)
	real*8 facsig(2,nfac,nt1,10)

	real*8 data(nbin,nm*nf),Qsqnew(nm*nf),dyda(nm*nf,np*nf)
	real*8 para(np*nf),covar(np*nf,np*nf),da(np*nf)
	real*8 chisq,chisqT,chisqold,dchisq,alamda

	real*8 dataaver(nm*nf),invcov(np*nf,np*nf)
	real*8 covmat(nm*nf,nm*nf),invmat(nm*nf,nm*nf)
	integer j,m,n,ia,ja
	real*8 errbar(7,1000)

	integer iter,iflag,istart(10)

         open(12, file='res_amp.dat', form='unformatted', &
                  access='direct', recl=16*nt1*27*3*4*7*nbin)

        read(12,rec=1)                                                 &
               ((((((tco(itq,im,i_i,i_q,i_f,ic),itq=1,nt1)          &
               ,im=1,27),i_i=1,3),i_q=1,4),i_f=1,7),ic=1,nbin)
       close(12)

!-------in is angular momentum of initial state
        in=1
!-------out is angular momentum of final state
        out=2
!-------dp is product of P patity of initial and final states.
        dp=-1
!-------if flag=1,initial state is still. if flag=2, final one is still.
        flag=2
!-------ns is the space lat number
        ns=8

!-------mass of initial state
        m1=0.696
!-------mass of final state
        m2=0.535

        call FacCG(in,out,dp,m1,m2,ns,flag,indall,ressort,resind)

!----------------------------------
	atm1=4.44

      open(66,file='new_result.2pp.G')

	p=0.0d0
        write(66,"(20A)") '----------------------'
        write(66,"(A,f9.5)") '   Q2=',-(m1-m2)**2

	iq=1
	ind(iq)=1
	tconew=0.0d0
	tcoavernew=0.0d0
	sigmanew=0.0d0
	p=0.0d0
	do i=1,3
	imom(i,1)=0
	end do
	p(4,1)=m1
	p(4,2)=m2
        do i=1,4
        pq(i)=p(i,1)-p(i,2)
        end do
	qsq(1)=-(m1-m2)**2	

	do ii=1,indall
        if(ii.gt.1) then
        if(resind(7,ii).gt.resind(7,ii-1)) then
        p(1,1)=0.4d0*dsin(3.14159265358979323d0*(resind(4,ii)-1)/ns)
        p(2,1)=0.4d0*dsin(3.14159265358979323d0*(resind(5,ii)-1)/ns)
        p(3,1)=0.4d0*dsin(3.14159265358979323d0*(resind(6,ii)-1)/ns)
        p(4,1)=sqrt(m1**2+p(1,1)**2+p(2,1)**2+p(3,1)**2)
	p(4,2)=m2
        do i=1,4
        pq(i)=p(i,1)-p(i,2)
        end do
        write(66,"(20A)") '----------------------'
        write(66,"(A,f9.5)") '   Q2=',  &
                (-pq(4)**2+pq(1)**2+pq(2)**2+pq(3)**2)
	ind(iq)=ind(iq)-1
	iq=iq+1
	ind(iq)=1

	do i=1,3
	imom(i,iq)=resind(3+i,ii)-1
	end do
	qsq(iq)=-pq(4)**2+pq(1)**2+pq(2)**2+pq(3)**2

        end if
        end if

	im=resind(4,ii)+3*(resind(5,ii)-1)+9*(resind(6,ii)-1)

	if(resind(2,ii).ne.4) then

	do it=1,nt1
	do ic=1,nbin
	tconew(ic,it,ind(iq),iq)= tconew(ic,it,ind(iq),iq)+   &
		tco(it,im,resind(1,ii),resind(2,ii),resind(3,ii)+2,ic)	&
		/resind(9,ii)
	end do
	end do

        write(66,"(6I3,,)") &
                resind(4,ii)-1,resind(5,ii)-1,resind(6,ii)-1,   &
                resind(1,ii),resind(2,ii),resind(3,ii)

        if(resind(8,ii).ne.resind(8,ii+1)) then

        write(66,"(A,f9.5)") '   Q2=',  &
                -pq(4)**2+pq(1)**2+pq(2)**2+pq(3)**2

        do it=1,nt1

        do ic=1,nbin
        tcoavernew(it,ind(iq),iq)=tcoavernew(it,ind(iq),iq)           &
                +tconew(ic,it,ind(iq),iq)/nbin
        end do


        do ic=1,nbin
          sigmanew(1,it,ind(iq),iq)=sigmanew(1,it,ind(iq),iq)         &
        	+dreal(tconew(ic,it,ind(iq),iq))**2		&
		-dreal(tcoavernew(it,ind(iq),iq))**2
          sigmanew(2,it,ind(iq),iq)=sigmanew(2,it,ind(iq),iq)         &
        	+dimag(tconew(ic,it,ind(iq),iq))**2		&
		-dimag(tcoavernew(it,ind(iq),iq))**2
        end do
            sigmanew(1,it,ind(iq),iq)=sqrt(sigmanew(1,it,ind(iq),iq)	&
             *(nbin-1)/nbin)
            sigmanew(2,it,ind(iq),iq)=sqrt(sigmanew(2,it,ind(iq),iq)	&
             *(nbin-1)/nbin)

        end do

	itmp=ii-abs(resind(9,ii))+1
	do i=1,nfac
		coef(i,ind(iq),iq)=ressort(i,itmp)
	end do

       write(66,"(3I3,5f10.5,4e13.5)")(         &
		iq,ind(iq),it,		&
                real(coef(1,ind(iq),iq)),real(coef(2,ind(iq),iq)),&
                real(coef(3,ind(iq),iq)),real(coef(4,ind(iq),iq)),&
                real(coef(5,ind(iq),iq)),		&
		real(tcoavernew(it,ind(iq),iq)),sigmanew(1,it,ind(iq),iq),&
                aimag(tcoavernew(it,ind(iq),iq)),sigmanew(2,it,ind(iq),iq),&
                it=1,nt1)

	ind(iq)=ind(iq)+1

	end if

	end if

	end do

	ind(iq)=ind(iq)-1

	close(66)
	
	if(iq.eq.0) then

	open(67,file='fac_result.2pp.G')

	mfac=5

	do i=1,ind(1)
	coef(1,i,1)=coef(1,i,1)+coef(4,i,1)
	end do

	resfac=0.0d0

	do iq=1,10

	if(iq.gt.1) then
	call qsqmin(coef(1,1,iq),tconew(1,1,1,iq),tcoavernew(1,1,iq),	&
		ind(iq),mfac,resfac(1,1,1,iq))
	else
        call qsqmin(coef(1,1,iq),tconew(1,1,1,iq),tcoavernew(1,1,iq),   &
                ind(iq),1,resfac(1,1,1,iq))
	end if

        do it=1,nt1
	do i=1,mfac

	facaver(i,it,iq)=0.0d0

        do ic=1,nbin
        facaver(i,it,iq)=facaver(i,it,iq)           &
                +resfac(i,ic,it,iq)/nbin
        end do

        do ic=1,nbin
          facsig(1,i,it,iq)=facsig(1,i,it,iq)         &
                +dreal(resfac(i,ic,it,iq))**2             &
                -dreal(facaver(i,it,iq))**2
          facsig(2,i,it,iq)=facsig(2,i,it,iq)         &
                +dimag(resfac(i,ic,it,iq))**2             &
                -dimag(facaver(i,it,iq))**2
        end do
            facsig(1,i,it,iq)=sqrt(facsig(1,i,it,iq)    &
             *(nbin-1)/nbin)
            facsig(2,i,it,iq)=sqrt(facsig(2,i,it,iq)    &
             *(nbin-1)/nbin)

        end do
	end do

              write(67,"(4I3,e13.5,A,2e13.5,A,2e13.5,A,2e13.5,A,2e13.5,A,&
			2e13.5,A,)")(     &
                imom(1,iq),imom(2,iq),imom(3,iq),it,qsq(iq)*atm1**2,'  ',&
                real(facaver(1,it,iq)),facsig(1,1,it,iq),'  ',&
                real(facaver(2,it,iq)),facsig(1,2,it,iq),'  ',&
                real(facaver(3,it,iq)),facsig(1,3,it,iq),'  ',&
                real(facaver(4,it,iq)),facsig(1,4,it,iq),'  ',&
                real(facaver(5,it,iq)),facsig(1,5,it,iq),'  ',&
                it=1,nt1)

	end do
	close(67)

         open(15, file='res_fac.dat', form='unformatted', &
                  access='direct', recl=16*nt1*5*10*nbin)

        write(15,rec=1)                                                 &
               ((((resfac(i_i,ic,itq,im),i_i=1,5),ic=1,nbin),itq=1,nt1) &
               ,im=1,10)
       close(15)

	end if	

         open(15, file='res_fac.dat', form='unformatted', &
                  access='direct', recl=16*5*nt1*10*nbin)

        read(15,rec=1)                                                 &
               ((((resfac(i_i,ic,itq,im),i_i=1,5),ic=1,nbin),itq=1,nt1)&
               ,im=1,10)
       close(15)

	data=0.0d0
	dataaver=0.0d0
	istart(1)=20
	istart(2)=20
        istart(3)=20
        istart(5)=16
        istart(6)=16
        istart(7)=16
        istart(8)=16
        istart(9)=19
	ii=1
	do iq=1,10
	if(iq.ne.4.and.iq.ne.10) then
	do i=1,nf
	itmp=nf*ii+i-nf
	qsqnew(itmp)=qsq(iq)*atm1**2

	do it=istart(iq),istart(iq)+4
	do ic=1,nbin
	data(ic,itmp)=data(ic,itmp)+dreal(resfac(i,ic,it,iq))/5
	dataaver(itmp)=dataaver(itmp)+dreal(resfac(i,ic,it,iq))/500
	end do
	end do

	end do

	ii=ii+1
	end if
	end do

!----------------------------method 1
        dyda=0.0d0
        do i=1,nm*nf
        j=np*mod(i-1,nf)+1
        dyda(i,j)=1
        dyda(i,j+1)=Qsqnew(i)
        dyda(i,j+2)=Qsqnew(i)*Qsqnew(i)
        end do

          call mkcovmat_r(data,dataaver,nm*nf,covmat)
	  covmat(2,2)=covmat(2,2)+0.0000001d0
	  covmat(3,3)=covmat(3,3)+0.0000001d0
          call inverse(covmat,invmat,nm*nf)

        ii=1
        do iq=1,10
        if(iq.ne.4.and.iq.ne.10) then
        do i=1,nf
        itmp=nf*ii+i-nf
	write(*,'(i3,3e13.5)') i,qsqnew(itmp),dataaver(itmp),	&
		sqrt(covmat(itmp,itmp))
        end do

        ii=ii+1
        end if
        end do

          covar=0.0d0
	  para =0.0d0
          do i=1,np*nf
          do j=1,np*nf
          do m=1,nm*nf
          do n=1,nm*nf
              covar(i,j)=covar(i,j)+invmat(m,n)*dyda(m,i)*dyda(n,j)
          end do
          end do
          end do
          end do

          call inverse(covar,invcov,np*nf)

          do j=1,np*nf
          do i=1,np*nf
          do m=1,nm*nf
          do n=1,nm*nf
              para(i)=para(i)+invcov(i,j)	&
                *invmat(m,n)*dataaver(n)*dyda(m,j)
          end do
          end do
          end do
          end do

	do i=1,nf
	j=3*i-2
	write(*,'(3f10.5)') para(j),para(j+1),para(j+2)
	end do

!----------------------------method 2
      	alamda = -1.0
	para  = 0.0d0
	covar = 0.0d0
     	chisq = 0.0d0
     	chisqT= 0.0d0

   	call mrqmin(qsqnew,data,nf*nm,para,covar,da,nf*np,chisq,chisqT,alamda)
      	iflag = 0
      	do iter=1,200
        	chisqold = chisq
		call mrqmin(qsqnew,data,nf*nm,para,covar,da,nf*np,chisq,&
			chisqT,alamda)
!		print*,chisq
        	dchisq = chisqold - chisq
        	if( dchisq .lt. 0.002)then
          		iflag = iflag + 1
        	else
          		iflag = 0
        	endif
        	if(iflag .ge. 4 .and. iter .ge. 10)go to 99
      	enddo
 99   	continue
!        call mrqmin(qsqnew,data,nf*nm,para,covar,da,nf*np,chisq,chisqT,alamda)

	write(*,*) chisq/(nf*nm-nf*np),alamda

          call inverse(covar,invcov,np*nf)

        do j=1,nf*np
        write(*,'(2f10.5)') para(j),sqrt(invcov(j,j))
        end do

	ii=1000
	do i=1,ii
	errbar(1,i)=-1.0d0+i*(2.5d0+1.0d0)/ii
	call func(errbar(1,i),para,invcov)
	end do

      	open(69,file='res_err.2pp.G')
              write(69,"(I5,e13.5,A,2e13.5,A,2e13.5,A,2e13.5)")( &
                i,errbar(1,i),' ',&
                errbar(2,i),errbar(3,i),' ',&
                errbar(4,i),errbar(5,i),' ',&
                errbar(6,i),errbar(7,i),&
                i=1,ii)
	
	close(69)

        end
!=======================================================================

        subroutine func(errbar,para,invcov)

        integer, parameter :: nf=3,np=3,nm=8
        real*8 qtmp,para(np*nf),invcov(np*nf,np*nf)
        real*8 errbar(7),dyda(3,3)

        integer i,j,k

	qtmp=errbar(1)

	do i=1,3
        errbar(2*i)=para(3*i-2)+para(3*i-1)*qtmp+para(3*i)*qtmp*qtmp
        dyda(i,1)=1
        dyda(i,2)=qtmp
        dyda(i,3)=qtmp*qtmp
        end do

	do k=1,3
	errbar(2*k+1)=0.0d0
	do j=1,3
	do i=1,3
	errbar(2*k+1)=errbar(2*k+1)+dyda(k,i)*dyda(k,j)	&
		*invcov(3*k-3+i,3*k-3+j)
	end do
	end do
	errbar(2*k+1)=sqrt(errbar(2*k+1))
	end do

!	 write(*,'(3f10.5)') errbar(1),errbar(2),errbar(3)

      return
	end 

!=======================================================================

        subroutine funcs(Qsq,nd,para,na,y,dyda)

        integer, parameter :: nf=3,np=3,nm=8
!	In fact, nd=nf*nm,na=nf*np

        integer nd,na
        real*8 Qsq(nd),y(nd)
        real*8 para(na),dyda(nd,na)
        real*8 chisq

	integer i,j,k

	dyda=0.0d0
	do i=1,nd
	j=np*mod(i-1,nf)+1
	y(i)=para(j)+para(j+1)*Qsq(i)+para(j+2)*Qsq(i)*Qsq(i)
	dyda(i,j)=1
	dyda(i,j+1)=Qsq(i)
	dyda(i,j+2)=Qsq(i)*Qsq(i)
!	y(i)=para(j)*(1+para(j+1)*Qsq(i))*exp(-para(j+2)*Qsq(i))
!	dyda(i,j)=(1+para(j+1)*Qsq(i))*exp(-para(j+2)*Qsq(i))
!	dyda(i,j+1)=para(j)*Qsq(i)*exp(-para(j+2)*Qsq(i))
!	dyda(i,j+2)=-para(j)*(1+para(j+1)*Qsq(i))*exp(-para(j+2)*Qsq(i))*Qsq(i)
	end do

      return
        end

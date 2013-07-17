!=======================================================================
!	Numerical recipes, C. 15.5.2
!	fit data use Levenberg-Marquardt method.
!	data is the sample of data used to be fitted.
!	ma is the modes overall of data
!	para is the parameters fitted by this subroutine
!	covar is the covariance matrix of fit parameter
!	na is the number of parameters
!	chisq is the final chi**2
!	funcs is the function calc fitting funciton and its derviation
!	alamda is the parameter to determin how to fit

      subroutine mrqmin(Qsq,data,nd,a,alpha,beta,na,		&
			chisq,chisqT,alamda)

	integer, parameter :: nbin=100

!        real*8, external :: funcs

	integer nd,na
	real*8 data(nbin,nd),Qsq(nd)
	real*8 a(na),beta(na),alpha(na,na)
	real*8 chisq,alamda

	real*8 aver(nd),covmat(nd,nd),invmat(nd,nd)
	real*8 alphaT(na,na),betaT(na),invalpha(na,na),atry(na)
	integer i,j

	do j=1,nd
	aver(j)=0.0d0
	do i=1,nbin
	  aver(j)=aver(j)+data(i,j)
	end do
	aver(j)=aver(j)/nbin
	end do
	
          call mkcovmat_r(data,aver,nd,covmat)
          covmat(2,2)=covmat(2,2)+0.0000001d0
          covmat(3,3)=covmat(3,3)+0.0000001d0
          call inverse(covmat,invmat,nd)

	if(alamda.lt.0.0d0) then
		call mrqcof(aver,Qsq,invmat,nd,a,na,alpha,beta,chisq)
		alamda=0.001d0
		chisqT=chisq
	end if

	do j=1,na
	do i=1,na
		alphaT(i,j)=alpha(i,j)
	end do
	alphaT(j,j)=alpha(j,j)*(1.0d0 + alamda)
	end do
!	if(alamda .eq. 0.0) return
        call inverse(alphaT,invalpha,na)
        do i=1,na
        atry(i)=a(i)
	do j=1,na
              atry(i)=atry(i)+invalpha(i,j)*beta(j)
        end do
        end do	

	call mrqcof(aver,Qsq,invmat,nd,atry,na,alphaT,betaT,chisq)
	if(chisq.lt.chisqT) then
	alamda=0.1*alamda
	chisqT=chisq
	alpha=alphaT
	beta=betaT
	a=atry
	else
	alamda=alamda*10
	chisq=chisqT
	end if

      return
        end


!=======================================================================
!       Numerical recipes, C. 15.5.2
!	function needed by mrqmin, calc alpha, beta and chisq
!	alpha(x,y) is dfdx(i).invmat(i,j).dfdy(j)
!	beta(x) is df(i).invmat(i,j).dfdx(j)

      subroutine mrqcof(aver,Qsq,invmat,nd,para,na,alpha,beta,chisq)

        integer, parameter :: nbin=100

	integer nd,na
        real*8 aver(nd),invmat(nd,nd),Qsq(nd)
	real*8 para(na),beta(na),alpha(na,na),dyda(nd,na)
	real*8 chisq

	integer i,j,ia,ja
	real*8 dy(nd),wt

	alpha=0.0d0
	beta=0.0d0
	chisq=0.0d0

	call funcs(Qsq,nd,para,na,dy,dyda)
	do i=1,nd
		dy(i)=dy(i)-aver(i)
	end do

	do j=1,nd
	do i=1,nd
	do ja=1,na
	wt=dyda(j,ja)*invmat(i,j)
	do ia=ja,na
		alpha(ia,ja)=alpha(ia,ja)+wt*dyda(i,ia)
	end do
	beta(ja)=beta(ja)-wt*dy(i)
	end do
	chisq=chisq+dy(i)*invmat(i,j)*dy(j)
        end do
        end do

        do ja=2,na
        do ia=1,ja-1
                alpha(ia,ja)=alpha(ja,ia)
        end do
        end do

      return
        end
		
!=======================================================================
!	Numerical recipes, C. 15.2 
!	fit data use inverse of derviation matrix of parameters.
!	coef is the coefficent list of form factors
!	tconew is the data to be fit
!	tcoavernew if the average of dat
!	ind(i) is the number of independence channel of monentum mode 'i'
!	mfac is the number of formfactors
!	resfac return formfactors calc here.

      subroutine qsqmin(coef,tconew,tcoavernew,ind,mfac,resfac)

        integer, parameter :: nbin=100,nt=96,nt1=nt/2+1,nfac=10

        complex*16 tconew(nbin,nt1,ind)
        complex*16 tcoavernew(nt1,ind)
        complex*16 coef(nfac,ind)
        complex*16 resfac(nfac,nbin,nt1)
	integer ind,mfac
	integer it,id,ic,i,j,k,m,n

	real*8 covmat(ind,ind,nt1),invmat(ind,ind,nt1)
	complex*16 data(nbin,ind,nt1),aver(ind,nt1)
	
	real*8 alpha(mfac,mfac),invalpha(mfac,mfac)

	print*,ind

	  resfac=0.0d0
	
	do it=1,nt1
	  do id=1,ind
	  do ic=1,nbin
	    data(ic,id,it)=tconew(ic,it,id)
	  end do
	  aver(id,it)=tcoavernew(it,id)
	  end do

	  call mkcovmat(data(1,1,it),aver(1,it),ind,covmat(1,1,it))
	  call inverse(covmat(1,1,it),invmat(1,1,it),ind)


	  alpha=0.0d0
	  do i=1,mfac
	  do j=1,mfac
	  do m=1,ind
	  do n=1,ind
	      alpha(i,j)=alpha(i,j)+	&
	  	dreal(coef(i,m)*invmat(m,n,it)*dconjg(coef(j,n)))
	  end do
	  end do
	  end do
	  end do

	  call inverse(alpha,invalpha,mfac)


	  do ic=1,nbin
          do i=1,mfac
          do j=1,mfac
          do m=1,ind
          do n=1,ind
	      resfac(i,ic,it)=resfac(i,ic,it)+invalpha(i,j)*coef(j,m) &
		*invmat(m,n,it)*data(ic,n,it)
          end do
          end do
          end do
	  end do
	  end do


	end do


      return
	end 

!=======================================================================
!       calc covariance matrix of data

        subroutine mkcovmat_r(data,aver,dim,covmat)

        integer, parameter :: nbin=100

        integer dim
        real*8 covmat(dim,dim)
        real*8 data(nbin,dim),aver(dim)
        integer i,j,ic

        covmat=0.0d0
        do i=1,dim
        do j=1,dim
        do ic=1,nbin
           covmat(i,j)=covmat(i,j)+     &
                (data(ic,i)-aver(i))*(data(ic,j)-aver(j))
        end do
        covmat(i,j)=covmat(i,j)*(nbin-1)/nbin
        end do
        end do

      return
        end


!=======================================================================
!	calc covariance matrix of data

        subroutine mkcovmat(data,aver,dim,covmat)

        integer, parameter :: nbin=100

	integer dim
        real*8 covmat(dim,dim)
        complex*16 data(nbin,dim),aver(dim)
	integer i,j,ic

	covmat=0.0d0
	do i=1,dim
	do j=1,dim
	do ic=1,nbin
	   covmat(i,j)=covmat(i,j)+	&
		dreal((data(ic,i)-aver(i))*dconjg(data(ic,j)-aver(j)))
	end do
	covmat(i,j)=covmat(i,j)*(nbin-1)/nbin
	end do
	end do

      return
        end

!=======================================================================
!	inverse matrix 'covmat' of rank 'dim' into 'invmat'

	subroutine inverse(covmat,invmat,dim)

        integer dim
        real*8 covmat(dim,dim),invmat(dim,dim)
	real*8 temp(dim,dim),w(dim),v(dim,dim)

	integer i,j,k,ncut
	real*8 wmax,wmin,epsl,wcut
	real*8 sum

	epsl=1.0d-8
	
	temp=covmat
	call dsvdcmp(temp,dim,dim,w,v)

	wmax= w(1)
	wmin= dabs(w(1))
      	do j=2,dim
        if(w(j) .gt. wmax)wmax = w(j)
        if(dabs(w(j)) .lt. wmin)wmin = dabs(w(j))
      	end do
      	wcut = wmax*eps1 
	ncut=0
      	do j=1,dim
        if(dabs(w(j)) .lt. wcut)w(j) = 0.0d0
        if(dabs(w(j)) .lt. wcut)ncut = ncut + 1
      	end do
	
      	do i=1,dim
        do j=1,dim
          invmat(i,j)=0.0d0
          do k=1,dim
            if(w(k) .ne. 0.0d0)then
              invmat(i,j) = invmat(i,j) + v(i,k)*temp(j,k)/w(k)
            endif
          end do
        end do
      	end do

!	sum=0.0d0
!        do i=1,dim
!        do j=1,dim
!              sum=sum+invmat(i,j)*covmat(j,i)
!        end do
!        end do
!
!	print*,sum/dim

      return
        end	

!=======================================================================
! Numerical recipes, P. 59
! Singular value decomposition of a matrix, used by inverse
!	a is matrix to be decomposition
!	m/n is the rank of a
!	w and v is the parameter of decomposition.
!	at last, a_ori-> v(i,k).w(k).u(j,k)(matrix 'u' storage in 'a')
!	inverse(a_ori)->v(i,k).inverse(w(k)).u(j,k)

      subroutine dsvdcmp(a,m,n,w,v)

        real*8 , external :: pythagd

      real*8 a(m,n),v(n,n),w(n)
      integer m,n

      integer i,j,k,l,its,jj,nm
      real*8 scale,anorm,s,f,g,h,c,x,y,z
      real*8 rv1(n)

      rv1 = 0.0d0
      v = 0.0d0
      w = 0.0d0

      g=0.0d0
      scale=0.0d0
      anorm=0.0d0

      do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i .le. m)then
           do k=i,m
             scale = scale + dabs(a(k,i))
           end do
           if(scale .ne. 0.0d0)then
             do k=i,m
               a(k,i) = a(k,i)/scale
               s = s +a(k,i)*a(k,i)
             end do
             f = a(i,i)
             g = -dsign(dsqrt(s),f)
             h = f*g - s
             a(i,i) = f - g
             do j=l,n
               s=0.0d0
               do k=i,m
                 s = s +a(k,i)*a(k,j)
               end do
               f = s/h
               do k=i,m
                 a(k,j) = a(k,j) + f*a(k,i)
               end do
             end do
             do k=i,m
               a(k,i) = scale*a(k,i)
             end do
          endif
        endif
        w(i) = scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0

        if( (i.le.m) .and. (i.ne.n))then
          do k=l,n
            scale = scale + dabs(a(i,k))
          enddo
          if(scale .ne. 0.0d0)then
            do k=l,n
              a(i,k) = a(i,k)/scale
              s = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -dsign(dsqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k=l,n
              rv1(k) = a(i,k)/h
            end do
            do j=l,m
              s=0.0d0
              do k=l,n
                s = s + a(j,k)*a(i,k)
              end do
              do k=l,n
                a(j,k) = a(j,k) + s*rv1(k)
              end do
            end do
            do k=l,n
              a(i,k) = scale*a(i,k)
            end do
          endif
        endif
        anorm = max(anorm, (dabs(w(i))+dabs(rv1(i))))
      end do

      do i=n,1,-1
        if(i .lt. n)then
          if(g .ne. 0.0d0)then
            do j=l,n
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j=l,n
              s=0.0d0
              do k=l,n
                s = s + a(i,k)*v(k,j)
              end do
              do k=l,n
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do
          endif
          do j = l,n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
          end do
        endif
        v(i,i) = 1.0d0
        g = rv1(i)
        l = i
      end do

      do i=min(m,n),1,-1
        l = i + 1
        g = w(i)
        do j=l,n
          a(i,j) = 0.0d0
        end do
        if(g .ne. 0.0d0)then
          g = 1.0d0/g
          do j=l,n
            s=0.0d0
            do k=l,m
              s = s + a(k,i)*a(k,j)
            end do
            f = (s/a(i,i))*g
            do k=i,m
              a(k,j) = a(k,j) + f*a(k,i)
            end do
          end do
          do j=i,m
            a(j,i) = a(j,i)*g
          end do
        else
          do j=i,m
            a(j,i) = 0.0d0
          end do
        endif
        a(i,i) = a(i,i) + 1.0d0
      end do

      do k=n,1,-1
        do its=1,50
          do l=k,1,-1
            nm = l - 1
            if((dabs(rv1(l))+anorm) .eq. anorm)goto 2
            if((dabs(w(nm))+anorm) .eq. anorm)goto 1
          end do
 1        c=0.0d0
          s=1.0d0
          do i=l,k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if((dabs(f)+anorm) .eq. anorm)goto 2
            g = w(i)
            h = pythagd(f,g)
            w(i) = h
            h = 1.0d0/h
            c = (g*h)
            s = -(f*h)
            do j=1,m
              y = a(j,nm)
              z = a(j,i)
              a(j,nm) =  (y*c) + (z*s)
              a(j,i) = -(y*s) + (z*c)
            end do
          end do

 2        z = w(k)

          if(l .eq. k)then
            if(z .lt. 0.0d0)then
              w(k) = -z
              do j=1,n
                v(j,k) = -v(j,k)
              end do
            endif
            goto 3
          endif

          if(its .eq. 50)then
	 print *, 'no convergence in dsvdcmp'
          return
          endif

          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.0d0*h*y)
          g = pythagd(f,1.0d0)
          f = ((x-z)*(x+z) + h*((y/(f+dsign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do j=l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = pythagd(f,h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj=1,n
              x = v(jj,j)
              z = v(jj,i)
              v(jj,j) = (x*c) + (z*s)
              v(jj,i) = -(x*s) + (z*c)
            end do
            z = pythagd(f,h)
            w(j) = z
            if(z .ne. 0.0d0)then
              z = 1.0d0/z
              c = f*z
              s = h*z
            endif
            f = (c*g) + (s*y)
            x = -(s*g) + (c*y)
            do jj=1,m
              y = a(jj,j)
              z = a(jj,i)
              a(jj,j) = (y*c) + (z*s)
              a(jj,i) = -(y*s) + (z*c)
            end do
          end do
          rv1(l) = 0.0d0
          rv1(k) = f
          w(k) = x
        end do
 3      continue
      end do

      return
      end

!=======================================================================
!	function used by Singular value decomposition, dsvdcmp

      function pythagd(a,b)
      implicit real*8(a-h,o-z)
      real*8 a,b,pythagd
      real*8 absa,absb

      absa = dabs(a)
      absb = dabs(b)

      if(absa .gt. absb)then
        pythagd = absa*dsqrt(1.0d0 + (absb/absa)**2)
      else
        if(absb .eq. 0.0d0)then
          pythagd = 0.0d0
        else
          pythagd = absb*dsqrt(1.0d0 + (absa/absb)**2)
        endif
      endif

      return
      end
	

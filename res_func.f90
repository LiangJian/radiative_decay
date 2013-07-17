! ftn -fastsse -Mbyteswapio -o amplitude amplitude.f90
      
!!=======================================================================
!      program amplitude
!
!        parameter ( nt = 96 ,ncon0 =100, nboot = 100,rmax=5)
!        parameter ( nq=1,nin=2,nout=2,nmax=10,nom=3)
!
!
!	integer in,out,dp,flag,i,ns,indall
!	real*8 p(4,2),m1,m2,pq(4)
!        complex*16 ressort(nmax,(2*nin+1)*(2*nq+2)*(2*nout+1)*nom*nom*nom)
!        integer resind(9,(2*nin+1)*(2*nq+2)*(2*nout+1)*nom*nom*nom)
!
!
!!-------in is angular momentum of initial state
!	in=2
!!-------out is angular momentum of final state
!	out=1
!!-------dp is product of P patity of initial and final states.
!	dp=-1
!!-------if flag=1,initial state is still. if flag=2, final one is still.
!	flag=1
!!-------ns is the space lat number
!	ns=8
!
!!-------mass of initial state
!	m1=0.79
!!-------mass of final state
!        m2=0.696
!
!	print*, 'please input J of initial, J of final, P patity product,&
!		 mom flag, and ns'
!	read(*,*) in,out,dp,flag,ns
!	print*, 'plese input mass of initial state and final state'
!	read(*,*) m1,m2
!
!	call FacCG(in,out,dp,m1,m2,ns,flag,indall,ressort,resind) 
!
!        open(12,file='res_sortform')
!
!        write(12,'(220A)') title
!
!        write(12,"(A,f9.5)") '   Q2=',        &
!                -(m1-m2)**2
!
!        write(*,*) indall
!
!        do ii=1,indall
!
!        if(ii.gt.1) then
!        if(resind(7,ii).gt.resind(7,ii-1)) then
!        if(flag.eq.2) then
!        p(1,1)=0.4d0*dsin(3.14159265358979323d0*(resind(4,ii)-1)/ns)
!        p(2,1)=0.4d0*dsin(3.14159265358979323d0*(resind(5,ii)-1)/ns)
!        p(3,1)=0.4d0*dsin(3.14159265358979323d0*(resind(6,ii)-1)/ns)
!        end if
!        p(4,1)=sqrt(m1**2+p(1,1)**2+p(2,1)**2+p(3,1)**2)
!        if(flag.eq.1) then
!        p(1,2)=-0.4d0*dsin(3.14159265358979323d0*(resind(4,ii)-1)/ns)
!        p(2,2)=-0.4d0*dsin(3.14159265358979323d0*(resind(5,ii)-1)/ns)
!        p(3,2)=-0.4d0*dsin(3.14159265358979323d0*(resind(6,ii)-1)/ns)
!        end if
!        p(4,2)=sqrt(m2**2+p(1,2)**2+p(2,2)**2+p(3,2)**2)
!        do i=1,4
!        pq(i)=p(i,1)-p(i,2)
!        end do
!        write(12,"(A,f9.5)") '   Q2=',  &
!                -pq(4)**2+pq(1)**2+pq(2)**2+pq(3)**2
!        end if
!        end if
!
!        write(12,"(6I3,3I4,,2f9.5,A,2f9.5,A,2f9.5,A,2f9.5,A,2f9.5,A     &
!                        ,2f9.5,A,2f9.5,A)") &
!                        resind(4,ii)-1,resind(5,ii)-1,resind(6,ii)-1,   &
!                        resind(1,ii),resind(2,ii),resind(3,ii),         &
!                        resind(7,ii),resind(8,ii),resind(9,ii),         &
!                        ressort(1,ii),"  ",    &
!                        ressort(2,ii),"  ",    &
!                        ressort(3,ii),"  ",    &
!                        ressort(4,ii),"  ",    &
!                        ressort(5,ii),"  ",    &
!                        ressort(6,ii),"  ",    &
!                        ressort(7,ii),"  "
!
!        if(ii.lt.indall) then
!        if(resind(8,ii).ne.resind(8,ii+1)) then
!        write(12,"(20A)") '----------------------'
!        end if
!        end if
!
!        end do
!
!        close(12)
!
!	end 
!
!=======================================================================
!!-------in is angular momentum of initial state
!!-------out is angular momentum of final state
!!-------dp is product of P patity of initial and final states.
!!-------if flag=1,initial state is still. if flag=2, final one is still.
!!-------ns is the space lat number
!!-------m1 is mass of initial state
!!-------m2 is mass of final state
!	ressort and resind return sorted coefficent of formfactor and 
!		sorted mode index
!	indall is the overall index of non-zero channel

        subroutine FacCG(in,out,dp,m1,m2,ns,flag,indall,ressort,resind)

        parameter ( nq=1,num=5,nmax=10,emax=4*4*5,nom=3)
!---------------------
!	nq is the angular momentum of jet.
!	num is max number of given kind of form factor, J_i=J_f=2,
!  		num=J_i+J_f+1=5
!	nmax is max number of kinds of form factors of given channel.
!	emax is max number of the dimention of epsl tensor. 
!		use to switch between different angular momentum of paticles
!	nom is max number of momentum of given dimention plus one.
!----------------------

        integer in,out,dp,flag,ns
        integer i,j,k,l,m,i1,j1,k1,ii,icount
	integer mom1,mom2,mom3
	real*8 m1,m2
        real*8 p(4,2),theta(2),pq(4),phi(3)
        real*8 res0(2*in+1,2*nq+1,2*out+1,num,3)
        complex*16 res(2*in+1,2*nq+2,2*out+1,num,3)
	complex*16 resnew(2*in+1,2*nq+2,2*out+1,nom,nom,nom,nmax)

        complex*16 ressort(nmax,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)
        integer resind(9,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)
	integer indall

	complex*16 e1(emax,num),met1(emax,num)
	complex*16 e2(emax,num),met2(emax,num)
	complex*16 e0(4*3)
	complex*16 eq(4,3),metq1(4,4)
	real*8 check(num,3)
	character form(3)
	character title(220)

	form(1)='E'
	form(2)='M'
	form(3)='C'
	do i=1,220
	title(i)=' '
	end do
	title(3)='x'
	title(6)='y'
	title(9)='z'
	title(12)='i'
	title(15)='q'
	title(18)='f'

	resnew=0.0d0

	do mom3=1,nom
	do mom2=1,nom
       	do mom1=1,nom

	p=0.0d0
	if(flag.eq.2) then
        p(1,1)=0.4d0*dsin(3.14159265358979323d0*(mom1-1)/ns)
        p(2,1)=0.4d0*dsin(3.14159265358979323d0*(mom2-1)/ns)
        p(3,1)=0.4d0*dsin(3.14159265358979323d0*(mom3-1)/ns)
	end if
        p(4,1)=sqrt(m1**2+p(1,1)**2+p(2,1)**2+p(3,1)**2)
	if(flag.eq.1) then
        p(1,2)=-0.4d0*dsin(3.14159265358979323d0*(mom1-1)/ns)
        p(2,2)=-0.4d0*dsin(3.14159265358979323d0*(mom2-1)/ns)
        p(3,2)=-0.4d0*dsin(3.14159265358979323d0*(mom3-1)/ns)
	end if
        p(4,2)=sqrt(m2**2+p(1,2)**2+p(2,2)**2+p(3,2)**2)

        do i=1,4
        pq(i)=p(i,1)-p(i,2)
        end do

	call lorentzgen(p(1,1),phi)

	call lorentzboost(pq,phi,-1)

        call mktheta(pq,theta)

        call mkepsl(pq,theta,e0)
	do j=1,3
	do i=1,4
        eq(i,j)=e0(i+j*4-4)
	end do
	call lorentzboost(eq(1,j),phi,1)
	end do
        call mkmetQ(eq,metq1)

	met2(1,1)=1.0d0
        call mkepsl(p(1,2),theta,e0)
        do j=1,3
	do i=1,4
	k=j*4+i-4
        e2(k,2)=e0(k)
        end do
	call lorentzboost(e2(j*4-3,2),phi,1)
	end do
	call mke2(e2(1,2),e2(1,3))
        call mkmet(e2(1,2),e2(1,3),met2(1,2),met2(1,3))

        met1(1,1)=1.0d0
        call mkepsl(p(1,1),theta,e0)
        do j=1,3
        do i=1,4
        k=j*4+i-4
        e1(k,2)=dconjg(e0(k))
        end do
        call lorentzboost(e1(j*4-3,2),phi,1)
        end do
        call mke2(e1(1,2),e1(1,3))
        call mkmet(e1(1,2),e1(1,3),met1(1,2),met1(1,3))

	res0=0.0d0
	call Amp(in,out,dp,res0)

!        do k=1,2*out+1
!        do j=1,2*nq+1
!        do i=1,2*in+1
!        write(*,"(3I4,10f10.5)") i,j,k,res0(i,j,k,2,1),res0(i,j,k,4,1),&
!                        res0(i,j,k,3,2),res0(i,j,k,2,3),res0(i,j,k,4,3)
!        end do
!        end do
!        end do

	res=0.0d0
	check=0.0d0
	icount=0

	do m=1,3
	do l=1,num

	  do k=1,2*out+1
	  do j=1,2*nq+2
	  do i=1,2*in+1

	    do k1=1,2*out+1
	    do j1=1,2*nq+1
	    do i1=1,2*in+1
	    	res(i,j,k,l,m)=res(i,j,k,l,m)+		&
	    		met1(i+(i1-1)*(2*in+1),in+1)*	&
	    		met2(k+(k1-1)*(2*out+1),out+1)*	&
	    		metq1(j,j1)*res0(i1,j1,k1,l,m)
	    end do 
	    end do
	    end do

	    check(l,m)=check(l,m)+cdabs(res(i,j,k,l,m))
	  end do
	  end do
	  end do

	  if(check(l,m).gt.0.0001) then
	    icount = icount+1
            do k=1,2*out+1
            do j=1,2*nq+2
            do i=1,2*in+1	
	      resnew(i,j,k,mom1,mom2,mom3,icount)=res(i,j,k,l,m)
            end do
            end do
            end do

	    if(mom1+mom2+mom3.eq.3) then
	      title(icount*20+14)=form(m)
	      title(icount*20+15)=char(48+l-1)
	    end if

	  end if

	end do
	end do

        end do
        end do
        end do


	open(11,file='res_formfacor')

	write(11,'(220A)') title
        do mom3=1,nom
        do mom2=1,nom
        do mom1=1,nom

	write(11,"(20A)") '----------------------'

        do k=1,2*out+1
        do j=1,2*nq+2
        do i=1,2*in+1
	write(11,"(6I3,2f9.5,A,2f9.5,A,2f9.5,A,2f9.5,A,2f9.5,A)")  &
			mom1-1,mom2-1,mom3-1,i,j,k,		&
			resnew(i,j,k,mom1,mom2,mom3,1),"  ",	&
                        resnew(i,j,k,mom1,mom2,mom3,2),"  ",    &
                        resnew(i,j,k,mom1,mom2,mom3,3),"  ",    &
                        resnew(i,j,k,mom1,mom2,mom3,4),"  ",    &
                        resnew(i,j,k,mom1,mom2,mom3,5),"  "
        end do
        end do
        end do


	end do
	end do
	end do

	close(11)

	call sortFAC(resnew,icount,in,out,ressort,resind,indall)

         open(13, file='sort_index.dat', form='unformatted', &
                  access='direct', recl=9*indall*4)

       write(13,rec=1)                                                 &
               ((resind(i,ii),i=1,9),ii=1,indall)
       close(13)
         open(14, file='sort_form.dat', form='unformatted', &
                  access='direct', recl=nmax*indall*16)

       write(14,rec=1)                                                 &
               ((ressort(i,ii),i=1,nmax),ii=1,indall)
       close(14)

        end subroutine FacCG
!	end  
!=======================================================================
!       generate lorentz boost parameter make initial state rest,
!       arctanh(v/c), theta(1),theta(2)
!       

        subroutine lorentzgen(p,phi)

        real*8 p(4),pa,phi(3)

        pa=dsqrt(p(1)**2+p(2)**2+p(3)**2)

	phi(1)=0.5d0*dlog((1+pa/p(4))/(1-pa/p(4)))

        if(pa.lt.0.000001) then
                phi(2)=0
                phi(3)=0
        else
                phi(2)=dacos(p(3)/pa)
                if(abs(pa-dabs(p(3))).lt.0.00001) then
                phi(3)=0
                else
                phi(3)=dacos(p(1)/sqrt(p(1)**2+p(2)**2))
                if(p(2).lt.0) phi(3)=2*dacos(-1.0d0)-phi(3)
                end if
        end if

!	phi(2)=dacos(-1.0d0)-phi(2)
!	phi(3)=dacos(-1.0d0)+phi(3)

        print*,phi(1),phi(2),phi(3)


        end subroutine lorentzgen

!=======================================================================
!       calc polarization vector of paticle

        subroutine lorentzboost(p,phi,flag)

        real*8 p(4),phi(3)
	real*8 pt(4),phit(3),lb(4,4)
        integer i,j,flag

	do i=1,3
	phit(i)=phi(i)
	end do
	if(flag.gt.0) then
	phit(1)=-phi(1)
	end if

	lb(4,4)=dcosh(phi(1))
	lb(1,4)=dsin(phit(2))*dcos(phit(3))
        lb(2,4)=dsin(phit(2))*dsin(phit(3))
        lb(3,4)=dcos(phit(2))
        lb(4,1)=lb(1,4)
        lb(4,2)=lb(2,4)
        lb(4,3)=lb(3,4)
	do j=1,3
	do i=1,3
	lb(i,j)=lb(i,4)*lb(4,j)*(dcosh(phit(1))-1)
	end do
	end do
	do j=1,3
	lb(j,4)=-lb(j,4)*dsinh(phit(1))
        lb(4,j)=-lb(4,j)*dsinh(phit(1))
	lb(j,j)=lb(j,j)+1
	end do

!        write(*,'(A,4f12.8)'),'lb0,',lb(4,4),lb(1,4),lb(2,4),lb(3,4)
!        write(*,'(A,4f12.8)'),'lb1,',lb(4,1),lb(1,1),lb(2,1),lb(3,1)
!        write(*,'(A,4f12.8)'),'lb2,',lb(4,2),lb(1,2),lb(2,2),lb(3,2)
!        write(*,'(A,4f12.8)'),'lb3,',lb(4,3),lb(1,3),lb(2,3),lb(3,3)

	do i=1,4
	pt(i)=0.0d0
	do j=1,4
	pt(i)=pt(i)+lb(j,i)*p(j)
	end do
	end do
	do i=1,4
	p(i)=pt(i)
	end do

        end subroutine lorentzboost

!=======================================================================

!=======================================================================
!	function to sort coefficent of channel by monemtum mode and 
!		independdence channels

        subroutine sortFAC(resnew,ncount,in,out,ressort,resind,indall)

        parameter ( nq=1,num=5,nmax=10,nom=3)

	integer in,out,icout
	integer i,j,k,ii,mom1,mom2,mom3,index
        complex*16 resnew(2*in+1,2*nq+2,2*out+1,nom,nom,nom,nmax)

	integer indall,weigth,start,ind1,ind2,flag
	real*8 tmp,tmpold
	complex*16 tmp2,tmp2old
	real*8 sort((2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)
	integer ind((2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)

	complex*16 sumtmp
        complex*16 ressort(nmax,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)
	integer resind(9,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)
	integer temp(6)

        do mom3=1,nom
        do mom2=1,nom
	do mom1=1,nom

        do k=1,2*out+1
        do j=1,2*nq+2
        do i=1,2*in+1
		index=  i+(j-1)*(2*in+1)+(k-1)*(2*in+1)*(2*nq+2)	&
			+(mom1-1)*(2*in+1)*(2*nq+2)*(2*out+1)		&
			+(mom2-1)*(2*in+1)*(2*nq+2)*(2*out+1)*nom	&
			+(mom3-1)*(2*in+1)*(2*nq+2)*(2*out+1)*nom**2
		ind(index)=index
		sumtmp=0.0d0
		do ii=1,ncount
		sumtmp=sumtmp+resnew(i,j,k,mom1,mom2,mom3,ii)*7**ii
		end do
                sort(index)=cdabs(sumtmp)+			&
			100000.d0*((mom1-1)**2+(mom2-1)**2+(mom3-1)**2)
        end do
        end do
        end do

        end do
        end do
        end do

	call quicksort(sort,ind,1,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom, &
			(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom)

!	write(*,"(I5,f30.5)") (ind(ii),sort(ii),ii=1,(2*in+1)*(2*nq+2)*(2*out+1)*nom**3)

	indall=0
	ind1=1
	ind2=1
	weight=0
	flag=1
	tmpold=0.0d0
	tmp2old=0.0d0
	start=1
        do ii=1,(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom
                index=ind(ii)-1
                temp(1)=mod(index,2*in+1)+1
                index=index/(2*in+1)
                temp(2)=mod(index,2*nq+2)+1
                index=index/(2*nq+2)
                temp(3)=mod(index,2*out+1)+1
                index=index/(2*out+1)
                temp(4)=mod(index,nom)+1
                index=index/nom
                temp(5)=mod(index,nom)+1
                index=index/nom
                temp(6)=mod(index,nom)+1

           itmp=sort(ii)/100000
           tmp=sort(ii)-100000.d0*itmp
           if(tmp.gt.0.00001) then
                indall=indall+1
                do i=1,6
                   resind(i,indall)=temp(i)
                end do
		tmp2=0.0d0
                do i=1,ncount
                   ressort(i,indall)=resnew(temp(1),temp(2),temp(3),temp(4),&
                        temp(5),temp(6),i)
		   tmp2=tmp2+ressort(i,indall)
                end do
                if(dabs(tmp-tmpold).gt.0.00001.and.tmpold.gt.0.00001) then
                   do i=start,indall-1
                       resind(9,i)=weight*resind(9,i)
		   end do
		   ind2=ind2+1
		   start =indall
	 	   weight=0
		   flag=1
		end if
		if(cdabs(tmp2old+tmp2).lt.0.00001) flag=-flag
		tmp2old=tmp2
		tmpold=tmp
		resind(7,indall)=ind1
		resind(8,indall)=ind2
		resind(9,indall)=flag
                weight=weight+1
		if(ii.lt.(2*in+1)*(2*nq+2)*(2*out+1)*nom*nom*nom) then
		if(itmp+1-sort(ii+1)/100000.lt.0.00001) then
		   ind1=ind1+1
		   ind2=0
		end if
		end if
	   end if

!	write(*,'(3I3,f20.5)') temp(4),temp(5),temp(6),sort(ii)

        end do

        do i=start,indall
           resind(9,i)=weight*resind(9,i)
        end do


	end subroutine sortFAC

!=======================================================================
!	standard quick sort function.
!	sort original channel index 'ind' to the oder corresponding 
!	the order of index for sort named as 'sort'

        recursive subroutine quicksort(sort,ind,start,end,nmax)

	integer ind(nmax),start,end,nmax
	real*8 sort(nmax)
	integer i,j,itmp
	real*8 tmp

	tmp=sort(start)
	itmp=ind(start)
	i=start
	j=end
	do
	  if(i.ge.j) goto 300
	  do 
	    if(sort(j).lt.tmp.or.i.ge.j) goto 100
	    j=j-1
	  end do
100       sort(i)=sort(j)
	  ind(i)=ind(j)
          do

            if(sort(i).gt.tmp.or.i.ge.j) goto 200
            i=i+1
          end do
200       sort(j)=sort(i)
          ind(j)=ind(i)
	end do

300	sort(i)=tmp
	ind(i)=itmp
	if(i-1.gt.start) call quicksort(sort,ind,start,i-1,nmax)
	if(end.gt.i+1) call quicksort(sort,ind,i+1,end,nmax)

        end subroutine quicksort

!=======================================================================
!	calc tranform matrix from CG reduce to amplitude for jet

        subroutine mkmetQ(eq,met1)
        complex*16 eq(4,3),met1(4,4)
        real*8 ex1(4,4),matric(4,4)
        integer i,j,i1,j1

        ex1=0.0d0
        matric=0.0d0
        matric(4,4)=1.0d0
        do i=1,4
        ex1(i,i)=1
        end do
        do i=1,3
        matric(i,i)=-1.0d0
        end do

        met1=0.0d0
        met2=0.0d0

        do j=1,3
        do i=1,4
        do j1=1,4
        met1(i,j)=met1(i,j)+eq(j1,4-j)*matric(j1,j1)*ex1(j1,i)*(-1)**(j-1)
        end do
        end do
        end do

        end subroutine mkmetQ

!=======================================================================
!       calc tranform matrix from CG reduce to amplitude for in/out
        subroutine mkmet(eq,eq2,met1,met2)
	complex*16 eq(4,3),eq2(4,4,5),met1(3,3),met2(5,5)
        real*8 ex1(4,4),ex2(4,4,5),matric(4,4)
	integer i,j,i1,j1

        ex1=0.0d0
        ex2=0.0d0
	matric=0.0d0
	matric(4,4)=1.0d0
        do i=1,4
        ex1(i,i)=1
        end do
        do i=1,3
        j=mod(i,3)+1
        ex2(j,6-i-j,i)=sqrt(0.5d0)
        ex2(6-i-j,j,i)=sqrt(0.5d0)
	matric(i,i)=-1.0d0
        end do
        ex2(1,1,4)=sqrt(0.5d0)
        ex2(2,2,4)=-sqrt(0.5d0)
        ex2(1,1,5)=-sqrt(1.0d0/6)
        ex2(2,2,5)=-sqrt(1.0d0/6)
        ex2(3,3,5)=sqrt(2.0d0/3)

	met1=0.0d0
	met2=0.0d0

	do j=1,3
	do i=1,3
	do j1=1,4
	met1(i,j)=met1(i,j)+eq(j1,4-j)*matric(j1,j1)*ex1(j1,i)*(-1)**(j-1)
	end do
	end do
	end do

        do j=1,5
        do i=1,5
	do i1=1,4
        do j1=1,4
        met2(i,j)= met2(i,j)+			&
		eq2(i1,j1,6-j)*matric(j1,j1)*ex2(j1,i1,i)*matric(i1,i1)	&
		*(-1)**(j-1)
        end do
        end do
	end do
        end do

!        do j=1,5 
!	write(*,"(i3,10f10.5)") j,met2(1,j),met2(2,j),met2(3,j),met2(4,j),met2(5,j)
!        end do

	end subroutine mkmet

!=======================================================================
!	calc polarization tensor of paticle	

        subroutine mke2(e1,e2)
        complex*16 e1(4,3),e2(4,4,5)
	integer i,j

	do j=1,4
	do i=1,4
	e2(i,j,1)=e1(i,1)*e1(j,1)
	e2(i,j,2)=(e1(i,1)*e1(j,2)+e1(i,2)*e1(j,1))*sqrt(0.5d0)
	e2(i,j,3)=sqrt(2.0d0/3)*(e1(i,1)*e1(j,3)*0.5d0+e1(i,3)*e1(j,1)*0.5d0	&
			+e1(i,2)*e1(j,2))
        e2(i,j,4)=(e1(i,3)*e1(j,2)+e1(i,2)*e1(j,3))*sqrt(0.5d0)
        e2(i,j,5)=e1(i,3)*e1(j,3)
	end do
	end do

!	do j=1,5
!	do i=1,4
!	write(*,"(2i3,8f10.5)") i,j,e2(1,i,j),e2(2,i,j),e2(3,i,j),e2(4,i,j)
!	end do
!	end do

        end subroutine mke2

!=======================================================================
!	calc polarization angle of paticle

        subroutine mktheta(p,theta)
        real*8 p(4),pa,theta(2)

        pa=sqrt(p(1)**2+p(2)**2+p(3)**2)

        if(pa.lt.0.000001) then
                theta(1)=0
                theta(2)=0
        else
                theta(1)=dacos(p(3)/pa)
		if(abs(pa-dabs(p(3))).lt.0.00001) then
		theta(2)=0
		else
		theta(2)=dacos(p(1)/sqrt(p(1)**2+p(2)**2))
		if(p(2).lt.0) theta(2)=2*dacos(-1.0d0)-theta(2)
		end if
	end if

	print*,theta(1),theta(2)

	end subroutine mktheta

!=======================================================================
!	calc polarization vector of paticle 

        subroutine mkepsl(p,theta,epsl)
	
	real*8 p(4),pa,mass,theta(2),sc(4)
	complex*16 epsl(4,3)
	integer i,j

	pa=sqrt(p(1)**2+p(2)**2+p(3)**2)
	if((p(1)+p(2)+p(3)).lt.0.0d0) pa=-pa
	if(abs(p(4)).lt.abs(pa)) then
	mass=sqrt(-p(4)**2+pa**2)
	else
	mass=sqrt(p(4)**2-pa**2)
	end if

	sc(1)=dcos(theta(1))
	sc(2)=dsin(theta(1))
        sc(3)=dcos(theta(2))
        sc(4)=dsin(theta(2))
        do i=1,4
        if(abs(sc(i)).lt.0.000001) sc(i)=0.0d0
        end do

	epsl(1,2)=-sc(2)*sc(3)*p(4)/mass
	epsl(2,2)=-sc(2)*sc(4)*p(4)/mass
	epsl(3,2)=-sc(1)*p(4)/mass
	epsl(4,2)=-pa/mass
	epsl(1,3)=dcmplx(sc(1)*sc(3),-sc(4))*sqrt(0.5d0)
        epsl(1,1)=dcmplx(-sc(1)*sc(3),-sc(4))*sqrt(0.5d0)
	epsl(2,3)=dcmplx(sc(1)*sc(4),sc(3))*sqrt(0.5d0)
        epsl(2,1)=dcmplx(-sc(1)*sc(4),sc(3))*sqrt(0.5d0)
	epsl(3,3)=-sc(2)*sqrt(0.5d0)
	epsl(3,1)=sc(2)*sqrt(0.5d0)
	epsl(4,3)=0.0d0
	epsl(4,1)=0.0d0

!        do j=1,3
!        write(*,"(i3,f10.5,8e10.3)") j,mass,epsl(1,j),epsl(2,j),epsl(3,j),epsl(4,j)
!        end do

	end subroutine mkepsl

!=======================================================================
!	calc reditive decay amplitude generate by CG coefficent

        subroutine Amp(in,out,dp,res)

        parameter ( nq=1,num=5)
	integer in,out,dp
	integer i,j,k,l,i0,j0,k0,l0
	real*8 res(2*in+1,2*nq+1,2*out+1,num,3),result1,result2

	res=0.0d0

	do l=1,num
	do k=1,2*out+1
	do i=1,2*in+1
	i0=i-in-1
	k0=k-out-1
	l0=l-1

	if(l0.ge.1.and.l0.le.(in+out)) then
        call CG(l0,-1,out,k0,in,i0,result1)
        call CG(l0,1,out,k0,in,i0,result2)	
!-----------E
	res(i,1,k,l,1)=			&
		sqrt((2*l0+1)/dble(2*in+1))*(1+dp*(-1)**l0)*0.5d0*result1
	res(i,3,k,l,1)=			&
                sqrt((2*l0+1)/dble(2*in+1))*(1+dp*(-1)**l0)*0.5d0*result2
!-----------M
        res(i,1,k,l,2)= -               &
                sqrt((2*l0+1)/dble(2*in+1))*(1-dp*(-1)**l0)*0.5d0*result1
        res(i,3,k,l,2)=                 &
                sqrt((2*l0+1)/dble(2*in+1))*(1-dp*(-1)**l0)*0.5d0*result2
	end if

	if(l0.ge.abs(in-out).and.l0.le.(in+out)) then
	call CG(l0,0,out,k0,in,i0,result1)
!-----------C
        res(i,2,k,l,3)=			&
		sqrt((2*l0+1)/dble(2*in+1))*(1+dp*(-1)**l0)*0.5d0*result1
	end if

	end do
	end do
	end do
	
        end subroutine Amp

!=======================================================================
!	calc CG coefficent

        subroutine CG(j1,m1,j2,m2,j,m,res)

        integer, external :: fac
	real*8 , external :: facm1

        integer j1,m1,j2,m2,j,m,v
	real*8 res

	res=0.0d0
	if(m.eq.(m1+m2)) then
	do v=-2*j1-2*j2,2*j1+2*j2
	res=res+((-1)**v)*		&
		facm1(j1+j2-j-v)*facm1(j1-m1-v)*facm1(j2+m2-v)*	&
		facm1(j-j2+m1+v)*facm1(j-j1-m2+v)*facm1(v)
	end do
	res=res*dsqrt((2*j+1)*fac(j1+j2-j)*fac(j+j1-j2)*fac(j2+j-j1)	&
		*facm1(j1+j2+j+1)		&
		*fac(j+m)*fac(j-m)		&
		*fac(j1+m1)*fac(j1-m1)		&
		*fac(j2+m2)*fac(j2-m2))
	end if

        end subroutine CG

      
!=======================================================================
!	inverse of factorial of n

        function facm1(n)

	integer, external :: fac

        integer n
	real*8 facm1

        if(n.lt.0) then
	facm1=0
	else
	facm1=1.0d0/fac(n)
	end if

        return
        end

!=======================================================================
!	factorial of n

        function fac(n)

	integer i,n,fac

	fac=1 
	if(n.lt.0) fac=10**10 
	if(n.gt.0) then 
	do i=1,n 
	fac=fac*i 
	end do
	end if	

	return
        end 


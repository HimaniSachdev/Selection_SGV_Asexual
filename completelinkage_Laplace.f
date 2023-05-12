      program introgression
      integer(KIND=8) tsteps,histories,N,i,j,c1,c2,k,h1,h2,tinit
      integer(KIND=8) j6,j5,parent,rec_flag,j3,j4,j8,j9,L
      integer(KIND=8) c3,c4,num_mark,rec_points,j8a,l1,l2,j8b,j8c
      integer(KIND=8) lplus,lminus,k3
      real(kind=8) y0,theta,temp1,temp2,y3,y4,x1,x2,x3,x4,V0,sd
      real(kind=8) temp,test1,test2,test3,tempg,error,mut,initf
      integer(KIND=8),allocatable, dimension(:)::nseq,nseq2
      real(KIND=8),allocatable, dimension(:)::fitness,relfit,varpop
      real(KIND=8),allocatable, dimension(:,:)::zfinal
      real(KIND=8),allocatable,dimension(:,:)::mom
      REAL(KIND=8), PARAMETER :: Pi = 3.1415926535897932





      N=25600       ! size of adapting population
      sd=9.5d+0/N
      V0=sd*sd       ! initial fitness variance

! v1:10^5 histories, v2: 10^6 histories      

      tsteps=1000000
      histories=1000 ! number of replicate populations
      theta=1.0D+0/histories


      allocate(fitness(N),nseq(N),relfit(N),mom(tsteps,7))
      allocate(nseq2(N),varpop(tsteps),zfinal(histories,3))


      do i=1,tsteps
      do j=1,7
      mom(i,j)=0.0d+0
      enddo
      varpop(i)=0.0d+0
      enddo


      call init_random_seed
 



      do k3=1,histories

       
       tempg=0.0d+0
       temp=sd/sqrt(2.0d+0)
       do i=1,N
       call random_number(x1)       

       if (x1.le. 0.5d+0)then
       fitness(i)=temp*log(2.0d+0*x1)
       else 
       fitness(i)=-temp*log(2.0d+0*(1.0d+0-x1))
       endif

       if (fitness(i)>tempg)then
       tempg=fitness(i)
       endif  
       enddo   
       tempg=tempg/sd

       do i=1,N
       nseq(i)=i
       enddo   


! start time evolution
       do i=1,tsteps

! calculate relative fitness of each individual
       temp1=tiny(double)

       do j=1,N

       relfit(j)=exp(fitness(nseq(j))) 
       if (relfit(j)>temp1)then
       temp1=relfit(j)
       endif
       enddo

       do j=1,N
       relfit(j)=relfit(j)/temp1
       enddo
        
 
! create offspring by randomly sampling a parent in proportion to rel fitness

        do j=1,N                     

 52     call random_number(x1)
        call random_number(x2)
        j6=aint(1+N*x2)
        if (x1>relfit(j6))then
         goto 52
        endif       
        nseq2(j)=nseq(j6)
        enddo

        do j=1,N                     
        nseq(j)=nseq2(j)
        enddo


! measure moments of trait value distribution
        if (mod(i,100)==0)then
        x1=0.0d+0
        x2=0.0d+0
        x3=0.0d+0
        x4=0.0d+0
        do j=1,N
        x1=x1+fitness(nseq(j))
        x2=x2+(fitness(nseq(j))**2.0d+0)
        x3=x3+(fitness(nseq(j))**3.0d+0)
        x4=x4+(fitness(nseq(j))**4.0d+0)
        enddo
        x1=x1/(N*sd)
        x2=x2/(N*(sd**2.0d+0))
        x3=x3/(N*(sd**3.0d+0))
        x4=x4/(N*(sd**4.0d+0))
        y2=x2-x1*x1
        y3=x3-3*x1*x2+2*x1*x1*x1
        y4=x4-4*x3*x1+6*x2*x1*x1-3*x1*x1*x1*x1
        mom(i/100,1)=mom(i/100,1)+theta*x1
        mom(i/100,2)=mom(i/100,2)+theta*y2
        mom(i/100,3)=mom(i/100,3)+theta*y3
        mom(i/100,4)=mom(i/100,4)+theta*y4
        mom(i/100,7)=mom(i/100,7)+theta*x1*x1

        if (y2>0.0000000002d+0)then
        varpop(i)=varpop(i)+theta
        y3=y3/(y2*sqrt(y2))
        y4=y4/(y2*y2)
        mom(i,5)=mom(i,5)+theta*y3
        mom(i,6)=mom(i,6)+theta*y4
        else
        do j=i+100,tsteps,100
        mom(j/100,1)=mom(j/100,1)+theta*x1
        mom(j/100,2)=mom(j/100,2)+theta*y2
        mom(j/100,3)=mom(j/100,3)+theta*y3
        mom(j/100,4)=mom(j/100,4)+theta*y4
        mom(j/100,7)=mom(j/100,7)+theta*x1*x1
        enddo
        go to 22
        endif
        endif

        enddo  

 22     zfinal(k3,1)=x1
        zfinal(k3,2)=tempg
        print*,k3,i
        enddo   ! all replicates completed.



       open(4,file='N025600Nsd9.5B.txt')
       do i=1,tsteps
       write(4,94) i, (mom(i,j),j=1,4),mom(i,7),(mom(i,j)/varpop(i)
     c ,j=5,6),varpop(i)
 94    Format(I8,1X,F11.8,1X,F11.8,1X,F11.8,1X,F11.8,1X,F14.8,1X,F11.8,
     c 1X,F11.8,1X,F11.8)
       enddo
       close(4)

       open(4,file='N025600Nsd9.5zfinalB.txt')
       do i=1,histories
       write(4,95) i,zfinal(i,1),zfinal(i,2)
 95    Format(I8,1X,F9.6,1X,F9.6)
       enddo
       close(4)

       open(4,file='temp.txt')
       do i=1,N
       write(4,96) i,nseq(i)
 96    Format(I8,1X,I8)
       enddo
       close(4)


         end 


      SUBROUTINE init_random_seed()
      INTEGER :: z5, m, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed, seedold
      
      CALL RANDOM_SEED(size = m)
      ALLOCATE(seed(m), seedold(m))
      
      CALL SYSTEM_CLOCK(COUNT=clock)
      
      seed = clock + 37 * (/ (z5 - 1, z5 = 1, m) /)
      CALL RANDOM_SEED(PUT = seed)
      CALL RANDOM_SEED(GET = seedold)
      
      DEALLOCATE(seed)
      END SUBROUTINE init_random_seed



c
c
c     ################################
c     ##  Added by Steve Rick 2015  ##
c     ################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine newcrg  --  get charge after charge transfer   ##
c     ##                                                            ##
c     ################################################################
c
c
c    DCT subroutine newcrg
c       set up for water-water interations only
c
      subroutine newcrg
      use sizes
      use limits
      use atoms
      use boxes
      use mpole
      use kmulti
      use neigh
      use potent
      implicit none
      integer i,j,kk,ii,l,m
      real*8 xr,yr,zr,rr,pi,sum1,sum2
      real*8, allocatable :: nacti(:)
      real*8, allocatable :: ndcti(:)

      allocate (nacti(n))
      allocate (ndcti(n))

c     write(*,*) "newcrg subroutine"
c MES : this subroutine is called

c  charge transfer matrix
c  zdqt(i,j) charge transfered to atom i due to hydrogen
c  bond of type j : j=1, H1 is donor, j=2, H2 is donor, j=3, O accepts
c  atom i :  1 = O, atom 2 = H1, atom 3 = H2

      pi = 4.d0*datan(1.d0)

      do i=1,n
         nacti(i) = 0.d0
         ndcti(i) = 0.d0
c        dedci(i) = 0.d0
      enddo

c get hydrogen bonds: nacti = number accepted
c                     ndcti = number donated
c MES : This is correct

c do this for oxygen atoms only
      do i=1,n-1
         if(type(i).eq.1) then
            do kk = 1, nelst(i)
               j=elst(kk,i)
               if(type(j).eq.1) then
                  xr = x(j) - x(i)
                  yr = y(j) - y(i)
                  zr = z(j) - z(i)
                  call image (xr,yr,zr)
                  rr = xr*xr + yr*yr + zr*zr

                  if(rr.lt.(rct2+1.200)**2) then
c loop over hydrogen atoms
                     do l=1,2
                        xr = x(j+l) - x(i)
                        yr = y(j+l) - y(i)
                        zr = z(j+l) - z(i)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)

                        if(rr.lt.rct1) then
                           nacti(i) = nacti(i) + 1.d0
                           ndcti(j+l) = ndcti(j+l) + 1.d0
                        endif
                        if(rr.ge.rct1.and.rr.le.rct2) then
                           nacti(i) = nacti(i)
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                           ndcti(j+l) = ndcti(j+l)
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                        endif

                        xr = x(j) - x(i+l)
                        yr = y(j) - y(i+l)
                        zr = z(j) - z(i+l)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)

                        if(rr.lt.rct1) then
                           nacti(j) = nacti(j) + 1.d0
                           ndcti(i+l) = ndcti(i+l) + 1.d0
                        endif
                        if(rr.ge.rct1.and.rr.le.rct2) then
                           nacti(j) = nacti(j)
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                           ndcti(i+l) = ndcti(i+l)
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                        endif

                     enddo
c                    end loop over H
                  endif
c                 end if rr within cutoff
               endif
c              end if neighbor is oxygen
            enddo
c           end loop over nelst
         endif
c        end if O_w 
      enddo
c     end loop over O_w


c transfer charge
c MES : this is correct
      do i=1,n
         if(type(i).eq.1) then
c oxygen charge
            rpole(1,i) = rpole0(i)
     & + zdqt(1,1)*ndcti(i+1) + zdqt(1,2)*ndcti(i+2)
     & + zdqt(1,3)*nacti(i)
c hydrogen hydrogen charges
            rpole(1,i+1) = rpole0(i+1)
     & + zdqt(2,1)*ndcti(i+1) + zdqt(2,2)*ndcti(i+2)
     & + zdqt(2,3)*nacti(i)
            rpole(1,i+2) = rpole0(i+2)
     & + zdqt(3,1)*ndcti(i+1) + zdqt(3,2)*ndcti(i+2)
     & + zdqt(3,3)*nacti(i)
         endif
      enddo


      do i=1,n
       pole(1,i) = rpole(1,i)
      enddo

      sum1 = 0.d0
      sum2 = 0.d0

c total charge - should equal charge of system (usually zero)
      do i=1,n
       sum1 = sum1 + rpole0(i)
       sum2 = sum2 + rpole(1,i)
      enddo

c     write(*,*) "Total system charge : "
c     write(*,*) sum1, sum2

      deallocate (nacti)
      deallocate (ndcti)

 993  format(3f14.6)
      return
      end
c end of newcrg subroutine

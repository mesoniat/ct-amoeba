c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient (energy,derivs)
      use sizes
      use atoms
      use bound
      use couple
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
c     use mpole
      implicit none
      integer i,j,kk
      real*8 energy,cutoff
      real*8 derivs(3,*)
      real*8 dr,delta,q0
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      eat = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(desum)) then
         if (size(desum) .lt. 3*n) then
            deallocate (desum)
            deallocate (deb)
            deallocate (dea)
            deallocate (deba)
            deallocate (deub)
            deallocate (deaa)
            deallocate (deopb)
            deallocate (deopd)
            deallocate (deid)
            deallocate (deit)
            deallocate (det)
            deallocate (dept)
            deallocate (debt)
            deallocate (deat)
            deallocate (dett)
            deallocate (dev)
            deallocate (dec)
            deallocate (decd)
            deallocate (ded)
            deallocate (dem)
            deallocate (dep)
            deallocate (der)
            deallocate (des)
            deallocate (delf)
            deallocate (deg)
            deallocate (dex)
         end if
      end if
      if (.not. allocated(desum)) then
         allocate (desum(3,n))
         allocate (deb(3,n))
         allocate (dea(3,n))
         allocate (deba(3,n))
         allocate (deub(3,n))
         allocate (deaa(3,n))
         allocate (deopb(3,n))
         allocate (deopd(3,n))
         allocate (deid(3,n))
         allocate (deit(3,n))
         allocate (det(3,n))
         allocate (dept(3,n))
         allocate (debt(3,n))
         allocate (deat(3,n))
         allocate (dett(3,n))
         allocate (dev(3,n))
         allocate (dec(3,n))
         allocate (decd(3,n))
         allocate (ded(3,n))
         allocate (dem(3,n))
         allocate (dep(3,n))
         allocate (der(3,n))
         allocate (des(3,n))
         allocate (delf(3,n))
         allocate (deg(3,n))
         allocate (dex(3,n))
      end if
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            deat(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
         end do
      end do
c
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 0.0d0
         end do
      end do
      einter = 0.0d0
c
c     DCT variables zeroed out elsewhere
c     dedci set to zero in newcrg subroutine
c     ect set to zero in ect1 subroutine
c     This is ok because dedci and ect are only used with CT routines
c         Their values are added to dem and em for later use
c MES debug
c     dr = 0.000d0
c     x(4) = x(4) + dr
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c MES : This is being called with mpole-list key
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_angtor)  call eangtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy and gradient routines
c
c DCT gets new charges from current geometry
      if (use_crgtr) call newcrg

c FD test
c     q0 = rpole(1,4)
c     delta = 0.0001
c     do kk = -1, 1
c       rpole(1,4) = q0 + real(kk)*delta
c       write(*,*) "Loop ",kk," q(4) = ",rpole(1,4)

      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar) then
         call empole1
      endif
      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1

c DCT gets charge transfer contribution to forces
c        and adds charge transfer energy to em
      if (use_crgtr) call ect1
c
c     sum up to get the total energy and first derivatives
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + eat + ett + ev + ec + ecd + ed
     &          + em + ep + er + es + elf + eg + ex
      energy = esum
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
     &                      + det(j,i) + dept(j,i) + debt(j,i)
     &                      + deat(j,i) + dett(j,i) + dev(j,i)
     &                      + dec(j,i) + decd(j,i) + ded(j,i)
     &                      + dem(j,i) 
     &                      + dep(j,i) + der(j,i)
     &                      + des(j,i) + delf(j,i) + deg(j,i)
     &                      + dex(j,i)
            derivs(j,i) = desum(j,i)
         end do
      end do

c     end do
c     end of FD test

c
c     check for an illegal value for the total energy
c
c     if (isnan(esum)) then
      if (esum .ne. esum) then
         write (iout,10)
   10    format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if

c     call gtest(derivs)
c     call qtest(derivs)
c     call ptest

      return
      end
c end of gradient subroutine

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
c    Calculates number of hydrogen bonds accepted and donated
c    by each water. Based on number and directions of HB, the new
c    charges are calculated.
c
c
      subroutine newcrg
      use sizes
      use limits
      use atoms
      use boxes
      use deriv
      use energi
      use mpole
      use potent
      use neigh
      use kmulti
      implicit none
      integer i,j,kk,ii,l,m
      real*8 xr,yr,zr,rr,pi,sum1,sum2,dq
      real*8, allocatable :: nacti(:)
      real*8, allocatable :: ndcti(:)

      allocate (nacti(n))
      allocate (ndcti(n))

c     write(*,*) "newcrg subroutine"

c  charge transfer matrix zdqt
c  zdqt(i,j) charge transfered to atom i due to hydrogen
c  bond of type j : j=1, H1 is donor, j=2, H2 is donor, j=3, O accepts
c  atom i :  1 = O, atom 2 = H1, atom 3 = H2
c     write(*,*) "CT matrix"
c     do i=1,3
c       do j=1,3
c         write(*,*) i,j,zdqt(i,j)
c       end do
c     end do

      pi = 4.d0*datan(1.d0)

      do i=1,n
         nacti(i) = 0.d0
         ndcti(i) = 0.d0
         dedci(i) = 0.0d0
         dedciX(i) = 0.0d0
         depdci(i) = 0.0d0
      enddo

c get hydrogen bonds: nacti = number accepted
c                     ndcti = number donated

c do this for oxygen atoms only
c for each oxygen, get neighborlist
      do i=1,n-1
         if(type(i).eq.1) then
            do kk = 1, nelst(i)
               j=elst(kk,i)
c if neighbor is an oxygen, calculate distance
               if(type(j).eq.1) then
                  xr = x(j) - x(i)
                  yr = y(j) - y(i)
                  zr = z(j) - z(i)
                  call image (xr,yr,zr)
                  rr = xr*xr + yr*yr + zr*zr
c if distance less than CT cut-off
                  if(rr.lt.(rct2+1.200)**2) then
c loop over hydrogen atoms of neighbor
                     do l=1,2
                        xr = x(j+l) - x(i)
                        yr = y(j+l) - y(i)
                        zr = z(j+l) - z(i)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)
        
c if O-H distance less than Rct1
                        if(rr.lt.rct1) then
                           nacti(i) = nacti(i) + 1.d0
                           ndcti(j+l) = ndcti(j+l) + 1.d0
                        endif
c else if O-H distance between Rct1 and Rct 2
                        if(rr.ge.rct1.and.rr.le.rct2) then
                           nacti(i) = nacti(i) 
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                           ndcti(j+l) = ndcti(j+l) 
     & + 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                        endif
    
c repeat loop for H of original oxygen
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
      do i=1,n
         if(type(i).eq.1) then
c oxygen charge
            rpole(1,i) = rpole0(i) 
     & + zdqt(1,1)*ndcti(i+1) + zdqt(1,2)*ndcti(i+2) 
     & + zdqt(1,3)*nacti(i)
c hydrogen charges
            rpole(1,i+1) = rpole0(i+1) 
     & + zdqt(2,1)*ndcti(i+1) + zdqt(2,2)*ndcti(i+2) 
     & + zdqt(2,3)*nacti(i)
            rpole(1,i+2) = rpole0(i+2) 
     & + zdqt(3,1)*ndcti(i+1) + zdqt(3,2)*ndcti(i+2) 
     & + zdqt(3,3)*nacti(i)
         endif
      enddo
 
c for testing dedci
      dq = 0.0000d0
      rpole(1,4) = rpole(1,4) + dq
      write(*,*) "dq = ",dq,"q(4) = ",rpole(1,4)

c copy new charge to other reference frame
      do i=1,n
        pole(1,i) = rpole(1,i)
      enddo

      sum1 = 0.d0
      sum2 = 0.d0
c total charge - should equal charge of system (usually zero)
c     do i=1,n
c       sum1 = sum1 + rpole0(i)
c       sum2 = sum2 + rpole(1,i)
c     enddo
c     write(*,*) "Total system charge : ",sum1,sum2

      deallocate (nacti)
      deallocate (ndcti)

 993  format(3f14.6)
      return
      end
c end of newcrg subroutine


c
c
c     ################################
c     ##  Added by Steve Rick 2015  ##
c     ################################
c
c   ###################################################################
c    ##                                                             ##
c     ##  subroutine ect1  --  gets charge transfer contribition   ## 
c     ##      to force and charge transfter energy                 ## 
c    ##                                                             ##
c   ###################################################################
c
c
c    DCT subroutine ect1
c       set up for water-water interations only
c    calculates Ect 
c    calculates dEct/dr and dq/dr
c    multiplies dedci by dq/dr and updates dem with the additional
c    force due to charge transfer
c
c
      subroutine ect1
      use sizes
      use deriv
      use energi
      use limits
      use atoms
      use boxes
      use mpole
      use potent
      use neigh
      use kmulti
      use virial
      implicit none
      integer i,j,kk,ii,l,m
      real*8 xr,yr,zr,rr,pi,sum1,sum2,dqt,ect
      real*8 xna,dxna
c xna ==> % of Qct or % of HB ; dxna is its derivative wrt rr
      real*8 tmp1,tmp2,tmp3
      real*8 xr0,yr0,zr0

c dedci already calculated, does not change within this routine
c     write(*,*) "ect1 subroutine"

c     write(*,*) "Pre-CT: dem_x, dem_y, dem_z"
c     do i=1,n
c       write(*,*) dedci(i)
c       dedci(i) = 0.0d0
c       write(*,*) dem(1,i),dem(2,i),dem(3,i) 
c     enddo

c     write(*,*) "CT variables :"
c     write(*,*) "xna            dxna            zr            rr "

c check virial
c     write(*,*) "Virial pre-CT"
c     do i=1,3
c       do j=1,3
c         write(*,*) vir(i,j)
c       enddo
c     enddo

c add forces to dem(i,j=1,2,3) (x,y,z)

      pi = 4.d0*datan(1.d0)

c i.e. total CT between waters = dqt
      dqt=zdqt(1,1)+zdqt(2,1)+zdqt(3,1)
      ect = 0.d0

c do this for oxygen atoms only
c for each oxygen, look at neighbor list
      do i=1,n-1
         if(type(i).eq.1) then
            do kk = 1, nelst(i)
               j=elst(kk,i)
               if(type(j).eq.1) then
c if neighbor is also an oxygen, calculate distance
                  xr = x(j) - x(i)
                  yr = y(j) - y(i)
                  zr = z(j) - z(i)
                  call image (xr,yr,zr)
                  rr = xr*xr + yr*yr + zr*zr
c                 write(*,*) i,j,"rr = ",rr

                  if(rr.lt.(rct2+1.200)**2) then
c check if water with oxygen atom j donates hydrogen bond to i
                     do l=1,2
                        xr = x(j+l) - x(i)
                        yr = y(j+l) - y(i)
                        zr = z(j+l) - z(i)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)
c                       write(*,*) i,j+l,"r = ",rr

                        if(rr.ge.rct1.and.rr.le.rct2) then
                           xna = 
     & 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                           dxna = 
     & -0.5d0*dsin(pi*(rr-rct1)/(rct2-rct1))
     & *pi/(rct2-rct1)

c                          write(*,*) xna,dxna,rr

c dE/dq * dq/dr
                           dem(1,i) = dem(1,i)
     & +(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*xr/rr
                           dem(2,i) = dem(2,i)
     & +(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*yr/rr
                           dem(3,i) = dem(3,i)
     & +(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*zr/rr
                           dem(1,j+l) = dem(1,j+l) 
     & -(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*xr/rr
                           dem(2,j+l) = dem(2,j+l) 
     & -(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*yr/rr
                           dem(3,j+l) = dem(3,j+l) 
     & -(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*zr/rr

                           vir(1,1) = vir(1,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*xr/rr
                           vir(1,2) = vir(1,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*xr/rr
                           vir(1,3) = vir(1,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*xr/rr
                           vir(2,1) = vir(2,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*yr/rr
                           vir(2,2) = vir(2,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*yr/rr
                           vir(2,3) = vir(2,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*yr/rr
                           vir(3,1) = vir(3,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*zr/rr
                           vir(3,2) = vir(3,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*zr/rr
                           vir(3,3) = vir(3,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3))*dxna*zr/rr
                        endif
c                       end if between rct1 and rct2

c molecule with oxygen atom i donates to hydrogen bond
                        xr = x(i+l) - x(j)
                        yr = y(i+l) - y(j)
                        zr = z(i+l) - z(j)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)
c                       write(*,*) i+l,j,"r = ",rr

                        if(rr.le.rct1) then
c add to charge transfer energy, no force (energy is constant)
                           ect=ect+muct*dqt+0.5d0*etact*dqt**2
                        endif

                        if(rr.gt.rct1.and.rr.le.rct2) then

                           xna = 
     & 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
                           dxna = 
     & -0.5d0*dsin(pi*(rr-rct1)/(rct2-rct1))
     & *pi/(rct2-rct1)
c dxna = N' in Eqn A9 of Lee 2011 JCP

c                          write(*,*) xna,dxna,rr

c ect here only (not above) to avoid double counting
                           ect=ect+muct*(xna*dqt)
     & +0.5d0*etact*(xna*dqt)**2

c dE/dq * dq/dr + dE_ct/dr
c MES debug -- I think this is the section I changed at one point
c    It is but changing l+1 to l in zdqt doesn't change dqdr
             dem(1,i+l) = dem(1,i+l)
     &+(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*xr/rr
             dem(2,i+l) = dem(2,i+l)
     &+(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*yr/rr
             dem(3,i+l) = dem(3,i+l)
     &+(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*zr/rr
             dem(1,j) = dem(1,j)
     &-(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*xr/rr
             dem(2,j) = dem(2,j)
     &-(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*yr/rr
             dem(3,j) = dem(3,j)
     &-(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     &)*dxna*zr/rr

c                          write(*,*) "dqdr in x analytic = "
c                          write(*,*) zdqt(1,l+1)*dxna*xr/rr,
c    & zdqt(2,l+1)*dxna*xr/rr,
c    & zdqt(3,l+1)*dxna*xr/rr,l+1
c                          write(*,*) "dqdr in y analytic = ",
c    & (zdqt(1,3)+zdqt(2,3)+zdqt(3,3))*dxna*yr/rr
c                          write(*,*) "dqdr in z analytic = ",
c    & (zdqt(1,3)+zdqt(2,3)+zdqt(3,3))*dxna*zr/rr
c                          write(*,*) "dEct/dr in x = ",
c    & ((muct+etact*xna*dqt)*dqt)*dxna*xr/rr
c                          write(*,*) "dEct/dr in y = ",
c    & ((muct+etact*xna*dqt)*dqt)*dxna*yr/rr
c                          write(*,*) "dEct/dr in z = ",
c    & ((muct+etact*xna*dqt)*dqt)*dxna*zr/rr

                           vir(1,1) = vir(1,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*xr/rr
                           vir(1,2) = vir(1,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*xr/rr
                           vir(1,3) = vir(1,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*xr/rr
                           vir(2,1) = vir(2,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*yr/rr
                           vir(2,2) = vir(2,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*yr/rr
                           vir(2,3) = vir(2,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*yr/rr
                           vir(3,1) = vir(3,1)
     & -xr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*zr/rr
                           vir(3,2) = vir(3,2)
     & -yr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*zr/rr
                           vir(3,3) = vir(3,3)
     & -zr*(dedci(i)*zdqt(1,3)+dedci(i+1)*zdqt(2,3)
     & +dedci(i+2)*zdqt(3,3)
     & +(muct+etact*xna*dqt)*dqt )*dxna*zr/rr
                        endif

                     enddo
c                    end loop over H
                  endif
c                 end if oxygens within cutoff
               endif
c              end if also O_w
            enddo
c           end loop over nblist
         endif
c        end if O_w
      enddo
c     end of loop over all O_w
  
c check virial
c     write(*,*) "Virial post-CT"
c     do i=1,3
c       do j=1,3
c         write(*,*) vir(i,j)
c       enddo
c     enddo

      em = em + ect
c     write(*,*) "Ect = ",ect,"    total em = ",em

c     write(*,*) "Post-CT: dem_x, dem_y, dem_z"
c     do i=1,n
c       write(*,*) dem(1,i),dem(2,i),dem(3,i)
c     enddo

      return
      end
c end of ect1 subroutine


c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gtest  --  find forces via finite-difference  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "gtest" compares dE/dr and dq/dr from analytical and 
c     finite difference methods.
c
c     call from gradient subroutine
c
c     Added by Marielle Soniat 2015 July, Nov
c
c
      subroutine gtest(derivs)
      use sizes
      use atoms
      use bound
      use couple
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use mpole
      use kmulti
      implicit none
      integer i,j,k
      integer atom,dir
      real*8 energy
      real*8 delta
      real*8 fxk,fyk,fzk,e0
      real*8 x0,y0,z0,q0
      real*8 xneg,yneg,zneg
      real*8 xpos,ypos,zpos
      real*8 exneg,eyneg,ezneg
      real*8 expos,eypos,ezpos
      real*8 qneg,qpos,dqdr_fd
      real*8 dedx_fd,dedy_fd,dedz_fd,dedx_an
      real*8 derivs(3,*)

      atom = 4
      dir = 1
      write(*,*) "Entering gtest for atom ",atom," in direction ",dir
      write(*,*) "Testing ect"

c     get initial forces
c     derivs(j,i) = derivatives, j=direction, i=atom
c     dedx_an = derivs(dir,atom)
      dedx_an = dep(dir,atom)

c     save initial position
      x0 = x(atom)
      y0 = y(atom)
      z0 = z(atom)
      q0 = rpole(1,atom)
      write(*,*) "(x,y,z) and q of ",atom," : ",x0,y0,z0,q0

c     change for position
      delta = 0.0001d0

c     for first atom, forces in x, y, and z
      xneg = x0 - delta
      x(atom) = xneg
      call empole1
      exneg = ep
      qneg = rpole(1,atom)

      xpos = x0 + delta
      x(atom) = xpos
      call empole1
      expos = ep
      qpos = rpole(1,atom)

      dedx_fd = ( expos - exneg ) / ( 2.0 * delta )
      write(*,*) "Force via FD   = ",dedx_fd
c     write(*,*) "Force analytic = ",dedx_an

      dqdr_fd = ( qpos - qneg ) / ( 2.0 * delta )
      write(*,*) "dqdr FD = ",dqdr_fd 

c     restore original position after calculation
      x(atom) = x0
      rpole(1,atom) = q0
      call empole1

      write(*,*) "Finished with gtest"

      return
      end
c end of gtest subroutine


c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine qtest  --  find forces via finite-difference  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "qtest" compares dE/dq from analytical and finite difference
c     methods.
c
c     Added by Marielle Soniat 2015 Nov
c
c still not working
c
      subroutine qtest(derivs)
      use sizes
      use atoms
      use bound
      use couple
      use deriv
      use energi
      use inter
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      use virial
      use mpole
      use kmulti
      implicit none
      integer i,j,k,dir
      real*8 energy
      real*8 delta
      real*8 et0,ep0,em0,q0
      real*8 etP,epP,emP
      real*8 etN,epN,emN
      real*8 fep,fem,ftot
      real*8 detdq_fd,detdq_an
      real*8 depdq_fd,depdq_an
      real*8 demdq_fd,demdq_an
      real*8 ldedci,ldepdci
      real*8 derivs(3,*)

      k = 4
      dir = 1
      delta = 0.0001
      write(*,*) "Entering qtest for atom ",k," for dir ",dir

c     get initial values
c     derivs(j,i) = derivatives, dir=direction, k=atom
      fep = dep(dir,k)
      fem = dem(dir,k)
      ftot = derivs(dir,k)
      ep0 = ep
      em0 = em
      et0 = esum
      ldedci = dedci(k)
      ldepdci = depdci(k)
      q0 = rpole(1,k)

      write(*,*) "Positive delta"
      rpole(1,k) = q0 + delta
      call empole1
      epP = ep
      emP = em
      etP = esum
      
      write(*,*) "Negative delta"
      rpole(1,k) = q0 - delta
      call empole1
      epN = ep
      emN = em
      etN = esum

      detdq_fd = ( etP - etN ) / ( 2.0d0 * delta )
      depdq_fd = ( epP - epN ) / ( 2.0d0 * delta )
      demdq_fd = ( emP - emN ) / ( 2.0d0 * delta )
      write(*,*) "Force on k in dir via FD   : "
      write(*,*) " total energy = ",detdq_fd
      write(*,*) " polarization = ",depdq_fd
      write(*,*) " perm m.p.    = ",demdq_fd

      write(*,*) "dedci on k analytic  = ",ldedci
      write(*,*) "depdci on k analytic = ",ldepdci

      rpole(1,k) = q0
      call empole1
c     ep = ep0
c     em = em0
c     esum = et0

      write(*,*) "Finished with qtest"

      return
      end
c end of qtest subroutine


c     ################################
c     ##  Added by Steve Rick 2015  ##
c     ################################
c
c   ###################################################################
c    ##                                                             ##
c     ##  subroutine ectE  --  gets charge transfer contribition   ## 
c     ##       and charge transfter energy                 ## 
c    ##                                                             ##
c   ###################################################################
c
c
c    DCT subroutine ectE 
c       set up for water-water interations only
c    Ect only, no derivatives
c       for use with energy subroutine and ptest
c
c
      subroutine ectE 
      use sizes
      use deriv
      use energi
      use limits
      use atoms
      use boxes
      use mpole
      use potent
      use neigh
      use kmulti
      use virial
      implicit none
      integer i,j,kk,ii,l,m
      real*8 xr,yr,zr,rr,pi,sum1,sum2,dqt,ect
      real*8 xna,dxna
      real*8 tmp1,tmp2,tmp3
      real*8 xr0,yr0,zr0


      write(*,*) "ectE subroutine"

      pi = 4.d0*datan(1.d0)

      dqt=zdqt(1,1)+zdqt(2,1)+zdqt(3,1)
c     dqt = -0.020d0
c 	i.e. total CT
      ect = 0.d0
c     write(*,*) "Reset Ect = ",ect
c     write(*,*) "Position of 1 : ",x(1),y(1),z(1)

c do this for oxygen atoms only
c for each oxygen, look at neighbor list
      do i=1,n-1
         if(type(i).eq.1) then
            do kk = 1, nelst(i)
               j=elst(kk,i)
               if(type(j).eq.1) then
c if neighbor is also an oxygen, calculate distance
                  xr = x(j) - x(i)
                  yr = y(j) - y(i)
                  zr = z(j) - z(i)
                  call image (xr,yr,zr)
                  rr = xr*xr + yr*yr + zr*zr
c                 write(*,*) "r(O-O) = ",dsqrt(rr)

                  if(rr.lt.(rct2+1.200)**2) then
c check if water with oxygen atom j donates hydrogen bond to i
                     do l=1,2
                        xr = x(j+l) - x(i)
                        yr = y(j+l) - y(i)
                        zr = z(j+l) - z(i)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)

                        if(rr.ge.rct1.and.rr.le.rct2) then
                           xna = 
     & 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
c                          write(*,*) i,j+l,rr,xna

                        endif
c                       end if between rct1 and rct2

c molecule with oxygen atom i donates to hydrogen bond
                        xr = x(i+l) - x(j)
                        yr = y(i+l) - y(j)
                        zr = z(i+l) - z(j)
                        call image (xr,yr,zr)
                        rr  = dsqrt(xr**2+yr**2+zr**2)

                        if(rr.le.rct1) then
c add to charge transfer energy, no force (energy is constant)
                           ect=ect+muct*dqt+0.5d0*etact*dqt**2
                        endif

                        if(rr.gt.rct1.and.rr.le.rct2) then

                           xna = 
     & 0.5d0*(1.d0+dcos(pi*(rr-rct1)/(rct2-rct1)))
c                          write(*,*) i,j+l,rr,xna

c ect here only (not above) to avoid double counting
                           ect=ect+muct*(xna*dqt)
     & +0.5d0*etact*(xna*dqt)**2


                        endif

                     enddo
c                    end loop over H
                  endif
c                 end if oxygens within cutoff
               endif
c              end if also O_w
            enddo
c           end loop over nblist
         endif
c        end if O_w
      enddo
c     end of loop over all O_w
  
c     write(*,*) "Ect = ",ect
      em = em + ect


      return
      end
c end of ectE subroutine


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
      implicit none
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)

      write(*,*) "Start of gradient subrout."
c
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
c

c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
c     if (use_list)  call nblist
      if (use_list) then
         call nblist
c        write(*,*) "using neighbor list"
c MES : This is being called with mpole-list key
      end if
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
      write(*,*) "Calling energy and gradient subrout."
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
  
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar)  call empole1
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
      write(*,*) "Summing energy and deriv. contributions"
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
     &                      + dem(j,i) + dep(j,i) + der(j,i)
     &                      + des(j,i) + delf(j,i) + deg(j,i)
     &                      + dex(j,i)
            derivs(j,i) = desum(j,i)
         end do
      end do
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
      write(*,*) "End of gradient subrout."
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
c
      subroutine newcrg
      use sizes
      use limits
      use atoms
      use boxes
      use mpole
      use potent
      use neigh
      use kmulti
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
         dedci(i) = 0.d0
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

c     sum1 = 0.d0
c     sum2 = 0.d0
c total charge - should equal charge of system (usually zero)
c     do i=1,n
c      sum1 = sum1 + rpole0(i)
c      sum2 = sum2 + rpole(1,i)
c     enddo
c     write(*,*) "Total system charge : "
c     write(*,*) sum1, sum2

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
      real*8, allocatable :: nacti(:)
      real*8, allocatable :: ndcti(:)

      allocate (nacti(n))
      allocate (ndcti(n))

c     write(*,*) "ect1 subroutine"
c MES : This subroutine is called.

c     write(*,*) "CT variables :"
c     write(*,*) "xna            dxna            zr            rr "

c check virial
      write(*,*) "Virial pre-CT"
      do i=1,3
        do j=1,3
          write(*,*) vir(i,j)
        enddo
      enddo

c add forces to dem(i,j=1,2,3) (x,y,z)

      pi = 4.d0*datan(1.d0)
      do i=1,n
         nacti(i) = 0.d0
         ndcti(i) = 0.d0
      enddo

      dqt=zdqt(1,1)+zdqt(2,1)+zdqt(3,1)
c     dqt = -0.20d0
c 	i.e. total CT
      ect = 0.d0

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
                           dxna = 
     & -0.5d0*dsin(pi*(rr-rct1)/(rct2-rct1))
     & *pi/(rct2-rct1)

c MES : this set of dem passes gtest
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

c OR ect once to avoid double-counting?
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

      write(*,*) xna," ",dxna," ",zr," ",rr

c MES : this is correct for ect
                           ect=ect+muct*(xna*dqt)
     & +0.5d0*etact*(xna*dqt)**2

                           dem(1,i+l) = dem(1,i+l)
     & +(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt 
     & )*dxna*xr/rr
                           dem(2,i+l) = dem(2,i+l)
     & +(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt 
     & )*dxna*yr/rr
                           dem(3,i+l) = dem(3,i+l)
     & +(
     &   dedci(i)*zdqt(1,l+1)
     & +dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     & )*dxna*zr/rr
                           dem(1,j) = dem(1,j)
     & -(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt 
     & )*dxna*xr/rr
                           dem(2,j) = dem(2,j)
     & -(dedci(i)*zdqt(1,l+1)+dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt 
     & )*dxna*yr/rr
                           dem(3,j) = dem(3,j)
     & -(
     &  dedci(i)*zdqt(1,l+1)
     & +dedci(i+1)*zdqt(2,l+1)
     & +dedci(i+2)*zdqt(3,l+1)
     & +(muct+etact*xna*dqt)*dqt
     & )*dxna*zr/rr

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
      write(*,*) "Virial post-CT"
      do i=1,3
        do j=1,3
          write(*,*) vir(i,j)
        enddo
      enddo

      write(*,*) ect
      em = em + ect
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
c     "gtest" compares dE/dr from analytical and finite difference
c     methods.
c
c     Added by Marielle Soniat 2015 July
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
      implicit none
      integer i,j
      real*8 energy
      real*8 delta
      real*8 fx1,fy1,fz1,e0
      real*8 fxn,fyn,fzn
      real*8 fxneg,fyneg,fzneg
      real*8 fxpos,fypos,fzpos
      real*8 x0,y0,z0
      real*8 xneg,yneg,zneg
      real*8 xpos,ypos,zpos
      real*8 exneg,eyneg,ezneg
      real*8 expos,eypos,ezpos
      real*8 dedx_fd,dedy_fd,dedz_fd
      real*8 derivs(3,*)

c MES : crashes if more than one option is used at a time.

      write(*,*) "Entering gtest"

c     get initial energy and forces
      e0 = energy ()
c     derivs(j,i) = derivatives, j=direction, i=atom
      fx1 = derivs(1,1)
      fy1 = derivs(2,1)
      fz1 = derivs(3,1)
      fxn = derivs(1,n)
      fyn = derivs(2,n)
      fzn = derivs(3,n)
c MES : need to recalculate from within this test,
c	not use from previous md step
c 	store forces for i=1 and i=n
      

c     save initial position
      x0 = x(1)
      y0 = y(1)
      z0 = z(1)
c     write(*,*) "(x,y,z) of 1 : ",x0," ",y0," ",z0

c     change for position
      delta = 0.000001d0

c     for first atom, forces in x, y, and z
c     restore original position after calculation
c     xneg = x0 - delta
c     x(1) = xneg
c     exneg = energy ()

c     xpos = x0 + delta
c     x(1) = xpos
c     expos = energy ()

c     dedx_fd = ( expos - exneg ) / ( 2.0 * delta )
c     write(*,*) "Force in x on 1 via FD  = ",dedx_fd
c     write(*,*) "Force in x on 1 analtic = ",fx1
c     x(1) = x0

c     yneg = y0 - delta
c     y(1) = yneg
c     eyneg = energy ()

c     ypos = y0 + delta
c     y(1) = ypos
c     eypos = energy ()

c     dedy_fd = ( eypos - eyneg ) / ( 2.0 * delta )
c     write(*,*) "Force in y on 1 via FD  = ",dedy_fd
c     write(*,*) "Force in y on 1 analtic = ",fy1
c     y(1) = y0

c     write(*,*) "(x,y,z) of 1 : ",x(1)," ",y(1)," ",z(1)

c     zneg = z0 - delta
c     z(1) = zneg
c     ezneg = energy ()

c     zpos = z0 + delta
c     z(1) = zpos
c     ezpos = energy ()

c     dedz_fd = ( ezpos - ezneg ) / ( 2.0 * delta )
c     write(*,*) "Force in z on 1 via FD  = ",dedz_fd
c     write(*,*) "Force in z on 1 analtic = ",fz1
c     z(1) = z0

c     write(*,*) "(x,y,z) of 1 : ",x(1)," ",y(1)," ",z(1)

c     for final atom, forces in x, y, and z
      x0 = x(n)
      y0 = y(n)
      z0 = z(n)
      write(*,*) "(x,y,z) of n : ",x(n)," ",y(n)," ",z(n)

c     xneg = x0 - delta
c     x(n) = xneg
c     exneg = energy ()

c     xpos = x0 + delta
c     x(n) = xpos
c     expos = energy ()

c     dedx_fd = ( expos - exneg ) / ( 2.0 * delta )
c     write(*,*) "Force in x on n via FD  = ",dedx_fd
c     write(*,*) "Force in x on n analtic = ",fxn
c     x(n) = x0

c     yneg = y0 - delta
c     y(n) = yneg
c     eyneg = energy ()

c     ypos = y0 + delta
c     y(n) = ypos
c     eypos = energy ()

c     dedy_fd = ( eypos - eyneg ) / ( 2.0 * delta )
c     write(*,*) "Force in y on n via FD  = ",dedy_fd
c     write(*,*) "Force in y on n analtic = ",fyn
c     y(n) = y0

c     use_crgtr = .false.

      zneg = z0 - delta
      z(n) = zneg
      ezneg = energy ()

      zpos = z0 + delta
      z(n) = zpos
      ezpos = energy ()

      dedz_fd = ( ezpos - ezneg ) / ( 2.0 * delta )
      write(*,*) "Force in z on n via FD   = ",dedz_fd
      write(*,*) "Force in z on n analytic = ",fzn
      z(n) = z0

c     use_crgtr = .true.

      write(*,*) "Finished with gtest"

      return
      end
c end of gtest subroutine



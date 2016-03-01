c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2010 by T. Darden, D. Gohara & Jay W. Ponder  ##
c     ##                      All Rights Reserved                     ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below implement various B-spline and coordinate  ##
c     ##  manipulations for particle mesh Ewald summation; spatial  ##
c     ##  grid assignment by David Gohara; modified from original   ##
c     ##  PME code by Thomas Darden, NIEHS, Research Triangle, NC   ##
c     ##                                                            ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine bspline_fill  --  get PME B-spline coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bspline_fill" finds B-spline coefficients and derivatives
c     for PME atomic sites along the fractional coordinate axes
c
c
      subroutine bspline_fill
      use sizes
      use atoms
      use boxes
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(igrid))  allocate (igrid(3,n))
      end if
c
c     offset used to shift sites off exact lattice bounds
c
      eps = 1.0d-8
c
      write(*,*) "get the B-spline coefficients for each atomic site"
c
      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bsplgen (w,thetai1(1,1,i))
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bsplgen (w,thetai2(1,1,i))
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(3,i) = ifr - bsorder
         call bsplgen (w,thetai3(1,1,i))
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bsplgen" gets B-spline coefficients and derivatives for
c     a single PME atomic site along a particular direction
c
c
      subroutine bsplgen (w,thetai)
      use pme
      use potent
      implicit none
      integer i,j,k
      integer level
      real*8 w,denom
      real*8 thetai(4,*)
c
c
c     set B-spline depth for partial charges or multipoles
c
      level = 2
      if (use_mpole .or. use_polar)  level = 4
c
c     initialization to get to 2nd order recursion
c
      bsbuild(2,2) = w
      bsbuild(2,1) = 1.0d0 - w
c
c     perform one pass to get to 3rd order recursion
c
      bsbuild(3,3) = 0.5d0 * w * bsbuild(2,2)
      bsbuild(3,2) = 0.5d0 * ((1.0d0+w)*bsbuild(2,1)
     &                       +(2.0d0-w)*bsbuild(2,2))
      bsbuild(3,1) = 0.5d0 * (1.0d0-w) * bsbuild(2,1)
c
c     compute standard B-spline recursion to desired order
c
      do i = 4, bsorder
         k = i - 1
         denom = 1.0d0 / dble(k)
         bsbuild(i,i) = denom * w * bsbuild(k,k)
         do j = 1, i-2
            bsbuild(i,i-j) = denom * ((w+dble(j))*bsbuild(k,i-j-1)
     &                               +(dble(i-j)-w)*bsbuild(k,i-j))
         end do
         bsbuild(i,1) = denom * (1.0d0-w) * bsbuild(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
      do i = bsorder-1, 2, -1
         bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
      end do
      bsbuild(k,1) = -bsbuild(k,1)
c
c     get coefficients for the B-spline second derivative
c
      if (level .eq. 4) then
         k = bsorder - 2
         bsbuild(k,bsorder-1) = bsbuild(k,bsorder-2)
         do i = bsorder-2, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
         do i = bsorder-1, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
c
c     get coefficients for the B-spline third derivative
c
         k = bsorder - 3
         bsbuild(k,bsorder-2) = bsbuild(k,bsorder-3)
         do i = bsorder-3, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder-1) = bsbuild(k,bsorder-2)
         do i = bsorder-2, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
         do i = bsorder-1, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
      end if
c
c     copy coefficients from temporary to permanent storage
c
      do i = 1, bsorder
         do j = 1, level
            thetai(j,i) = bsbuild(bsorder-j+1,i)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine table_fill  --  spatial chunks for each site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "table_fill" constructs an array which stores the spatial
c     regions of the particle mesh Ewald grid with contributions
c     from each electrostatic site
c
c
      subroutine table_fill
      use sizes
      use atoms
      use chunks
      use pme
      implicit none
      integer i,k
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      logical negx,negy,negz
      logical posx,posy,posz
      logical midx,midy,midz
c
c     write(*,*) "table_fill"
c
c     zero out the PME table marking chunks per site
c
      do k = 1, nchunk
         do i = 1, n
            pmetable(i,k) = 0
         end do
      end do
c
c     loop over sites to find the spatial chunks for each
c
      do i = 1, n
         nearpt(1) = igrid(1,i) + grdoff
         nearpt(2) = igrid(2,i) + grdoff
         nearpt(3) = igrid(3,i) + grdoff
         if (nearpt(1) .lt. 1) then
            nearpt(1) = mod(nearpt(1),nfft1) + nfft1
         else if (nearpt(1) .gt. nfft1) then
            nearpt(1) = mod(nearpt(1),nfft1)
         end if
         if (nearpt(2) .lt. 1) then
            nearpt(2) = mod(nearpt(2),nfft2) + nfft2
         else if (nearpt(2) .gt. nfft2) then
            nearpt(2) = mod(nearpt(2),nfft2)
         end if
         if (nearpt(3) .lt. 1) then
            nearpt(3) = mod(nearpt(3),nfft3) + nfft3
         else if (nearpt(3) .gt. nfft3) then
            nearpt(3) = mod(nearpt(3),nfft3)
         end if
         abound(1) = nearpt(1) - nlpts
         abound(2) = nearpt(1) + nrpts
         abound(3) = nearpt(2) - nlpts
         abound(4) = nearpt(2) + nrpts
         abound(5) = nearpt(3) - nlpts
         abound(6) = nearpt(3) + nrpts
         cid(1) = (nearpt(1)-1)/ngrd1 + 1
         cid(2) = (nearpt(2)-1)/ngrd2 + 1
         cid(3) = (nearpt(3)-1)/ngrd3 + 1
         cbound(1) = (cid(1)-1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = (cid(2)-1)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = (cid(3)-1)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
c
c     set and store central chunk where the site is located
c
         k = (cid(3)-1)*nchk1*nchk2 + (cid(2)-1)*nchk1 + cid(1)
         pmetable(i,k) = 1
c
c     flags for atom bounds to left or right of central chunk
c
         negx = (abound(1) .lt. cbound(1))
         negy = (abound(3) .lt. cbound(3))
         negz = (abound(5) .lt. cbound(5))
         posx = (abound(2) .gt. cbound(2))
         posy = (abound(4) .gt. cbound(4))
         posz = (abound(6) .gt. cbound(6))
c
c     flags for atom bounds fully inside the central chunk
c
         midx = (.not.negx .and. .not.posx)
         midy = (.not.negy .and. .not.posy)
         midz = (.not.negz .and. .not.posz)
         if (midx .and. midy .and. midz)  goto 10
c
c     flags for atom bounds that overlap the central chunk
c
         midx = (.not.negx .or. .not.posx)
         midy = (.not.negy .or. .not.posy)
         midz = (.not.negz .or. .not.posz)
c
c     check for overlap with any of the neighboring chunks
c
         if (midx .and. midy .and. negz)  call setchunk (i,cid,0,0,-1)
         if (midx .and. midy .and. posz)  call setchunk (i,cid,0,0,1)
         if (midx .and. negy .and. midz)  call setchunk (i,cid,0,-1,0)
         if (midx .and. posy .and. midz)  call setchunk (i,cid,0,1,0)
         if (negx .and. midy .and. midz)  call setchunk (i,cid,-1,0,0)
         if (posx .and. midy .and. midz)  call setchunk (i,cid,1,0,0)
         if (midx .and. negy .and. negz)  call setchunk (i,cid,0,-1,-1)
         if (midx .and. negy .and. posz)  call setchunk (i,cid,0,-1,1)
         if (midx .and. posy .and. negz)  call setchunk (i,cid,0,1,-1)
         if (midx .and. posy .and. posz)  call setchunk (i,cid,0,1,1)
         if (negx .and. midy .and. negz)  call setchunk (i,cid,-1,0,-1)
         if (negx .and. midy .and. posz)  call setchunk (i,cid,-1,0,1)
         if (posx .and. midy .and. negz)  call setchunk (i,cid,1,0,-1)
         if (posx .and. midy .and. posz)  call setchunk (i,cid,1,0,1)
         if (negx .and. negy .and. midz)  call setchunk (i,cid,-1,-1,0)
         if (negx .and. posy .and. midz)  call setchunk (i,cid,-1,1,0)
         if (posx .and. negy .and. midz)  call setchunk (i,cid,1,-1,0)
         if (posx .and. posy .and. midz)  call setchunk (i,cid,1,1,0)
         if (negx .and. negy .and. negz)  call setchunk (i,cid,-1,-1,-1)
         if (negx .and. negy .and. posz)  call setchunk (i,cid,-1,-1,1)
         if (negx .and. posy .and. negz)  call setchunk (i,cid,-1,1,-1)
         if (posx .and. negy .and. negz)  call setchunk (i,cid,1,-1,-1)
         if (negx .and. posy .and. posz)  call setchunk (i,cid,-1,1,1)
         if (posx .and. negy .and. posz)  call setchunk (i,cid,1,-1,1)
         if (posx .and. posy .and. negz)  call setchunk (i,cid,1,1,-1)
         if (posx .and. posy .and. posz)  call setchunk (i,cid,1,1,1)
   10    continue
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setchunk  --  site overlaps neighboring chunk  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setchunk" marks a chunk in the PME spatial table which is
c     overlapped by the B-splines for an electrostatic site
c
c
      subroutine setchunk (i,cid,off1,off2,off3)
      use sizes
      use chunks
      use pme
      implicit none
      integer i,k
      integer off1,off2,off3
      integer cid(3),temp(3)
c
c
c     mark neighboring chunk overlapped by an electrostatic site
c
      temp(1) = cid(1) + off1
      if (temp(1) .lt. 1)  temp(1) = nchk1
      if (temp(1) .gt. nchk1)  temp(1) = 1
      temp(2) = cid(2) + off2
      if (temp(2) .lt. 1)  temp(2) = nchk2
      if (temp(2) .gt. nchk2)  temp(2) = 1
      temp(3) = cid(3) + off3
      if (temp(3) .lt. 1)  temp(3) = nchk3
      if (temp(3) .gt. nchk3)  temp(3) = 1
      k = (temp(3)-1)*nchk1*nchk2 + (temp(2)-1)*nchk1 + temp(1)
      pmetable(i,k) = 1
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_pchg  --  put partial charges on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_pchg" places the fractional atomic partial charges onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_pchg
      use sizes
      use atoms
      use charge
      use chunks
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 term
c
      write(*,*) "Entered grid_pchg"
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,jsite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,u0,term,t0)
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, nion
            iatm = iion(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm) * pchg(isite)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     term = v0 * u0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term*t0
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      write(*,*) "Done with grid_pchg"
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine grid_mpole  --  put multipoles on PME grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "grid_mpole" places the fractional atomic multipoles onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_mpole (fmp)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,*)
c
c     write(*,*) "Entered grid_mpole"
c     write(*,*) "    nchunk = ",nchunk
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
               do isite = 1,npole 
                 dqgrdci(1,i,j,k,isite) = 0.0d0
                 dqgrdci(2,i,j,k,isite) = 0.0d0
               end do
            end do
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,v2,u0,u1,u2,term0,term1,term2,t0,t1,t2)
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, npole
            iatm = ipole(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  v2 = thetai3(3,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     u2 = thetai2(3,m,iatm)
                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
     &                          + fmp(9,isite)*u0*v1
                     term2 = fmp(5,isite) * u0 * v0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        t2 = thetai1(3,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
     &                                      + term1*t1 + term2*t2
c DCT    dfmpdci = 1
c    dqgrdci(1,i,j,k,isite) is the change in qgrid(1,i,j,k) 
c    due to a change in charge at isite
                        dqgrdci(1,i,j,k,isite) = dqgrdci(1,i,j,k,isite) 
     & + u0*v0*t0
c                       if ( i.eq.11 .and. j.eq.8 .and. k.eq.7 ) then
c                         write(*,*) i,j,k,isite
c                         write(*,*) qgrid(1,i,j,k),
c    & dqgrdci2(1,i,j,k), dqgrdci(1,i,j,k,isite)
c                       endif
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

c     write(*,*) "outside OMP parallel"
c     do k = 1, nfft3
c        do j = 1, nfft2
c           do i = 1, nfft1
c              if(qgrid(1,i,j,k).ne.0.000000d0) then
c                write(*,*) i,j,k,dqgrdci(1,i,j,k)
c                write(*,*) i,j,k,qgrid(1,i,j,k)
c              end if
c           end do
c        end do
c     end do
c     write(*,*) "11,8,7",qgrid(1,11,8,7),dqgrdci(1,11,8,7)


c     write(*,*) "Done with grid_mpole"
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_uind  --  put induced dipoles on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_uind" places the fractional induced dipoles onto the
c     particle mesh Ewald grid
c
c
c     subroutine grid_uind (fuind,fuinp)
      subroutine grid_uind (fuind,fuinp,dfuinddci,dfuinpdci)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,jsite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8, allocatable :: dterm01dci(:)
      real*8, allocatable :: dterm11dci(:)
      real*8, allocatable :: dterm02dci(:)
      real*8, allocatable :: dterm12dci(:)
c     real*8, allocatable :: dqgrdciX(:,:,:,:,:,:)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 term01,term11
      real*8 term02,term12
      real*8 fuind(3,*)
      real*8 fuinp(3,*)
      real*8 dfuinddci(3,npole,*)
      real*8 dfuinpdci(3,npole,*)

      allocate(dterm01dci(npole))
      allocate(dterm11dci(npole))
      allocate(dterm02dci(npole))
      allocate(dterm12dci(npole))
c     allocate(dqgrdciX(2,nfft1,nfft2,nfft3,npole,npole))
c
      write(*,*) "Entered grid_uind"
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
               do isite = 1, npole
                 dqgrdci(1,i,j,k,isite) = 0.0d0
                 dqgrdci(2,i,j,k,isite) = 0.0d0
                 do jsite = 1, npole
                   dqgrdciX(1,i,j,k,isite,jsite) = 0.0d0
                   dqgrdciX(2,i,j,k,isite,jsite) = 0.0d0
                 end do
               end do
            end do
         end do
      end do
c     write(*,*) dqgrdci(1,10,8,7,3),dqgrdciX(1,10,8,7,3,4)
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,u0,u1,term01,term11,term02,term12,t0,t1)
!$OMP DO
c
c     put the induced dipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, npole
            iatm = ipole(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     term01 = fuind(2,isite)*u1*v0
     &                           + fuind(3,isite)*u0*v1
                     term11 = fuind(1,isite)*u0*v0
                     term02 = fuinp(2,isite)*u1*v0
     &                           + fuinp(3,isite)*u0*v1
                     term12 = fuinp(1,isite)*u0*v0
c DCT
                     do jsite = 1, npole
                        dterm01dci(jsite) = 
     &                             dfuinddci(2,isite,jsite)*u1*v0
     &                           + dfuinddci(3,isite,jsite)*u0*v1
                        dterm11dci(jsite) = 
     &                             dfuinddci(1,isite,jsite)*u0*v0
                        dterm02dci(jsite) = 
     &                             dfuinpdci(2,isite,jsite)*u1*v0
     &                           + dfuinpdci(3,isite,jsite)*u0*v1
                        dterm12dci(jsite) = 
     &                             dfuinpdci(1,isite,jsite)*u0*v0
                     end do
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
c                       if (i.eq.10 .and. j.eq.8 .and. k.eq.7) then
c                         write(*,*) "Pre ",i,j,k,isite
c                         write(*,*) qgrid(1,i,j,k),
c    & dqgrdci(1,i,j,k,isite),dqgrdciX(1,i,j,k,isite,4)
c                       end if 
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term01*t0
     &                                      + term11*t1
                        qgrid(2,i,j,k) = qgrid(2,i,j,k) + term02*t0
     &                                      + term12*t1
c DCT
c MES debug : here dqgrdci is actually the charge at qgrid due to isite
c    Note: this is DIFFERENT from use of dqgrdci in mpole pme rouitnes.
                           dqgrdci(1,i,j,k,isite) = 
     & dqgrdci(1,i,j,k,isite) + term01*t0 + term11*t1
                           dqgrdci(2,i,j,k,isite) =
     & dqgrdci(2,i,j,k,isite) + term02*t0 + term12*t1

                        do jsite = 1, npole
                           dqgrdciX(1,i,j,k,isite,jsite) = 
     & + dterm01dci(jsite)*t0 + dterm11dci(jsite)*t1
                           dqgrdciX(2,i,j,k,isite,jsite) = 
     & + dterm02dci(jsite)*t0 + dterm12dci(jsite)*t1
                        end do
c                       if (i.eq.10 .and. j.eq.8 .and. k.eq.7) then
c                         write(*,*) "Post ",i,j,k,isite
c                         write(*,*) term01,dterm01dci(4)
c                         write(*,*) qgrid(1,i,j,k),
c    & dqgrdci(1,i,j,k,isite),dqgrdciX(1,i,j,k,isite,4)
c                       end if
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

      deallocate(dterm01dci)
      deallocate(dterm11dci)
      deallocate(dterm02dci)
      deallocate(dterm12dci)
c     deallocate(dqgrdciX)

      write(*,*) "Done with grid_uind"
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine adjust  --  alter site bounds for the PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "adjust" modifies site bounds on the PME grid and returns
c     an offset into the B-spline coefficient arrays
c
c
      subroutine adjust (offset,nfft,nchk,amin,amax,cmin,cmax)
      implicit none
      integer offset
      integer nfft,nchk
      integer amin,amax
      integer cmin,cmax
c
c
c     modify grid offset and bounds for site at edge of chunk
c
      offset = 0
      if (nchk .ne. 1) then
         if (amin.lt.cmin .or. amax.gt.cmax) then
            if (amin.lt.1 .or. amax.gt.nfft) then
               if (cmin .eq. 1) then
                  offset = 1 - amin
                  amin = 1
               else if (cmax .eq. nfft) then
                  amax = nfft
                  amin = amin + nfft
               end if
            else
               if (cmin .gt. amin) then
                  offset = cmin - amin
                  amin = cmin
               else
                  amax = cmax
               end if
            end if
         end if
      end if
      offset = offset + 1 - amin
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpole  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpole" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpole (fphi)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 dtuv000dqg,dt0,dtu00, dtuv000
      real*8 fphi(20,*)
c
      write(*,*) "Entered fphi_mpole"
c     write(*,*) "bsorder = ",bsorder
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,fphi)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
c        location of isite on grid
c        write(*,*) "isite ",isite,"   iatm ",iatm,"   igrd0 ",igrd0
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,iatm)
                  t1 = t1 + tq*thetai1(2,it1,iatm)
                  t2 = t2 + tq*thetai1(3,it1,iatm)
                  t3 = t3 + tq*thetai1(4,it1,iatm)
               end do
c              end of it1
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
c           end of it2
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
c        end of it3
c
c        mpole of atom/isite calc from grid
         fphi(1,isite) = tuv000
         fphi(2,isite) = tuv100
         fphi(3,isite) = tuv010
         fphi(4,isite) = tuv001
         fphi(5,isite) = tuv200
         fphi(6,isite) = tuv020
         fphi(7,isite) = tuv002
         fphi(8,isite) = tuv110
         fphi(9,isite) = tuv101
         fphi(10,isite) = tuv011
         fphi(11,isite) = tuv300
         fphi(12,isite) = tuv030
         fphi(13,isite) = tuv003
         fphi(14,isite) = tuv210
         fphi(15,isite) = tuv201
         fphi(16,isite) = tuv120
         fphi(17,isite) = tuv021
         fphi(18,isite) = tuv102
         fphi(19,isite) = tuv012
         fphi(20,isite) = tuv111
      end do
c     end of isite
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

      write(*,*) "end of fphi_mpole"

      return
      end
c     end of fphi_mpole
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fphi_uind  --  induced potential from grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fphi_uind" extracts the induced dipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3
      real*8 t0_1,t0_2,t1_1,t1_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu30,tu21,tu12,tu03
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fdip_phi1(10,*)
      real*8 fdip_phi2(10,*)
      real*8 fdip_sum_phi(20,*)
c
      write(*,*) "Entered fphi_uind"
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,
!$OMP& igrid,bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv100_1 = 0.0d0
         tuv010_1 = 0.0d0
         tuv001_1 = 0.0d0
         tuv200_1 = 0.0d0
         tuv020_1 = 0.0d0
         tuv002_1 = 0.0d0
         tuv110_1 = 0.0d0
         tuv101_1 = 0.0d0
         tuv011_1 = 0.0d0
         tuv100_2 = 0.0d0
         tuv010_2 = 0.0d0
         tuv001_2 = 0.0d0
         tuv200_2 = 0.0d0
         tuv020_2 = 0.0d0
         tuv002_2 = 0.0d0
         tuv110_2 = 0.0d0
         tuv101_2 = 0.0d0
         tuv011_2 = 0.0d0
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00_1 = 0.0d0
            tu01_1 = 0.0d0
            tu10_1 = 0.0d0
            tu20_1 = 0.0d0
            tu11_1 = 0.0d0
            tu02_1 = 0.0d0
            tu00_2 = 0.0d0
            tu01_2 = 0.0d0
            tu10_2 = 0.0d0
            tu20_2 = 0.0d0
            tu11_2 = 0.0d0
            tu02_2 = 0.0d0
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0_1 = 0.0d0
               t1_1 = 0.0d0
               t2_1 = 0.0d0
               t0_2 = 0.0d0
               t1_2 = 0.0d0
               t2_2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq_1 = qgrid(1,i,j,k)
                  tq_2 = qgrid(2,i,j,k)
                  t0_1 = t0_1 + tq_1*thetai1(1,it1,iatm)
                  t1_1 = t1_1 + tq_1*thetai1(2,it1,iatm)
                  t2_1 = t2_1 + tq_1*thetai1(3,it1,iatm)
                  t0_2 = t0_2 + tq_2*thetai1(1,it1,iatm)
                  t1_2 = t1_2 + tq_2*thetai1(2,it1,iatm)
                  t2_2 = t2_2 + tq_2*thetai1(3,it1,iatm)
                  t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,iatm)
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
               t0 = t0_1 + t0_2
               t1 = t1_1 + t1_2
               t2 = t2_1 + t2_2
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fdip_phi1(2,isite) = tuv100_1
         fdip_phi1(3,isite) = tuv010_1
         fdip_phi1(4,isite) = tuv001_1
         fdip_phi1(5,isite) = tuv200_1
         fdip_phi1(6,isite) = tuv020_1
         fdip_phi1(7,isite) = tuv002_1
         fdip_phi1(8,isite) = tuv110_1
         fdip_phi1(9,isite) = tuv101_1
         fdip_phi1(10,isite) = tuv011_1
         fdip_phi2(2,isite) = tuv100_2
         fdip_phi2(3,isite) = tuv010_2
         fdip_phi2(4,isite) = tuv001_2
         fdip_phi2(5,isite) = tuv200_2
         fdip_phi2(6,isite) = tuv020_2
         fdip_phi2(7,isite) = tuv002_2
         fdip_phi2(8,isite) = tuv110_2
         fdip_phi2(9,isite) = tuv101_2
         fdip_phi2(10,isite) = tuv011_2
         fdip_sum_phi(1,isite) = tuv000
         fdip_sum_phi(2,isite) = tuv100
         fdip_sum_phi(3,isite) = tuv010
         fdip_sum_phi(4,isite) = tuv001
         fdip_sum_phi(5,isite) = tuv200
         fdip_sum_phi(6,isite) = tuv020
         fdip_sum_phi(7,isite) = tuv002
         fdip_sum_phi(8,isite) = tuv110
         fdip_sum_phi(9,isite) = tuv101
         fdip_sum_phi(10,isite) = tuv011
         fdip_sum_phi(11,isite) = tuv300
         fdip_sum_phi(12,isite) = tuv030
         fdip_sum_phi(13,isite) = tuv003
         fdip_sum_phi(14,isite) = tuv210
         fdip_sum_phi(15,isite) = tuv201
         fdip_sum_phi(16,isite) = tuv120
         fdip_sum_phi(17,isite) = tuv021
         fdip_sum_phi(18,isite) = tuv102
         fdip_sum_phi(19,isite) = tuv012
         fdip_sum_phi(20,isite) = tuv111
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      write(*,*) "Done with fphi_uind"
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cmp_to_fmp" transforms the atomic multipoles from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp (cmp,fmp)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ctf(10,10)
      real*8 cmp(10,*)
      real*8 fmp(10,*)
c
c     write(*,*) "Entered cmp_to_fmp"
c
c     find the matrix to convert Cartesian to fractional
c
c     write(*,*) "calling cart_to_frac"
      call cart_to_frac (ctf)
c     write(*,*) "    ctf(1,1) = ",ctf(1,1)
c
c     apply the transformation to get the fractional multipoles
c
      do i = 1, npole
         fmp(1,i) = ctf(1,1) * cmp(1,i)
c        write(*,*) i,fmp(1,i),cmp(1,i)
         do j = 2, 4
            fmp(j,i) = 0.0d0
            do k = 2, 4
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
         do j = 5, 10
            fmp(j,i) = 0.0d0
            do k = 5, 10
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
      end do
c     write(*,*) "End of cmp_to_fmp"
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine cart_to_frac  --  Cartesian to fractional  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "cart_to_frac" computes a transformation matrix to convert
c     a multipole object in Cartesian coordinates to fractional
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine cart_to_frac (ctf)
      use sizes
      use boxes
      use pme
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ctf(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
c
c     get the Cartesian to fractional conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ctf(j,i) = 0.0d0
         end do
      end do
      ctf(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ctf(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphi  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_to_cphi" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c
c
      subroutine fphi_to_cphi (fphi,cphi)
      use sizes
      use mpole
      implicit none
      integer i,j,k,ii
      real*8 ftc(10,10)
      real*8 cphi(10,*)
      real*8 fphi(20,*)
c     real*8 dfphidciX(20,npole,*)
c     real*8 dcphidciX(20,npole,*)
c
c     write(*,*) "Entered fphi_to_cphi"
c
c     find the matrix to convert fractional to Cartesian
c
c     write(*,*) "Calling frac_to_cart"
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      do i = 1, npole
         cphi(1,i) = ftc(1,1) * fphi(1,i)
         do j = 2, 4
            cphi(j,i) = 0.0d0
            do k = 2, 4
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
c           end of k
         end do
c        end of j

         do j = 5, 10
            cphi(j,i) = 0.0d0
            do k = 5, 10
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
c           end of k
         end do
c        end of j

      end do
c     end of i

c     write(*,*) "End of fphi_to_cphi"
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine frac_to_cart  --  fractional to Cartesian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "frac_to_cart" computes a transformation matrix to convert
c     a multipole object in fraction coordinates to Cartesian
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine frac_to_cart (ftc)
      use sizes
      use boxes
      use pme
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ftc(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
c
c     get the fractional to Cartesian conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ftc(j,i) = 0.0d0
         end do
      end do
      ftc(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ftc(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(k,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = 2.0d0 * a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(m,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
         end do
      end do
      return
      end



c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpoleCT  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpoleCT" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpoleCT (fphi,dfphidciX)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm,jsite
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111

      real*8, allocatable :: xdt0(:)
      real*8, allocatable :: xdt1(:)
      real*8, allocatable :: xdt2(:)
      real*8, allocatable :: xdt3(:)
      real*8, allocatable :: xdtu00(:)
      real*8, allocatable :: xdtu10(:)
      real*8, allocatable :: xdtu01(:)
      real*8, allocatable :: xdtu20(:)
      real*8, allocatable :: xdtu11(:)
      real*8, allocatable :: xdtu02(:)
      real*8, allocatable :: xdtu21(:)
      real*8, allocatable :: xdtu12(:)
      real*8, allocatable :: xdtu30(:)
      real*8, allocatable :: xdtu03(:)
      real*8, allocatable :: xdtuv000(:)
      real*8, allocatable :: xdtuv100(:)
      real*8, allocatable :: xdtuv010(:)
      real*8, allocatable :: xdtuv001(:)
      real*8, allocatable :: xdtuv200(:)
      real*8, allocatable :: xdtuv020(:)
      real*8, allocatable :: xdtuv002(:)
      real*8, allocatable :: xdtuv110(:)
      real*8, allocatable :: xdtuv101(:)
      real*8, allocatable :: xdtuv011(:)
      real*8, allocatable :: xdtuv300(:)
      real*8, allocatable :: xdtuv030(:)
      real*8, allocatable :: xdtuv003(:)
      real*8, allocatable :: xdtuv210(:)
      real*8, allocatable :: xdtuv201(:)
      real*8, allocatable :: xdtuv120(:)
      real*8, allocatable :: xdtuv021(:)
      real*8, allocatable :: xdtuv102(:)
      real*8, allocatable :: xdtuv012(:)
      real*8, allocatable :: xdtuv111(:)

      real*8 fphi(20,*)
      real*8 dfphidciX(20,npole,npole)
c
      write(*,*) "Entered fphi_mpoleCT"
c     write(*,*) "bsorder = ",bsorder
c     write(*,*) ""

      allocate (xdt0(npole))
      allocate (xdt1(npole))
      allocate (xdt2(npole))
      allocate (xdt3(npole))
      allocate (xdtu00(npole))
      allocate (xdtu10(npole))
      allocate (xdtu01(npole))
      allocate (xdtu20(npole))
      allocate (xdtu11(npole))
      allocate (xdtu02(npole))
      allocate (xdtu21(npole))
      allocate (xdtu12(npole))
      allocate (xdtu30(npole))
      allocate (xdtu03(npole))
      allocate (xdtuv000(npole))
      allocate (xdtuv100(npole))
      allocate (xdtuv010(npole))
      allocate (xdtuv001(npole))
      allocate (xdtuv200(npole))
      allocate (xdtuv020(npole))
      allocate (xdtuv002(npole))
      allocate (xdtuv110(npole))
      allocate (xdtuv101(npole))
      allocate (xdtuv011(npole))
      allocate (xdtuv300(npole))
      allocate (xdtuv030(npole))
      allocate (xdtuv003(npole))
      allocate (xdtuv210(npole))
      allocate (xdtuv201(npole))
      allocate (xdtuv120(npole))
      allocate (xdtuv021(npole))
      allocate (xdtuv102(npole))
      allocate (xdtuv012(npole))
      allocate (xdtuv111(npole))

c
c     set OpenMP directives for the major loop structure
c
c MES added dqgrdci and dfphidci to OMP shared var.
!$OMP PARALLEL default(private) shared(npole,ipole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,fphi,
!$OMP& dqgrdci,dfphidciX)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
c        location of isite on grid
c        write(*,*) "isite ",isite,"   iatm ",iatm,"   igrd0 ",igrd0

         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0

c xdtuv000 is the change in tuv000 at isite (implied) 
c    due to change in charge at jsite
         do jsite = 1, npole
           xdtuv000(jsite) = 0.0d0
           xdtuv001(jsite) = 0.0d0
           xdtuv010(jsite) = 0.0d0
           xdtuv100(jsite) = 0.0d0
           xdtuv200(jsite) = 0.0d0
           xdtuv020(jsite) = 0.0d0
           xdtuv002(jsite) = 0.0d0
           xdtuv110(jsite) = 0.0d0
           xdtuv101(jsite) = 0.0d0
           xdtuv011(jsite) = 0.0d0
           xdtuv300(jsite) = 0.0d0
           xdtuv030(jsite) = 0.0d0
           xdtuv003(jsite) = 0.0d0
           xdtuv210(jsite) = 0.0d0
           xdtuv201(jsite) = 0.0d0
           xdtuv120(jsite) = 0.0d0
           xdtuv021(jsite) = 0.0d0
           xdtuv102(jsite) = 0.0d0
           xdtuv012(jsite) = 0.0d0
           xdtuv111(jsite) = 0.0d0
         end do

         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)

            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0

            do jsite = 1, npole
              xdtu00(jsite) = 0.0d0
              xdtu10(jsite) = 0.0d0
              xdtu01(jsite) = 0.0d0
              xdtu20(jsite) = 0.0d0
              xdtu11(jsite) = 0.0d0
              xdtu02(jsite) = 0.0d0
              xdtu30(jsite) = 0.0d0
              xdtu21(jsite) = 0.0d0
              xdtu12(jsite) = 0.0d0
              xdtu03(jsite) = 0.0d0
            end do

            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)

               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0

               do jsite = 1, npole
                 xdt0(jsite) = 0.0d0
                 xdt1(jsite) = 0.0d0
                 xdt2(jsite) = 0.0d0
                 xdt3(jsite) = 0.0d0
               end do

               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2

c                 if ( i.eq.11 .and. j.eq.8 .and. k.eq.7 ) then
c                   write(*,*) i,j,k,isite
c                   write(*,*) qgrid(1,i,j,k),dqgrdci2(1,i,j,k)
c                   do jsite = 1 , npole
c                     write(*,*) dqgrdci(1,i,j,k,jsite)
c                   end do
c                 endif

                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,iatm)
                  t1 = t1 + tq*thetai1(2,it1,iatm)
                  t2 = t2 + tq*thetai1(3,it1,iatm)
                  t3 = t3 + tq*thetai1(4,it1,iatm)

                  do jsite = 1, npole
                    xdt0(jsite) = xdt0(jsite)
     &                + thetai1(1,it1,iatm)*dqgrdci(1,i,j,k,jsite)
                    xdt1(jsite) = xdt1(jsite)
     &                + thetai1(2,it1,iatm)*dqgrdci(1,i,j,k,jsite)
                    xdt2(jsite) = xdt2(jsite)
     &                + thetai1(3,it1,iatm)*dqgrdci(1,i,j,k,jsite)
                    xdt3(jsite) = xdt3(jsite)
     &                + thetai1(4,it1,iatm)*dqgrdci(1,i,j,k,jsite)
                  end do

               end do
c              end of it1

               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3

               do jsite = 1, npole
                 xdtu00(jsite) = xdtu00(jsite) + xdt0(jsite)*u0
                 xdtu10(jsite) = xdtu10(jsite) + xdt1(jsite)*u0
                 xdtu01(jsite) = xdtu01(jsite) + xdt0(jsite)*u1
                 xdtu20(jsite) = xdtu20(jsite) + xdt2(jsite)*u0
                 xdtu11(jsite) = xdtu11(jsite) + xdt1(jsite)*u1
                 xdtu02(jsite) = xdtu02(jsite) + xdt0(jsite)*u2
                 xdtu30(jsite) = xdtu30(jsite) + xdt3(jsite)*u0
                 xdtu21(jsite) = xdtu21(jsite) + xdt2(jsite)*u1
                 xdtu12(jsite) = xdtu12(jsite) + xdt1(jsite)*u2
                 xdtu03(jsite) = xdtu03(jsite) + xdt0(jsite)*u3
               end do

            end do
c           end of it2

            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1

            do jsite = 1, npole
              xdtuv000(jsite) = xdtuv000(jsite) + xdtu00(jsite)*v0
              xdtuv100(jsite) = xdtuv100(jsite) + xdtu10(jsite)*v0
              xdtuv010(jsite) = xdtuv010(jsite) + xdtu01(jsite)*v0
              xdtuv001(jsite) = xdtuv001(jsite) + xdtu00(jsite)*v1
              xdtuv200(jsite) = xdtuv200(jsite) + xdtu20(jsite)*v0
              xdtuv020(jsite) = xdtuv020(jsite) + xdtu02(jsite)*v0
              xdtuv002(jsite) = xdtuv002(jsite) + xdtu00(jsite)*v2
              xdtuv110(jsite) = xdtuv110(jsite) + xdtu11(jsite)*v0
              xdtuv101(jsite) = xdtuv101(jsite) + xdtu10(jsite)*v1
              xdtuv011(jsite) = xdtuv011(jsite) + xdtu01(jsite)*v1
              xdtuv300(jsite) = xdtuv300(jsite) + xdtu30(jsite)*v0
              xdtuv030(jsite) = xdtuv030(jsite) + xdtu03(jsite)*v0
              xdtuv003(jsite) = xdtuv003(jsite) + xdtu00(jsite)*v3
              xdtuv210(jsite) = xdtuv210(jsite) + xdtu21(jsite)*v0
              xdtuv201(jsite) = xdtuv201(jsite) + xdtu20(jsite)*v1
              xdtuv120(jsite) = xdtuv120(jsite) + xdtu12(jsite)*v0
              xdtuv021(jsite) = xdtuv021(jsite) + xdtu02(jsite)*v1
              xdtuv102(jsite) = xdtuv102(jsite) + xdtu10(jsite)*v2
              xdtuv012(jsite) = xdtuv012(jsite) + xdtu01(jsite)*v2
              xdtuv111(jsite) = xdtuv111(jsite) + xdtu11(jsite)*v1
            end do

         end do
c        end of it3
c
c        mpole of atom/isite calc from grid
         fphi(1,isite) = tuv000
         fphi(2,isite) = tuv100
         fphi(3,isite) = tuv010
         fphi(4,isite) = tuv001
         fphi(5,isite) = tuv200
         fphi(6,isite) = tuv020
         fphi(7,isite) = tuv002
         fphi(8,isite) = tuv110
         fphi(9,isite) = tuv101
         fphi(10,isite) = tuv011
         fphi(11,isite) = tuv300
         fphi(12,isite) = tuv030
         fphi(13,isite) = tuv003
         fphi(14,isite) = tuv210
         fphi(15,isite) = tuv201
         fphi(16,isite) = tuv120
         fphi(17,isite) = tuv021
         fphi(18,isite) = tuv102
         fphi(19,isite) = tuv012
         fphi(20,isite) = tuv111

         do jsite = 1, npole
           dfphidciX(1,isite,jsite) = xdtuv000(jsite)
           dfphidciX(2,isite,jsite) = xdtuv100(jsite)
           dfphidciX(3,isite,jsite) = xdtuv010(jsite)
           dfphidciX(4,isite,jsite) = xdtuv001(jsite)
           dfphidciX(5,isite,jsite) = xdtuv200(jsite)
           dfphidciX(6,isite,jsite) = xdtuv020(jsite)
           dfphidciX(7,isite,jsite) = xdtuv002(jsite)
           dfphidciX(8,isite,jsite) = xdtuv110(jsite)
           dfphidciX(9,isite,jsite) = xdtuv101(jsite)
           dfphidciX(10,isite,jsite) = xdtuv011(jsite)
           dfphidciX(11,isite,jsite) = xdtuv300(jsite)
           dfphidciX(12,isite,jsite) = xdtuv030(jsite)
           dfphidciX(13,isite,jsite) = xdtuv003(jsite)
           dfphidciX(14,isite,jsite) = xdtuv210(jsite)
           dfphidciX(15,isite,jsite) = xdtuv201(jsite)
           dfphidciX(16,isite,jsite) = xdtuv120(jsite)
           dfphidciX(17,isite,jsite) = xdtuv021(jsite)
           dfphidciX(18,isite,jsite) = xdtuv102(jsite)
           dfphidciX(19,isite,jsite) = xdtuv012(jsite)
           dfphidciX(20,isite,jsite) = xdtuv111(jsite)
         end do


      end do
c     end of isite
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

      deallocate (xdt0)
      deallocate (xdt1)
      deallocate (xdt2)
      deallocate (xdt3)
      deallocate (xdtu00)
      deallocate (xdtu10)
      deallocate (xdtu01)
      deallocate (xdtu20)
      deallocate (xdtu11)
      deallocate (xdtu02)
      deallocate (xdtu21)
      deallocate (xdtu12)
      deallocate (xdtu30)
      deallocate (xdtu03)
      deallocate (xdtuv000)
      deallocate (xdtuv100)
      deallocate (xdtuv010)
      deallocate (xdtuv001)
      deallocate (xdtuv200)
      deallocate (xdtuv020)
      deallocate (xdtuv002)
      deallocate (xdtuv110)
      deallocate (xdtuv101)
      deallocate (xdtuv011)
      deallocate (xdtuv300)
      deallocate (xdtuv030)
      deallocate (xdtuv003)
      deallocate (xdtuv210)
      deallocate (xdtuv201)
      deallocate (xdtuv120)
      deallocate (xdtuv021)
      deallocate (xdtuv102)
      deallocate (xdtuv012)
      deallocate (xdtuv111)

      write(*,*) "End of fphi_mpoleCT"
      return
      end
c     end of fphi_mpoleCT

c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphiCT  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_to_cphiCT" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c     routine for CT adds conversion of dfphidciX to dcphidciX
c
c

      subroutine fphi_to_cphiCT (fphi,cphi,dfphidciX,dcphidciX)
      use sizes
      use mpole
      implicit none
      integer i,j,k,ii
      real*8 ftc(10,10)
      real*8 cphi(10,*)
      real*8 fphi(20,*)
      real*8 dfphidciX(20,npole,*)
      real*8 dcphidciX(20,npole,*)
c
c     write(*,*) "Entered fphi_to_cphi"
c
c     find the matrix to convert fractional to Cartesian
c
c     write(*,*) "Calling frac_to_cart"
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      do i = 1, npole
         cphi(1,i) = ftc(1,1) * fphi(1,i)
         do ii = 1, npole
           dcphidciX(1,i,ii) = ftc(1,1)*dfphidciX(1,i,ii)
         end do

         do j = 2, 4
            cphi(j,i) = 0.0d0

            do k = 2, 4
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               do ii = 1, npole
                 dcphidciX(j,i,ii) = dcphidciX(j,i,ii)
     & + ftc(j,k)*dfphidciX(k,i,ii)
               end do
            end do
c           end of k

         end do
c        end of j

         do j = 5, 10
            cphi(j,i) = 0.0d0

            do k = 5, 10
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
               do ii = 1, npole
                 dcphidciX(j,i,ii) = dcphidciX(j,i,ii)
     & + ftc(j,k)*dfphidciX(k,i,ii)
               end do
            end do
c           end of k

         end do
c        end of j

      end do
c     end of i

c     write(*,*) "End of fphi_to_cphiCT"
      return
      end
c     end of fphi_to_cphiCT

c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fphi_uindCT  --  induced potential from grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fphi_uindCT" extracts the induced dipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uindCT (fdip_phi1,fdip_phi2,fdip_sum_phi,
     & dfdip_phi1,dfdip_phi2,dfdip_sum_phi,fdip_phi1a,fdip_phi2a,
     & fdip_sum_phia)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm,jsite
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3
      real*8 t0_1,t0_2,t1_1,t1_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu30,tu21,tu12,tu03
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111

      real*8 t0a,t1a,t2a,t3a
      real*8 t0_1a,t0_2a,t1_1a,t1_2a
      real*8 t2_1a,t2_2a,tq_1a,tq_2a
      real*8 tu00a,tu10a,tu01a,tu20a,tu11a
      real*8 tu02a,tu30a,tu21a,tu12a,tu03a
      real*8 tu00_1a,tu01_1a,tu10_1a
      real*8 tu00_2a,tu01_2a,tu10_2a
      real*8 tu20_1a,tu11_1a,tu02_1a
      real*8 tu20_2a,tu11_2a,tu02_2a
      real*8 tuv100_1a,tuv010_1a,tuv001_1a
      real*8 tuv100_2a,tuv010_2a,tuv001_2a
      real*8 tuv200_1a,tuv020_1a,tuv002_1a
      real*8 tuv110_1a,tuv101_1a,tuv011_1a
      real*8 tuv200_2a,tuv020_2a,tuv002_2a
      real*8 tuv110_2a,tuv101_2a,tuv011_2a
      real*8 tuv000a,tuv100a,tuv010a,tuv001a
      real*8 tuv200a,tuv020a,tuv002a,tuv110a
      real*8 tuv101a,tuv011a,tuv300a,tuv030a
      real*8 tuv003a,tuv210a,tuv201a,tuv120a
      real*8 tuv021a,tuv102a,tuv012a,tuv111a

      real*8 fdip_phi1(10,*)
      real*8 fdip_phi2(10,*)
      real*8 fdip_sum_phi(20,*)

      real*8 fdip_phi1a(10,*)
      real*8 fdip_phi2a(10,*)
      real*8 fdip_sum_phia(20,*)
      real*8 dfdip_phi1(10,npole,npole)
      real*8 dfdip_phi2(10,npole,npole)
      real*8 dfdip_sum_phi(20,npole,npole)
      
      real*8, allocatable :: dtq_1(:)
      real*8, allocatable :: dtq_2(:)
      real*8, allocatable :: dt0(:)
      real*8, allocatable :: dt1(:)
      real*8, allocatable :: dt2(:)
      real*8, allocatable :: dt3(:)
      real*8, allocatable :: dt0_1(:)
      real*8, allocatable :: dt1_1(:)
      real*8, allocatable :: dt2_1(:)
      real*8, allocatable :: dt3_1(:)
      real*8, allocatable :: dt0_2(:)
      real*8, allocatable :: dt1_2(:)
      real*8, allocatable :: dt2_2(:)
      real*8, allocatable :: dt3_2(:)
      real*8, allocatable :: dtu00(:)
      real*8, allocatable :: dtu10(:)
      real*8, allocatable :: dtu01(:)
      real*8, allocatable :: dtu20(:)
      real*8, allocatable :: dtu11(:)
      real*8, allocatable :: dtu02(:)
      real*8, allocatable :: dtu30(:)
      real*8, allocatable :: dtu21(:)
      real*8, allocatable :: dtu12(:)
      real*8, allocatable :: dtu03(:)
      real*8, allocatable :: dtu00_1(:)
      real*8, allocatable :: dtu10_1(:)
      real*8, allocatable :: dtu01_1(:)
      real*8, allocatable :: dtu00_2(:)
      real*8, allocatable :: dtu10_2(:)
      real*8, allocatable :: dtu01_2(:)
      real*8, allocatable :: dtu20_1(:)
      real*8, allocatable :: dtu11_1(:)
      real*8, allocatable :: dtu02_1(:)
      real*8, allocatable :: dtu20_2(:)
      real*8, allocatable :: dtu11_2(:)
      real*8, allocatable :: dtu02_2(:)
      real*8, allocatable :: dtuv000(:)
      real*8, allocatable :: dtuv100(:)
      real*8, allocatable :: dtuv010(:)
      real*8, allocatable :: dtuv001(:)
      real*8, allocatable :: dtuv200(:)
      real*8, allocatable :: dtuv020(:)
      real*8, allocatable :: dtuv002(:)
      real*8, allocatable :: dtuv110(:)
      real*8, allocatable :: dtuv101(:)
      real*8, allocatable :: dtuv011(:)
      real*8, allocatable :: dtuv300(:)
      real*8, allocatable :: dtuv030(:)
      real*8, allocatable :: dtuv003(:)
      real*8, allocatable :: dtuv210(:)
      real*8, allocatable :: dtuv201(:)
      real*8, allocatable :: dtuv120(:)
      real*8, allocatable :: dtuv021(:)
      real*8, allocatable :: dtuv102(:)
      real*8, allocatable :: dtuv012(:)
      real*8, allocatable :: dtuv111(:)
      real*8, allocatable :: dtuv100_1(:)
      real*8, allocatable :: dtuv010_1(:)
      real*8, allocatable :: dtuv001_1(:)
      real*8, allocatable :: dtuv100_2(:)
      real*8, allocatable :: dtuv010_2(:)
      real*8, allocatable :: dtuv001_2(:)
      real*8, allocatable :: dtuv200_1(:)
      real*8, allocatable :: dtuv020_1(:)
      real*8, allocatable :: dtuv002_1(:)
      real*8, allocatable :: dtuv200_2(:)
      real*8, allocatable :: dtuv020_2(:)
      real*8, allocatable :: dtuv002_2(:)
      real*8, allocatable :: dtuv110_1(:)
      real*8, allocatable :: dtuv101_1(:)
      real*8, allocatable :: dtuv011_1(:)
      real*8, allocatable :: dtuv110_2(:)
      real*8, allocatable :: dtuv101_2(:)
      real*8, allocatable :: dtuv011_2(:)
c
      allocate (dtq_1(npole))
      allocate (dtq_2(npole))
      allocate (dt0(npole))
      allocate (dt1(npole))
      allocate (dt2(npole))
      allocate (dt3(npole))
      allocate (dt0_1(npole))
      allocate (dt1_1(npole))
      allocate (dt2_1(npole))
      allocate (dt3_1(npole))
      allocate (dt0_2(npole))
      allocate (dt1_2(npole))
      allocate (dt2_2(npole))
      allocate (dt3_2(npole))
      allocate (dtu00(npole))
      allocate (dtu10(npole))
      allocate (dtu01(npole))
      allocate (dtu20(npole))
      allocate (dtu11(npole))
      allocate (dtu02(npole))
      allocate (dtu30(npole))
      allocate (dtu21(npole))
      allocate (dtu12(npole))
      allocate (dtu03(npole))
      allocate (dtu00_1(npole))
      allocate (dtu10_1(npole))
      allocate (dtu01_1(npole))
      allocate (dtu00_2(npole))
      allocate (dtu10_2(npole))
      allocate (dtu01_2(npole))
      allocate (dtu20_1(npole))
      allocate (dtu11_1(npole))
      allocate (dtu02_1(npole))
      allocate (dtu20_2(npole))
      allocate (dtu11_2(npole))
      allocate (dtu02_2(npole))
      allocate (dtuv000(npole))
      allocate (dtuv100(npole))
      allocate (dtuv010(npole))
      allocate (dtuv001(npole))
      allocate (dtuv200(npole))
      allocate (dtuv020(npole))
      allocate (dtuv002(npole))
      allocate (dtuv110(npole))
      allocate (dtuv101(npole))
      allocate (dtuv011(npole))
      allocate (dtuv300(npole))
      allocate (dtuv030(npole))
      allocate (dtuv003(npole))
      allocate (dtuv210(npole))
      allocate (dtuv201(npole))
      allocate (dtuv120(npole))
      allocate (dtuv021(npole))
      allocate (dtuv102(npole))
      allocate (dtuv012(npole))
      allocate (dtuv111(npole))
      allocate (dtuv100_1(npole))
      allocate (dtuv010_1(npole))
      allocate (dtuv001_1(npole))
      allocate (dtuv100_2(npole))
      allocate (dtuv010_2(npole))
      allocate (dtuv001_2(npole))
      allocate (dtuv200_1(npole))
      allocate (dtuv020_1(npole))
      allocate (dtuv002_1(npole))
      allocate (dtuv200_2(npole))
      allocate (dtuv020_2(npole))
      allocate (dtuv002_2(npole))
      allocate (dtuv110_1(npole))
      allocate (dtuv101_1(npole))
      allocate (dtuv011_1(npole))
      allocate (dtuv110_2(npole))
      allocate (dtuv101_2(npole))
      allocate (dtuv011_2(npole))
c
      write(*,*) "Entered fphi_uindCT"
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,
!$OMP& igrid,bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid,fdip_phi1,fdip_phi2,fdip_sum_phi,
!$OMP& dqgrdciX,dfdip_phi1,dfdip_phi2,dfdip_sum_phi)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)

         tuv100_1 = 0.0d0
         tuv010_1 = 0.0d0
         tuv001_1 = 0.0d0
         tuv200_1 = 0.0d0
         tuv020_1 = 0.0d0
         tuv002_1 = 0.0d0
         tuv110_1 = 0.0d0
         tuv101_1 = 0.0d0
         tuv011_1 = 0.0d0

         tuv100_2 = 0.0d0
         tuv010_2 = 0.0d0
         tuv001_2 = 0.0d0
         tuv200_2 = 0.0d0
         tuv020_2 = 0.0d0
         tuv002_2 = 0.0d0
         tuv110_2 = 0.0d0
         tuv101_2 = 0.0d0
         tuv011_2 = 0.0d0

         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0

         tuv100_1a = 0.0d0
         tuv010_1a = 0.0d0
         tuv001_1a = 0.0d0
         tuv200_1a = 0.0d0
         tuv020_1a = 0.0d0
         tuv002_1a = 0.0d0
         tuv110_1a = 0.0d0
         tuv101_1a = 0.0d0
         tuv011_1a = 0.0d0

         tuv100_2a = 0.0d0
         tuv010_2a = 0.0d0
         tuv001_2a = 0.0d0
         tuv200_2a = 0.0d0
         tuv020_2a = 0.0d0
         tuv002_2a = 0.0d0
         tuv110_2a = 0.0d0
         tuv101_2a = 0.0d0
         tuv011_2a = 0.0d0

         tuv000a = 0.0d0
         tuv001a = 0.0d0
         tuv010a = 0.0d0
         tuv100a = 0.0d0
         tuv200a = 0.0d0
         tuv020a = 0.0d0
         tuv002a = 0.0d0
         tuv110a = 0.0d0
         tuv101a = 0.0d0
         tuv011a = 0.0d0
         tuv300a = 0.0d0
         tuv030a = 0.0d0
         tuv003a = 0.0d0
         tuv210a = 0.0d0
         tuv201a = 0.0d0
         tuv120a = 0.0d0
         tuv021a = 0.0d0
         tuv102a = 0.0d0
         tuv012a = 0.0d0
         tuv111a = 0.0d0

         do jsite = 1, npole
            dtuv100_1(jsite) = 0.0d0
            dtuv010_1(jsite) = 0.0d0
            dtuv001_1(jsite) = 0.0d0
            dtuv200_1(jsite) = 0.0d0
            dtuv020_1(jsite) = 0.0d0
            dtuv002_1(jsite) = 0.0d0
            dtuv110_1(jsite) = 0.0d0
            dtuv101_1(jsite) = 0.0d0
            dtuv011_1(jsite) = 0.0d0

            dtuv100_2(jsite) = 0.0d0
            dtuv010_2(jsite) = 0.0d0
            dtuv001_2(jsite) = 0.0d0
            dtuv200_2(jsite) = 0.0d0
            dtuv020_2(jsite) = 0.0d0
            dtuv002_2(jsite) = 0.0d0
            dtuv110_2(jsite) = 0.0d0
            dtuv101_2(jsite) = 0.0d0
            dtuv011_2(jsite) = 0.0d0

            dtuv000(jsite) = 0.0d0
            dtuv001(jsite) = 0.0d0
            dtuv010(jsite) = 0.0d0
            dtuv100(jsite) = 0.0d0
            dtuv200(jsite) = 0.0d0
            dtuv020(jsite) = 0.0d0
            dtuv002(jsite) = 0.0d0
            dtuv110(jsite) = 0.0d0
            dtuv101(jsite) = 0.0d0
            dtuv011(jsite) = 0.0d0
            dtuv300(jsite) = 0.0d0
            dtuv030(jsite) = 0.0d0
            dtuv003(jsite) = 0.0d0
            dtuv210(jsite) = 0.0d0
            dtuv201(jsite) = 0.0d0
            dtuv120(jsite) = 0.0d0
            dtuv021(jsite) = 0.0d0
            dtuv102(jsite) = 0.0d0
            dtuv012(jsite) = 0.0d0
            dtuv111(jsite) = 0.0d0
         end do

         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)

            tu00_1 = 0.0d0
            tu01_1 = 0.0d0
            tu10_1 = 0.0d0
            tu20_1 = 0.0d0
            tu11_1 = 0.0d0
            tu02_1 = 0.0d0

            tu00_2 = 0.0d0
            tu01_2 = 0.0d0
            tu10_2 = 0.0d0
            tu20_2 = 0.0d0
            tu11_2 = 0.0d0
            tu02_2 = 0.0d0

            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0

            tu00_1a = 0.0d0
            tu01_1a = 0.0d0
            tu10_1a = 0.0d0
            tu20_1a = 0.0d0
            tu11_1a = 0.0d0
            tu02_1a = 0.0d0

            tu00_2a = 0.0d0
            tu01_2a = 0.0d0
            tu10_2a = 0.0d0
            tu20_2a = 0.0d0
            tu11_2a = 0.0d0
            tu02_2a = 0.0d0

            tu00a = 0.0d0
            tu10a = 0.0d0
            tu01a = 0.0d0
            tu20a = 0.0d0
            tu11a = 0.0d0
            tu02a = 0.0d0
            tu30a = 0.0d0
            tu21a = 0.0d0
            tu12a = 0.0d0
            tu03a = 0.0d0

            do jsite = 1, npole
               dtu00_1(jsite) = 0.0d0
               dtu01_1(jsite) = 0.0d0
               dtu10_1(jsite) = 0.0d0
               dtu20_1(jsite) = 0.0d0
               dtu11_1(jsite) = 0.0d0
               dtu02_1(jsite) = 0.0d0

               dtu00_2(jsite) = 0.0d0
               dtu01_2(jsite) = 0.0d0
               dtu10_2(jsite) = 0.0d0
               dtu20_2(jsite) = 0.0d0
               dtu11_2(jsite) = 0.0d0
               dtu02_2(jsite) = 0.0d0

               dtu00(jsite) = 0.0d0
               dtu10(jsite) = 0.0d0
               dtu01(jsite) = 0.0d0
               dtu20(jsite) = 0.0d0
               dtu11(jsite) = 0.0d0
               dtu02(jsite) = 0.0d0
               dtu30(jsite) = 0.0d0
               dtu21(jsite) = 0.0d0
               dtu12(jsite) = 0.0d0
               dtu03(jsite) = 0.0d0
            end do

            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)

               t0_1 = 0.0d0
               t1_1 = 0.0d0
               t2_1 = 0.0d0
               t0_2 = 0.0d0
               t1_2 = 0.0d0
               t2_2 = 0.0d0
               t3 = 0.0d0

               t0_1a = 0.0d0
               t1_1a = 0.0d0
               t2_1a = 0.0d0
               t0_2a = 0.0d0
               t1_2a = 0.0d0
               t2_2a = 0.0d0
               t3a = 0.0d0

               do jsite = 1, npole
                  dt0_1(jsite) = 0.0d0
                  dt1_1(jsite) = 0.0d0
                  dt2_1(jsite) = 0.0d0
                  dt0_2(jsite) = 0.0d0
                  dt1_2(jsite) = 0.0d0
                  dt2_2(jsite) = 0.0d0
                  dt3(jsite) = 0.0d0
               end do

               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq_1 = qgrid(1,i,j,k)
                  tq_2 = qgrid(2,i,j,k)
                  t0_1 = t0_1 + tq_1*thetai1(1,it1,iatm)
                  t1_1 = t1_1 + tq_1*thetai1(2,it1,iatm)
                  t2_1 = t2_1 + tq_1*thetai1(3,it1,iatm)
                  t0_2 = t0_2 + tq_2*thetai1(1,it1,iatm)
                  t1_2 = t1_2 + tq_2*thetai1(2,it1,iatm)
                  t2_2 = t2_2 + tq_2*thetai1(3,it1,iatm)
                  t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,iatm)

                  tq_1a = dqgrdci(1,i,j,k,isite)
                  tq_2a = dqgrdci(2,i,j,k,isite)
                  t0_1a = t0_1a + tq_1a*thetai1(1,it1,iatm)
                  t1_1a = t1_1a + tq_1a*thetai1(2,it1,iatm)
                  t2_1a = t2_1a + tq_1a*thetai1(3,it1,iatm)
                  t0_2a = t0_2a + tq_2a*thetai1(1,it1,iatm)
                  t1_2a = t1_2a + tq_2a*thetai1(2,it1,iatm)
                  t2_2a = t2_2a + tq_2a*thetai1(3,it1,iatm)
                  t3a = t3a + (tq_1a+tq_2a)*thetai1(4,it1,iatm)

                  do jsite = 1, npole
                     dtq_1(jsite) = dqgrdciX(1,i,j,k,isite,jsite)
                     dtq_2(jsite) = dqgrdciX(2,i,j,k,isite,jsite)
                     dt0_1(jsite) = dt0_1(jsite) 
     &                     + dtq_1(jsite)*thetai1(1,it1,iatm)
                     dt1_1(jsite) = dt1_1(jsite)
     &                     + dtq_1(jsite)*thetai1(2,it1,iatm)
                     dt2_1(jsite) = dt2_1(jsite)
     &                     + dtq_1(jsite)*thetai1(3,it1,iatm)
                     dt0_2(jsite) = dt0_2(jsite)
     &                     + dtq_2(jsite)*thetai1(1,it1,iatm)
                     dt1_2(jsite) = dt1_2(jsite)
     &                     + dtq_2(jsite)*thetai1(2,it1,iatm)
                     dt2_2(jsite) = dt2_2(jsite)
     &                     + dtq_2(jsite)*thetai1(3,it1,iatm)
                     dt3(jsite) = dt3(jsite) 
     & + ( dtq_1(jsite) + dtq_2(jsite) )*thetai1(4,it1,iatm)
                  end do
c                 end of jsite

c                 if ( i.eq.10 .and. j.eq.8 .and. k.eq.7 ) then
c                   write(*,*) i,j,k,isite
c                   write(*,*) qgrid(1,i,j,k),
c    & dqgrdciX(1,i,j,k,isite,4)
c                 endif

               end do
c              end of it1 

               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2

               t0 = t0_1 + t0_2
               t1 = t1_1 + t1_2
               t2 = t2_1 + t2_2

               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3

               tu00_1a = tu00_1a + t0_1a*u0
               tu10_1a = tu10_1a + t1_1a*u0
               tu01_1a = tu01_1a + t0_1a*u1
               tu20_1a = tu20_1a + t2_1a*u0
               tu11_1a = tu11_1a + t1_1a*u1
               tu02_1a = tu02_1a + t0_1a*u2
               tu00_2a = tu00_2a + t0_2a*u0
               tu10_2a = tu10_2a + t1_2a*u0
               tu01_2a = tu01_2a + t0_2a*u1
               tu20_2a = tu20_2a + t2_2a*u0
               tu11_2a = tu11_2a + t1_2a*u1
               tu02_2a = tu02_2a + t0_2a*u2

               t0a = t0_1a + t0_2a
               t1a = t1_1a + t1_2a
               t2a = t2_1a + t2_2a

               tu00a = tu00a + t0a*u0
               tu10a = tu10a + t1a*u0
               tu01a = tu01a + t0a*u1
               tu20a = tu20a + t2a*u0
               tu11a = tu11a + t1a*u1
               tu02a = tu02a + t0a*u2
               tu30a = tu30a + t3a*u0
               tu21a = tu21a + t2a*u1
               tu12a = tu12a + t1a*u2
               tu03a = tu03a + t0a*u3

               do jsite = 1, npole
                  dtu00_1(jsite) = dtu00_1(jsite) 
     &                   + dt0_1(jsite)*u0
                  dtu10_1(jsite) = dtu10_1(jsite)
     &                   + dt1_1(jsite)*u0
                  dtu01_1(jsite) = dtu01_1(jsite)
     &                   + dt0_1(jsite)*u1
                  dtu20_1(jsite) = dtu20_1(jsite)
     &                   + dt2_1(jsite)*u0
                  dtu11_1(jsite) = dtu11_1(jsite)
     &                   + dt1_1(jsite)*u1
                  dtu02_1(jsite) = dtu02_1(jsite)
     &                   + dt0_1(jsite)*u2
                  dtu00_2(jsite) = dtu00_2(jsite)
     &                   + dt0_2(jsite)*u0
                  dtu10_2(jsite) = dtu10_2(jsite)
     &                   + dt1_2(jsite)*u0
                  dtu01_2(jsite) = dtu01_2(jsite)
     &                   + dt0_2(jsite)*u1
                  dtu20_2(jsite) = dtu20_2(jsite)
     &                   + dt2_2(jsite)*u0
                  dtu11_2(jsite) = dtu11_2(jsite)
     &                   + dt1_2(jsite)*u1
                  dtu02_2(jsite) = dtu02_2(jsite)
     &                   + dt0_2(jsite)*u2

                  dt0(jsite) = dt0_1(jsite) 
     &                             + dt0_2(jsite)
                  dt1(jsite) = dt1_1(jsite) 
     &                             + dt1_2(jsite)
                  dt2(jsite) = dt2_1(jsite) 
     &                             + dt2_2(jsite)

                  dtu00(jsite) = dtu00(jsite) 
     &                  + dt0(jsite)*u0
                  dtu10(jsite) = dtu10(jsite) 
     &                  + dt1(jsite)*u0
                  dtu01(jsite) = dtu01(jsite) 
     &                  + dt0(jsite)*u1
                  dtu20(jsite) = dtu20(jsite) 
     &                  + dt2(jsite)*u0
                  dtu11(jsite) = dtu11(jsite) 
     &                  + dt1(jsite)*u1
                  dtu02(jsite) = dtu02(jsite) 
     &                  + dt0(jsite)*u2
                  dtu30(jsite) = dtu30(jsite) 
     &                  + dt3(jsite)*u0
                  dtu21(jsite) = dtu21(jsite) 
     &                  + dt2(jsite)*u1
                  dtu12(jsite) = dtu12(jsite) 
     &                  + dt1(jsite)*u2
                  dtu03(jsite) = dtu03(jsite) 
     &                  + dt0(jsite)*u3
               end do
c              end  of jsite

            end do
c           end of it2

            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1

            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1

            tuv100_1a = tuv100_1a + tu10_1a*v0
            tuv010_1a = tuv010_1a + tu01_1a*v0
            tuv001_1a = tuv001_1a + tu00_1a*v1
            tuv200_1a = tuv200_1a + tu20_1a*v0
            tuv020_1a = tuv020_1a + tu02_1a*v0
            tuv002_1a = tuv002_1a + tu00_1a*v2
            tuv110_1a = tuv110_1a + tu11_1a*v0
            tuv101_1a = tuv101_1a + tu10_1a*v1
            tuv011_1a = tuv011_1a + tu01_1a*v1
            tuv100_2a = tuv100_2a + tu10_2a*v0
            tuv010_2a = tuv010_2a + tu01_2a*v0
            tuv001_2a = tuv001_2a + tu00_2a*v1
            tuv200_2a = tuv200_2a + tu20_2a*v0
            tuv020_2a = tuv020_2a + tu02_2a*v0
            tuv002_2a = tuv002_2a + tu00_2a*v2
            tuv110_2a = tuv110_2a + tu11_2a*v0
            tuv101_2a = tuv101_2a + tu10_2a*v1
            tuv011_2a = tuv011_2a + tu01_2a*v1

            tuv000a = tuv000a + tu00a*v0
            tuv100a = tuv100a + tu10a*v0
            tuv010a = tuv010a + tu01a*v0
            tuv001a = tuv001a + tu00a*v1
            tuv200a = tuv200a + tu20a*v0
            tuv020a = tuv020a + tu02a*v0
            tuv002a = tuv002a + tu00a*v2
            tuv110a = tuv110a + tu11a*v0
            tuv101a = tuv101a + tu10a*v1
            tuv011a = tuv011a + tu01a*v1
            tuv300a = tuv300a + tu30a*v0
            tuv030a = tuv030a + tu03a*v0
            tuv003a = tuv003a + tu00a*v3
            tuv210a = tuv210a + tu21a*v0
            tuv201a = tuv201a + tu20a*v1
            tuv120a = tuv120a + tu12a*v0
            tuv021a = tuv021a + tu02a*v1
            tuv102a = tuv102a + tu10a*v2
            tuv012a = tuv012a + tu01a*v2
            tuv111a = tuv111a + tu11a*v1

            do jsite = 1, npole
               dtuv100_1(jsite) = dtuv100_1(jsite) 
     &                        + dtu10_1(jsite)*v0
               dtuv010_1(jsite) = dtuv010_1(jsite) 
     &                        + dtu01_1(jsite)*v0
               dtuv001_1(jsite) = dtuv001_1(jsite) 
     &                        + dtu00_1(jsite)*v1
               dtuv200_1(jsite) = dtuv200_1(jsite) 
     &                        + dtu20_1(jsite)*v0
               dtuv020_1(jsite) = dtuv020_1(jsite) 
     &                        + dtu02_1(jsite)*v0
               dtuv002_1(jsite) = dtuv002_1(jsite) 
     &                        + dtu00_1(jsite)*v2
               dtuv110_1(jsite) = dtuv110_1(jsite) 
     &                        + dtu11_1(jsite)*v0
               dtuv101_1(jsite) = dtuv101_1(jsite) 
     &                        + dtu10_1(jsite)*v1
               dtuv011_1(jsite) = dtuv011_1(jsite) 
     &                        + dtu01_1(jsite)*v1
               dtuv100_2(jsite) = dtuv100_2(jsite) 
     &                        + dtu10_2(jsite)*v0
               dtuv010_2(jsite) = dtuv010_2(jsite) 
     &                        + dtu01_2(jsite)*v0
               dtuv001_2(jsite) = dtuv001_2(jsite) 
     &                        + dtu00_2(jsite)*v1
               dtuv200_2(jsite) = dtuv200_2(jsite) 
     &                        + dtu20_2(jsite)*v0
               dtuv020_2(jsite) = dtuv020_2(jsite) 
     &                        + dtu02_2(jsite)*v0
               dtuv002_2(jsite) = dtuv002_2(jsite) 
     &                        + dtu00_2(jsite)*v2
               dtuv110_2(jsite) = dtuv110_2(jsite) 
     &                        + dtu11_2(jsite)*v0
               dtuv101_2(jsite) = dtuv101_2(jsite) 
     &                        + dtu10_2(jsite)*v1
               dtuv011_2(jsite) = dtuv011_2(jsite) 
     &                        + dtu01_2(jsite)*v1

               dtuv000(jsite) = dtuv000(jsite) 
     &                      + dtu00(jsite)*v0
               dtuv100(jsite) = dtuv100(jsite) 
     &                      + dtu10(jsite)*v0
               dtuv010(jsite) = dtuv010(jsite) 
     &                      + dtu01(jsite)*v0
               dtuv001(jsite) = dtuv001(jsite) 
     &                      + dtu00(jsite)*v1
               dtuv200(jsite) = dtuv200(jsite) 
     &                      + dtu20(jsite)*v0
               dtuv020(jsite) = dtuv020(jsite) 
     &                      + dtu02(jsite)*v0
               dtuv002(jsite) = dtuv002(jsite) 
     &                      + dtu00(jsite)*v2
               dtuv110(jsite) = dtuv110(jsite) 
     &                      + dtu11(jsite)*v0
               dtuv101(jsite) = dtuv101(jsite) 
     &                      + dtu10(jsite)*v1
               dtuv011(jsite) = dtuv011(jsite) 
     &                      + dtu01(jsite)*v1
               dtuv300(jsite) = dtuv300(jsite) 
     &                      + dtu30(jsite)*v0
               dtuv030(jsite) = dtuv030(jsite) 
     &                      + dtu03(jsite)*v0
               dtuv003(jsite) = dtuv003(jsite) 
     &                      + dtu00(jsite)*v3
               dtuv210(jsite) = dtuv210(jsite) 
     &                      + dtu21(jsite)*v0
               dtuv201(jsite) = dtuv201(jsite) 
     &                      + dtu20(jsite)*v1
               dtuv120(jsite) = dtuv120(jsite) 
     &                      + dtu12(jsite)*v0
               dtuv021(jsite) = dtuv021(jsite) 
     &                      + dtu02(jsite)*v1
               dtuv102(jsite) = dtuv102(jsite) 
     &                      + dtu10(jsite)*v2
               dtuv012(jsite) = dtuv012(jsite) 
     &                      + dtu01(jsite)*v2
               dtuv111(jsite) = dtuv111(jsite) 
     &                      + dtu11(jsite)*v1
            end do
c           end of jsite

         end do
c        end of it3

c        if (isite.eq.3) then
c           write(*,*) isite,tuv100a,dtuv100(4)
c        end if

         fdip_phi1(2,isite) = tuv100_1
         fdip_phi1(3,isite) = tuv010_1
         fdip_phi1(4,isite) = tuv001_1
         fdip_phi1(5,isite) = tuv200_1
         fdip_phi1(6,isite) = tuv020_1
         fdip_phi1(7,isite) = tuv002_1
         fdip_phi1(8,isite) = tuv110_1
         fdip_phi1(9,isite) = tuv101_1
         fdip_phi1(10,isite) = tuv011_1

         fdip_phi2(2,isite) = tuv100_2
         fdip_phi2(3,isite) = tuv010_2
         fdip_phi2(4,isite) = tuv001_2
         fdip_phi2(5,isite) = tuv200_2
         fdip_phi2(6,isite) = tuv020_2
         fdip_phi2(7,isite) = tuv002_2
         fdip_phi2(8,isite) = tuv110_2
         fdip_phi2(9,isite) = tuv101_2
         fdip_phi2(10,isite) = tuv011_2

         fdip_sum_phi(1,isite) = tuv000
         fdip_sum_phi(2,isite) = tuv100
         fdip_sum_phi(3,isite) = tuv010
         fdip_sum_phi(4,isite) = tuv001
         fdip_sum_phi(5,isite) = tuv200
         fdip_sum_phi(6,isite) = tuv020
         fdip_sum_phi(7,isite) = tuv002
         fdip_sum_phi(8,isite) = tuv110
         fdip_sum_phi(9,isite) = tuv101
         fdip_sum_phi(10,isite) = tuv011
         fdip_sum_phi(11,isite) = tuv300
         fdip_sum_phi(12,isite) = tuv030
         fdip_sum_phi(13,isite) = tuv003
         fdip_sum_phi(14,isite) = tuv210
         fdip_sum_phi(15,isite) = tuv201
         fdip_sum_phi(16,isite) = tuv120
         fdip_sum_phi(17,isite) = tuv021
         fdip_sum_phi(18,isite) = tuv102
         fdip_sum_phi(19,isite) = tuv012
         fdip_sum_phi(20,isite) = tuv111

         fdip_phi1a(2,isite) = tuv100_1a
         fdip_phi1a(3,isite) = tuv010_1a
         fdip_phi1a(4,isite) = tuv001_1a
         fdip_phi1a(5,isite) = tuv200_1a
         fdip_phi1a(6,isite) = tuv020_1a
         fdip_phi1a(7,isite) = tuv002_1a
         fdip_phi1a(8,isite) = tuv110_1a
         fdip_phi1a(9,isite) = tuv101_1a
         fdip_phi1a(10,isite) = tuv011_1a

         fdip_phi2a(2,isite) = tuv100_2a
         fdip_phi2a(3,isite) = tuv010_2a
         fdip_phi2a(4,isite) = tuv001_2a
         fdip_phi2a(5,isite) = tuv200_2a
         fdip_phi2a(6,isite) = tuv020_2a
         fdip_phi2a(7,isite) = tuv002_2a
         fdip_phi2a(8,isite) = tuv110_2a
         fdip_phi2a(9,isite) = tuv101_2a
         fdip_phi2a(10,isite) = tuv011_2a

         fdip_sum_phia(1,isite) = tuv000a
         fdip_sum_phia(2,isite) = tuv100a
         fdip_sum_phia(3,isite) = tuv010a
         fdip_sum_phia(4,isite) = tuv001a
         fdip_sum_phia(5,isite) = tuv200a
         fdip_sum_phia(6,isite) = tuv020a
         fdip_sum_phia(7,isite) = tuv002a
         fdip_sum_phia(8,isite) = tuv110a
         fdip_sum_phia(9,isite) = tuv101a
         fdip_sum_phia(10,isite) = tuv011a
         fdip_sum_phia(11,isite) = tuv300a
         fdip_sum_phia(12,isite) = tuv030a
         fdip_sum_phia(13,isite) = tuv003a
         fdip_sum_phia(14,isite) = tuv210a
         fdip_sum_phia(15,isite) = tuv201a
         fdip_sum_phia(16,isite) = tuv120a
         fdip_sum_phia(17,isite) = tuv021a
         fdip_sum_phia(18,isite) = tuv102a
         fdip_sum_phia(19,isite) = tuv012a
         fdip_sum_phia(20,isite) = tuv111a

         do jsite = 1, npole
            dfdip_phi1(2,isite,jsite) = dtuv100_1(jsite)
            dfdip_phi1(3,isite,jsite) = dtuv010_1(jsite)
            dfdip_phi1(4,isite,jsite) = dtuv001_1(jsite)
            dfdip_phi1(5,isite,jsite) = dtuv200_1(jsite)
            dfdip_phi1(6,isite,jsite) = dtuv020_1(jsite)
            dfdip_phi1(7,isite,jsite) = dtuv002_1(jsite)
            dfdip_phi1(8,isite,jsite) = dtuv110_1(jsite)
            dfdip_phi1(9,isite,jsite) = dtuv101_1(jsite)
            dfdip_phi1(10,isite,jsite) = dtuv011_1(jsite)

            dfdip_phi2(2,isite,jsite) = dtuv100_2(jsite)
            dfdip_phi2(3,isite,jsite) = dtuv010_2(jsite)
            dfdip_phi2(4,isite,jsite) = dtuv001_2(jsite)
            dfdip_phi2(5,isite,jsite) = dtuv200_2(jsite)
            dfdip_phi2(6,isite,jsite) = dtuv020_2(jsite)
            dfdip_phi2(7,isite,jsite) = dtuv002_2(jsite)
            dfdip_phi2(8,isite,jsite) = dtuv110_2(jsite)
            dfdip_phi2(9,isite,jsite) = dtuv101_2(jsite)
            dfdip_phi2(10,isite,jsite) = dtuv011_2(jsite)
         
            dfdip_sum_phi(1,isite,jsite) = dtuv000(jsite)
            dfdip_sum_phi(2,isite,jsite) = dtuv100(jsite)
            dfdip_sum_phi(3,isite,jsite) = dtuv010(jsite)
            dfdip_sum_phi(4,isite,jsite) = dtuv001(jsite)
            dfdip_sum_phi(5,isite,jsite) = dtuv200(jsite)
            dfdip_sum_phi(6,isite,jsite) = dtuv020(jsite)
            dfdip_sum_phi(7,isite,jsite) = dtuv002(jsite)
            dfdip_sum_phi(8,isite,jsite) = dtuv110(jsite)
            dfdip_sum_phi(9,isite,jsite) = dtuv101(jsite)
            dfdip_sum_phi(10,isite,jsite) = dtuv011(jsite)
            dfdip_sum_phi(11,isite,jsite) = dtuv300(jsite)
            dfdip_sum_phi(12,isite,jsite) = dtuv030(jsite)
            dfdip_sum_phi(13,isite,jsite) = dtuv003(jsite)
            dfdip_sum_phi(14,isite,jsite) = dtuv210(jsite)
            dfdip_sum_phi(15,isite,jsite) = dtuv201(jsite)
            dfdip_sum_phi(16,isite,jsite) = dtuv120(jsite)
            dfdip_sum_phi(17,isite,jsite) = dtuv021(jsite)
            dfdip_sum_phi(18,isite,jsite) = dtuv102(jsite)
            dfdip_sum_phi(19,isite,jsite) = dtuv012(jsite)
            dfdip_sum_phi(20,isite,jsite) = dtuv111(jsite)
         end do
c        end of jsite

      end do
c     end of isite
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL

      deallocate (dtq_1)
      deallocate (dtq_2)
      deallocate (dt0)
      deallocate (dt1)
      deallocate (dt2)
      deallocate (dt3)
      deallocate (dt0_1)
      deallocate (dt1_1)
      deallocate (dt2_1)
      deallocate (dt3_1)
      deallocate (dt0_2)
      deallocate (dt1_2)
      deallocate (dt2_2)
      deallocate (dt3_2)
      deallocate (dtu00)
      deallocate (dtu10)
      deallocate (dtu01)
      deallocate (dtu20)
      deallocate (dtu11)
      deallocate (dtu02)
      deallocate (dtu30)
      deallocate (dtu21)
      deallocate (dtu12)
      deallocate (dtu03)
      deallocate (dtu00_1)
      deallocate (dtu10_1)
      deallocate (dtu01_1)
      deallocate (dtu00_2)
      deallocate (dtu10_2)
      deallocate (dtu01_2)
      deallocate (dtu20_1)
      deallocate (dtu11_1)
      deallocate (dtu02_1)
      deallocate (dtu20_2)
      deallocate (dtu11_2)
      deallocate (dtu02_2)
      deallocate (dtuv000)
      deallocate (dtuv100)
      deallocate (dtuv010)
      deallocate (dtuv001)
      deallocate (dtuv200)
      deallocate (dtuv020)
      deallocate (dtuv002)
      deallocate (dtuv110)
      deallocate (dtuv101)
      deallocate (dtuv011)
      deallocate (dtuv300)
      deallocate (dtuv030)
      deallocate (dtuv003)
      deallocate (dtuv210)
      deallocate (dtuv201)
      deallocate (dtuv120)
      deallocate (dtuv021)
      deallocate (dtuv102)
      deallocate (dtuv012)
      deallocate (dtuv111)
      deallocate (dtuv100_1)
      deallocate (dtuv010_1)
      deallocate (dtuv001_1)
      deallocate (dtuv100_2)
      deallocate (dtuv010_2)
      deallocate (dtuv001_2)
      deallocate (dtuv200_1)
      deallocate (dtuv020_1)
      deallocate (dtuv002_1)
      deallocate (dtuv200_2)
      deallocate (dtuv020_2)
      deallocate (dtuv002_2)
      deallocate (dtuv110_1)
      deallocate (dtuv101_1)
      deallocate (dtuv011_1)
      deallocate (dtuv110_2)
      deallocate (dtuv101_2)
      deallocate (dtuv011_2)

c     write(*,*) "fdip_phi1"
c     write(*,*) fdip_phi1(2,3),dfdip_phi1(2,3,4)
c     write(*,*) "fdip_phi2"
c     write(*,*) fdip_phi2(2,3),dfdip_phi2(2,3,4)
c     write(*,*) "fdip_sum_phi"
c     write(*,*) fdip_sum_phi(2,3),dfdip_sum_phi(2,3,4)

      write(*,*) "Done with fphi_uindCT"

      return
      end
c     end of NOfphi_uindCT

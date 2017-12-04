c
      program first_order
      IMPLICIT NONE
c
      real mscat,windis,xfp,xptar,yfp,ytar,yptar,delta
      real cham_pos_res,mscat1,mscat2,mirdis,mscat3

c
c SHMS 1st order optics
c xptar (mr) = 0.26 xfp ( in mm) - 1.38 xpfp ( in mr)
c delta (%)  = 0.06 xfp  (in mm)  - 0.0012 xpfp ( in mr)
c ytar  (mm) = -.61 yfp (in mm) - 0.04ypfp (in mr)
c yptar  (mr)= 0.27 yfp (in mm) - 1.6 ypfp (in mr)
c
      write(*,'(8(a8,2x))') " type","exit dis","xpfp","xfp"
     > ,"xptar","yptar","ytar","delta"
      write(*,'(8(a8,2x))') "     ","   mm     "," mr "," mm"
     > ,"mr","mr"," mm "," % "
      cham_pos_res = 0.200
      windis=0. ! distance in mm
      mscat=0.0002 ! just look at wire chamber resolution 
      xfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
      yfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
c calculate RMS for target quantities
      xptar = sqrt( (0.26*xfp)**2 + (1.38*mscat*1000)**2)
      delta = sqrt( (0.06*xfp)**2 + (0.0012*mscat*1000)**2)
      yptar = sqrt( (0.27*yfp)**2 + (1.6*mscat*1000)**2)
      ytar =  sqrt( (0.61*yfp)**2 + (0.04*mscat*1000)**2)
      write(*,'(a8,7f10.5)') ' wc  '
     >,windis,mscat,xfp,xptar,yptar,ytar,delta
      windis=61.*10. ! distance in mm
      call multscat('e',2000.,8.89,0.0508,mscat)
      xfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
      yfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
c calculate RMS for target quantities
      xptar = sqrt( (0.26*xfp)**2 + (1.38*mscat*1000)**2)
      delta = sqrt( (0.06*xfp)**2 + (0.0012*mscat*1000)**2)
      yptar = sqrt( (0.27*yfp)**2 + (1.6*mscat*1000)**2)
      ytar =  sqrt( (0.61*yfp)**2 + (0.04*mscat*1000)**2)
      write(*,'(a8,7f10.5)') 
     >'vac   ',windis,mscat,xfp,xptar,yptar,ytar,delta
c
      windis=291.*10. ! distance in mm
      call multscat('e',2000.,8.89,0.0508,mscat)
      xfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
      yfp=sqrt((windis*mscat)**2+cham_pos_res*cham_pos_res)
c calculate RMS for target quantities
      xptar = sqrt( (0.26*xfp)**2 + (1.38*mscat*1000)**2)
      delta = sqrt( (0.06*xfp)**2 + (0.0012*mscat*1000)**2)
      yptar = sqrt( (0.27*yfp)**2 + (1.6*mscat*1000)**2)
      ytar =  sqrt( (0.61*yfp)**2 + (0.04*mscat*1000)**2)
      write(*,'(a8,7f10.5)') 'vac   '
     >,windis,mscat,xfp,xptar,yptar,ytar,delta
c Cerenkov
      windis=291.*10. ! distance in mm
      call multscat('e',2000.,8.89,0.0508,mscat1)
      mirdis=61.*10. ! distance in mm
      call multscat('e',2000.,12.29,0.3,mscat2)
      xfp=sqrt((windis*mscat1)**2+(mirdis*mscat2)**2
     >+cham_pos_res*cham_pos_res)
      yfp=sqrt((windis*mscat1)**2+(mirdis*mscat2)**2
     >+cham_pos_res*cham_pos_res)
c calculate RMS for target quantities
      mscat = sqrt(mscat1**2+mscat2**2)
      xptar = sqrt( (0.26*xfp)**2 + (1.38*mscat*1000)**2)
      delta = sqrt( (0.06*xfp)**2 + (0.0012*mscat*1000)**2)
      yptar = sqrt( (0.27*yfp)**2 + (1.6*mscat*1000)**2)
      ytar =  sqrt( (0.61*yfp)**2 + (0.04*mscat*1000)**2)
      write(*,'(a8,7f10.5)') 'cer  '
     >,windis,mscat,xfp,xptar,yptar,ytar,delta
c
      end
c
      subroutine MULTSCAT(ch,pmom,radlen,len,mscat)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Program to calculate the multiple scattering
c       angle. Formula from V. L. Highland, NIM 129,p497 (1975). 5% accurate
c       for elements with z>20 and 10-3<L/Lr<10, 10-20% accurate for z<20 
c       and/or low velocity. Set to compute the planar angle.
C-
C-   Inputs  : 
C-   Outputs : 
C-   Controls: 
C-
C-   Created  17-MAR-1993   M. K. Jones
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      real length,density,pmom
      real len,den,pm
      real radlen,mp,ep,beta,zinc,lenrad,temp,mscat
      character*1 ch
c
c
      if (ch .eq. 'p') then
      mp = 938.279
      else
      mp=.511
      endif
      ep = sqrt(pmom*pmom + mp*mp)
      beta = pmom/ep
      zinc = 1
       temp = 1 + .038*log(len/radlen)
      mscat = 13.6/pmom/beta*zinc*temp*sqrt(len/radlen)
      return
      end
C----------------------------------------------------------------------
c

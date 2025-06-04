! *******************************************************************
! *   INTERPOLATION PARAMETER'S LIBRARY READING                     *
! *******************************************************************
      SUBROUTINE YMNCFALI(fullfilename)
 
      implicit      none
      include       'YMNC_DAT.fh'
      integer*4     libuni0       !CURRENT UNIT
      character*(*) fullfilename  !FILE NAME
      character*5   version       !VERSION IDENTIFIER
      integer*4     narg          !TOTAL NUMBER OF ARGUMENTS
      integer*4     nnf           !TOTAL NUMBER OF INTERPOLATION PARAM.
      integer*4     nsort         !TOTAL NUMBER OF SORTS
      integer*4     ksdimtotal    !COEFF.ARRAY DIM. FOR ONE SORT
      integer*4     nextmult      !INTERP.COEF.NUMBER OF EXT.MULT.FACTOR
      INTEGER*4     multint       !INTERNAL MULT.FACT.EXISTS FLAG VALUE
 
!____ ERROR FLAG FOR OBJECT YMNC
      YMNCerr = 0
 
!____ OPENING LIBRARY BINARY FILE FOR READING
      libuni0 = 1
      open(libuni0,file=fullfilename, access='stream', &
        form='unformatted',status='old')
 
!____ READING LIBRARY VERSION IDENTIFIER:
      read(libuni0) version
      if(version.ne.'NF201') YMNCerr = 1
 
!____ READING TOTAL NUMBER OF ARGUMENTS:
      read(libuni0) narg
      if(narg.ne.YMNCnarg) YMNCerr = 2
 
!____ READING TOTAL NUMBER OF INTERPOLATION PARAMETERS:
      read(libuni0) nnf
      if(nnf.ne.YMNCnfnum) YMNCerr = 3
 
!____ READING TOTAL NUMBER OF PARAM'S SORTS:
      read(libuni0) nsort
      if(nsort.ne.YMNCnsort) YMNCerr = 4
 
!____ READING COEFFICIENTS ARRAY TOTAL DIMENSION FOR ONE SORT:
      read(libuni0) ksdimtotal
      if(ksdimtotal.ne.YMNCksdimtotal) YMNCerr = 5
 
!____ READING INTERP. COEF. NUMBER OF EXT. MULT. FACTOR:
      read(libuni0) nextmult
 
!____ READING INTERP. POINTS NUMBER FOR EACH ARGUMENT:
      read(libuni0) YMNCnnodd
 
!____ READING PARAMETER'S ARGUMENTS NUMBER SUBJECT TO REDUCING:
      read(libuni0) YMNCredarg
 
!____ READING INPUT REDUCING VALUES(0-NO REDUCING;!=0-REDUCING PRESENCE):
      read(libuni0) YMNCreduce
 
!____ READING (ARGUMENT POINTS NUMBER-1) OR (1) IF ARGUMENT IS REDUCED:
      read(libuni0) YMNCnnodd1
 
!____ READING INTERP. PARAMETER'S OFFSETS ARRAY IN THE COEFF.ARRAY:
      read(libuni0) YMNCksoff
 
!____ READING COEFFICIENTS DIMEENSION ARRAY FOR EACH INTERP. PARAM.
      read(libuni0) YMNCksdimarr
 
!____ READING KOEFF. "2**ARG.NUM." FOR EACH INTERP. PARAM.
      read(libuni0) YMNCksnargarr
 
!____ READING SORT OFFSETS ARRAY
      read(libuni0) YMNCsortoff
 
!____ READING SORT INDEXES ARRAY
      read(libuni0) YMNCsortind
 
!____ READING 1./GRID STEP FOR EACH ARGUMENTS:
      read(libuni0) YMNCostep
 
!____ READING MINIMUM ARGUMENTS VALUES:
      read(libuni0) YMNCfval
 
!____ READING MAXIMUM ARGUMENTS VALUES:
      read(libuni0) YMNClval
 
!____ READING POINTS VALUES
      read(libuni0) YMNCuzlval
 
!____ READING NUMBER OF DERIVATIVES FOR EACH PARAMETER
      read(libuni0) YMNCdfnum
 
!____ READING DERIVATIVE'S POSITION FOR EACH PARAMETER
      read(libuni0) YMNCdflag
 
!____ READING EXTERNAL MULT. FACT. EXISTS FLAG VALUE:
      read(libuni0) YMNCemfl
 
!____ READING INTERNAL MULT. FACT. EXISTS FLAG VALUE:
      read(libuni0) multint
 
!____ READING STEP FLAG
      read(libuni0) YMNCstepfl
 
!____ READING TOTAL COEFFICIENTS ARRAY FOR ALL SORTS
      read(libuni0) YMNCksarr
 
!____ READING EXTERNAL MULT.FACTORS IF THEY EXISTS:
      if(YMNCemfl.ne.0) then
        read(libuni0) YMNCextmultarr
      end if
 
!____ READING INTERNAL MULT.FACTORS IF THEY EXISTS:
      if(multint.ne.0) then
        read(libuni0) YMNCintmultarr
      end if
 
!____ CLOSING LIBRARY BINARY FILE
      close(libuni0)
 
      RETURN
      END
! *******************************************************************
! *   PARAMETER'S INTERPOLATION FOR CURRENT PARAMETER'S SORT        *
! *******************************************************************
      SUBROUTINE YMNCCALC(X,NSORT,RESUL,DRESUL,KSDFA)
 
      implicit    none
      include     'YMNC_DAT.fh'
      integer*4   nsort                !CURRENT PARAM.'S SORT
      real*4      resul    (YMNCnfnum) !INTERP.FUNCT.RETURNING VALUE
      real*4      dresul   (YMNCnfnum*YMNCnarg)
      real*4      x        (YMNCnarg)  !CURRENT ARGUMENT VALUES
      integer*4   i,j,k,jb,ju,js,ii
      integer*4   ksdfa    (YMNCnfnum) !OFFSET FOR DERIVITIES ARRAY
      integer*4   nfint    (YMNCnarg)  !ARG.INTERV.INDEX FOR CURR. PARAM.
      real*4      nfxof    (YMNCnarg)  !ARG.OFFSETS FOR CURR. PARAM.
      integer*4   npnt_1   (YMNCnarg)  !INTERP.PNT.NUM.-1 FOR CURR. PARAM.
      real*4      nostep   (YMNCnarg)  !1./STEP FOR CURR.PARAM.
      real*4      ostepv   (YMNCnarg)  !1./STEP FOR CURR.PARAM.IF FLSTEP=1
      real*4      vv
      real*4      p        (10)
      real*4      y        (YMNCnarg)
 
!____ CURRENT INTERVAL ARRAY CALCULATION
      if(YMNCstepfl.eq.0) then
        do i=YMNCnarg,1,-1
          y(i) = x(i)
          if(y(i).gt.YMNClval(i)) y(i) = YMNClval(i)
          if(y(i).lt.YMNCfval(i)) y(i) = YMNCfval(i)
          YMNCinterv(i) = (y(i) - YMNCfval(i))*YMNCostep(i)
          if(YMNCinterv(i).lt.0) YMNCinterv(i)=0
          if((YMNCinterv(i).ge.YMNCnnodd(i)-1).and.&
            (YMNCnnodd(i).gt.1)) YMNCinterv(i)=YMNCnnodd(i)-2
          YMNCxof(i)=(x(i) - YMNCfval(i))*YMNCostep(i)-YMNCinterv(i)
        end do
      else
        ju=1
        do i=1,YMNCnarg
          if(YMNCnnodd(i).gt.1) then
            vv=x(i)
            js=ju
            jb=ju+YMNCnnodd(i)-1
            do while(js.lt.jb-1)
              j=(js+jb)/2
              if(YMNCuzlval(j).le.vv) then
                js=j
              else
                jb=j
              end if
            end do
            YMNCinterv(i)=js-ju
            ostepv(i) = 1./(YMNCuzlval(js+1)-YMNCuzlval(js))
            YMNCxof(i)=(vv-YMNCuzlval(js))*ostepv(i)
          else
            YMNCinterv(i)= 0
            YMNCxof(i)   = 0.
            ostepv(i)    = 0.
          end if
          ju=ju+YMNCnnodd(i)
        end do
      end if
 
      dresul = 0.
 
!____ PROCCESSING FOR ALL INTERPOLATION PARAMETERS:
      do i = 1,YMNCnfnum
 
        ksdfa(i) = (i-1)*YMNCnarg
 
!____ FORMING INTERVAL ARR.(NFINT),ARG. OFFSETS ARR.(NFXOF)
!____ AND INTERP.PNT.NUM.-1 ARR.(NPNT_1) FOR CURRENT INTERP. PARAMETER
        k=0
        do j = 1,YMNCnarg
          if(YMNCreduce((i-1)*YMNCnarg+j).eq.0) then
            k = k + 1
            nfint(k)  = YMNCinterv(j)
            nfxof(k)  = YMNCxof(j)
            npnt_1(k) = YMNCnnodd1((i-1)*YMNCnarg+j)
          end if
          if(YMNCstepfl.eq.0) then
            nostep(j) = YMNCostep(j)
          else
            nostep(j) = ostepv(j)
          end if
        end do
 
        if(YMNCredarg(i).lt.1) then
!____ CALCULATION RESULT FOR ARGUMENTS NUMBER == 0
          k = YMNCsortoff(YMNCsortind(nsort)) + YMNCksoff(i) + 1
          resul(i) = YMNCksarr(k)
        else
!____ DEFINE INTERPOLATION FUNCTION AND RECEIVE RESULT
!____ FINDING OFFSET IN THE K_ARR SUBJECT TO INTERVAL ARRAY
          k  = nfint(1)
          if(YMNCredarg(i).gt.1) then
            do ii = 2,YMNCredarg(i)
              k = k * npnt_1(ii) + nfint(ii)
            end do
          end if
          k  = k * YMNCksnargarr(i)+&
               YMNCsortoff(YMNCsortind(nsort)) + YMNCksoff(i)
 
!____ CALLING INTERPOLATION FUNCTION
          select case (YMNCredarg(i))
            case (1)
              resul(i) = YMNCksarr(k+1) +&
                        YMNCksarr(k+2)*nfxof(1)
              if(btest(YMNCdflag(i),1))&
                dresul(ksdfa(i)+1) = YMNCksarr(k+2)
            case (2)
              dresul(ksdfa(i)+1) = YMNCksarr(k+3) +&
                                  YMNCksarr(k+4) * nfxof(2)
              resul(i) = YMNCksarr(k+1) +&
                        YMNCksarr(k+2) * nfxof(2) +&
                        dresul(ksdfa(i)+1) * nfxof(1)
              if(btest(YMNCdflag(i),2))&
               dresul(ksdfa(i)+2) = YMNCksarr(k+2) +&
                                    YMNCksarr(k+4) * nfxof(1)
            case (3)
              p(1) = YMNCksarr(k+3)+YMNCksarr(k+4)*nfxof(3)
              p(2) = YMNCksarr(k+7)+YMNCksarr(k+8)*nfxof(3)
              dresul(ksdfa(i)+1) = YMNCksarr(k+5) +&
                                  YMNCksarr(k+6)*nfxof(3) +&
                                  p(2)*nfxof(2)
              resul(i) = YMNCksarr(k+1) +&
                        YMNCksarr(k+2)*nfxof(3) +&
                        p (1)*nfxof(2) +&
                        dresul(ksdfa(i)+1)*nfxof(1)
              if(btest(YMNCdflag(i),2))&
               dresul(ksdfa(i)+2) = p(1) + p(2)*nfxof(1)
              if(btest(YMNCdflag(i),3))&
               dresul(ksdfa(i)+3) = YMNCksarr(k+2) +&
                                    YMNCksarr(k+4)*nfxof(2) +&
                                   (YMNCksarr(k+6) +&
                                    YMNCksarr(k+8)*nfxof(2))*&
                                    nfxof(1)
            case (4)
              p(5) = YMNCksarr(k+15)+YMNCksarr(k+16)*nfxof(4)
              p(1) = YMNCksarr(k+13)+YMNCksarr(k+14)*nfxof(4)+&
                    p ( 5)*nfxof( 3)
              p(6) = YMNCksarr(k+7)+YMNCksarr(k+8)*nfxof(4)
              p(2) = YMNCksarr(k+5)+YMNCksarr(k+6)*nfxof(4) +&
                    p ( 6)*nfxof( 3)
              p(3) = YMNCksarr(k+3)+YMNCksarr(k+4)*nfxof(4)
              p(4) = YMNCksarr(k+11)+YMNCksarr(k+12)*nfxof(4)
              dresul(ksdfa(i)+1) = YMNCksarr(k+9)+&
                                  YMNCksarr(k+10)*nfxof(4)+&
                                  p ( 4)*nfxof( 3)+&
                                  p ( 1)*nfxof( 2)
              resul(i) = YMNCksarr(k+1)+YMNCksarr(k+2)*nfxof(4)+&
                        p ( 3)*nfxof( 3) +&
                        p ( 2)*nfxof( 2) +&
                        dresul(ksdfa(i)+1)*nfxof(1)
              if(btest(YMNCdflag(i),2))&
               dresul(ksdfa(i)+2) = p(2) + p(1)*nfxof(1)
              if(btest(YMNCdflag(i),3))&
               dresul(ksdfa(i)+3) = p(3) +&
                                    p(6)*nfxof(2) +&
                                   (p(4)+p(5)*nfxof(2))*nfxof(1)
              if(btest(YMNCdflag(i),4))&
               dresul(ksdfa(i)+4) = YMNCksarr(k+2)+&
                                    YMNCksarr(k+4)*nfxof(3) +&
                                   (YMNCksarr(k+6)+&
                                    YMNCksarr(k+8)*nfxof(3))*&
                                    nfxof(2) +&
                                   (YMNCksarr(k+10)+&
                                    YMNCksarr(k+12)*nfxof(3) +&
                                   (YMNCksarr(k+14)+&
                                    YMNCksarr(k+16)*nfxof(3))*&
                                    nfxof(2))*nfxof(1)
            case (5)
              p(4) = YMNCksarr(k+13)+YMNCksarr(k+14)*nfxof(5)+&
                   (YMNCksarr(k+15)+YMNCksarr(k+16)*nfxof(5))*&
                    nfxof(4)
              p(1) = YMNCksarr(k+9)+YMNCksarr(k+10)*nfxof(5)+&
                   (YMNCksarr(k+11)+YMNCksarr(k+12)*nfxof(5))*&
                    nfxof(4) + p(4)*nfxof(3)
              p(2) = YMNCksarr(k+5)+YMNCksarr(k+6)*nfxof(5)+&
                   (YMNCksarr(k+7)+YMNCksarr(k+8)*nfxof(5))*&
                    nfxof(4)
              p(6) = YMNCksarr(k+29)+YMNCksarr(k+30)*nfxof(5)+&
                   (YMNCksarr(k+31)+YMNCksarr(k+32)*nfxof(5))*&
                    nfxof(4)
              p(3) = YMNCksarr(k+25)+YMNCksarr(k+26)*nfxof(5)+&
                   (YMNCksarr(k+27)+YMNCksarr(k+28)*nfxof(5))*&
                    nfxof(4)+p(6)*nfxof(3)
              p(5) = YMNCksarr(k+21)+YMNCksarr(k+22)*nfxof(5)+&
                   (YMNCksarr(k+23)+YMNCksarr(k+24)*nfxof(5))*&
                    nfxof(4)
              p(7) = YMNCksarr(k+3)+YMNCksarr(k+4)*nfxof(5)
              dresul(ksdfa(i)+1) = YMNCksarr(k+17)+&
                                  YMNCksarr(k+18)*nfxof(5)+&
                                 (YMNCksarr(k+19)+&
                                  YMNCksarr(k+20)*nfxof(5))*&
                                  nfxof(4)+p(5)*nfxof(3)+&
                                  p(3)*nfxof(2)
              resul(i) = YMNCksarr(k+1)+YMNCksarr(k+2)*nfxof(5)+&
                        p(7)*nfxof(4)+&
                        p(2)*nfxof(3)+&
                        p(1)*nfxof(2)+&
                        dresul(ksdfa(i)+1)*nfxof(1)
              if(btest(YMNCdflag(i),2))&
               dresul(ksdfa(i)+2) = p(1)+p(3)*nfxof(1)
              if(btest(YMNCdflag(i),3))&
               dresul(ksdfa(i)+3) = p(2)+p(4)*nfxof(2) +&
                                   (p(5)+p(6)*nfxof(2))*nfxof(1)
              if(btest(YMNCdflag(i),4)) dresul(ksdfa(i)+4) = p(7) +&
               (YMNCksarr(k+7)+YMNCksarr(k+8)*nfxof(5))*&
                nfxof(3)+&
               (YMNCksarr(k+11)+YMNCksarr(k+12)*nfxof(5) +&
               (YMNCksarr(k+15)+YMNCksarr(k+16)*nfxof(5))*&
                nfxof(3))*nfxof(2)+&
               (YMNCksarr(k+19)+YMNCksarr(k+20)*nfxof(5) +&
               (YMNCksarr(k+23)+YMNCksarr(k+24)*nfxof(5))*&
                nfxof(3) +&
               (YMNCksarr(k+27)+YMNCksarr(k+28)*nfxof(5) +&
               (YMNCksarr(k+31)+YMNCksarr(k+32)*nfxof(5))*&
                nfxof(3))*nfxof(2))*nfxof(1)
              if(btest(YMNCdflag(i),5)) dresul(ksdfa(i)+5) =&
                YMNCksarr(k+2)+YMNCksarr(k+4)*nfxof(4)+&
               (YMNCksarr(k+6)+YMNCksarr(k+8)*nfxof(4))*&
                nfxof(3)+&
               (YMNCksarr(k+10)+YMNCksarr(k+12)*nfxof(4) +&
               (YMNCksarr(k+14)+YMNCksarr(k+16)*nfxof(4))*&
                nfxof(3))*nfxof(2)+&
               (YMNCksarr(k+18)+YMNCksarr(k+20)*nfxof(4)+&
               (YMNCksarr(k+22)+YMNCksarr(k+24)*nfxof(4))*&
                nfxof(3)+&
               (YMNCksarr(k+26)+YMNCksarr(k+28)*nfxof(4) +&
               (YMNCksarr(k+30)+YMNCksarr(k+31)*nfxof(4))*&
                nfxof(3))*nfxof(2))*nfxof(1)
            case (6)
              dresul(ksdfa(i)+1) = YMNCksarr(k+33)+&
               YMNCksarr(k+34)*nfxof(6)+(YMNCksarr(k+35)+&
               YMNCksarr(k+36)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+37)+YMNCksarr(k+38)*nfxof(6)+&
              (YMNCksarr(k+39)+YMNCksarr(k+40)*nfxof(6))*&
               nfxof(5))*nfxof(4)&
             +(YMNCksarr(k+41)+YMNCksarr(k+42)*nfxof(6)+&
              (YMNCksarr(k+43)+YMNCksarr(k+44)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+45)+YMNCksarr(k+46)*nfxof(6)&
             +(YMNCksarr(k+47)+YMNCksarr(k+48)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3)&
             +(YMNCksarr(k+49)+YMNCksarr(k+50)*nfxof(6)+&
              (YMNCksarr(k+51)+YMNCksarr(k+52)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+53)+YMNCksarr(k+54)*nfxof(6)+&
              (YMNCksarr(k+55)+YMNCksarr(k+56)*nfxof(6))*&
               nfxof(5))*nfxof(4)&
             +(YMNCksarr(k+57)+YMNCksarr(k+58)*nfxof(6)+&
              (YMNCksarr(k+59)+YMNCksarr(k+60)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+61)+YMNCksarr(k+62)*nfxof(6)&
             +(YMNCksarr(k+63)+YMNCksarr(k+64)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3))*nfxof(2)
              resul(i) = YMNCksarr(k+1)+YMNCksarr(k+2)*nfxof(6)+&
              (YMNCksarr(k +3)+YMNCksarr(k+ 4)*nfxof(6))*&
               nfxof(5)+&
              (YMNCksarr(k +5)+YMNCksarr(k+ 6)*nfxof(6)+&
              (YMNCksarr(k+7)+YMNCksarr(k+8)*nfxof(6))*&
               nfxof(5))*nfxof(4) +&
              (YMNCksarr(k +9)+YMNCksarr(k+10)*nfxof(6)+&
              (YMNCksarr(k+11)+YMNCksarr(k+12)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+13)+YMNCksarr(k+14)*nfxof(6)&
             +(YMNCksarr(k+15)+YMNCksarr(k+16)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3)+&
              (YMNCksarr(k+17)+YMNCksarr(k+18)*nfxof(6)+&
              (YMNCksarr(k+19)+YMNCksarr(k+20)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+21)+YMNCksarr(k+22)*nfxof(6)+&
              (YMNCksarr(k+23)+YMNCksarr(k+24)*nfxof(6))*&
               nfxof(5))*nfxof(4)&
             +(YMNCksarr(k+25)+YMNCksarr(k+26)*nfxof(6)+&
              (YMNCksarr(k+27)+YMNCksarr(k+28)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+29)+YMNCksarr(k+30)*nfxof(6)&
             +(YMNCksarr(k+31)+YMNCksarr(k+32)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3))*nfxof(2)+&
               dresul(ksdfa(i)+1)*nfxof(1)
              if(btest(YMNCdflag(i),2)) dresul(ksdfa(i)+2) =&
               YMNCksarr(k+17)+YMNCksarr(k+18)*nfxof(6)+&
              (YMNCksarr(k+19)+YMNCksarr(k+20)*nfxof(6))*&
               nfxof(5)+&
              (YMNCksarr(k+21)+YMNCksarr(k+22)*nfxof(6)+&
              (YMNCksarr(k+23)+YMNCksarr(k+24)*nfxof(6))*&
               nfxof(5))*nfxof(4)+&
              (YMNCksarr(k+25)+YMNCksarr(k+26)*nfxof(6)+&
              (YMNCksarr(k+27)+YMNCksarr(k+28)*nfxof(6))*&
               nfxof(5) +&
              (YMNCksarr(k+29)+YMNCksarr(k+30)*nfxof(6)+&
              (YMNCksarr(k+31)+YMNCksarr(k+32)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3)+&
              (YMNCksarr(k+49)+YMNCksarr(k+50)*nfxof(6)+&
              (YMNCksarr(k+51)+YMNCksarr(k+52)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+53)+YMNCksarr(k+54)*nfxof(6)+&
              (YMNCksarr(k+55)+YMNCksarr(k+56)*nfxof(6))*&
               nfxof(5))*nfxof(4)&
             +(YMNCksarr(k+57)+YMNCksarr(k+58)*nfxof(6)+&
              (YMNCksarr(k+59)+YMNCksarr(k+60)*nfxof(6))*nfxof(5)&
             +(YMNCksarr(k+61)+YMNCksarr(k+62)*nfxof(6)&
             +(YMNCksarr(k+63)+YMNCksarr(k+64)*nfxof(6))*&
               nfxof(5))*nfxof(4))*nfxof(3))*nfxof(1)
              if(btest(YMNCdflag(i),3))&
               dresul(ksdfa(i)+3) =  0.0
              if(btest(YMNCdflag(i),4))&
               dresul(ksdfa(i)+4) =  0.0
              if(btest(YMNCdflag(i),5))&
               dresul(ksdfa(i)+5) =  0.0
              if(btest(YMNCdflag(i),6))&
               dresul(ksdfa(i)+6) =  0.0
          end select
        end if
 
        if(YMNCdfnum(i).gt.0) then
          do j = 1,YMNCnarg
            if(btest(YMNCdflag(i),j)) then
              dresul(ksdfa(i)+j) = dresul(ksdfa(i)+j) * nostep(j)
            end if
          end do
        end if
 
      end do
 
      RETURN
      END

      subroutine YNREFsection (  ARGUM, M, CONST)
      implicit none

      integer M

      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST ( YMNCNFNUM )

      real NFconst(7,5)

      data NFconst/&
! tip 1 - main reflector cell
             1.262,  & !D1
             0.3237,  & !D2
             0.002,& !Sa1
             0.05077,& !Sa2
             1.813e-2,& !S1-2
!             2.6e-2,& !S1-2
             1.,   & !v1
             1.,   & !v2
! tip 2 - neighbor to corner
             1.155,  & !D1
             0.299,  & !D2
             0.0022,& !Sa1
             0.0335,& !Sa2
             1.032e-2,& !S1-2
!             2.532e-2,& !S1-2
             1.,   & !v1
             1.,   & !v2
! tip 3 - corner
             1.505,  & !D1
             0.291,  & !D2
             0.00133,& !Sa1
             0.03035,& !Sa2
             2.778e-2,& !S1-2
             1.,   & !v1
             1.,   & !v2
! axial reflector
             1.4,  & !D1
             0.3,  & !D2
             0.003,& !Sa1
             5.e-2,& !Sa2
             3.e-2,& !S1-2
             1.,   & !v1
             1.,   & !v2
! axial reflector
             1.4,  & !D1
             0.3,  & !D2
             0.003,& !Sa1
             5.e-2,& !Sa2
             3.e-2,& !S1-2
             1.,   & !v1
             1./   !v2

      M=abs(M)
      const(:)=0.
      const(1)=NFconst(1,M)
      const(9)=NFconst(2,M)
      const(2)=NFconst(3,M)
      const(10)=NFconst(4,M)
      const(5)=NFconst(5,M)

      return
      end

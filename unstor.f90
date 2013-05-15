!**********************************************************************
! UNSTOR.FOR
! RAD FEB. 2009, BASED ON MAESUS.FOR
!
! This file contains all the subroutines related to understorey radiation
! calculations as part of MAESPA. Based on MAESUS, but with a number of
! changes. For example, understorey no longer consists of 'clumps', but
! LAI of the understorey is instead provided as input.
! Also contains the BEWDY model (see Medlyn et al. 2000, Can.J.For.Res.).
!
! INPUTUSSTR - Read understorey structure pars.
! INPUTUSPHY - Read understorey physiological pars.
! OUTPUTUS - Output understorey results (point-wise).
! INTERPUS - Interpolate understorey dimensions.
! BEWDYPARMS - Get parameters for BEWDY.
! BEWDY - The BEWDY model, gives APAR, PS and ET for understorey.
! PSMOSS - Moss photosynthesis.
!
!**********************************************************************

        
!**********************************************************************
SUBROUTINE INPUTUSSTR(NOUSPOINTS,X0,Y0,GRDAREAI,XLU,YLU,ZLU,USLAI,NOFUDATES,DATESFU,HTUS,&
                        NOHUDATES,DATESHU,FOLNUS,NONUDATES,DATESNU,EXTKUS)
!**********************************************************************
    USE maestcom
    IMPLICIT NONE
    INTEGER I,J,NOUSPOINTS,IOERROR,NOX,NOHUDATES
    INTEGER NONUDATES,NOFUDATES,NALPHA
    INTEGER DATESFU(maxdate)
    INTEGER DATESHU(maxdate),DATESNU(maxdate)
    REAL XLU(MAXP),YLU(MAXP),ZLU(MAXP)
    REAL COORDS(MAXP*2), XMAX, YMAX,X0,Y0
    REAL XWIDTH,YWIDTH,GRDAREAI,SPACING,EXTKUS
    REAL USLAI(maxdate,MAXT)
    REAL HTUS(maxdate,MAXT),FOLNUS(maxdate,MAXT)
    REAL ALPHA(MAXANG),FALPHA(MAXANG),BEXTANG(MAXANG)
    NAMELIST /CONTROL/ NOUSPOINTS,XMAX,YMAX,X0,Y0
 
    ! Read control flags: no of points and type of input
    READ (USTOREYI, CONTROL, IOSTAT = IOERROR)
    IF ((IOERROR.NE.0).OR.(NOUSPOINTS.EQ.0)) CALL SUBERROR('ERROR: MISSING CONTROL INFO IN USTOREY FILE',IFATAL,IOERROR)
    IF (NOUSPOINTS.GT.MAXT) THEN
        CALL SUBERROR('WARNING: TOO MANY USTOREY POINTS SPECIFIED',IWARN,IOERROR)
        NOUSPOINTS = MAXP
    END IF

    ! Place understorey points in regular grid (a square).
    XWIDTH = XMAX - X0 ! Width of square for understorey points.
    YWIDTH = YMAX - Y0
    GRDAREAI = XWIDTH * YWIDTH
      
    XLU = 0.0
    YLU = 0.0
    ZLU = 0.0

    SPACING = SQRT(XWIDTH*YWIDTH/NOUSPOINTS)
    NOX = FLOOR(SQRT(REAL(NOUSPOINTS)))

    DO I = 1,NOUSPOINTS
        XLU(I) = X0 + (I-1)/NOX * SPACING
        YLU(I) = Y0 + MOD(I-1,NOX) * SPACING
    END DO

    ! Read in dimensions of plants
    CALL READTREEARRAY(USTOREYI,3,NOUSPOINTS,NOHUDATES,DATESHU,HTUS)
    !CALL READTREEARRAY(USTOREYI,5,NOUSPOINTS,NOFUDATES,DATESFU,FOLUS)
    !CALL READTREEARRAY(USTOREYI,6,NOUSPOINTS,NODUDATES,DATESDU,DIAMUS)
    CALL READTREEARRAY(USTOREYI,7,NOUSPOINTS,NONUDATES,DATESNU,FOLNUS)
    CALL READTREEARRAY(USTOREYI,8,NOUSPOINTS,NOFUDATES,DATESFU,USLAI)
    ! Calculate foliar extinction coefficient - is in vertical direction only
!    CALL READLIA(USTOREYI, NALPHA, ALPHA, FALPHA)              ! a remettre !!!!!!!!!!! ATTENTION
    CALL EXBEAM(NALPHA,ALPHA,FALPHA,1.0,0.0,EXTKUS,BEXTANG)

    RETURN
END SUBROUTINE INPUTUSSTR

!**********************************************************************
SUBROUTINE INPUTUSPHY(JMAXN25I,IECOUI,EAVJI,EDVJI,DELSJI,TVJUPI,TVJDNI,     &
                        VCMAXN25I,EAVCI,EDVCI,DELSCI,UNMINI,AJQI,ABSRPI,    &
                        GSBG0U,GSBG1U,CICARATI,RD0I,RDKI,RDTI,SLAI,EFFYI,   &
                        MOSS,JMAX25MI,VCMAX25MI,THETAMI)
! Get understorey physiology parameters.
! 16/8/00 Change to take Ball-Berry params for gs only.
! 22/12/03 Add option to specify Ci:Ca ratio - helpful for moss
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    
    INTEGER UFILE,IECO,MOSS,IOERROR,IECOUI
    REAL JMAXN25,JMAXN25I,JMAX25M,JMAX25MI
    REAL EAVJ,EDVJ,DELSJ,AJQ,EAVC,EDVC,DELSC,TVJUP,TVJDN
    REAL VCMAXN25,UNMIN,ABSRP,JMAX25,VCMAX25,THETA,CICARAT
    REAL G0,G1,RD,RDK,RDT,SLA,EFFY,VCMAX25MI
    REAL THETAMI,VCMAXN25I,UNMINI,ABSRPI,EAVCI
    REAL EDVCI,EAVJI,EDVJI,DELSJI,AJQI,DELSCI,TVJUPI
    REAL TVJDNI,CICARATI,GSBG0U,GSBG1U,RD0I,RDKI
    REAL RDTI,SLAI,EFFYI

    NAMELIST /JMAXPARS/ IECO,EAVJ,EDVJ,DELSJ,AJQ
    NAMELIST /VCMAXPARS/ EAVC,EDVC,DELSC,TVJUP,TVJDN
    NAMELIST /BEWDYPARS/ JMAXN25,VCMAXN25,UNMIN,ABSRP
    NAMELIST /MOSSPARS/ JMAX25,VCMAX25,THETA
    NAMELIST /CICA/ CICARAT
    NAMELIST /BBGS/ G0, G1
    NAMELIST /RDPARS/ RD,RDK,RDT,SLA,EFFY

    UFILE = USTOREYI

    ! Read in moss params if applicable
    REWIND (UFILE)
    READ (UFILE,MOSSPARS,IOSTAT = IOERROR)
    IF (IOERROR.EQ.0) THEN
        JMAX25MI = JMAX25
        VCMAX25MI = VCMAX25
        THETAMI = THETA
        MOSS = 1
    ELSE
        MOSS = 0
	      ! Otherwise Read in BEWDY-specific params
    REWIND (UFILE)
    READ (UFILE,BEWDYPARS,IOSTAT = IOERROR)
    IF (IOERROR.NE.0) CALL SUBERROR('INPUT ERROR: MISSING UNDERSTOREY BEWDY PARAMS',IFATAL,IOERROR)
        JMAXN25I = JMAXN25
        VCMAXN25I = VCMAXN25
        UNMINI = UNMIN
        ABSRPI = ABSRP
    END IF

    ! Read in T-response parameters
    REWIND (UFILE)
    AJQ = ALPHAQ !Default values
    IECO = 1 ! Ecocraft formulation of T-deps of Km and Gamma. 
             ! For Montpied formulation, put 0. 
    EDVC = 0.0
    DELSC = 0.0
    TVJUP = -100.0
    TVJDN = -100.0

    READ (UFILE,JMAXPARS,IOSTAT = IOERROR)
    IF (IOERROR.NE.0)CALL SUBERROR('INPUT ERROR: MISSING UNDERSTOREY JMAXPARS',IFATAL,IOERROR)
    REWIND (UFILE)
    READ (UFILE,VCMAXPARS,IOSTAT = IOERROR)
    IF (IOERROR.NE.0)CALL SUBERROR('INPUT ERROR: MISSING UNDERSTOREY VCMAXPARS',IFATAL,IOERROR)

    EAVCI = EAVC
    EDVCI = EDVC
    DELSCI = DELSC
    EAVJI = EAVJ
    EDVJI = EDVJ
    DELSJI = DELSJ
    AJQI = AJQ
    IECOUI = IECO
    TVJUPI = TVJUP
    TVJDNI = TVJDN

    ! Read in stomatal conductance params (Ball-Berry model only or Ci:Ca ratio)
    REWIND (UFILE)
    CICARATI = 0.0
    READ (UFILE, CICA,IOSTAT = IOERROR)
    IF (IOERROR.EQ.0) THEN
        CICARATI = CICARAT
        GSBG0U = 0.0
        GSBG1U = 0.0
    ELSE
        REWIND (UFILE)
        READ (UFILE, BBGS,IOSTAT = IOERROR)
        IF (IOERROR.NE.0)CALL SUBERROR('INPUT ERROR: MISSING UNDERSTOREY GS PARS',IFATAL,IOERROR)
        GSBG0U = G0
        GSBG1U = G1
        IF (G1.LT.0.0) CALL SUBERROR('ERROR IN GS PARAMETERS: G1 MUST BE > 0',IFATAL,0)
    END IF
    
    ! Read in respiration parameters
    REWIND(UFILE)
    RD = 0.0
    RDK = 0.0
    RDT = 0.0
    SLA = -1.0
    EFFY = 0.0
    READ (UFILE, RDPARS, IOSTAT = IOERROR)
    IF (IOERROR.NE.0) CALL SUBERROR('WARNING: MISSING UNDERSTOREY RD PARS', IWARN,IOERROR)
    RD0I = RD
    RDKI = RDK
    RDTI = RDT
    SLAI = SLA
    EFFYI = EFFY

    RETURN
END SUBROUTINE InputUSPhy

		

!**********************************************************************
SUBROUTINE OUTPUTUS(IDAY,NOUSPOINTS,XLU,YLU,ZLU,UIBEAM,UIDIFF,PARUS,APARUS,PSUS,ETUS)

! Output understorey results.
! RAD, January 2009.
!**********************************************************************

    USE maestcom
    USE switches
    IMPLICIT NONE
    INTEGER I,J,NOUSPOINTS,IDAY
    
    REAL UIBEAM(MAXHRS,MAXP),UIDIFF(MAXHRS,MAXP)
    REAL APARUS(MAXHRS,MAXP),PARUS(MAXHRS,MAXP)
    REAL PSUS(MAXHRS,MAXP),ETUS(MAXHRS,MAXP)
    REAL XLU(MAXP),YLU(MAXP),ZLU(MAXP)
    
    IF (IOFORMAT .EQ. 0) THEN 
        DO I = 1,KHRS
            DO J = 1,NOUSPOINTS
                WRITE(UPARUS,997)IDAY,I,J,XLU(J),YLU(J),ZLU(J),UIBEAM(I,J),UIDIFF(I,J),&
                        PARUS(I,J),APARUS(I,J),PSUS(I,J),ETUS(I,J)
            END DO
        END DO    
        997 FORMAT (I3,1X,2(I4,1X),10(F7.2, 1X))
    ELSE IF (IOFORMAT .EQ. 0) THEN
        DO I = 1,KHRS
            DO J = 1,NOUSPOINTS
                WRITE(UPARUS) REAL(IDAY), REAL(I),REAL(J),XLU(J),YLU(J),ZLU(J),&
                      UIBEAM(I,J),UIDIFF(I,J),PARUS(I,J),APARUS(I,J),PSUS(I,J),ETUS(I,J)
            END DO
        END DO  
    END IF
    RETURN
END SUBROUTINE OUTPUTUS


!**********************************************************************
SUBROUTINE INTERPUS(IDAY,ISTART,NUMPNT,FNMIN,EXTK,GRDAREAI,DATESFU,NOFUDATES,USLAITAB,&
                    USLAI,DATESHU,NOHUDATES,HTUS,ZLU,DATESNU,NONUDATES,FOLNUS,FN0US,AREAUS)
! Interpolate values of understorey parameters
!**********************************************************************
    USE maestcom
    IMPLICIT NONE
    INTEGER IDAY,ISTART,NUMPNT,NOHUDATES,NONUDATES,NOFUDATES,IPT

    ! Dates for tree dimensions
    INTEGER DATESFU(maxdate),DATESHU(maxdate)
    INTEGER DATESNU(maxdate)
	
    ! Understorey plant dimensions
    REAL USLAITAB(maxdate,MAXT)
    REAL HTUS(maxdate,MAXT),FOLNUS(maxdate,MAXT)
    REAL USLAI(MAXP),AREAUS(MAXP),FN0US(MAXP)
    REAL FNUS(MAXP),LAI,ZLU(MAXP)
    REAL FNMIN,EXTK,GRDAREAI
    
    CALL TREEINTERP(IDAY,ISTART,NOFUDATES,DATESFU,USLAITAB,NUMPNT,USLAI)
    CALL TREEINTERP(IDAY,ISTART,NOHUDATES,DATESHU,HTUS,NUMPNT,ZLU)
	CALL TREEINTERP(IDAY,ISTART,NONUDATES,DATESNU,FOLNUS,NUMPNT,FNUS)
    
    DO IPT = 1,NUMPNT
        AREAUS(IPT) = GRDAREAI / REAL(NUMPNT) ! Probably not needed anymore !Was: PI*(DUS(IPT)**2)/4.
	    ! Assign nitrogen contents.
        LAI = USLAI(IPT)
        FN0US(IPT) = (FNUS(IPT)-FNMIN)*LAI*EXTK/(1.-EXP(-LAI*EXTK))+ FNMIN
    END DO

    RETURN
END SUBROUTINE INTERPUS

!**********************************************************************
SUBROUTINE BEWDYPARMS(IHOUR,TAIR,RH,CA,JMAXN25,IECOU,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN,&
                        VCMAXN25,EAVC,EDVC,DELSC,AJQ,G0,G1,CICARAT,BALPHA,BLAMBDA)
! Calculates parameters for BEWDY model (Alpha, Lambda) from the Farquhar/Leuning models.
!**********************************************************************
    IMPLICIT NONE
    INTEGER IHOUR,IECOU
    
    REAL TAIR,CA,RH
    REAL JMAXN25,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN
    REAL VCMAXN25,EAVC,EDVC,DELSC
    REAL AJQ,G0,G1
    REAL BALPHA,BLAMBDA,CICARAT,DRAWDOWN
    REAL GAMMASTAR,KM,JMAXN,VCMAXN,CI,AJMAX,ACMAX
    REAL, EXTERNAL :: GAMMAFN
    REAL, EXTERNAL :: KMFN
    REAL, EXTERNAL :: JMAXTFN
    REAL, EXTERNAL :: VCMAXTFN
    
    ! Calculate photosynthetic parameters from leaf temperature - here assumed same as air T.
    GAMMASTAR = GAMMAFN(TAIR,IECOU)                   ! CO2 compensation point, umol mol-1
    KM = KMFN(TAIR,IECOU)                             ! Michaelis-Menten for Rubisco, umol mol-1
    JMAXN = JMAXTFN(JMAXN25,TAIR,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Potential electron transport rate, umol m-2 s-1
    VCMAXN = VCMAXTFN(VCMAXN25,TAIR,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Maximum Rubisco activity, umol m-2 s-1

    ! Calculate ACMAX
    !      GSDIVA = G1 /1.6 * RH / (CS - GAMMASTAR)
    !      A = G0 + GSDIVA * VCMAXN
    !      B = (1. - CS*GSDIVA) * VCMAXN + G0 * (KM - CS)
    !     +      - GSDIVA * VCMAXN*GAMMASTAR 
    !      ! = -(1. - CS*GSDIVA) * VCMAXN*GAMMASTAR - G0*KM*CS
    !      CIC = QUADP(A,B,!,IQERROR)

    !      IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
    !        ACMAX = 0.0
    !      ELSE
    !        ACMAX = VCMAXN * (CIC - GAMMASTAR) / (CIC + KM)
    !      END IF

    ! Calculate AJMAX
    !      A = G0 + GSDIVA * JMAXN
    !      B = (1. - CS*GSDIVA) * JMAXN + G0 * (2.*GAMMASTAR - CS)
    !     +      - GSDIVA * JMAXN*GAMMASTAR
    !      ! = -(1. - CS*GSDIVA) * GAMMASTAR * JMAXN 
    !     +      - G0*2.*GAMMASTAR*CS
    !      CIJ = QUADP(A,B,!,IQERROR)

    !      IF ((IQERROR.EQ.1).OR.(CIJ.LE.0.0).OR.(CIJ.GT.CS)) THEN
    !        AJMAX = 0.0
    !      ELSE
    !        AJMAX = JMAXN/4. * (CIJ - GAMMASTAR) / (CIJ + 2*GAMMASTAR)
    !      END IF

    ! Calculate BALPHA
    !      A = G0 + GSDIVA * AJQ
    !      B = (1. - CS*GSDIVA) * AJQ + G0 * (2.*GAMMASTAR - CS)
    !     +      - GSDIVA * AJQ*GAMMASTAR
    !      ! = -(1. - CS*GSDIVA) * GAMMASTAR * AJQ 
    !     +      - G0*2.*GAMMASTAR*CS
    !      CIJ = QUADP(A,B,!,IQERROR)

    !      IF ((IQERROR.EQ.1).OR.(CIJ.LE.0.0).OR.(CIJ.GT.CS)) THEN
    !        BALPHA = 0.0
    !      ELSE
    !        BALPHA = AJQ/4. * (CIJ - GAMMASTAR) / (CIJ + 2*GAMMASTAR)
    !      END IF

    ! Calculate Ci from Ball-Berry function
    !      CI = CA - (CA-GAMMASTAR)*(1+VPD/GSLD0)*1.6/GSLG1
    IF (CICARAT.GT.0.0) THEN
        CI = CA*CICARAT
    ELSE 
        drawdown = (CA-GAMMASTAR)/RH*1.6/G1
        CI = CA - drawdown
        IF(CI.LT.0)CI = 50
    END IF

    ! Calculate Bewdy parameters
    BALPHA = AJQ/4.0*(CI-GAMMASTAR)/(CI+2.*GAMMASTAR)  
    AJMAX = JMAXN/4.0*(CI-GAMMASTAR)/(CI+2.*GAMMASTAR) 
    ACMAX = VCMAXN*(CI-GAMMASTAR)/(CI+KM) 

    IF (AJMAX.GT.ACMAX) THEN
        BLAMBDA = ACMAX
    ELSE 
        BLAMBDA = AJMAX
	END IF
    RETURN
END SUBROUTINE BEWDYPARMS



!**********************************************************************
 SUBROUTINE BEWDY(IHOUR,BEAM,DIFF,SUNLA,BALPHA,BLAMBDA,LEAFN0,LEAFNMIN,&
                    KEXT,ABSRP,LAI,APAR,PS,PARUNDER)
     
! Calculates ! assimilation using the BEWDY model. 
! See Medlyn et al. 2000 (CJFR, Appendix A).
!**********************************************************************
    IMPLICIT NONE
    INTEGER IHOUR
    REAL L, LAI, KEXT, LEAFN0, LEAFNMIN
    REAL APAR, PARUNDER, B, D, PS, DIFF, SUNLA, BEAM, ABSRP
    REAL BALPHA, BLAMBDA
    
    APAR = (DIFF+SUNLA*BEAM)*ABSRP*(1.-EXP(-KEXT*LAI))
    PARUNDER = (DIFF+SUNLA*BEAM) - APAR  ! RAD
    B = BALPHA * KEXT * BEAM * ABSRP
    D = BALPHA * KEXT * DIFF * ABSRP
    L = BLAMBDA * (LEAFN0 - LEAFNMIN)
    PS = 1./KEXT*(1.-EXP(-KEXT*LAI))*(L*D*(L+D)+B*L*L)/(L+D)/(L+D)
    PS = PS + 1./KEXT*(B*B*L*L)/((L+D)**3)*LOG(((L+D)*EXP(-KEXT*LAI)+B)/(L+D+B))
    PS = SUNLA*PS + (1-SUNLA)*1./KEXT*(1.-EXP(-KEXT*LAI))*(L*D)/(L+D)
    RETURN
END SUBROUTINE BEWDY

!**********************************************************************
SUBROUTINE PSMOSS(APAR,TLEAF,RH,CA,JMAX25,IECO,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN,&
                    VCMAX25,EAVC,EDVC,DELSC,AJQ,THETA,PS)
! Calculates assimilation by moss. 
!**********************************************************************
    IMPLICIT NONE
    INTEGER IECO, IQERROR 
    REAL APAR,TLEAF,RH,CA
    REAL JMAX25,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN
    REAL VCMAX25,EAVC,EDVC,DELSC
    REAL AJQ,THETA,PS
    REAL GAMMASTAR,JMAX,VCMAX,KM,J,VJ,AJ,AC
    REAL, EXTERNAL :: GAMMAFN, KMFN, JMAXTFN, VCMAXTFN, QUADM
    ! I presume this is 0 as it certainly isn't set anywhere
    ! Seems strange that it is passed but neve set, why bother
    ! passing it?!
    
    ! Calculate photosynthetic parameters from leaf temperature.
    GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensation point, umol mol-1
    KM = KMFN(TLEAF,IECO)                             ! Michaelis-Menten for Rubisco, umol mol-1
    JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Potential electron transport rate, umol m-2 s-1
    VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Maximum Rubisco activity, umol m-2 s-1
    J = QUADM(THETA,-(AJQ*APAR+JMAX),AJQ*APAR*JMAX,IQERROR) ! Actual e- transport rate, umol m-2 s-1
    VJ = J/4.0                                        ! RuBP-regen rate, umol m-2 s-1

    ! Farquhar model
    AJ = VJ*(CA-GAMMASTAR)/(CA+2*GAMMASTAR)
    AC = VCMAX*(CA-GAMMASTAR)/(CA+KM)
    PS = MIN(AJ,AC)

    RETURN
END SUBROUTINE PSMOSS

************************************************************************
*
*  subroutines for the computation of insolation
*  J. Laskar, F. Joutel, F. Boudin
* (c) Astronomie et Systemes Dynamiques/ Institut de Mecanique Celeste (1993-2010)
*
* 05/11/2001 Mickael : modification de l'appel de polint dans qromb
*                      (mauvais type fourni)
* 22/09/2003 Mickael : utilisation de e et pibarh au lieu de ak,ah, psi
* 05/11/2003 Mickael : changement du common date en ddate d
* 26/05/2004 Mickael : correction qromb suivant numerical recipes 2nd ed.
* 12/01/2010 JxL : correction de vraimoy et moyvrai
* 15/01/2010 Mickael : correction dans les routines des "." en ".D0"
*************************************************************************
       subroutine vraimoy(hl,e,pibar,hlm)
c***********************************************************************
c  calcule la longitude moyenne (hlm) en fonction de la longitude      *
c  vraie (hl), de l'excentricite (e) et de la longitude du perihelie   *
c  (pibar)                                                             *
c  (c) ASD/BDL (1/10/92)  
c  correction  27/1/93                                                 *
c  correction d'une erreur et extension a l'ordre 5 (JxL 12/1/10)
c***********************************************************************
      implicit none
      real(8), intent(in) :: hl,e,pibar
      real(8), intent(out) :: hlm
      real(8)  eV
      eV = hl-pibar
      hlm = hl 
     & -2*e*sin(eV)
     & +(3*e**2/4+e**4/8)*sin(2*eV)
     & -(e**3/3+e**5/8)*sin(3*eV)
     & +5*e**4/32*sin(4*eV) -3*e**5/40*sin(5*eV)
      return
      end

      subroutine moyvrai(hlm,e,pibar,hl)
c***********************************************************************
c  calcule la longitude vraie (hl) en fonction de la longitude         *
c  moyenne (hlm), de l'excentricite (e) et de la longitude du          *
c  perihelie (pibar)                                                   *
c  (c) ASD/BDL (1/10/92)
c  etendu a l'ordre 5 (JxL 12/1/10)
c***********************************************************************
      implicit none
      real(8), intent(in) :: hlm,e,pibar
      real(8), intent(out) :: hl
      real(8) eM
      eM = hlm-pibar
      hl = hlm 
     &+(2*e-e**3/4 + 5*e**5/96)*sin(eM)
     &+(5*e**2/4 - 11*e**4/24)*sin(2*eM)
     &+(13*e**3/12 - 43*e**5/64)*sin(3*eM)
     &+ 103*e**4/96*sin(4*eM) + 1097*e**5/960*sin(5*eM)
      return
      end

       SUBROUTINE cwj(wd,e,pibar,eps,phi,w)
c***********************************************************************
c                                                                      *
c      INSOLATION JOURNALIERE D'UN POINT DE LATITUDE DONNEE            *
c                                                                      *
c A L'ENTREE:                                                          *
c   wd : la longitude vraie du soleil (en radian)  compte a partir     *
c        de l'equinoxe de la date                                      *
c   e : excentricite                                                   *
c   pibar: longitude du perihelie comptee a partir de l'equinoxe de la *
c   date + pi (repere geocentrique)                                    *
c   eps: l'obliquite                                                   *
c   phi: latitude du point sur la terre                                *
c                                                                      *
c EN SORTIE:                                                           *
c   w:l'insolation journaliere d'un point de latitude donnee           *
c  (c) ASD/BDL (1/10/92)                                              *
c***********************************************************************
! 15/01/2010 Mickael : ajout des .D0 
      implicit none
      real(8), intent(in) :: e, eps, pibar, phi, wd
      real(8), intent(out) :: w
      real(8) ::  so,pi,v,cosv,aux,rho,sinusdelta,delta,cho,ho,a1,a2
      parameter(pi=3.1415926535897932d0)
      common/cdo/so
c
c-----------------------------------------------------------------------
c     CALCUL D'INSOLATION
c-----------------------------------------------------------------------
c
c anomalie vraie (v):
      v = wd-pibar
c
c Distance terre-soleil (rho):
      cosv= cos(v)
      aux=1.D0+e*cosv
      rho=(1.D0-e**2)/aux
c
c Declinaison du soleil (delta):
      sinusdelta=sin(eps)*sin(wd)
      delta=asin(sinusdelta)
c
c
c LATITUDES DE LEVER ET COUCHER DE SOLEIL:
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      aux=pi/2.d0-abs(delta)
      if ((-aux).lt.phi.and.phi.lt.aux) then
c |Angle horaire| du lever et coucher du soleil (ho):
          cho=-tan(phi)*tan(delta)
          ho=acos(cho)
c Insolation (w)
          w=(ho*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(ho))
          w=w*so/(pi*rho**2)
          goto 200
      endif
c
c LATITUDES SANS COUCHER:
c^^^^^^^^^^^^^^^^^^^^^^^^
      a1=pi/2.d0-delta
      a2=pi/2.d0+delta
      if (phi.ge.a1.or.phi.le.(-a2)) then
          w=so*sin(phi)*sin(delta)/rho**2
          goto 200
      endif
c
c LATITUDES SANS LEVER:
c^^^^^^^^^^^^^^^^^^^^^^
      if (phi.le.(-a1).or.phi.ge.a2) w=0.d0
c
 200  return
      end


      SUBROUTINE wmcal(mois,e,eps,pibarh,phi,w)
c***********************************************************************
c
c      INSOLATION MENSUELLE D'UN POINT DE LATITUDE DONNEE
c
c A L'ENTREE:
c   mois : le numero du mois  (l'annee est divisee en 12 mois de 30
c          degres et le numero 3 correspond a  fin fevrier- fin mars)
c   e,eps,pibarh: les elements orbitaux pour la date tps
c   e : excentricite 
c   eps:l'obliquite
c   pibarh: longitude du perihelie comptee a partir de l'equinoxe de la
c   date
c   phi: latitude du point sur la terre
c
c EN SORTIE:
c   w:l'insolation mensuelle d'un point de latitude donnee
c (c) ASD/IMC (22/09/2003)
c***********************************************************************
      implicit none
      real(8), intent(in) :: e, eps, pibarh, phi
      integer, intent(in) :: mois
      real(8), intent(out) :: w
      real(8) :: edu, pibar,phidu,epsdu, so, pi, hlm0, hlm1, hlm2
      parameter(pi=3.1415926535897932d0)
      common/cdo/so
      common/cfunc/edu,pibar,phidu,epsdu
      external F
c variables 'dummy'
      phidu=phi
      epsdu=eps
      edu = e
c
c-----------------------------------------------------------------------
c     CALCUL D'INSOLATION
c-----------------------------------------------------------------------
c
c longitude du perihelie par rapport a l'equinoxe de reference(pitild),
c par rapport a l'equinoxe de la date(pitbar), excentricite(e):
      pibar=pibarh+pi
c
c Longitude moyenne au 21 mars
c
      call vraimoy(0.D0,e,pibar,hlm0)
c
c  longitude moyenne (hlm1) et vraie (hl1) au debut du mois
      hlm1 = hlm0+(mois-4)*pi*30.D0/180D0
c longitude moyenne (hlm2) et vraie (hl1) a la fin du mois
      hlm2 = hlm1+30.D0*pi/180D0
c calcul de l'insolation moyenne a l'aide de la
c methode de Romdberg
      call qromb(F,hlm1,hlm2,w)
      w=w/30.D0/pi*180.D0
      return
      end

      double precision function F(hlm)
      implicit double precision (a-h,p-y)
      common/cfunc/e,pibar,phi,eps
      call moyvrai(hlm,e,pibar,hl)
      call cwj(hl,e,pibar,eps,phi,w)
      F=w
      return
      end


      SUBROUTINE wjour(date,e,eps,pibarh,phi,w)
c***********************************************************************
c
c      INSOLATION JOURNALIERE D'UN POINT DE LATITUDE DONNEE
c
c A L'ENTREE:
c   date : la longitude vraie du soleil (en degres)  comptee a partir
c          de l'equinoxe vrai
c   e,eps,pibarh: les elements orbitaux pour la date tps
c   e : excentricite 
c   eps:l'obliquite
c   pibarh: longitude du perihelie comptee a partir de l'equinoxe de la
c   date
c   phi: latitude du point sur la terre
c
c EN SORTIE:
c   w:l'insolation journaliere d'un point de latitude donnee
c (c) ASD/IMC (22/09/2003)
c***********************************************************************
      implicit none
      real(8), intent(in) :: e, eps, pibarh, phi, date
      real(8), intent(out) :: w
      real(8) ::  pibar, so, pi, hl
      parameter(pi=3.1415926535897932d0)
      common/cdo/so
c
c-----------------------------------------------------------------------
c     CALCUL D'INSOLATION
c-----------------------------------------------------------------------
c
c longitude du perihelie par rapport a l'equinoxe de reference(pitild),
c par rapport a l'equinoxe de la date(pibar), excentricite(e):
      pibar=pibarh+pi
c
c Longitude vraie comptee a partir de l'equinoxe vrai et
c anomalie vraie (wd,v):
      hl = date*pi/180.D0
      call cwj(hl,e,pibar,eps,phi,w)
      end


      SUBROUTINE wjcal(datecal,e,eps,pibarh,phi,w)
c***********************************************************************
c
c      INSOLATION JOURNALIERE D'UN POINT DE LATITUDE DONNEE
c
c A L'ENTREE:
c   datecal: la longitude moyenne (en degres) du soleil comptee a partir
c          de l'equinoxe de la date
c   e,eps,pibarh: les elements orbitaux pour la date tps
c   e : excentricite 
c   eps:l'obliquite
c   pibarh: longitude du perihelie comptee a partir de l'equinoxe de la
c   date
c   phi: latitude du point sur la terre
c
c EN SORTIE:
c   w:l'insolation journaliere d'un point de latitude donnee
c (c) ASD/BDL (1/10/92)
c***********************************************************************
      implicit none
      real(8), intent(in) :: e, eps, pibarh, phi, datecal
      real(8), intent(out) :: w
      real(8) ::  so, pi, hlm, hlm0, pibar, wd
      parameter(pi=3.1415926535897932d0)
      common/cdo/so
c
c-----------------------------------------------------------------------
c     CALCUL D'INSOLATION
c-----------------------------------------------------------------------
c
c longitude du perihelie par rapport a l'equinoxe de reference(pitild),
c par rapport a l'equinoxe de la date(pibar), excentricite(e):
      pibar=pibarh+pi
c calcul de la longitude vraie
      hlm=datecal*pi/180.D0
c longitude moyenne au 21 mars
      call vraimoy(0.D0,e,pibar,hlm0)
c longitude moyenne a la date
      hlm=hlm0+hlm
c longitude vraie a la date
      call moyvrai(hlm,e,pibar,wd)
      call cwj(wd,e,pibar,eps,phi,w)
      end

      SUBROUTINE Telinsol(tps,Tel)
c********************************************************************
c
c           LES ELEMENTS CLIMATIQUES k,h,eps,phi
c           POUR UNE DATE DONNEE tps (en nombre entier de pas).
c
c A L'ENTREE:
c  tps : la date pour laquelle on dsire les elements (en annees)
c
c EN SORTIE:
c  Tel(nbf) : tableau des lments orbitaux (k,h,eps,phi) pour la
c             date tps
c  (c) ASD/BDL (1/10/92)
c********************************************************************
      implicit double precision (a-h,o-z)
      parameter (nbel1=5,nbf=4)
      parameter (nacd=100)
      dimension Tel(nbf),Taux(nbel1,nacd)
      data ncour /-1/
      save ncour
      save Taux
      common/ddate/pas,itdebut
c nrec est le numero d'ordre de l'enregistrement qui contient tps
      tt=tps/abs(pas)
c itt est l'entier le plus proche  de tt
      itt=nint(tt)
      nrec=(itt-itdebut)/nacd+1
      if (nrec.ne.ncour) then
         read(10,rec=nrec) ((Taux(i,j),i=1,nbel1),j=1,nacd)
         ncour=nrec
      endif
      itdeb=itdebut+(nrec-1)*nacd
      iplace=itt-itdeb+1
      do i=2,nbel1
       Tel(i-1)=Taux(i,iplace)
      end do
      end

      SUBROUTINE wam(e,w)
c***********************************************************************
c
c                  INSOLATION  ANNUELLE MOYENNEE
c
c A L"ENTREE:
c  e:excentricite pour une date donnee
c
c EN SORTIE:
c  w:insolation annuelle moyenne
c  (c) ASD/BDL (1/10/92)
c***********************************************************************
! 15/01/2010 Mickael : ajout des .D0 
      implicit none
      real(8), intent(in) :: e
      real(8), intent(out) :: w
      real(8) rho, so, pi
      parameter(pi=3.1415926535897932d0)
      common/cdo/so
c-----------------------------------------------------------------------
c     CALCUL D'INSOLATION
c-----------------------------------------------------------------------
c Carre de la distance terre-soleil moyennee sur une revolution (rho):
      rho=dsqrt(1.D0-e**2)
c Calcul d"insolation (w):
      w=(so)/(4.D0*rho)
      return
      end


      SUBROUTINE soresul(ifich,x,y)
c***********************************************************************
c
c                     SORTIE DES RESULTATS
c
c  PARAMETRES D'ENTREE:
c  ifich : le numero de connection du fichier de sortie
c  x     :la date
c  y     : la valeur de la fonction de la date x
c (c) ASD/BDL (1/10/92)
c***********************************************************************
c
      implicit double precision (a-h,o-z)
      logical ecriture
      common/fichier/ecriture
c
1000  format(1x,f10.0,d24.16)
      if (ecriture) then
        write(ifich,1000) x*1D-3,y
      else
        write(*,1000) x*1D-3,y
      endif
      return
      end


      SUBROUTINE QROMB(FUNC,A,B,SS)
c******************************************************************************
c   NUMERICAL RECIPES                                                         *
c   ROMBERG INTEGRATION                                                       *
c******************************************************************************
c 05/11/2001 Mickael : ajout de 0.D0 (au lieu de 0.) pour l'appel de polint
c 26/05/2004 Mickael : correction du test IF (ABS(DSS).LT. en .LE.
c correction effectuee dans "numerical recipes, second edition"
! 15/01/2010 Mickael : ajout des .D0 + replacement des "0.25*" en "/4.D0"
      IMPLICIT DOUBLE PRECISION (A-H,P-Y)
      EXTERNAL FUNC
      PARAMETER(EPS=1.D-6,JMAX=20,JMAXP=JMAX+1,K=5,KM=4)
      DIMENSION S(JMAXP),H(JMAXP)
      H(1)=1.D0
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          L=J-KM
          CALL POLINT(H(L),S(L),K,0.D0,SS,DSS)
C          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/4.D0
11    CONTINUE
      PAUSE 'Too many steps.'
      END


      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H,P-Y)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.D0)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


      SUBROUTINE TRAPZD(FUNC,A,B,S,N)
      IMPLICIT DOUBLE PRECISION (A-H,P-Y)
      EXTERNAL FUNC
      SAVE IT
      IF (N.EQ.1) THEN
        S=0.5D0*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.D0
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

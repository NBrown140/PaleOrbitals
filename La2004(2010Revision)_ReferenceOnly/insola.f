       program insola
*******************************************************************************
*                                                                             *
*              Computation of the various insolation functions                *
*                                                                             *
*              this program uses                                              *
*                   insolsub.f   : subroutines for the computation            *
*                                  of insolation                              *
*                   prepinsol.f  : subroutines for the binary temporary file  *
*                   prepsub.f    : subroutines for computing elements at      *
*                                  specified time                             *
*                   insola.par   : parameters for insola.f                    *
*                   nomclimaneg  ,                                            *
*                   nomclimapos  : ASCII files supplied                       *
*                                  which contains the necessary               *
*                                  orbital and precession quantities          *
*                                                                             *
*                                                                             *
*                                                                             *
*   the following parameters are taken from the NAMELIST insola.par           *
*           nomclimapos : ASCII file for the elements t,e,eps, psibar on      *
*                         positive time.                                      *
*           nomclimaneg : ASCII file for the elements t,e,eps, psibar on      *
*                         negative time.                                      *
*           nominsolbin : temporary binary file for the elements t,e,eps,     *
*                         cos(psibar), sin(psibar)                            *
*                         (variables interne = nomfich)                       *
*           pas         : stepsize (years)                                    *
*           datefin     : final time (million of years)                       *
*           datedebut   : starting time  (million of years)                   *
*           so          : solar constant   W*m^-2                             *
*                                                                             *
*   other parameters:                                                         *
*                                                                             *
*           nbf         : 4 (number of variables)                             *
*           nbel1       : 4 + 1                                               *
*           nacd        : numbers of n-uplets in a record                     *
*                                                                             *
*  (c) Astronomie et Systemes Dynamiques/ institut de Mecanique Celeste (1993-2010)*
*  J. Laskar, F. Joutel, M. Gastineau                                         *
*                                                                             *
*  version 0.82 (7/7/1993)                                                    *
*  version 0.83 (27/1/93)                                                     *
*  version la2001b4 (26/06/2001) Mickael : changement du fichier de parametres*
*  version la2001b9 (12/02/2002) Mickael : ajout dans le fichier de parametres*
*  version la2003   (04/03/2003) Mickael : fichier de parametres              *
*  version la2003   (04/09/2003) Mickael : creation du fichier binaires       *
*                   (05/11/2003) Mickael : changement du common date en ddate *
*  version la2004   (26/05/2004) Mickael : correction qromb (insolsub.f)      *
*  version la2004b2 (12/01/2010) JxL : correction de vraimoy et moyvrai       *
*                   (15/01/2010) Mickael : correction . en .D0                *
*******************************************************************************

      implicit none
      integer nacd, nbel, nbf
      real(8), parameter :: pi=3.1415926535897932d0
      parameter(nacd=100,nbel=5,nbf=4)
      real(8) :: Tel(nbel)
      character(LEN=50) :: nomfich,fichout
      real(8) pasbin, so, aux, phi,e ,eps, pibarh,w, tps, pas, debut
      real(8) datejour, datecal  
      real(8) datedebut,datefin
      integer itdebut, nrecl, npt, i, ICHOI, moiscal
      common/cdo/so
      common/ddate/pasbin,itdebut
c
      character(LEN=50) nominsolbin, nomclimapos,nomclimaneg
      NAMELIST/NAMSTD/ datedebut,datefin,pas, 
     &nominsolbin, nomclimapos,nomclimaneg,so
     
      open(8,file='insola.par',form='formatted',status='old')
      read(8,NAMSTD)
      close(8)
      
c
c     verification du pas
c
      if (pas.lt.1.D0) then
      write (*,*) 'Step size (pas) must be greater than 1 year'
      write (*,*) 'INSOLA stopped'
      goto 7000
      endif

      nomfich = nominsolbin
      

c
c Creation du fichier binaire et recupere le pas des donnees 
c
      call transbin(datedebut, datefin, nominsolbin, nomclimapos, 
     &              nomclimaneg, nbf, nbel, pasbin) 
c
c Les fichiers ont une longueur d'enregistrement correspondant
      nrecl=8*nbel*nacd

      open(11,file=nomfich,status='old',access='direct',recl=nrecl)
c
      read(11,rec=1) debut
      itdebut=int(debut*1D3/pasbin)

c Nombre de lignes lues:
      datedebut=datedebut*1.D6
      datefin=datefin*1.D6
      aux=dabs(datefin-datedebut)
      npt=int((aux/dabs(pas))+1)
*-----------------------------------------------------------------------
*   principal menu
*-----------------------------------------------------------------------
      write (*,*)
      write (*,*)
      write (*,*) '***************************************************'
      write (*,*) '*                                                 *'
      write (*,*) '*                                                 *'
      write (*,*) '*                                                 *'
      write (*,*) '*                  INSOLA                         *'
      write (*,*) '*                                                 *'
      write (*,*) '*   computation of various insolation quantities  *'
      write (*,*) '*                                                 *'
      write (*,*) '*              (c) ASD/IMC (1993-2010)            *'
      write (*,*) '*                                                 *'
      write (*,*) '*                (revision 2010/01/18)            *'
      write (*,*) '*                                                 *'
      write (*,*) '*                                                 *'
      write (*,*) '*                                                 *'
      write (*,*) '*                                                 *'
      write (*,*) '***************************************************'
999   continue
      write (*,*)
      write (*,*)
      write (*,*) '    1 : mean daily    insolation/true longitude '
      write (*,*) '    2 : mean daily    insolation/mean longitude '
      write (*,*) '    3 : mean monthly  insolation '
      write (*,*) '    4 : mean annual   insolation '
      write (*,*) '    5 : HELP'
      write (*,*) '    6 : QUIT'
      write (*,*)
      write (*,*) ' your choice ? '
      read  (*,*) ICHOI
      GOTO (1000,2000,3000,4000,5000,7000), ICHOI
      write (*,*) ' incorrect choice '
      GOTO 999

*---------------------------------- mean daily/ true longitude
1000  CONTINUE
      write (*,*)' Latitude on the Earth (in degrees) ? '
      read (*,*)  phi
      write (*,*) ' latitude = ',phi
      write (*,*)' True longitude (in degrees) ?'
      read (*,*) datejour
      write (*,*) 'true longitude = ',datejour
      phi=phi*pi/180.d0
      GOTO 6000


*---------------------------------- mean daily/ mean longitude
2000  CONTINUE
      write (*,*)' Latitude on the Earth (in degrees) ? '
      read (*,*)  phi
      write (*,*) ' latitude = ',phi
      write (*,*)
     $'approximate conventional dates for a given mean longitude'
      write (*,*)
      write (*,*)'  0: 21 march    '
      write (*,*)' 30: 21 april    '
      write (*,*)' 60: 21 may      '
      write (*,*)' 90: 21 june     '
      write (*,*)'120: 21 july     '
      write (*,*)'150: 21 august   '
      write (*,*)'180: 21 september'
      write (*,*)'210: 21 october  '
      write (*,*)'240: 21 november '
      write (*,*)'270: 21 december  '
      write (*,*)'300: 21 january   '
      write (*,*)'330: 21 february  '
      write (*,*)
      write (*,*)' MEAN longitude (in degrees)  (0-360)?'
      read (*,*) datecal
      write (*,*) 'MEANlongitude = ',datecal
      phi=phi*pi/180.d0
      GOTO 6000

*---------------------------------- mean monthly
3000  CONTINUE
      write (*,*)' Latitude on the Earth (in degrees) ? '
      read (*,*)  phi
      write (*,*) ' latitude = ',phi
      write (*,*)
      write (*,*)'approximate conventional dates of the months'
      write (*,*)
      write (*,*)' 1: 21 december  - 20 january'
      write (*,*)' 2: 21 january   - 20 february'
      write (*,*)' 3: 21 february  - 20 march'
      write (*,*)' 4: 21 march     - 20 april'
      write (*,*)' 5: 21 april     - 20 may'
      write (*,*)' 6: 21 may       - 20 june'
      write (*,*)' 7: 21 june      - 20 july'
      write (*,*)' 8: 21 july      - 20 august'
      write (*,*)' 9: 21 august    - 20 september'
      write (*,*)'10: 21 september - 20 october'
      write (*,*)'11: 21 october   - 20 november'
      write (*,*)'12: 21 november  - 20 december'
      write (*,*)
      write (*,*) 'Your choice ? (1-12)'
      read (*,*) moiscal
      write (*,*) 'month = ',moiscal
      phi=phi*pi/180.d0
      GOTO 6000
*---------------------------------- mean annual
4000  CONTINUE
      GOTO 6000


*--------------------------- Computation of insolation ------------
6000  continue
      write (*,*)
      write (*,*) 'Name of the output file'
      read(*,1010) fichout
1010  format (A)
      open (22,file = fichout,form='formatted')
      tps=datedebut
      do i=1,npt

         call Telor(tps,Tel)
         e=Tel(1)
         eps=Tel(2)
         pibarh=atan2(Tel(4),Tel(3))
         
         select case(ICHOI)
          case (1) 
            call wjour(datejour,e,eps,pibarh,phi,w)
            
          case (2) 
            call wjcal(datecal,e,eps,pibarh,phi,w)
            
          case (3) 
            call wmcal(moiscal,e,eps,pibarh,phi,w)
            
          case (4) 
            call wam(e,w)
            
          case default 
            write(*,*) 'Incorrect choice'
         end select

         write(22,1020) tps*1D-3,w 
1020  format(1x,f10.3,d20.12)
         tps=tps+pas
      end do
      close (22)
      GOTO 999
*--------------------------- HELP ------------------------------------
5000  write (*,*)
      write (*,*)' 1,2 :'
      write (*,*)'The daily insolation is computed for a given latitude'
      write (*,*)'of the Earth, and for a given position of the Earth '
      write (*,*)'on its orbit. This position is sometime given by a '
      write (*,*)'conventional date where the origine is at the moving'
      write (*,*)'equinox with conventional date 21 march'
      write (*,*)'As this terminology is confusing, in the present '
      write (*,*)'program, the position of the Earth on its orbit is'
      write (*,*)'given either by its true longitude (1), or by the'
      write (*,*)'mean longitude (2), both expressed in degrees'
      write (*,*)'with origine at the moving equinox'
      write (*,*)
      write (*,*)'The true longitude is the true position angle of the'
      write (*,*)'Sun-Earth direction with respect to the equinox. '
      write (*,*)'The mean longitude is a fictif angle which is '
      write (*,*)'described by the Earth at constant velocity.'
      write (*,*)'The mean longitude is thus proportional to the time,'
      write (*,*)'and  should be the normal choice in many cases.'
      write (*,*)
      write (*,*)
      write (*,*) ' press RETURN for MORE '
      PAUSE
      write (*,*)' 3 :'
      write (*,*)' The monthly insolation correspond to the mean daily'
      write (*,*)'insolation over 1/12 of a year, that is 30 degrees'
      write (*,*)'of mean longitude, with origine at the equinox.'
      write (*,*)
      write (*,*)' 4 :'
      write (*,*)'The mean annual insolation is computed over a full'
      write (*,*)'revolution of the Earth around the Sun.'
      write (*,*)
      write (*,*)
      GOTO 999
7000  CONTINUE
      end


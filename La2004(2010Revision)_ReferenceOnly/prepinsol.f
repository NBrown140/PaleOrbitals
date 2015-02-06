      subroutine TRANSBIN(datedebut,datefin,nomfich,
     &                    nomascpos, nomascneg, nbel, nbel1,pasbin)
c****************************************************************************
c  (c) Astronomie et Systemes Dynamiques/ institut de Mecanique Celeste (2003)
c  J. Laskar, F. Joutel, M. Gastineau                                         
c  Ce sous programme lit les fichiers sequentiels ASCII et les reecrit en un
c   fichier a acces direct binaire (nbel1 colonnes) , avec une longueur 
c   d'enregistrement correspondant aux donnees (nbel colonnes) pour 100 000 ans.
c   De nbel+2 a nbel1 les donnees sont mises a 0.
c   En sortie : pasbin : contient le pas de temps
c****************************************************************************
      implicit none
      character(50), intent(in) :: nomfich,nomascpos,nomascneg
      real(8), intent(in) :: datedebut,datefin
      integer, intent(in) :: nbel, nbel1
      real(8), intent(out) :: pasbin
      integer nacd, nrecl, ios, nrec
      integer idebut, ifin, ndebut, nfin, pfin, pdebut, i, n, j, ndate
      parameter(nacd=100)
      real(8) Taux(nbel1,nacd)
      character(20) statut
      
      Taux(nbel+1:nbel1,:) = 0

      nrecl=8*nbel1*nacd

c
c
c Determination des donnees a lire :
      idebut=int(datedebut)*1000.d0-nacd
      ifin=int(datefin)*1000.d0+nacd

      if (idebut.lt.0.and.ifin.lt.0) then
         ndebut=idebut/nacd
         nfin=ifin/nacd-1
         pdebut=0
         pfin=-1
         ndate = (nfin-ndebut)+2
      endif
      if (idebut.gt.0.and.ifin.gt.0) then
         ndebut=1
         nfin=0
         pdebut=idebut/nacd
         pfin=ifin/nacd-1
         ndate =(pfin-pdebut)+2
      endif
      if (idebut.lt.0.and.ifin.gt.0) then
         ndebut=idebut/nacd
         nfin=-1
         pdebut=0
         pfin=ifin/nacd-1
         ndate = (pfin-ndebut)+2
      endif

c
c----------------------------------------------------------------------
c verifie si le fichier binaire existe et est compatible avec les dates
c----------------------------------------------------------------------
c
      open(10,file=nomfich,status='old',access='direct',
     *recl=nrecl,err=101,iostat=ios)

c lecture du dernier enregistrement 
c qui contient la date de debut et la date de fin et du pas
      read(10,rec=ndate, err=103) (Taux(j,1),j=1,3)
      pasbin = Taux(3,1)

      if ((Taux(1,1).ge.datedebut).and.(Taux(2,1).le.datefin)) then
c     fichier compatible => reecriture non necessaire
       close(10)
       goto 102
      endif

103   close(10)

c
c----------------------------------------------------------------------
c Lecture et reecriture necessaire
c----------------------------------------------------------------------
c

c
c Les fichiers ASCII:
c--------------------
101   open(11,file=nomascneg,status='old')
      open(12,file=nomascpos,status='old')
c
c Le fichier a acces direct:
c---------------------------
c      write(*,*) 'nrecl=', nrecl
      open(10,file=nomfich,status='unknown',access='direct',
     *recl=nrecl,err=99,iostat=ios)
c
c Pour les temps negatifs:
c-------------------------
      read(11,*)
      do i=-1,ifin,-1
         read(11,*)
      enddo
      do n=nfin,ndebut,-1
         do i=nacd,1,-1
            read(11,*) (Taux(j,i),j=1,nbel)
            Taux(nbel+1,i)=sin(Taux(nbel,i))
            Taux(nbel,i)=cos(Taux(nbel,i))
          enddo
          pasbin = dabs(Taux(1,2)-Taux(1,1))*1D3
          nrec=(n-ndebut)+1          
          write(10,rec=nrec)((Taux(j,i),j=1,nbel1),i=1,nacd)
      enddo
c
c Pour les temps positifs:
c-------------------------
      do i=0,idebut-1
         read(12,*)
      enddo
      do n=pdebut,pfin
         do i=1,nacd
            read(12,*) (Taux(j,i),j=1,nbel)
            Taux(nbel+1,i)=sin(Taux(nbel,i))
            Taux(nbel,i)=cos(Taux(nbel,i))
         enddo
         pasbin = dabs(Taux(1,2)-Taux(1,1))*1D3
         nrec=(n-ndebut)+1
         if (idebut.gt.0) nrec=n-pdebut+1
         write(10,rec=nrec)((Taux(j,i),j=1,nbel1),i=1,nacd)
      enddo

c ecriture du dernier enregistrement 
c qui contient la date de debut et la date de fin et du pas
c-------------------------
      Taux(1,1) = datedebut
      Taux(2,1) = datefin
      Taux(3,1) = pasbin
      write(10,rec=ndate) (Taux(j,1),j=1,3)

c
c
      goto 100
 99   write(*,*) 'error ',ios,' in the file  :',nomfich
      stop
100   write(*,*) 'The file  ', nomfich, ' is created.'
      close(10)
      close(11)
      close(12)
102   return
      end

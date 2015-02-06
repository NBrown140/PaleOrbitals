      SUBROUTINE Telor(tps,Tel)
c********************************************************************
c
c           LES ELEMENTS ORBITAUX k,h,q,p
c                 POUR UNE DATE DONNEE tps.
c
c A L'ENTREE:
c  tps:la date pour laquelle on desire les elements (en annees)
c
c EN SORTIE:
c  Tel(nbf):tableau des elements orbitaux(k,h,q,p) pour la date tps
c 
c  (c) ASD/IMC 26/06/2001
c   (vf)
c  05/11/2003 Mickael : changement du common date en ddate
c********************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (ni=8,nbf1=4,nbel1=5)
      parameter (nbf=4,nbel=5,nacd=100)
      dimension ta(ni),ya(nbf,ni),y(nbf)
      dimension Tel(nbf),Taux(2,nbel1,nacd)
      data ittprev/-2000000/
      data ttprev/-2.D9/
      data nreccour1/-1/
      data nreccour2/-1/
      save ittprev,ncour1,ncour2,ttprev
      save Taux,nreccour1,nreccour2
      save ta,ya,y
      common/ddate/pasout,itdebut
c
c--------------------------------------------------------------------
c Numero d'ordre de l'enregistrement qui contient tps
c--------------------------------------------------------------------
c
      tt=tps/pasout
      if (tt.eq.ttprev) then
       go to 201
      else
       ttprev=tt
      endif
c
c itt est la partie entiere de tt
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      itt=int(tt)
      if ((tt.lt.0.D0).and.(tt.lt.itt)) then
       itt=itt-1
      endif
      if (itt.eq.ittprev) then
        go to 200
      endif
      ittprev=itt
      nrec=(itt-itdebut)/nacd+1
      itdeb=itdebut+(nrec-1)*nacd
      iplace=itt-itdeb+1
      if (iplace.lt.ni/2) then
        idebord=-1
      else
        if (iplace.gt.nacd-ni/2) then
          idebord=+1
        else
          idebord=0
        endif
      endif
c
c  Pas de debordement
c  ^^^^^^^^^^^^^^^^^^^
      if (idebord.eq.0) then
        m=ni
        nl1=iplace-ni/2+1
        if (nrec.eq.nreccour1) then
          n1=ncour1
        else
          if (nrec.eq.nreccour2) then
            n1=ncour2
          else
            n1=1
            ncour1=1
            nreccour1=nrec
            call Listableau(n1,nrec,Taux)
          endif
        endif
      else
        nrecp=nrec+idebord
        if (idebord.lt.0) then
          nl1=nacd+(iplace-ni/2)+1
          m=ni/2-iplace
        else
          nl1=iplace-ni/2+1
          m=nacd-nl1+1
        endif
        nl2=1
c
c   Chargement eventuel des tableaux en m‚moire
c   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (nreccour1.eq.nrec) then
          if (nreccour2.ne.nrecp) then
            call calncour(ncour1,ncour2)
            nreccour2=nrecp
            call Listableau(ncour2,nrecp,Taux)
          endif
        else
          if (nreccour2.eq.nrec) then
            if (nreccour1.ne.nrecp) then
              call calncour(ncour2,ncour1)
              nreccour1=nrecp
              call Listableau(ncour1,nrecp,Taux)
             endif
          else
            if (nreccour1.eq.nrecp) then
              call calncour(ncour1,ncour2)
              nreccour2=nrec
              call Listableau(ncour2,nrec,Taux)
            else
              if (nreccour2.eq.nrecp) then
                call calncour(ncour2,ncour1)
                nreccour1=nrec
                call Listableau(ncour1,nrec,Taux)
              else
                ncour1=1
                ncour2=2
                nreccour1=nrec
                nreccour2=nrecp
                call Listableau(ncour1,nrec,Taux)
                call Listableau(ncour2,nrecp,Taux)
              endif
            endif
          endif
        endif    
c
c   recherche du premier tableau 
c   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (nreccour1.eq.nrec) then
          if (idebord.lt.0) then
            n1=ncour2
            n2=ncour1
          else
            n1=ncour1
            n2=ncour2
          endif
         else
           if (idebord.lt.0) then
             n1=ncour1
             n2=ncour2
           else
             n1=ncour2
             n2=ncour1
            endif
         endif
       endif
c
c  chargement du tableau d'interpolation
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
100    call Remplistableau(n1,nl1,m,n2,nl2,taux,ta,ya)
        
c
c Quand les tableaux necessaires sont en memoire:
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 200  call pint(ta,ya,ni,tps,nbf,y)
 201  do i=1,nbf
         Tel(i)=y(i)
      enddo
      return
      end

        subroutine Listableau(n,nrec,Taux)
        implicit double precision (a-h,p-y)
        parameter (nacd=100,nbel1=5)
        dimension Taux(2,nbel1,nacd)
c
        read(11,rec=nrec) ((Taux(n,i,j),i=1,nbel1),j=1,nacd)
        return
        end

        subroutine Remplistableau(n1,nl1,m,n2,nl2,Taux,ta,ya)
c
c  Remplit les tableaux ta et ya a l'aide de Taux(n1, , ) sur une longueur m
c  a partir de l'indice nl1 puis de Taux(n2, , ) sur une longueur ni-m a partir
c  de l'indice nl2
c  
        implicit double precision (a-h,o-z)
        parameter (ni=8,nbel1=5,nbf=4,nacd=100)
        dimension ta(ni),ya(nbf,ni)
        dimension Taux(2,nbel1,nacd)
c
        nl1p=nl1-1
        nl2p=nl2-1
        do i=1,m
          ta(i)=Taux(n1,1,nl1p+i)*1.D3  
          do j=1,nbel1-1
            ya(j,i)=Taux(n1,j+1,nl1p+i)
          end do
        end do
        if (m.lt.ni) then
          do i=1,ni-m
            ta(m+i)=Taux(n2,1,nl2p+i)*1.D3
            do j=1,nbel1-1
              ya(j,m+i)=Taux(n2,j+1,nl2p+i)
            end do
          end do
        endif
        return
        end

        subroutine calncour(n,m)
          implicit double precision (a-h,p-y)
          if (n.eq.1) then
            m=2
          else
            m=1
          endif
          return
          end


      SUBROUTINE PINT(xa,ya,n,x,nbf,y)
c     ********************************
c------------------------------------------------------------------------
c
c A partir de fonctions tabulees (xa,ya),ce sous-programme construit
c un polynome de degre (n-1) permettant de calculer la valeur de la
c fonction, y, au point, x.
c (cf algorithme de Neville, Numerical recipes)
c
c A L'ENTREE:
c  xa(n):tableau de dates
c  ya(nbf,n):valeurs des fonctions aux dates xa(n)
c  n:nombre de points necessaires a l'interpolation
c  x:date pour laquelle on effectue l'interpolation
c  nbf:nombre de fonctions a interpoler
c
c EN SORTIE:
c  y(nbf):tableau des valeurs interpolees a la date x
c  
c  (c) ASD/IMC (2001)
c------------------------------------------------------------------------
c
c Nombre de fonctions interpolees
c
      implicit double precision (a-h,o-z) 
      parameter (nmax=100,nbfmax=20)
      dimension xa(n),ya(nbf,n),c(nmax),d(nmax)
      dimension y(nbf),dy(nbfmax)
c
      do 100 nvar=1,nbf
         ns=1
         dif=dabs(x-xa(1))
         do  10 i=1,n
             dift=dabs(x-xa(i))
             if (dift.lt.dif) then
                 ns=i
                 dif=dift
             endif
             c(i)=ya(nvar,i)
             d(i)=ya(nvar,i)
 10      continue
         y(nvar)=ya(nvar,ns)
         ns=ns-1
         do 20 m=1,n-1
            do 30 i=1,n-m
               ho=xa(i)-x
               hp=xa(i+m)-x
               w=c(i+1)-d(i)
               den=ho-hp
               if (den.eq.0) then
                  write(*,*)' pause pint'
                  pause
               endif
               den=w/den
               d(i)=hp*den
               c(i)=ho*den
 30         continue
            if (2*ns.lt.n-m) then
               dy(nvar)=c(ns+1)
            else
               dy(nvar)=d(ns)
               ns=ns-1
            endif
            y(nvar)=y(nvar)+dy(nvar)
 20      continue
 100  continue
      return
      end

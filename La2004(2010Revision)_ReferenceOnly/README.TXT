===========================                         ============================
                         ASTRONOMIE ET SYSTEMES DYNAMIQUES
                         
                         INSTITUT DE MECANIQUE CELESTE

                                    la2004

                                2010, January 18 
===========================                         ============================



The program insola is derived from the same program from
the Earth paleoclimate solution La93 (Laskar, Joutel, Boudin, 1993) with 
several differences.

The utilisation has been greatly simplified, as everything is now 
in a single application, insola, that allows to compute all 
insolation quantities for any time interval between the limits 
of the input files defined in insola.par (nomclimapos, nomclimaneg).

Once the application is generated, the only file to edit is 
insola.par.


We have provided solution La2004 from -50 to + 21 Myr.




La2004    : 
==========

Main reference :
Laskar, J., Gastineau, M., Joutel, F., Robutel, P., Levrard, B., Correia, A.,: 2004,
A long term numerical solution for the insolation quantities of Earth.
{\it in preparation}


********************************************************************************
*  Authors: J. Laskar, M. Gastineau, F. Joutel                                 *
*  (c) Astronomie et Systemes Dynamiques, Institut de Mecanique Celeste,       *
*      Paris (2004)                                                            *
*                                                                              *
*                               Jacques Laskar                                 *
*                               Astronomie et Systemes Dynamiques,             *
*                               Institut de Mecanique Celeste                  *
*                               77 av. Denfert-Rochereau                       *
*                               75014 Paris                                    *
*                       email:  laskar@imcce.fr                                *
*                                                                              *
********************************************************************************




CONTAINS:
=========



File Summary:
--------------------------------------------------------------------------------
 FileName                 Explanations
--------------------------------------------------------------------------------
README                    This file
INSOLN.LA2004.BTL.ASC    Nominal solution La2004,
                                    before present years (-51Myr to 0)
INSOLP.LA2004.BTL.ASC    Nominal solution La2004
                                    after present years (0 to +21Myr)
INSOLN.LA2004.BTL.100.ASC Nominal solution La2004,
                                    before present years (-101Myr to 0)
INSOLN.LA2004.BTL.250.ASC Nominal solution La2004,
                                    before present years (-249Myr to 0)


Makefile                  to compile and generate derived files
insola.par                Parameter file for insola.
insola.f                  FORTRAN program for computations of various
                                    insolation quantities
insolsub.f                subroutines for insola.f
prepsub.f                 FORTRAN subroutine for insola.f
prepinsol.f               FORTRAN subroutine for preparation step for insola.f
--------------------------------------------------------------------------------


Byte-per-byte Description of file: INSOLN.LA2004.BTL.ASC
Byte-per-byte Description of file: INSOLP.LA2004.BTL.ASC
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1-14   F13.3  1000yr  t        Time from J2000  in 1000 years
  18-39   D25.16 ---     e        eccentricity
  43-64   D25.16 rad     eps      obliquity (radians)
  68-89   D25.16 rad     pibar    longitude of perihelion from moving equinox
                                  (radians)
--------------------------------------------------------------------------------

Byte-per-byte Description of file: INSOLN.LA2004.BTL.100.ASC
Byte-per-byte Description of file: INSOLN.LA2004.BTL.250.ASC
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1-9    F8.0   1000yr  t        Time from J2000  in 1000 years
  10-18   D9.6   ---     e        eccentricity
  19-27   D9.6   rad     eps      obliquity (radians)
  28-37   D10.6  rad     pibar    longitude of perihelion from moving equinox
                                  (radians)
--------------------------------------------------------------------------------

Description:

  La2004 is the nominal solution for precessional quantities and orbital quantities of 
  the Earth. The solution La2004 is provided with fortran subroutine in order to 
  compute the insolation quantities of Earth.
  


BIBLIOGRAPHY:
============



  La90: Laskar, J.: 1990, The chaotic motion of the solar system.
                  A numerical estimate of the chaotic zones
                  Icarus, 88, 266

  La93: Laskar, J., Joutel, F., Boudin, F.: 1993, Orbital, precessional
                  and insolation quantities for the Earth
                  from -20 Myr to + 10Myr
                  Astron. Astrophys. 270, 522

  La2004 :Laskar, J., Gastineau, M., Joutel, F., Robutel, P., Levrard, B., 
                  Correia, A.,
                  A long term numerical solution for the insolation quantities 
                  of Earth. {\it in preparation}



================================================================================


================================================================================
Compilation of the Program
--------------------------

  Edit the makefile to choose a fortran-90 compiler and run "make".




Detailed Description of the Program
-----------------------------------

  insola.par : THE USER should edit the file insola.par to set
          all parameters.

  insola : interactive insolation computations


--------------------------------------------------------------------------------
INSOLA.PAR

The user should edit the parameter file insola.par.

 datefin   : ending   time (Myr)
 datedebut : starting time (Myr)
 pas       : sampling step (in years). 
             The default value is 1000. The minimum step is 1 year
             (it is usually not necessary for pas to be smaller than 100 (yr)). 
 nominsolbin   : Binary intermediate file for insolation
 nomclimapos : ASCII file for climatic variables for positive time  
 nomclimaneg : ASCII file for climatic variables for negative time
 so          : Solar constant expressed in W/m2 
               solar constant at 1 AU so0 =1368 W/m2

 

Starting and ending date must verify :   -50 <= datedebut <= 0 <= datefin <= +20 
The user may restrict himself to e.g. -2Myr to 0 by putting instead
     datedebut=-2.D0
     datefin=0.D0



--------------------------------------------------------------------------------
INSOLA

  The routines insola are designed to compute all
  necessary insolation quantities derived from the orbital and
  precessional quantities computed above.

  They are given on the form of FORTRAN code, so the user can check
  if they correspond to his needs. He can also design his own insolation
  routines.

  Insola will construct a binary intermediate file 
  when this file is not present before computations.

  insola is a self documented program
  For more details, the user can refer to  (Laskar et al, 1993,2004)

================================================================================
History
-------
 2010/01/15 : fix bug in the function vraimoy and truncate the developpement at the order 5

================================================================================
User feed-back is encouraged. Unless otherwise specified, send comments and bug 
reports to:                    E-mail     : laskar@imcce.fr
                               Fax        : (33) 1 40 51 20 55
                               Postal mail: Institut de Mecanique Celeste
                                            77 avenue Denfert Rochereau
                                            F-75014 PARIS
================================================================================

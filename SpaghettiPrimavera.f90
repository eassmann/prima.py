!! SpaghettiPrimavera.f90 Version 20-03-2008
!!                                Wien2k
!!
!! Author  Maurits W. Haverkort
!!         Written at Max Planck Institute Stuttgart
!!
!!         (C) 2008

!! Modified for inclusion in Python (prima.py):
!!         Elias Assmann <elias.assmann@gmail.com> (2013-2015)

!! prima.py version 0.3
!!
!! $Id: SpaghettiPrimavera.f90 618 2015-09-14 15:58:11Z assmann $

!! Spaghetti Primavera is a program that creates a band-structure plot (eps)
!! given a list of k-points and a list of eigen-energies and eigen-functions
!! or characters.
!!
!! Sperate functions can specify the line thickness and color of the bands
!! as a function of the character.
!!
!! The Wien2k version expects a klist file and a qtl file

module sppv_data
  character(512) :: version

  integer, parameter :: MaxAtoms=1000, MaxOrbitals=50
  character(len=512), parameter :: progname = "sppv.so"

  character(512) :: qtlName, bandName, cmdLine, fontName="Times-Roman"
  integer :: MinorTicks=4, nkpoints
  real(8) :: BaseThickness=.1, AxesThickness=.5
  real(8) :: Xsize=500, Ysize=700, MajorTicks=1.0, TextSize=12
  real(8) :: Emin=-5, Emax=5, Efermi=0
  real(8), allocatable :: legcolors(:,:) ! (3, #legend)
  real(8), allocatable :: orbcolors(:,:) ! (3, #selected)
  real(8), allocatable :: orbthicks(:)   ! (#selected)
  integer, allocatable :: atomsorbs(:,:) ! (2, #selected)

  character, allocatable :: legend(:)
  logical :: writeLegend=.false., writeKLines=.true., writeKLabels=.true.
  logical :: writeFermiLine=.true.
  integer, allocatable :: legentries(:)
end module sppv_data

module sppv_formats
  ! f2py does not like to have ‘len=*’ here
  character(len=256), parameter :: &
       fmt_line  = "(2(F0.6,1X), 'lineto')",       &
       fmt_move  = "(2(F0.6,1X), 'moveto')",       &
       fmt_width = "(F0.3,1X,    'setlinewidth')", &
       fmt_color = "(3(F0.4,1X), 'setrgbcolor')",  &
       fmt_scale = "(F0.6,1X,    'scalefont')",    &
       fmt_bgcol = "(3(F0.4,1X), 'setrgbcolor clippath fill % added Sep 2014&
       & to make current ‘evince’ versions display white bg')"
end module sppv_formats

module sppv_colors
  real(8), parameter, dimension(3) :: &
       WHITE=(/1,1,1/), BLACK=(/0,0,0/), &
       RED=(/1,0,0/), GREEN=(/1,0,0/), BLUE=(/1,0,0/)

  real(8), parameter, dimension(3) :: &
       col_bg=WHITE, col_axes=BLACK, col_ticks=BLACK, col_leg=BLACK, &
       col_err=RED, col_klab=BLACK
end module sppv_colors

subroutine SpaghettiPrimavera
  use sppv_data

  implicit none

  call WritePSHead(Xsize, Ysize, version, cmdLine, qtlName)
  call WritePSKpoints(bandName, Nkpoints, Xsize, Ysize, AxesThickness, &
       &              TextSize, fontName, writeKLines, writeKLabels)
  if (writeLegend) then
     call WritePSLegend(TextSize, fontName, XSize, YSize, &
          &             legcolors, legend, legentries)
  end if
  call WritePSAxesA(Xsize, Ysize, Emin, Emax, AxesThickness, writeFermiLine)
  call WritePSBands(qtlName, nkpoints, Xsize, Ysize, &
       &            Efermi, TextSize, Emin, Emax, &
       &            MaxAtoms, MaxOrbitals, fontName)
  call WritePSAxesB(Xsize, Ysize, Emin, Emax, MajorTicks, MinorTicks, &
       &            AxesThickness, TextSize, fontName)
  call WritePSTail()

contains

!!! Excerpted from util_w2w.F
subroutine croak(message, status)
  use iso_fortran_env, only: ERROR_UNIT
  use sppv_data,       only: Ysize, progname, fontName, TextSize
  use sppv_colors,     only: col_err
  use sppv_formats,    only: fmt_color, fmt_move, fmt_scale

  character(len=*), intent(in), optional   :: message
  integer,          intent(in), optional   :: status

  integer            :: s
  s=1

  if (present(status)) s=status

  if (present(message)) then
     write(ERROR_UNIT, '(A, ": ", A)') trim(progname), trim(message)
     write(6,'(A)')"newpath"
     write(6,'(A)')"/" // trim(fontName) // " findfont"
     write(6,fmt_scale) TextSize
     write(6,'(A)')"setfont"
     write(6,fmt_color) col_err
     write(6,fmt_move)50.,Ysize/2+50
     write(6,'("(",A,") show")') message
  end if

  call WritePSTail()
  close(6)
  call exit(s)
end subroutine croak

subroutine carp(message)
  use iso_fortran_env, only: ERROR_UNIT
  use sppv_data,       only: progname

  character(len=*), intent(in) :: message

  write(ERROR_UNIT, '(A, ": ", A)') trim(progname), message
end subroutine carp

subroutine cluck(message)
  use iso_fortran_env, only: OUTPUT_UNIT
  use sppv_data,       only: progname

  character(len=*), intent(in) :: message

  write(OUTPUT_UNIT, '(A, ": ", A)') trim(progname), message
end subroutine cluck
!!! End util_w2w.F excerpt

function CharacterToRGBColor(OrbCharacter) result(RGBColor)
  use sppv_data
  real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
  real(8)              :: RGBColor(3)

  integer :: i

  RGBColor = 0

  do i = 1, size(atomsorbs,2)
     RGBColor = RGBColor + &
          OrbCharacter(atomsorbs(1, i), atomsorbs(2, i)) * orbcolors(:,i)
  end do
end function CharacterToRGBColor

function CharacterToThickness(OrbCharacter) result(Thickness)
  use sppv_data
  real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
  real(8)              :: Thickness

  integer :: i

  Thickness = BaseThickness

  do i = 1, size(atomsorbs,2)
     Thickness = Thickness &
          + OrbCharacter(atomsorbs(1, i), atomsorbs(2, i)) * orbthicks(i)
  end do
end function CharacterToThickness

subroutine WritePSHead(Xsize, Ysize, version, cmdLine, qtlName)
  !
  ! writes the header of the postscript to standard output
  !

  use sppv_formats
  use sppv_colors

  implicit none

  real(8),          intent(in) :: Xsize, Ysize
  character(len=*), intent(in) :: version, cmdLine, qtlName

  write(6,'(A)')"%!PS-Adobe-2.0 EPSF-2.0"
  write(6,'(A,2(1X,I0))')"%%BoundingBox: 0 0",int(Xsize+100),int(Ysize+100)
  write(6,'(A,2(1X,F0.6))')"%%HiResBoundingBox: 0.000000 0.000000",Xsize+100,Ysize+100
  write(6,'(A)')"%% Creator: Spaghetti Primavera by Maurits W. Haverkort"
  write(6,'(A)')"%%          and prima.py v. "//trim(version)
  write(6,'(A)')"%% command line: "//trim(cmdLine)
  write(6,'(A)')"%% using file ‘"  //trim(qtlName)//"’"
  write(6,'(A)')"%%EndComments"
  write(6,'(A)')"%%BeginProlog"
  write(6,fmt_bgcol) col_bg
  write(6,'(A)')"save"
  write(6,'(A)')"countdictstack"
  write(6,'(A)')"mark"
  write(6,'(A)')"newpath"
  write(6,'(A)')"/setpagedevice {pop} def"
  write(6,'(A)')"/L{ newpath setlinewidth setrgbcolor moveto lineto&
       & lineto 0 setlinecap 1 setlinejoin stroke } def"
  write(6,'(A)')"/L1{ newpath setlinewidth setrgbcolor moveto&
       & lineto 0 setlinecap 1 setlinejoin stroke } def"
  write(6,'(A)')"%%EndProlog"
  write(6,'(A)')"%%Page 1 1"
end subroutine WritePSHead

subroutine WritePSTail()
  !
  ! Writes the tail of the postscript to standard output
  !
  implicit none
  write(6,'(A)')"%%Trailer"
  write(6,'(A)')"cleartomark"
  write(6,'(A)')"countdictstack"
  write(6,'(A)')"exch sub { end } repeat"
  write(6,'(A)')"restore"
  write(6,'(A)')"%%EOF"
end subroutine WritePSTail

subroutine WritePSAxesA(Xsize, Ysize, Emin, Emax, AxesThickness, writeFermiLine)
  !
  ! Writes the part of the axes system to the standard output
  ! that should be behind the graph 
  !
  use sppv_formats
  use sppv_colors

  implicit none
  real(8), intent(in) :: Xsize, Ysize, Emin, Emax, AxesThickness
  logical, intent(in) :: writeFermiLine

  real(8) :: PtEfermi,PtYunit
  ! A line at the fermi energy
  PtEfermi=50-(Ysize*Emin)/(Emax-Emin)
  PtYunit=(Ysize)/(Emax-Emin)
  if (writeFermiLine .and. Emin<0 .and. Emax>0) then
     write(6,'(A)')"newpath"
     write(6,fmt_width) AxesThickness
     write(6,fmt_color) col_axes
     write(6,fmt_move) 50.,PtEfermi
     write(6,fmt_line) 50+Xsize,PtEfermi
     write(6,'(A)')"stroke"
  endif
end subroutine WritePSAxesA

subroutine WritePSAxesB(Xsize, Ysize, Emin, Emax, MajorTicks, MinorTicks, &
     &                  AxesThickness, TextSize, fontName)
  !
  ! Writes the part of the axes system that should be on top of the graph
  !

  use sppv_formats
  use sppv_colors

  implicit none

  integer,          intent(in) :: MinorTicks
  real(8),          intent(in) :: Xsize, Ysize, Emin, Emax
  real(8),          intent(in) :: MajorTicks, TextSize, AxesThickness
  character(len=*), intent(in) :: fontName

  real(8) :: PtEfermi,PtYunit,PtTick

  character(len=*), parameter :: fmt_rj = "(A,F5.1,A)"


  ! some positions in pts
  PtEfermi=50-(Ysize*Emin)/(Emax-Emin)
  PtYunit=(Ysize)/(Emax-Emin)

  ! write the ticks
  if((Emin+MajorTicks).LT.Emax) then
     write(6,'(A)')"newpath"
     write(6,fmt_width) AxesThickness
     write(6,fmt_color) col_ticks
     PtTick=PtEfermi
     do while (PtTick.LT.(Ysize+50))
        write(6,fmt_move)50.,PtTick
        write(6,fmt_line)50+Xsize/100.0,PtTick
        write(6,fmt_move)50+Xsize,PtTick
        write(6,fmt_line)50+Xsize-Xsize/100.0,PtTick
        PtTick=PtTick+MajorTicks*PtYunit
     enddo
     PtTick=PtEfermi
     PtTick=PtTick-MajorTicks*PtYunit
     do while (PtTick.GE.(50))
        write(6,fmt_move)50.,PtTick
        write(6,fmt_line)50+Xsize/100.0,PtTick
        write(6,fmt_move)50+Xsize,PtTick
        write(6,fmt_line)50+Xsize-Xsize/100.0,PtTick
        PtTick=PtTick-MajorTicks*PtYunit
     enddo
     write(6,'(A)')"stroke"
     write(6,'(A)')"newpath"
     write(6,fmt_width) AxesThickness
     write(6,fmt_color) col_axes
     PtTick=PtEfermi
     do while (PtTick.LT.(Ysize+50))
        write(6,fmt_move)50.,PtTick
        write(6,fmt_line)50+Xsize/150.0,PtTick
        write(6,fmt_move)50+Xsize,PtTick
        write(6,fmt_line)50+Xsize-Xsize/150.0,PtTick
        PtTick=PtTick+MajorTicks*PtYunit/(1.0+MinorTicks)
     enddo
     PtTick=PtEfermi
     PtTick=PtTick-MajorTicks*PtYunit/(1.0+MinorTicks)
     do while (PtTick.GE.(50))
        write(6,fmt_move)50.,PtTick
        write(6,fmt_line)50+Xsize/150.0,PtTick
        write(6,fmt_move)50+Xsize,PtTick
        write(6,fmt_line)50+Xsize-Xsize/150.0,PtTick
        PtTick=PtTick-MajorTicks*PtYunit/(1.0+MinorTicks)
     enddo
     write(6,'(A)')"stroke"
  endif

  ! definition for right-justifying tick numbers
  ! from http://www.cs.utsa.edu/~wagner/CS3723/postrecs/justify/just.html
  write(6, '(A)') "% procedure to right justify against x"
  write(6, '(A)') "/rj { % stack: string "
  write(6, '(A)') "  dup stringwidth pop % stack: str width"
!  write(6, '(A)') "  x y moveto % move to right end"
  write(6, '(A)') "  neg 0 rmoveto % move to left end"
  write(6, '(A)') "  show % print string"
!  write(6, '(A)') "  next % increment y"
  write(6, '(A)') "} def"

  ! write the numbers next to the ticks
  write(6,'(A)')"/" // trim(fontName) // " findfont"
  write(6,fmt_scale) TextSize
  write(6,'(A)')"setfont"
  write(6,fmt_color) col_axes
  PtTick=PtEfermi
  do while (PtTick.LT.(Ysize+50))
     write(6,fmt_move)45., PtTick-TextSize/4.0
     write(6,fmt_rj)"(",(PtTick-PtEfermi)/PtYunit,") rj"
     PtTick=PtTick+MajorTicks*PtYunit
  enddo
  PtTick=PtEfermi
  PtTick=PtTick-MajorTicks*PtYunit
  do while (PtTick.GE.(50))
     write(6,fmt_move)45., PtTick-TextSize/4.0
     write(6,fmt_rj)"(",(PtTick-PtEfermi)/PtYunit,") rj"
     PtTick=PtTick-MajorTicks*PtYunit
  enddo
  ! write a frame around the graph
  write(6,'(A)')"newpath"
  write(6,fmt_width) AxesThickness
  write(6,'(A)')"0 setlinejoin"
  write(6,fmt_color) col_axes
  write(6,fmt_move) 50., 50.
  write(6,fmt_line)Xsize+50,50.
  write(6,fmt_line)Xsize+50,Ysize+50
  write(6,fmt_line)50.,Ysize+50
  write(6,'(A)')"closepath"
  write(6,'(A)')"stroke"
end subroutine WritePSAxesB

subroutine WritePSKpoints(bandName, NKpoints, Xsize, Ysize, AxesThickness, &
     &                    TextSize, fontName, writeKLines, writeKLabels)
  !
  ! reads the file bandName and 
  ! makes vertical lines for each k-point that has a label in the 
  ! file band-name. These lines are labeled
  ! returns the number of points that should be between each k-point
  ! and the number of k-points
  !
  use sppv_formats
  use sppv_colors

  implicit none

  character(len=*), intent(in)  :: bandName, fontName
  real(8),          intent(in)  :: Xsize, Ysize, TextSize, AxesThickness
  integer,          intent(out) :: nkpoints
  logical,          intent(in)  :: writeKLines, writeKLabels

  integer   :: OK,nlabels,labelpoint(100),ilabel
  character :: RDWD*(128),label(100)*(2)

  open(unit=1, file=bandName, action='read', iostat=OK)
  if(OK.NE.0) call croak("Could not open klist file"//trim(bandName))

  ! read file "bandName" and filter the number of k-points, 
  ! the labels and their position
  nkpoints=0
  nlabels=0
  read(unit=1,fmt='(A)',end=201) RDWD
  do while (RDWD.NE."END")
     nkpoints=nkpoints+1
     if(RDWD(1:1).NE." ") then
        nlabels=nlabels+1
        label(nlabels)=RDWD(1:2)
        labelpoint(nlabels)=nkpoints
     endif
     read(unit=1,fmt='(A)',end=201) RDWD
  enddo

  if (writeKLines) then
     ! plot lines for specific k-points
     do ilabel = 1, nlabels
        write(6,'(A)')"newpath"
        write(6,fmt_width) AxesThickness
        write(6,fmt_color) col_klab
        write(6,fmt_move) 50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0),&
             &            50.
        write(6,fmt_line) 50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0),&
             &            50+Ysize
        write(6,'(A)')"stroke"
     enddo
  end if

  if (writeKLabels) then
     ! write the labels for specific k-points
     do ilabel = 1, nlabels
        if (label(ilabel)(2:2).EQ." ") then
           write(6,'(A)')"/" // trim(fontName) // " findfont"
        else
           write(6,'(A)')"/Symbol findfont"
        endif
        write(6,fmt_scale)TextSize
        write(6,'(A)')"setfont"
        write(6,fmt_color) col_klab
        write(6,fmt_move) 50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0)-&
             &            TextSize/2.0,50-TextSize
        write(6,'(A)')"(",label(ilabel)(1:1),") show"
     enddo
  end if

  close(1)
  return

201 call croak("Unexpected end of file while reading qtl file "//&
         trim(bandName))
end subroutine WritePSKpoints

subroutine WritePSLegend(TextSize, fontName, XSize, YSize, &
     &                   legcolors, legend, legentries)
  implicit none

  real(8),          intent(in) :: TextSize, XSize, YSize
  real(8),          intent(in) :: legcolors(:,:)
  character(len=*), intent(in) :: fontName
  character(len=*), intent(in) :: legend(:)
  integer,          intent(in) :: legentries(:)

  integer       :: i, s=0
  character(64) :: l
  real(8)       :: XPos=48, sep

  sep=.66*TextSize + 1.5*TextSize

  do i = 1, size(legentries)
     write(l, '(64A)'), legend(s+1 : legentries(i))
     call WritePSLegE(TextSize, fontName, &
          &           XPos, YSize+70, &
          &           legcolors(:, i), l)

     XPos = XPos + sep + .66*len_trim(l)*TextSize

     s = legentries(i)
  end do
end subroutine WritePSLegend

subroutine WritePSLegE(TextSize, fontName, XPos, YPos, Color, Name)
  use sppv_formats
  use sppv_colors

  implicit none

  real(8), intent(in) :: TextSize, XPos, YPos, Color(3)
  
  character(len=*) :: Name, fontName
  real(8) :: wd, ht, sep

  wd=.66*TextSize
  ht=wd
  sep=.5*TextSize

  write(6,'(A)')"/" // trim(fontName) // " findfont"
  write(6,fmt_scale)TextSize
  write(6,'(A)')"setfont"
  write(6,fmt_move) XPos+sep, YPos
  write(6,fmt_color) col_leg
  write(6,'(A)')"(", trim(Name) ,") show"
  write(6,'(A)')"newpath"
  write(6,fmt_move)XPos,    YPos+ht
  write(6,fmt_line)XPos+wd, YPos+ht
  write(6,fmt_line)XPos+wd, YPos
  write(6,fmt_line)XPos,    YPos
  write(6,'(A)')" closepath"
  write(6,fmt_color) Color
  write(6,'(A)')" fill"
  write(6,'(A)')" stroke"
end subroutine WritePSLegE

subroutine WritePSBands(qtlName, nkpoints, Xsize, Ysize, &
     &                  Efermi, TextSize, Emin, Emax, &
     &                  MaxAtoms, MaxOrbitals, fontName)

  use sppv_formats
  use sppv_colors

  implicit none

  character(len=*), intent(in) :: qtlName, fontName
  integer,          intent(in) :: nkpoints, MaxAtoms, MaxOrbitals
  real(8),          intent(in) :: Xsize, Ysize, Emin, Emax, Efermi, TextSize

  real(8), parameter :: RydToEv =13.6056981

  real(8)   :: OrbCharacter(MaxAtoms,MaxOrbitals), RGBColor(3),Thickness
  real(8)   :: Energy,Eprev,Enow,Enext, ptsPerE,ptsPerK
  integer   :: OK,nAtom,nOrb(MaxAtoms),i,j,k,junki,iband
  character :: RDWD*(256)

  character(len=*), parameter :: &
       fmt_L = '(10(F0.4,1X), "L")', fmt_L1 = '(8(F0.4,1X), "L1")'

  ptsPerE=Ysize/(Emax-Emin)
  ptsPerK=Xsize/(nkpoints-1)
  ! Open the qtl file
  open(unit=2, file=qtlName, action='read', iostat=OK)
  if(OK.NE.0) call croak("Could not open qtl file"//trim(qtlName))

  nOrb=0
  OrbCharacter=0

  ! find out how many atoms and orbitals there are
  nAtom=0
  read(unit=2,fmt='(A)',end=201) RDWD
  ! newer Wien2k versions seem to omit the colon in “BAND:”
  do while (index(RDWD,"BAND").EQ.0)
     if(index(RDWD,"JATOM").NE.0) then
        nAtom=nAtom+1

        if (nAtom > MaxAtoms) then
           write(0, '("sppv: maximum number of atoms (", I0, ") exceeded")')&
                MaxAtoms
           call exit(1)
        end if

        i=1
        do
           nOrb(nAtom)=nOrb(nAtom)+1
           if (nOrb(nAtom) > MaxOrbitals) then
              write(0, '("sppv: maximum number of orbitals (", I0, ") &
                   & exceeded for atom ", I0)') MaxOrbitals, nAtom
              call exit(1)
           end if

           j = index(RDWD(i:), ",")
           i = i + j
           if (j==0 .or. RDWD(i:)=='') then
              exit
           end if
        enddo
     endif
     read(unit=2,fmt='(A)',end=201) RDWD
  enddo
  ! add the intersitial character
  nAtom=nAtom+1
  nOrb(nAtom)=1
  ! read nkpoints k points then read the string BAND: XX if the 
  ! string BAND: XX can not be read while EOF has been reached stop
  iband=0
  do
     iband=iband+1
     do i = 1, nkpoints
        do j = 1, nAtom
           read(unit=2,fmt=*,end=201) Energy,junki,&
                &                     (OrbCharacter(j,k),k=1,nOrb(j))
        enddo
        Energy=(Energy-Efermi)*RydToEv

        Thickness = CharacterToThickness(OrbCharacter)
        RGBColor  = CharacterToRGBColor (OrbCharacter)

        if(i.EQ.1)then
           Eprev=Energy
        elseif(i.EQ.2)then
           Enow=Energy
           if((Eprev.GT.Emin).AND.(Eprev.LT.Emax).AND.&
                &    (((Eprev+Enow)/2.0).GT.Emin).AND.(((Eprev+Enow)/2.0).LT.Emax)) &
                &  then
              write(6,fmt_L1) &
                   50+(i-2.0)*ptsPerK, 50+(Eprev-Emin)           *ptsPerE, &
                   50+(i-1.5)*ptsPerK, 50+(0.5*(Enow+Eprev)-Emin)*ptsPerE, &
                   RGBColor,Thickness
           endif
        else
           Enext=Energy
           if((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Eprev).GT.Emin).AND.  &
                &     (0.5*(Enow+Eprev).LT.Emax).AND.  &
                &     (0.5*(Enow+Enext).GT.Emin).AND.  &
                &     (0.5*(Enow+Enext).LT.Emax))then
              write(6,fmt_L) &
                   50+(i-2.5)*ptsPerK,50+(0.5*(Eprev+Enow)-Emin)*ptsPerE, &
                   50+(i-2.0)*ptsPerK,50+(Enow-Emin)            *ptsPerE, &
                   50+(i-1.5)*ptsPerK,50+(0.5*(Enext+Enow)-Emin)*ptsPerE, &
                   RGBColor,Thickness
           elseif((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Eprev).GT.Emin).AND. &
                &     (0.5*(Enow+Eprev).LT.Emax))then
              write(6,fmt_L1) &
                   50+(i-2.5)*ptsPerK, 50+(0.5*(Eprev+Enow)-Emin)*ptsPerE, &
                   50+(i-2.0)*ptsPerK, 50+(Enow-Emin)            *ptsPerE, &
                   RGBColor,Thickness
           elseif((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Enext).GT.Emin).AND. &
                &     (0.5*(Enow+Enext).LT.Emax))then
              write(6,fmt_L1) &
                   50+(i-1.5)*ptsPerK, 50+(0.5*(Enext+Enow)-Emin)*ptsPerE, &
                   50+(i-2.0)*ptsPerK, 50+(Enow-Emin)            *ptsPerE, &
                   RGBColor,Thickness
           endif
           if(i.EQ.nkpoints)then
              if((Enext.GT.Emin).AND.(Enext.LT.Emax).AND. &
                   &    (((Enext+Enow)/2.0).GT.Emin).AND.(((Enext+Enow)/2.0).LT.Emax)) &
                   &   then
                 write(6,fmt_L1) &
                      50+(i-1.0)*ptsPerK, 50+(Enext-Emin)           *ptsPerE, &
                      50+(i-1.5)*ptsPerK, 50+(0.5*(Enow+Enext)-Emin)*ptsPerE, &
                      RGBColor, Thickness
              endif
           endif

           Eprev=Enow
           Enow=Enext
        endif
     enddo
     ! read the string BAND: XX or iff EOF finished
     read(unit=2,fmt='(A)',end=300) RDWD
  enddo
300 continue
  close(2)
  return

  ! Goto 201 if the file "qtlName" ended before the end was expected
201 call croak("Unexpected end of file while reading qtl file "//&
         qtlName)
end subroutine WritePSBands
end subroutine SpaghettiPrimavera

!! Time-stamp: <2015-09-14 16:05:19 assman@faepop23.tu-graz.ac.at>

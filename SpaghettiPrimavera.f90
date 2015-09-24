!! SpaghettiPrimavera.f90 Version 20-03-2008
!!                                Wien2k
!!
!! Author  Maurits W. Haverkort
!!         Written at Max Planck Institute Stuttgart
!!
!!         (C) 2008

!! Modified for inclusion in Python (prima.py):
!!         Elias Assmann <elias.assmann@gmail.com> (2013-2014)

!! prima.py version 0.2
!!
!! $Id: SpaghettiPrimavera.f90 354 2014-09-01 20:53:41Z assmann $

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

  character(512) :: qtlName, bandName, cmdLine, fontName="Times-Roman"
  integer :: MinorTicks=4, nkpoints
  real(8) :: Xsize=500, Ysize=700, MajorTicks=1.0, TextSize=12
  real(8) :: Emin=-5, Emax=5, Efermi=0
  real(8), allocatable :: colors(:, :)

  character, allocatable :: legend(:)
  logical :: writeLegend = .false.
  integer, allocatable :: legentries(:)
end module sppv_data

subroutine SpaghettiPrimavera(CharacterToRGBColor, CharacterToThickness)
  use sppv_data
  implicit none

  ! These are callbacks to be passed from Python.  Interestingly, when
  ! I give “intent(callback)” for them explicitly, the program dies of
  ! SIGSEGV.
  external CharacterToRGBColor
  external CharacterToThickness

  interface
     ! It seems CharacterToRGBColor has to be a subroutine to be able
     ! to return RGBColor as an array.  For consistency then,
     ! CharacterToThickness is also a subroutine.
     subroutine CharacterToRGBColor(OrbCharacter, RGBColor)
       use sppv_data
       real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
       real(8), intent(out) :: RGBColor(3)
     end subroutine CharacterToRGBColor

     subroutine CharacterToThickness(OrbCharacter, Thickness)
       use sppv_data
       real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
       real(8), intent(out) :: Thickness
     end subroutine CharacterToThickness
  end interface

  interface
     subroutine WritePSLegend(TextSize, fontName, XSize, YSize, &
          &                   colors, legend, legentries)
       real(8) :: TextSize, XSize, YSize, colors(:,:)
       character(*) :: fontName, legend(:)
       integer :: legentries(:)
     end subroutine WritePSLegend
  end interface

  call WritePSHead(Xsize, Ysize, version, cmdLine)
  call WritePSKpoints(bandName, Nkpoints, Xsize, Ysize, &
       &              TextSize, fontName)
  if (writeLegend) then
     call WritePSLegend(TextSize, fontName, XSize, YSize, &
          &             colors, legend, legentries)
  end if
  call WritePSAxesA(Xsize, Ysize, Emin, Emax)
  call WritePSBands(qtlName, nkpoints, Xsize, Ysize, &
       &            Efermi, TextSize, Emin, Emax, &
       &            CharacterToRGBColor, CharacterToThickness, &
       &            MaxAtoms, MaxOrbitals, fontName)
  call WritePSAxesB(Xsize, Ysize, Emin, Emax, MajorTicks, &
       &            MinorTicks, TextSize, fontName)
  call WritePSTail()
end subroutine SpaghettiPrimavera

subroutine WritePSHead(Xsize, Ysize, version, cmdLine)
  !
  ! writes the header of the postscript to standard output
  !

  implicit none

  real*8 Xsize,Ysize
  character(*) :: version, cmdLine

  write(6,100)"%!PS-Adobe-2.0 EPSF-2.0"
  write(6,101)"%%BoundingBox: 0 0",int(Xsize+100),int(Ysize+100)
  write(6,102)"%%HiResBoundingBox: 0.000000 0.000000",Xsize+100,Ysize+100
  write(6,100)"%%Creator: Spaghetti Primavera by Maurits W. Haverkort"
  write(6,100)"%%prima.py V. "//trim(version)//", command line: " // trim(cmdLine)
  write(6,100)"%%EndComments"
  write(6,100)"%%BeginProlog"
  write(6,100)"1 1 1 setrgbcolor clippath fill % added Sep 2014 to make current ‘evince’ versions display white bg"
  write(6,100)"save"
  write(6,100)"countdictstack"
  write(6,100)"mark"
  write(6,100)"newpath"
  write(6,100)"/setpagedevice {pop} def"
  write(6,100)"/L{ newpath setlinewidth setrgbcolor moveto lineto&
       & lineto 0 setlinecap 1 setlinejoin stroke } def"
  write(6,100)"/L1{ newpath setlinewidth setrgbcolor moveto&
       & lineto 0 setlinecap 1 setlinejoin stroke } def"
  write(6,100)"%%EndProlog"
  write(6,100)"%%Page 1 1"
100 FORMAT(A)
101 FORMAT(A,2I5)
102 FORMAT(A,2F12.6)
end subroutine WritePSHead

subroutine WritePSTail()
  !
  ! Writes the tail of the postscript to standard output
  !
  implicit none
  write(6,100)"%%Trailer"
  write(6,100)"cleartomark"
  write(6,100)"countdictstack"
  write(6,100)"exch sub { end } repeat"
  write(6,100)"restore"
  write(6,100)"%%EOF"
100 FORMAT(A)
end subroutine WritePSTail

subroutine WritePSAxesA(Xsize, Ysize, Emin, Emax)
  !
  ! Writes the part of the axes system to the standard output
  ! that should be behind the graph 
  !
  implicit none
  real*8 Xsize, Ysize, Emin, Emax

  real*8 PtEfermi,PtYunit
  ! A line at the fermi energy
  PtEfermi=50-(Ysize*Emin)/(Emax-Emin)
  PtYunit=(Ysize)/(Emax-Emin)
  if ((Emin.LT.0).AND.(Emax.GE.0)) then
     write(6,100)"newpath"
     write(6,100)"1 setlinewidth"
     write(6,100)"0.5 0.5 0.5 setrgbcolor"
     write(6,103)"50 ",PtEfermi," moveto"
     write(6,102)50+Xsize,PtEfermi," lineto"
     write(6,100)"stroke"
  endif

100 FORMAT(A)
102 FORMAT(2F12.6,A)
103 FORMAT(A,F12.6,A)
end subroutine WritePSAxesA

subroutine WritePSAxesB(Xsize, Ysize, Emin, Emax, MajorTicks, &
     &                  MinorTicks, TextSize, fontName)
  !
  ! Writes the part of the axes system that should be on top of the graph
  !

  implicit none

  integer MinorTicks
  real*8 Xsize, Ysize, Emin, Emax, MajorTicks, TextSize
  character(*) :: fontName

  real*8 PtEfermi,PtYunit,PtTick
  ! some positions in pts
  PtEfermi=50-(Ysize*Emin)/(Emax-Emin)
  PtYunit=(Ysize)/(Emax-Emin)
  ! write the ticks
  if((Emin+MajorTicks).LT.Emax) then
     write(6,100)"newpath"
     write(6,100)"1 setlinewidth"
     write(6,100)"0 0 0 setrgbcolor"
     PtTick=PtEfermi
     do while (PtTick.LT.(Ysize+50))
        write(6,103)"50 ",PtTick," moveto"
        write(6,102)50+Xsize/100.0,PtTick," lineto"
        write(6,102)50+Xsize,PtTick," moveto"
        write(6,102)50+Xsize-Xsize/100.0,PtTick," lineto"
        PtTick=PtTick+MajorTicks*PtYunit
     enddo
     PtTick=PtEfermi
     PtTick=PtTick-MajorTicks*PtYunit
     do while (PtTick.GE.(50))
        write(6,103)"50 ",PtTick," moveto"
        write(6,102)50+Xsize/100.0,PtTick," lineto"
        write(6,102)50+Xsize,PtTick," moveto"
        write(6,102)50+Xsize-Xsize/100.0,PtTick," lineto"
        PtTick=PtTick-MajorTicks*PtYunit
     enddo
     write(6,100)"stroke"
     write(6,100)"newpath"
     write(6,100)"0.75 setlinewidth"
     write(6,100)"0 0 0 setrgbcolor"
     PtTick=PtEfermi
     do while (PtTick.LT.(Ysize+50))
        write(6,103)"50 ",PtTick," moveto"
        write(6,102)50+Xsize/150.0,PtTick," lineto"
        write(6,102)50+Xsize,PtTick," moveto"
        write(6,102)50+Xsize-Xsize/150.0,PtTick," lineto"
        PtTick=PtTick+MajorTicks*PtYunit/(1.0+MinorTicks)
     enddo
     PtTick=PtEfermi
     PtTick=PtTick-MajorTicks*PtYunit/(1.0+MinorTicks)
     do while (PtTick.GE.(50))
        write(6,103)"50 ",PtTick," moveto"
        write(6,102)50+Xsize/150.0,PtTick," lineto"
        write(6,102)50+Xsize,PtTick," moveto"
        write(6,102)50+Xsize-Xsize/150.0,PtTick," lineto"
        PtTick=PtTick-MajorTicks*PtYunit/(1.0+MinorTicks)
     enddo
     write(6,100)"stroke"
  endif
  ! write the numbers next to the ticks
  write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,100)"0 0 0 setrgbcolor"
  PtTick=PtEfermi
  do while (PtTick.LT.(Ysize+50))
     write(6,102)50-2*TextSize,PtTick-TextSize/4.0," moveto"
     write(6,104)"(",(PtTick-PtEfermi)/PtYunit,") show"
     PtTick=PtTick+MajorTicks*PtYunit
  enddo
  PtTick=PtEfermi
  PtTick=PtTick-MajorTicks*PtYunit
  do while (PtTick.GE.(50))
     write(6,102)50-2*TextSize,PtTick-TextSize/4.0," moveto"
     write(6,104)"(",(PtTick-PtEfermi)/PtYunit,") show"
     PtTick=PtTick-MajorTicks*PtYunit
  enddo
  ! write a frame around the graph
  write(6,100)"newpath"
  write(6,100)"1 setlinewidth"
  write(6,100)"0 setlinejoin"
  write(6,100)"0 0 0 setrgbcolor"
  write(6,100)"50 50 moveto"
  write(6,101)Xsize+50," 50 lineto"
  write(6,102)Xsize+50,Ysize+50," lineto"
  write(6,103)"50",Ysize+50," lineto"
  write(6,100)"closepath"
  write(6,100)"stroke"

100 FORMAT(A)
101 FORMAT(F12.6,A)
102 FORMAT(2F12.6,A)
103 FORMAT(A,F12.6,A)
104 FORMAT(A,F5.1,A)
end subroutine WritePSAxesB

subroutine WritePSKpoints(bandName, nkpoints, Xsize, Ysize, &
     &                    TextSize, fontName)
  !
  ! reads the file bandName and 
  ! makes vertical lines for each k-point that has a label in the 
  ! file band-name. These lines are labeled
  ! returns the number of points that should be between each k-point
  ! and the number of k-points
  !
  implicit none

  character(len=*) bandName, fontName
  real*8 Xsize, Ysize, TextSize
  integer nkpoints

  integer OK,nlabels,labelpoint(100),ilabel
  character RDWD*(128),label(100)*(2)

  open(unit=1,file=bandName,status="old",iostat=OK)
  if(OK.NE.0) then
     write(6,100)"/" // trim(fontName) // " findfont"
     write(6,101)TextSize," scalefont"
     write(6,100)"setfont"
     write(6,100)"newpath"
     write(6,100)"1 0 0 setrgbcolor"
     write(6,103)"50 ",Ysize/2+50," moveto"
     write(6,100)"(Could not open file",bandName,") show"
     call WritePSTail()
     stop
  endif
  ! read file "bandName" and filter the number of k-points, 
  ! the labels and their position
  nkpoints=0
  nlabels=0
  read(unit=1,fmt=100,iostat=OK,err=200,end=201) RDWD
  do while (RDWD.NE."END")
     nkpoints=nkpoints+1
     if(RDWD(1:1).NE." ") then
        nlabels=nlabels+1
        label(nlabels)=RDWD(1:2)
        labelpoint(nlabels)=nkpoints
     endif
     read(unit=1,fmt=100,iostat=OK,err=200,end=201) RDWD
  enddo
  ! plot lines for specific k-points
  ilabel=1
  do while (ilabel.LT.nlabels-1)
     ilabel=ilabel+1
     write(6,100)"newpath"
     write(6,100)"1 setlinewidth"
     write(6,100)"0.5 0.5 0.5 setrgbcolor"
     write(6,101)50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0),&
          &              " 50 moveto"
     write(6,102)50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0),&
          &              50+Ysize," lineto"
     write(6,100)"stroke"
  enddo
  ! write the labels for specific k-points
  ilabel=0
  do while (ilabel.LT.nlabels)
     ilabel=ilabel+1
     if (label(ilabel)(2:2).EQ." ") then
        write(6,100)"/" // trim(fontName) // " findfont"
     else
        write(6,100)"/Symbol findfont"
     endif
     write(6,101)TextSize," scalefont"
     write(6,100)"setfont"
     write(6,100)"0 0 0 setrgbcolor"
     write(6,102)50+(labelpoint(ilabel)-1.0)*Xsize/(nkpoints-1.0)-&
          &              TextSize/2.0,50-TextSize," moveto"
     write(6,100)"(",label(ilabel)(1:1),") show"
  enddo
  close(1)
  return
  ! Goto 200 if the file "bandName" had an error upon reading
200 write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,100)"newpath"
  write(6,100)"1 0 0 setrgbcolor"
  write(6,103)"50 ",Ysize/2+50," moveto"
  write(6,100)"(Unspecified error while reading file ",&
       &              bandName,") show"
  call WritePSTail()
  close(1)
  stop
  ! Goto 201 if the file "bandName" ended before the string "END" was read
201 write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,100)"newpath"
  write(6,100)"1 0 0 setrgbcolor"
  write(6,103)"50 ",Ysize/2+50," moveto"
  write(6,100)"(Unexpected end of file while reading ",&
       &              bandName,") show"
  call WritePSTail()
  stop
  close(1)

100 FORMAT(A)
101 FORMAT(F12.6,A)
102 FORMAT(2F12.6,A)
103 FORMAT(A,F12.6,A)
end subroutine WritePSKpoints

subroutine WritePSLegend(TextSize, fontName, XSize, YSize, &
     &                   colors, legend, legentries)
  implicit none

  real*8 TextSize, XSize, YSize
  real(8) :: colors(:,:)
  character(*) :: fontName
  character    :: legend(:)
  integer      :: legentries(:)

  integer :: i, s=0
  character(64) l
  real(8) :: XPos=48, sep

  sep=.66*TextSize + 1.5*TextSize

  do i = 1, size(legentries)
     write(l, '(64A)'), legend(s+1 : legentries(i))
     call WritePSLegE(TextSize, fontName, &
          &           XPos, YSize+70, &
          &           colors(i, :), l)

     XPos = XPos + sep + .66*len_trim(l)*TextSize

     s = legentries(i)
  end do

  XSize=XSize                    ! avoid warning
end subroutine WritePSLegend

subroutine WritePSLegE(TextSize, fontName, XPos, YPos, Color, Name)
  implicit none
  real*8 TextSize, XPos, YPos, Color(3)
  
  character(*) Name, fontName
  real(8) :: wd, ht, sep

  wd=.66*TextSize
  ht=wd
  sep=.5*TextSize

  write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,102) XPos+sep, YPos, " moveto"
  write(6,105)0d0, 0d0, 0d0, " setrgbcolor"
  write(6,100)"(", trim(Name) ,") show"
  write(6,100)"newpath"
  write(6,102)XPos,    YPos+ht, " moveto"
  write(6,102)XPos+wd, YPos+ht, " lineto"
  write(6,102)XPos+wd, YPos,             " lineto"
  write(6,102)XPos,    YPos,             " lineto"
  write(6,100)" closepath"
  write(6,105)Color, " setrgbcolor"
  write(6,100)" fill"
  write(6,100)" stroke"

100 FORMAT(A)
101 FORMAT(F12.6,A)
102 FORMAT(2F12.6,A)
105 FORMAT(3F6.3,A)
end subroutine WritePSLegE

subroutine WritePSBands(qtlName, nkpoints ,Xsize, Ysize, &
     &                  Efermi, TextSize, Emin, Emax, &
     &                  CharacterToRGBColor, CharacterToThickness, &
     &                  MaxAtoms, MaxOrbitals, fontName)

  implicit none

  external CharacterToRGBColor, CharacterToThickness

  interface
     subroutine CharacterToRGBColor(OrbCharacter, RGBColor)
       use sppv_data
       real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
       real(8), intent(out) :: RGBColor(3)
     end subroutine CharacterToRGBColor

     subroutine CharacterToThickness(OrbCharacter, Thickness)
       use sppv_data
       real(8), intent(in)  :: OrbCharacter(MaxAtoms, MaxOrbitals)
       real(8), intent(out) :: Thickness
     end subroutine CharacterToThickness
  end interface

  character(len=*) :: qtlName, fontName
  integer          :: nkpoints, MaxAtoms, MaxOrbitals
  real*8           :: Xsize, Ysize, Emin, Emax, Efermi, TextSize

  real RydToEv
  parameter (RydToEv=13.6056981)

  real*8 OrbCharacter(MaxAtoms,MaxOrbitals),Thickness,RGBColor(3),&
       &       Energy,Eprev,Enow,Enext,ThicknessB,RGBColorB(3),&
       &       ptsPerE,ptsPerK
  integer OK,nAtom,nOrb(MaxAtoms),i,j,k,junki
  character RDWD*(256)

  ptsPerE=Ysize/(Emax-Emin)
  ptsPerK=Xsize/(nkpoints-1)
  ! Open the qtl file
  open(unit=2,file=qtlName,status="old",iostat=OK)
  if(OK.NE.0) then
     write(6,100)"/" // trim(fontName) // " findfont"
     write(6,101)TextSize,"12 scalefont"
     write(6,100)"setfont"
     write(6,100)"newpath"
     write(6,100)"1 0 0 setrgbcolor"
     write(6,103)"50 ",Ysize/2+50," moveto"
     write(6,100)"(Could not open file",qtlName,") show"
     call WritePSTail()
     stop
  endif

  nOrb=0
  OrbCharacter=0

  ! find out how many atoms and orbitals there are
  nAtom=0
  read(unit=2,fmt=100,iostat=OK,err=200,end=201) RDWD
  do while (index(RDWD,"BAND").EQ.0) ! newer Wien2k versions seem to omit the colon in “BAND:”
     if(index(RDWD,"JATOM").NE.0) then
        nAtom=nAtom+1

        if (nAtom > MaxAtoms) then
           write(0, '("sppv: maximum number of atoms (", I0, ") exceeded")') &
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
     read(unit=2,fmt=100,iostat=OK,err=200,end=201) RDWD
  enddo
  ! add the intersitial character
  nAtom=nAtom+1
  nOrb(nAtom)=1
  ! read nkpoints k points then read the string BAND: XX if the 
  ! string BAND: XX can not be read while EOF has been reached stop
  do
     do i = 1, nkpoints
        do j = 1, nAtom
           read(unit=2,fmt=*,iostat=OK,err=200,end=201) Energy,junki,&
                &                            (OrbCharacter(j,k),k=1,nOrb(j))
           Energy=(Energy-Efermi)*RydToEv
        enddo

        call CharacterToThickness(OrbCharacter, Thickness)
        call CharacterToRGBColor (OrbCharacter, RGBColor)

        if(i.EQ.1)then
           Eprev=Energy
           RGBColorB=RGBColor
           ThicknessB=Thickness
        elseif(i.EQ.2)then
           Enow=Energy
           if((Eprev.GT.Emin).AND.(Eprev.LT.Emax).AND.&
                &    (((Eprev+Enow)/2.0).GT.Emin).AND.(((Eprev+Enow)/2.0).LT.Emax)) &
                &  then
              write(6,105)50+(i-2.0)*ptsPerK,50+(Eprev-Emin)*ptsPerE, &
                   &                50+(i-1.5)*ptsPerK,50+(0.5*(Enow+Eprev)-Emin)* &
                   &                ptsPerE,RGBColorB,ThicknessB," L1"
           endif
           RGBColorB=RGBColor
           ThicknessB=Thickness
        else
           Enext=Energy
           if((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Eprev).GT.Emin).AND.  &
                &     (0.5*(Enow+Eprev).LT.Emax).AND.  &
                &     (0.5*(Enow+Enext).GT.Emin).AND.  &
                &     (0.5*(Enow+Enext).LT.Emax))then
              write(6,106)50+(i-2.5)*ptsPerK,50+(0.5*(Eprev+Enow)-Emin)* &
                   &                ptsPerE,50+(i-2.0)*ptsPerK,50+(Enow-Emin)* &
                   &                ptsPerE,50+(i-1.5)*ptsPerK,50+(0.5*(Enext+Enow)- &
                   &                Emin)*ptsPerE,RGBColorB,ThicknessB," L"
           elseif((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Eprev).GT.Emin).AND. &
                &     (0.5*(Enow+Eprev).LT.Emax))then
              write(6,105)50+(i-2.5)*ptsPerK,50+(0.5*(Eprev+Enow)-Emin)* &
                   &                ptsPerE,50+(i-2.0)*ptsPerK,50+(Enow-Emin)* &
                   &                ptsPerE,RGBColorB,ThicknessB," L1"
           elseif((Enow.GT.Emin).AND.(Enow.LT.Emax).AND. &
                &     (0.5*(Enow+Enext).GT.Emin).AND. &
                &     (0.5*(Enow+Enext).LT.Emax))then
              write(6,105)50+(i-1.5)*ptsPerK,50+(0.5*(Enext+Enow)-Emin)* &
                   &                ptsPerE,50+(i-2.0)*ptsPerK,50+(Enow-Emin)* &
                   &                ptsPerE,RGBColorB,ThicknessB," L1"
           endif
           RGBColorB=RGBColor
           ThicknessB=Thickness
           if(i.EQ.nkpoints)then
              if((Enext.GT.Emin).AND.(Enext.LT.Emax).AND. &
                   &    (((Enext+Enow)/2.0).GT.Emin).AND.(((Enext+Enow)/2.0).LT.Emax)) &
                   &   then
                 write(6,105)50+(i-1.0)*ptsPerK,50+(Enext-Emin)*ptsPerE, &
                      &                50+(i-1.5)*ptsPerK,50+(0.5*(Enow+Enext)-Emin)* &
                      &                ptsPerE,RGBColorB,ThicknessB," L1"
              endif
           endif
           Eprev=Enow
           Enow=Enext
        endif
     enddo
     ! read the string BAND: XX or iff EOF finished
     read(unit=2,fmt=100,iostat=OK,err=200,end=300) RDWD
  enddo
300 continue
  close(2)
  return
  ! Goto 200 if the file "qtlName" had an error upon reading
200 write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,100)"newpath"
  write(6,100)"1 0 0 setrgbcolor"
  write(6,103)"50 ",Ysize/2+50," moveto"
  write(6,100)"(Unspecified error while reading file ", &
       &              qtlName,") show"
  call WritePSTail()
  close(2)
  stop
  ! Goto 201 if the file "qtlName" ended before the end was expected
201 write(6,100)"/" // trim(fontName) // " findfont"
  write(6,101)TextSize," scalefont"
  write(6,100)"setfont"
  write(6,100)"newpath"
  write(6,100)"1 0 0 setrgbcolor"
  write(6,103)"50 ",Ysize/2+50," moveto"
  write(6,100)"(Unexpected end of file while reading ", &
       &              qtlName,") show"
  call WritePSTail()
  stop
  close(2)

100 FORMAT(A)
101 FORMAT(F12.6,A)
103 FORMAT(A,F12.6,A)
105 FORMAT(8F10.4,A)
106 FORMAT(10F10.4 A)
end subroutine WritePSBands

# Microsoft Developer Studio Project File - Name="trove" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=trove - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "trove.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "trove.mak" CFG="trove - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "trove - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "trove - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "trove - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "trove - Win32 Release"
# Name "trove - Win32 Debug"
# Begin Source File

SOURCE=.\accuracy.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\diag.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_DIAG_=\
	".\Debug\accuracy.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\dipole.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_DIPOL=\
	".\Debug\accuracy.mod"\
	".\Debug\fields.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	".\Debug\tran.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\fields.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_FIELD=\
	".\Debug\accuracy.mod"\
	".\Debug\input.mod"\
	".\Debug\lapack.mod"\
	".\Debug\me_bnd.mod"\
	".\Debug\me_numer.mod"\
	".\Debug\me_rot.mod"\
	".\Debug\me_str.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\input.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\lapack.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_LAPAC=\
	".\Debug\accuracy.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\me_bnd.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_ME_BN=\
	".\Debug\accuracy.mod"\
	".\Debug\me_numer.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\me_numer.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_ME_NU=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\me_rot.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_ME_RO=\
	".\Debug\accuracy.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\me_str.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_ME_ST=\
	".\Debug\accuracy.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_abcd.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_A=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_ch3oh.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_C=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_xy.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_X=\
	".\Debug\accuracy.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_xy2.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_XY=\
	".\Debug\accuracy.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_xy3.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_XY3=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_xy4.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_XY4=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mol_zxy2.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOL_Z=\
	".\Debug\accuracy.mod"\
	".\Debug\moltype.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\molecules.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOLEC=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	".\Debug\mol_abcd.mod"\
	".\Debug\mol_ch3oh.mod"\
	".\Debug\mol_xy.mod"\
	".\Debug\mol_xy2.mod"\
	".\Debug\mol_xy3.mod"\
	".\Debug\mol_xy4.mod"\
	".\Debug\mol_zxy2.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\moltype.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_MOLTY=\
	".\Debug\accuracy.mod"\
	".\Debug\lapack.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\perturbation.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_PERTU=\
	".\Debug\accuracy.mod"\
	".\Debug\diag.mod"\
	".\Debug\fields.mod"\
	".\Debug\lapack.mod"\
	".\Debug\me_numer.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\refinement.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_REFIN=\
	".\Debug\accuracy.mod"\
	".\Debug\fields.mod"\
	".\Debug\lapack.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	".\Debug\tran.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\symmetry.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_SYMME=\
	".\Debug\accuracy.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\timer.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_TIMER=\
	".\Debug\accuracy.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\tran.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_TRAN_=\
	".\Debug\accuracy.mod"\
	".\Debug\fields.mod"\
	".\Debug\me_numer.mod"\
	".\Debug\molecules.mod"\
	".\Debug\moltype.mod"\
	".\Debug\perturbation.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\trove.f90

!IF  "$(CFG)" == "trove - Win32 Release"

!ELSEIF  "$(CFG)" == "trove - Win32 Debug"

DEP_F90_TROVE=\
	".\Debug\accuracy.mod"\
	".\Debug\dipole.mod"\
	".\Debug\fields.mod"\
	".\Debug\moltype.mod"\
	".\Debug\perturbation.mod"\
	".\Debug\refinement.mod"\
	".\Debug\symmetry.mod"\
	".\Debug\timer.mod"\
	".\Debug\tran.mod"\
	

!ENDIF 

# End Source File
# End Target
# End Project

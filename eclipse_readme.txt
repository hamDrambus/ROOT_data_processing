Setups:
C/C++ Build->Settings->
    GCC C++ Compiler->
	->Dialect->Other dialect flags: -std=c++11
	->Preprocessor->-D: NDEBUG
			    R__HAVE_CONFIG
	->Includes: /home/Frolov/Documents/ROOT/root_v6.08.04_install/include/root
		    /home/Frolov/Documents/ElectronDriftInAr/OneDimensional/include
	->Miscellaneous: -fPIC -pthread
    GCC C++ Linker->
	->General: -pthread
	->Libraries: Gui
Core
RIO
Net
Hist
Graf
Graf3d
Gpad
Tree
Rint
Postscript
Matrix
Physics
MathCore
Thread
MultiProc
Geom
m
Spectrum
Thread
	/home/Frolov/Documents/ROOT/root_v6.08.04_install/lib/root
/home/Frolov/Documents/ElectronDriftInAr/OneDimensional/include
	->Miscellaneous->Linker flags: -m64 -fPIC -fsigned-char -pipe -std=c++11

C/C++ General->Paths and Symbols->
    ->Source Location
    

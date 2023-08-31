Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.
2) CERN Root. Do not forget to run '$source ../bin/thisroot.sh' and to add root's library paths to LD_LIBRARY_PATH. Best to add 'source ../bin/thisroot.sh' to startup script in /etc/profile.d/

Not necessary for compilation, but necessary for overall functioning:
3) gnuplot - for plotting data (example waveforms).

=================================================
To compile this code without eclipse, use CleanCompile.sh script.
Use Recompile.sh to recompile if there is no changes in file structure of the project (no need to call cmake).
Run.sh can be used to run the program from the 'correct' directory. 

Alternatively, run GenEclipseProject.sh script to create Eclipse project.

Open in eclipse as [File->Import->General->Existing project in folder]->browse to generated -build folder. Build via [Project Explorer->Build Targets]. Debug as 'C/C++ Attach to Application' or 'C/C++ Application', with set binary location and disabled auto build (if necessary). When any files or libraries are added to/removed from the project, it must be regenerated with GenEclipseProject.sh.

=================================================
How to run data processing:
	1) Install required and recommended dependencies.
	2) If eclipse is used then set up the project (script GenEclipseProject.sh is provided).
	3) Select manifests corresponding the the data to be analyzed in GlobalParameters.cpp::Init_globals
	4) Compile the code using CMake or Eclipse.
	5) Run the code. Make sure to run from correct directory (so that paths in manifests point to the data).

=================================================
Setting up the processing

Simulation parameters are hard-coded in Manifest files. Adding new data analysis requires adding new manifest file and recompiling the project, although it is a simple matter.

Take note that this program is a counterpart of Post_processing. Data_processing produces photoelectron peak data in binary format which is an input for Post_processing.



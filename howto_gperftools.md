Usage of gperftools profiler with GeantV on linux
=================================================

1. Download and install.

For extensive instructions, see: https://github.com/gperftools/gperftools/blob/master/INSTALL
Two options are possible:

a. Download the development package from your preferred distibution (e.g. libgoogle-perftools-dev on Ubuntu)
or:
b. Install manually from the github repository:
Make sure libunwind is installed on your system. This is critical to determine the call-chain of a program
git clone https://github.com/gperftools/gperftools.git
cd gperftools
./autogen.sh
./configure --enable-frame-pointers --with-tcmalloc-pagesize=32
sudo make -j install (in case you want libs in /usr/local/lib)

2. Compile  with gperftools support

Make sure your installation of gperftoold is findable by adding: -DCMAKE_PREFIX_PATH=/usr/local or whatever installation path you have chosen. Then just add -DGPERFTOOLS=ON and compile, watching that the tool detection went fine. For project dependencies it is recommended to compile with -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer -gdwarf" (done already in the project)

Compile with debugging info on (RelWithDebInfo or Debug)

3. Profiling can be enabled via API:

// Include the profiler header if Gperftools support is ON
#ifdef USE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

// Surrond the profiled function
#ifdef USE_GPERFTOOLS
  ProfilerStart("gperfprof.out");
#endif

  ProfiledFunction();

#ifdef USE_GPERFTOOLS
  ProfilerStop();
#endif

Run your executable like below. The default sampling frequency is 100, enough for long running tests, but too small for short ones:
  CPUPROFILE_FREQUENCY=500 MyExecutable <params>

4. Inspect the profile

Output the profile tree in output.pdf:
  pprof --pdf MyExecutable gperfprof.out > output.pdf

There are other options and profiling tweaks described here: https://github.com/gperftools/gperftools/blob/master/README

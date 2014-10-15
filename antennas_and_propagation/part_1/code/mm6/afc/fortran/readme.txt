The Fortran kernel needs to be compiled first by running the appropriate make batch. G95 Fortran compiler must be installed on your machine. Choose the batch that suits your target machine or make your own one:

:: make_opteron.bat  ... AMD Opteron
:: make_prescott.bat ... Intel Pentium 4, Prescott core (SSE3)
:: make_nocona.bat   ... Intel with 64-bit extensions (e.g. Core2)

:: make_fast ... fast speed of compiler, but no optimizations -- should work on most platforms

More info on the compiler options:

http://www.g95.org/docs.html#rung95
http://gcc.gnu.org/onlinedocs/gcc-4.1.1/gcc/
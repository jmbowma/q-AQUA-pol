Details and instructions of q-AQUA-pol water potential:

I. General information
(1) Folder "coef_qAQUA-pol" contains all fitted coefficients for 2b,3b,4b fits
(2) Folder "src" contains all source codes for the water potential together with examples

II. Steps for compiling the code
(1) Go to folder "src"
(2) Check the Makefile
In this example, we use IFORT compiler and OPENMP parallelization can be enabled with modifications.
User can change the compiler as needed, such as GFORTRAN

(3) One example of calling the potential and gradients is given in "getpot.f90"

(4) Compile the code through "Make getpot.x". Note, it may take some time to compile the 2b,3b, and 4b library files such as bemsa2b.f90. Just be patient :)

(5) A benchmark test is give with one water hexamer geometry "hexamer.xyz"
Energy: -45.7738923286911 kcal/mol
Gradient(a.u.):
   1     0.00394754    -0.00182171     0.00069273
   2    -0.00228976    -0.00142203     0.00400568
   3    -0.00075949    -0.00061204     0.00420559
   4    -0.00292673    -0.00130259    -0.00122316
   5    -0.00168868    -0.00239567     0.00269531
   6    -0.00187900    -0.00180394    -0.00272606
   7     0.00003051    -0.00049893    -0.00338875
   8     0.00253757     0.00353316    -0.00065256
   9    -0.00338173     0.00147785    -0.00057620
  10     0.00113714     0.00272610     0.00042893
  11    -0.00328271     0.00208605    -0.00000976
  12     0.00267610     0.00330605    -0.00039869
  13    -0.00129751     0.00319493    -0.00475530
  14     0.00379564     0.00193308    -0.00310067
  15     0.00382317     0.00406095     0.00004663
  16    -0.00307279    -0.00295002     0.00401148
  17     0.00213928    -0.00433354     0.00042017
  18     0.00049144    -0.00517769     0.00032463

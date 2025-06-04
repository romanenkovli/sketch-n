!==============================================================================!
!                                   SKETCH-N                                   !
!  Nodal Neutron Diffusion Code for Solving Steady-State and Kinetics Problems !
!              Version 1.0 (c) Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !
!      Date & Time of Calculation :  1 January  2002   5:04:42 pm              !
!==============================================================================!
 
!------------------------------------------------------------------------------!
!     INPUT/OUTPUT FILES:                                                      !
!------------------------------------------------------------------------------!
 
   Input General Data: Input/PWR_NEACRP_A1.DAT                                                                             
   Output Kinetics Restart File: Restart/pwr_neacrp_A1_0.dat                                                                         
 
 
!------------------------------------------------------------------------------!
!     SKETCH-N OPTIONS                                                         !
!------------------------------------------------------------------------------!
 
   Space Discretization            : Semi-Analytic Method (SANM)                                                     
   Transverse Leakage Approximation: QLA 
   Nonlinear Iteration Procedure   : Smith
 
!------------------------------------------------------------------------------!
!     ITERATIVE CONVERGENCE CRITERIA                                           !
!------------------------------------------------------------------------------!
 
   Number of Outer Iterations per Nonlinear (convergence) :            1
   Number of Outer Iterations per Nonlinear (calculation) :           10
   Eigenvalue Accuracy when Critical Boron Search Starts  :  0.10000E-01
   Eigenvalue Accuracy when Critical Boron Search Finished:  0.10000E-04
   Convergence Criteria for Eigenvalue                    :  0.10000E-04
   Convergence Criteria on Neutron Flux                   :  0.10000E-04
   Number of Inner Iterations per Outer                   :            2
   Convergence Criteria for Inner Iterations              :  0.10000E-07
 
 
!------------------------------------------------------------------------------!
!     ITERATIVE PARAMETERS                                                     !
!------------------------------------------------------------------------------!
 
   Estimate of Minimum Eigenvalue of the Iteration Matrix  :  0.00000E+00
   Estimate of Maximum Eigenvalue or Dominance Ratio
                                  of the Iteration Matrix  :  0.80000E+00
   Adaptivity Parameter of the Chebyshev Acceleration      :  0.00000E+00
   Outer Iteration Number when Chebyshev Acceleration Start:            5
   Inverse Value of the Wieland Shift                      :  0.00000E+00
 
!------------------------------------------------------------------------------!
!   3D NEACRP PWR Rod Ejection Benchmark, case A1                              !
!      (c) Slava March 1998 JAERI                                              !
!                                                                              !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!     PROBLEM DESCRIPTION                                                      !
!------------------------------------------------------------------------------!
 
 Reactor Type:   PWR
 Reactor Power [MWt]: 0.693750E-03
 Type of the Computing Problem: BoronSearch         
 Thermal-Hydraulics Model: SKAZKA              
 
 Geometry Parameters
    Number of Dimensions:                                      3
    Initial Spatial Mesh (x,y,z):                        9x 9x18
    Neuronics Spatial Mesh (x,y,z):                     17x17x18
    Number of Bundles (Assemblies) in x-y plane:              64
    Number of Nodes per Bundle (Assembly) in x-y plane:        4
    Number of Neutronics Nodes in x-y plane:                 221
    Number of Bundle (Assembly) Types:                         9
 Geometry Data :

    Numbering the Reactor Bundles (Assemblies):                                 

        1   2   3   4   5   6   7   8   9
      --- --- --- --- --- --- --- --- ---
  1:    1   2   3   4   5   6   7   8   9
  2:   10  11  12  13  14  15  16  17  18
  3:   19  20  21  22  23  24  25  26  27
  4:   28  29  30  31  32  33  34  35
  5:   36  37  38  39  40  41  42  43
  6:   44  45  46  47  48  49  50
  7:   51  52  53  54  55  56
  8:   57  58  59  60  61
  9:   62  63  64
 
    Initial Spatial Mesh X:
         10.803 21.606 21.606 21.606 21.606 21.606 21.606 21.606 21.606
    Initial Spatial Mesh Y:
         10.803 21.606 21.606 21.606 21.606 21.606 21.606 21.606 21.606
    Initial Spatial Mesh Z:
         30.000  7.700 11.000 15.000 30.000 30.000 30.000 30.000 30.000 30.000
         30.000 30.000 30.000 30.000 12.800 12.800  8.000 30.000
    Neutronics Spatial Mesh X:
         10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803
         10.803 10.803 10.803 10.803 10.803 10.803 10.803
    Neutronics Spatial Mesh Y:
         10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803 10.803
         10.803 10.803 10.803 10.803 10.803 10.803 10.803
    Neutronics Spatial Mesh Z:
         30.000  7.700 11.000 15.000 30.000 30.000 30.000 30.000 30.000 30.000
         30.000 30.000 30.000 30.000 12.800 12.800  8.000 30.000

    Reactor Loading by the Bundle Types                                         

        1   2   3   4   5   6   7   8   9
      --- --- --- --- --- --- --- --- ---
  1:    1   2   1   2   1   3   1   4   5
  2:    2   1   6   1   6   1   7   4   5
  3:    1   6   1   6   1   6   4   8   5
  4:    2   1   6   1   6   9   4   5
  5:    1   6   1   6   1   4   8   5
  6:    3   1   6   9   4   8   5
  7:    1   7   4   4   8   5
  8:    4   4   8   5   5
  9:    5   5   5
 
    Bundle Material Compositions:
          1:   1   4   4   4   4   4   4   4   4   4
               4   4   4   4   4   4   4   1
          2:   1   5   9   9   9   9   9   9   9   9
               9   9   9   9   9   9   9   1
          3:   1   5   7   7   7   7   7   7   7   7
               7   7   7   7   7   7   7   1
          4:   1   6   6   6   6   6   6   6   6   6
               6   6   6   6   6   6   6   1
          5:   2   2   2   2   2   2   2   2   2   2
               2   2   2   2   2   2   2   2
          6:   1   5   8   8   8   8   8   8   8   8
               8   8   8   8   8   8   8   1
          7:   1   6  11  11  11  11  11  11  11  11
              11  11  11  11  11  11  11   1
          8:   3   3   3   3   3   3   3   3   3   3
               3   3   3   3   3   3   3   3
          9:   1   6  10  10  10  10  10  10  10  10
              10  10  10  10  10  10  10   1
 
    Boundary Conditions Types (left, right, up, down, top, bottom) :
      "0" is zero flux, "1" is logarithmic derivative
         1 0 1 0 0 0
    Boundary Conditions Constants:
       left  x:    0.000  0.000
       right x:    0.000  0.000
          up y:    0.000  0.000
        down y:    0.000  0.000
         top z:    0.000  0.000
      bottom z:    0.000  0.000
!------------------------------------------------------------------------------!
 
    Reactor Core Parameters:
       Number of Bundles (Assemblies) in x-y plane:              47
       Number of Axial Layers                     :              16
       Reactor Core Height, [cm]                  :         367.300

    Reactor Core Geometry:                                                      

        1   2   3   4   5   6   7   8
      --- --- --- --- --- --- --- ---
  1:    1   2   3   4   5   6   7   8
  2:    9  10  11  12  13  14  15  16
  3:   17  18  19  20  21  22  23
  4:   24  25  26  27  28  29  30
  5:   31  32  33  34  35  36
  6:   37  38  39  40  41
  7:   42  43  44  45
  8:   46  47
 
!------------------------------------------------------------------------------!
 
 Control Rod Parameters:
    Number of Control Rods                                 :      15
    Number of Control Rod Types                            :       2
    Number of Control Rod Elements (Absorber, Driver, etc.):       2
    Number of Control Rod Material Compositions            :       3
    Number of Bundles Covered by Control Rod               :       1
 Control Rod Data :

    Location of the Control Rods in the Core  :                                 

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:     1    -    2    -    -    -    3    -
  2:     -    4    -    -    -    5    -    -
  3:     6    -    7    -    8    -    9
  4:     -    -    -   10    -    -    -
  5:     -    -   11    -   12    -
  6:     -   13    -    -    -
  7:    14    -   15    -
  8:     -    -
 

    Control Rod Types                         :                                 

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:     1    -    1    -    -    -    1    -
  2:     -    1    -    -    -    1    -    -
  3:     1    -    1    -    1    -    2
  4:     -    -    -    1    -    -    -
  5:     -    -    1    -    1    -
  6:     -    1    -    -    -
  7:     1    -    2    -
  8:     -    -
 
    Control Rod Material Compositions:
          1:   1   3
          2:   2   3
    Length of the Control Rod Elements (Absorber, Driver): 
           362.159  1000.000
!------------------------------------------------------------------------------!
 
   Cross Section parameters:
       Number of Material Compositions                 :      11
       Number of Neutron Energy Groups                 :       2
       Number of Delayed Neutron Precursor Groups      :       6
       Number of Feedbacks                             :       4
       Macro Cross Section Model                       : POLYNOM   
       Homogenization Method for Partially Roded Nodes : FLUX-WEIGHTING  
 Macro Cross Section Data :
    Diffusion Coefficient Flag in the Input Data :  0
    Basic Macro Cross Section Set:
  0.532058E-01  0.264554E-01  0.373279E-03  0.000000E+00  0.000000E+00
  0.386406E+00  0.177215E-01  0.000000E+00  0.000000E+00             1

  0.295609E+00  0.231613E-01  0.118782E-02  0.000000E+00  0.000000E+00
  0.245931E+01  0.252618E+00  0.000000E+00  0.000000E+00             2

  0.295609E+00  0.200808E-01  0.118782E-02  0.000000E+00  0.000000E+00
  0.245931E+01  0.252618E+00  0.000000E+00  0.000000E+00             3

  0.222117E+00  0.182498E-01  0.871774E-02  0.498277E-02  0.190224E-02
  0.803140E+00  0.652550E-01  0.839026E-01  0.343581E-01             4

  0.221914E+00  0.180040E-01  0.906133E-02  0.557659E-02  0.214498E-02
  0.795538E+00  0.723354E-01  0.998629E-01  0.408938E-01             5

  0.221715E+00  0.177670E-01  0.938496E-02  0.615047E-02  0.237972E-02
  0.789253E+00  0.789203E-01  0.114667E+00  0.469561E-01             6

  0.222039E+00  0.171381E-01  0.931692E-02  0.555010E-02  0.213629E-02
  0.776230E+00  0.796328E-01  0.985576E-01  0.403596E-01             7

  0.222083E+00  0.168501E-01  0.940032E-02  0.554083E-02  0.213318E-02
  0.769969E+00  0.821087E-01  0.980059E-01  0.401338E-01             8

  0.222127E+00  0.165626E-01  0.948286E-02  0.553137E-02  0.213003E-02
  0.763813E+00  0.845912E-01  0.974109E-01  0.398902E-01             9

  0.221836E+00  0.169043E-01  0.963720E-02  0.612382E-02  0.237097E-02
  0.770705E+00  0.861187E-01  0.113241E+00  0.463724E-01            10

  0.221878E+00  0.166175E-01  0.971937E-02  0.611444E-02  0.236781E-02
  0.764704E+00  0.885488E-01  0.112635E+00  0.461246E-01            11

    Polynomial Coefficients for Macro Cross Sections
           Material Composition number :   1
  0.611833E-07  0.791457E-09  0.187731E-06  0.000000E+00  0.000000E+00
  0.517535E-05  0.102635E-04  0.000000E+00  0.000000E+00  0.120020E+04
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.306600E+03
  0.745756E-01  0.371310E-01  0.207688E-03  0.000000E+00  0.000000E+00
  0.533634E+00  0.758421E-02  0.000000E+00  0.000000E+00  0.712500E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.891450E+03
           Material Composition number :   2
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.776184E-03  0.844695E-04  0.000000E+00  0.000000E+00  0.120020E+04
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.306600E+03
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.712500E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.891450E+03
           Material Composition number :   3
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.776184E-03  0.844695E-04  0.000000E+00  0.000000E+00  0.120020E+04
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.306600E+03
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.712500E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00
  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.891450E+03
           Material Composition number :   4
  0.347809E-07 -0.108590E-06  0.128505E-06 -0.112099E-08 -0.548361E-09
 -0.976510E-05  0.708807E-05 -0.243045E-05 -0.995273E-06  0.120020E+04
 -0.203310E-05  0.809676E-06  0.212191E-06  0.124709E-06  0.445177E-07
 -0.108674E-03 -0.315597E-04 -0.416439E-04 -0.170531E-04  0.306600E+03
  0.135665E+00  0.293195E-01  0.155185E-02  0.920694E-03  0.318680E-03
  0.992628E+00  0.252662E-01  0.247746E-01  0.101452E-01  0.712500E+00
 -0.309197E-04 -0.275536E-04  0.349709E-04  0.640134E-06  0.222662E-06
 -0.137292E-03 -0.371806E-04 -0.563037E-04 -0.230564E-04  0.891450E+03
           Material Composition number :   5
  0.353826E-07 -0.106951E-06  0.126709E-06 -0.167880E-08 -0.777980E-09
 -0.850169E-05  0.682311E-05 -0.272445E-05 -0.111566E-05  0.120020E+04
 -0.198080E-05  0.858474E-06  0.226000E-06  0.135145E-06  0.488317E-07
 -0.906150E-04 -0.321435E-04 -0.453102E-04 -0.185545E-04  0.306600E+03
  0.135748E+00  0.292696E-01  0.161491E-02  0.964160E-03  0.336573E-03
  0.981985E+00  0.286667E-01  0.314993E-01  0.128990E-01  0.712500E+00
 -0.308607E-04 -0.276766E-04  0.351798E-04  0.997431E-06  0.369389E-06
 -0.117481E-03 -0.377039E-04 -0.604155E-04 -0.247402E-04  0.891450E+03
           Material Composition number :   6
  0.359838E-07 -0.105374E-06  0.124986E-06 -0.221038E-08 -0.996653E-09
 -0.746251E-05  0.659798E-05 -0.295883E-05 -0.121164E-05  0.120020E+04
 -0.192434E-05  0.903494E-06  0.239939E-06  0.149084E-06  0.545978E-07
 -0.762786E-04 -0.323776E-04 -0.478475E-04 -0.195937E-04  0.306600E+03
  0.135827E+00  0.292154E-01  0.168015E-02  0.101410E-02  0.357207E-03
  0.972267E+00  0.319571E-01  0.381097E-01  0.156060E-01  0.712500E+00
 -0.309165E-04 -0.278390E-04  0.353841E-04  0.141847E-05  0.542384E-06
 -0.101337E-03 -0.377558E-04 -0.630960E-04 -0.258379E-04  0.891450E+03
           Material Composition number :   7
  0.337806E-07 -0.100873E-06  0.119869E-06 -0.171323E-08 -0.777980E-09
 -0.673744E-05  0.629310E-05 -0.255359E-05 -0.104561E-05  0.120020E+04
 -0.269634E-05  0.701311E-06  0.248530E-06  0.140773E-06  0.488317E-07
 -0.762435E-04 -0.300119E-04 -0.420202E-04 -0.172074E-04  0.306600E+03
  0.131033E+00  0.282489E-01  0.168397E-02  0.981951E-03  0.336573E-03
  0.934697E+00  0.314240E-01  0.351588E-01  0.144016E-01  0.712500E+00
 -0.313746E-04 -0.273550E-04  0.348699E-04  0.945431E-06  0.369389E-06
 -0.108271E-03 -0.372748E-04 -0.579662E-04 -0.237320E-04  0.891450E+03
           Material Composition number :   8
  0.332495E-07 -0.988578E-07  0.117585E-06 -0.172421E-08 -0.793325E-09
 -0.619725E-05  0.611904E-05 -0.248880E-05 -0.101904E-05  0.120020E+04
 -0.307905E-05  0.617380E-06  0.261854E-06  0.143235E-06  0.522555E-07
 -0.733397E-04 -0.291929E-04 -0.407701E-04 -0.166956E-04  0.306600E+03
  0.129379E+00  0.278895E-01  0.171972E-02  0.988437E-03  0.346474E-03
  0.918171E+00  0.324715E-01  0.363251E-01  0.148808E-01  0.712500E+00
 -0.315503E-04 -0.272381E-04  0.347274E-04  0.926078E-06  0.339044E-06
 -0.105521E-03 -0.371808E-04 -0.571108E-04 -0.233804E-04  0.891450E+03
           Material Composition number :   9
  0.327201E-07 -0.968489E-07  0.115319E-06 -0.173502E-08 -0.796917E-09
 -0.568220E-05  0.594711E-05 -0.242240E-05 -0.991817E-06  0.120020E+04
 -0.353877E-05  0.516547E-06  0.274313E-06  0.146019E-06  0.534283E-07
 -0.713711E-04 -0.283041E-04 -0.394319E-04 -0.161475E-04  0.306600E+03
  0.127682E+00  0.275202E-01  0.174989E-02  0.995175E-03  0.349236E-03
  0.901293E+00  0.335945E-01  0.374499E-01  0.153431E-01  0.712500E+00
 -0.317281E-04 -0.271169E-04  0.346026E-04  0.905802E-06  0.330425E-06
 -0.102525E-03 -0.370201E-04 -0.561543E-04 -0.229872E-04  0.891450E+03
           Material Composition number :  10
  0.343859E-07 -0.993312E-07  0.118186E-06 -0.224335E-08 -0.996653E-09
 -0.586898E-05  0.608443E-05 -0.277657E-05 -0.113696E-05  0.120020E+04
 -0.263907E-05  0.744320E-06  0.264289E-06  0.155858E-06  0.545978E-07
 -0.639554E-04 -0.303509E-04 -0.444431E-04 -0.181997E-04  0.306600E+03
  0.131116E+00  0.281877E-01  0.175528E-02  0.103522E-02  0.357207E-03
  0.924925E+00  0.349853E-01  0.420693E-01  0.172298E-01  0.712500E+00
 -0.314192E-04 -0.275049E-04  0.350637E-04  0.135642E-05  0.542384E-06
 -0.938886E-04 -0.371403E-04 -0.605052E-04 -0.247739E-04  0.891450E+03
           Material Composition number :  11
  0.338559E-07 -0.973291E-07  0.115917E-06 -0.225369E-08 -0.101112E-08
 -0.538345E-05  0.591697E-05 -0.270780E-05 -0.110878E-05  0.120020E+04
 -0.302147E-05  0.659521E-06  0.279060E-06  0.158814E-06  0.586765E-07
 -0.616984E-04 -0.295626E-04 -0.431588E-04 -0.176738E-04  0.306600E+03
  0.129463E+00  0.278259E-01  0.179499E-02  0.104291E-02  0.368921E-03
  0.908456E+00  0.361032E-01  0.433215E-01  0.177435E-01  0.712500E+00
 -0.315908E-04 -0.273835E-04  0.349119E-04  0.133336E-05  0.506596E-06
 -0.917126E-04 -0.369909E-04 -0.596284E-04 -0.244141E-04  0.891450E+03
    Initial Feedback Values                         
  0.100000E+04  0.286000E+03  0.754039E+00  0.559150E+03
    Differential Cross Sections for Control Rods    
  0.373220E-02 -0.319253E-02  0.247770E-02 -0.102786E-03 -0.377989E-04
 -0.219926E-01  0.255875E-01 -0.282319E-02 -0.115483E-02             1

  0.374092E-02 -0.314239E-02  0.242926E-02 -0.122634E-03 -0.459250E-04
 -0.167503E-01  0.256478E-01 -0.328086E-02 -0.134262E-02             2

  0.697102E-02 -0.119034E-02  0.879034E-04 -0.655496E-04 -0.197926E-04
 -0.113498E-01  0.170043E-02 -0.146252E-02 -0.599154E-03             3

    Prompt Neutron Spectrum:   
   1.00000   0.00000
    Energy Release per Fission:   
  0.321300E-10  0.320600E-10
    Part of the Power Directly Deposited into Coolant:   
  0.190000E-01
    Neutron Velocities:   
  0.280000E+08  0.440000E+06
    Delayed Neutron Decay Constants:   
 0.128000E-01 0.318000E-01 0.119000E+00 0.318100E+00 0.140270E+01 0.392860E+01
    Delayed Neutron Yields:   
 0.258400E-03 0.152000E-02 0.139080E-02 0.307040E-02 0.110200E-02 0.258400E-03
    Total Yield of Delayed Neutrons:   
 0.760000E-02
    Delayed Neutron Spectrum:   
   1.00000   0.00000
!------------------------------------------------------------------------------!
!==============================================================================!
!                                   SKETCH-N                                   !
!  Nodal Neutron Diffusion Code for Solving Steady-State and Kinetics Problems !
!              Version 1.0 (c) Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !
!      Date & Time of Calculation :  1 January  2002   5:04:42 pm              !
!==============================================================================!
                 INITIAL DATA SET
                 STEADY-STATE CALCULATIONS
!------------------------------------------------------------------------------!
     k_ef                                  :  1.000000
     Boron Concentration                   :  0.10000E+04 ppm
!------------------------------------------------------------------------------!

    Control Rod Positions (cm, PWR from the bottom, BWR from the top):          

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:   38.   -   38.   -    -    -   38.   - 
  2:    -  401.   -    -    -   38.   -    - 
  3:   38.   -   38.   -   38.   -   38.
  4:    -    -    -  401.   -    -    - 
  5:    -    -   38.   -   38.   - 
  6:    -   38.   -    -    - 
  7:   38.   -   38.   - 
  8:    -    - 
 

    Control Rod Positions (% of the Covered Reactor Core)           :           

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:   98.   -   98.   -    -    -   98.   - 
  2:    -    0.   -    -    -   98.   -    - 
  3:   98.   -   98.   -   98.   -   98.
  4:    -    -    -    0.   -    -    - 
  5:    -    -   98.   -   98.   - 
  6:    -   98.   -    -    - 
  7:   98.   -   98.   - 
  8:    -    - 
 
  Number of the CR in the Axial Layer:
          0.000  0.000 13.000 13.000 13.000 13.000 13.000 13.000 13.000 13.000
         13.000 13.000 13.000 13.000 13.000 13.000 13.000 14.741
!------------------------------------------------------------------------------!
     POWER DISTRIBUTION            
!------------------------------------------------------------------------------!
     Total Reactor Power                   :  0.69375E-03 MWt
     Average Power Density                 :  0.10308E-03 Wt/cm^3

     NODAL POWER PEAKING FACTORS                                                
     Scaling Factor                        :  0.97008E+04
     Average Value                         :  1.00000

     3D Maximum Value                      :  1.16788
     Reactor Location of 3D Maximum (X,Y,Z):  2  8 16
     Core    Location of 3D Maximum (X,Y,Z):  2  8 15

     3D Minimum Value                      :  0.83208
     Reactor Location of 3D minimum (X,Y,Z):  1  7 14
     Core    Location of 3D minimum (X,Y,Z):  1  7 13

     2D Maximum Value                      :  1.16788
     Reactor Location of 2D Maximum (X,Y)  :  8  1
     Core    Location of 2D Maximum (X,Y)  :  8  1

     2D Minimum Value                      :  0.83266
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  1.01462
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.99969
     Reactor Location of the 1D minimum (Z):       15
     Core    Location of the 1D minimum (Z):       14
!------------------------------------------------------------------------------!

     2D RADIAL ASSEMBLY-AVERAGED POWER DENSITY                                  

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.833 0.997 0.833 0.997 0.859 1.008 0.833 1.168
  2:  0.997 0.859 1.002 0.859 1.002 0.833 1.149 1.168
  3:  0.833 1.002 0.833 1.002 0.833 1.002 1.137
  4:  0.997 0.859 1.002 0.859 1.002 1.155 1.168
  5:  0.859 1.002 0.833 1.002 0.833 1.168
  6:  1.008 0.833 1.002 1.155 1.168
  7:  0.833 1.149 1.137 1.168
  8:  1.168 1.168
 
!------------------------------------------------------------------------------!
     1D AXIAL AVERAGE POWER DENSITY

   1:   1.0146
   2:   0.9997
   3:   0.9997
   4:   0.9997
   5:   0.9997
   6:   0.9997
   7:   0.9997
   8:   0.9997
   9:   0.9997
  10:   0.9997
  11:   0.9997
  12:   0.9997
  13:   0.9997
  14:   0.9997
  15:   0.9997
  16:   0.9997
!------------------------------------------------------------------------------!
            FEEDBACK DISTRIBUTION
!------------------------------------------------------------------------------!

     BORON CONCENTRATION, [ppm]                                                      
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  1.00000

     3D Maximum Value                      :  1.00000
     Reactor Location of 3D Maximum (X,Y,Z):  2  8 14
     Core    Location of 3D Maximum (X,Y,Z):  2  8 13

     3D Minimum Value                      :  1.00000
     Reactor Location of 3D minimum (X,Y,Z):  2  8  2
     Core    Location of 3D minimum (X,Y,Z):  2  8  1

     2D Maximum Value                      :  1.00000
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  1.00000
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  1.00000
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  1.00000
     Reactor Location of the 1D minimum (Z):        3
     Core    Location of the 1D minimum (Z):        2
!------------------------------------------------------------------------------!

     COOLANT TEMPERATURE, [C]                                                        
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.28600

     3D Maximum Value                      :  0.28600
     Reactor Location of 3D Maximum (X,Y,Z):  2  8 16
     Core    Location of 3D Maximum (X,Y,Z):  2  8 15

     3D Minimum Value                      :  0.28600
     Reactor Location of 3D minimum (X,Y,Z):  2  8  2
     Core    Location of 3D minimum (X,Y,Z):  2  8  1

     2D Maximum Value                      :  0.28600
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  0.28600
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.28600
     Reactor Location of the 1D Maximum (Z):        4
     Core    Location of the 1D Maximum (Z):        3

     1D Minimum                            :  0.28600
     Reactor Location of the 1D minimum (Z):        3
     Core    Location of the 1D minimum (Z):        2
!------------------------------------------------------------------------------!

     COOLANT DENSITY, [g/cm^3]                                                       
     Scaling Factor                        :  0.10000E+01
     Average Value                         :  0.75404

     3D Maximum Value                      :  0.75404
     Reactor Location of 3D Maximum (X,Y,Z):  2  8 14
     Core    Location of 3D Maximum (X,Y,Z):  2  8 13

     3D Minimum Value                      :  0.75404
     Reactor Location of 3D minimum (X,Y,Z):  2  8 17
     Core    Location of 3D minimum (X,Y,Z):  2  8 16

     2D Maximum Value                      :  0.75404
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  0.75404
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.75404
     Reactor Location of the 1D Maximum (Z):        4
     Core    Location of the 1D Maximum (Z):        3

     1D Minimum                            :  0.75404
     Reactor Location of the 1D minimum (Z):       15
     Core    Location of the 1D minimum (Z):       14
!------------------------------------------------------------------------------!

     DOPPLER FUEL TEMPERATURE, [K]                                                   
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.55915

     3D Maximum Value                      :  0.55915
     Reactor Location of 3D Maximum (X,Y,Z):  2  8 17
     Core    Location of 3D Maximum (X,Y,Z):  2  8 16

     3D Minimum Value                      :  0.55915
     Reactor Location of 3D minimum (X,Y,Z):  2  2  2
     Core    Location of 3D minimum (X,Y,Z):  2  2  1

     2D Maximum Value                      :  0.55915
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  0.55915
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.55915
     Reactor Location of the 1D Maximum (Z):        3
     Core    Location of the 1D Maximum (Z):        2

     1D Minimum                            :  0.55915
     Reactor Location of the 1D minimum (Z):        4
     Core    Location of the 1D minimum (Z):        3
!------------------------------------------------------------------------------!
            2D FEEDBACK DISTRIBUTION
!------------------------------------------------------------------------------!

     BORON CONCENTRATION, [ppm]                                                      

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  2:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  3:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  4:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  5:  1.000 1.000 1.000 1.000 1.000 1.000
  6:  1.000 1.000 1.000 1.000 1.000
  7:  1.000 1.000 1.000 1.000
  8:  1.000 1.000
 

     COOLANT TEMPERATURE, [C]                                                        

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.286 0.286 0.286 0.286 0.286 0.286 0.286 0.286
  2:  0.286 0.286 0.286 0.286 0.286 0.286 0.286 0.286
  3:  0.286 0.286 0.286 0.286 0.286 0.286 0.286
  4:  0.286 0.286 0.286 0.286 0.286 0.286 0.286
  5:  0.286 0.286 0.286 0.286 0.286 0.286
  6:  0.286 0.286 0.286 0.286 0.286
  7:  0.286 0.286 0.286 0.286
  8:  0.286 0.286
 

     COOLANT DENSITY, [g/cm^3]                                                       

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.754 0.754 0.754 0.754 0.754 0.754 0.754 0.754
  2:  0.754 0.754 0.754 0.754 0.754 0.754 0.754 0.754
  3:  0.754 0.754 0.754 0.754 0.754 0.754 0.754
  4:  0.754 0.754 0.754 0.754 0.754 0.754 0.754
  5:  0.754 0.754 0.754 0.754 0.754 0.754
  6:  0.754 0.754 0.754 0.754 0.754
  7:  0.754 0.754 0.754 0.754
  8:  0.754 0.754
 

     DOPPLER FUEL TEMPERATURE, [K]                                                   

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  2:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  3:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  4:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  5:  0.559 0.559 0.559 0.559 0.559 0.559
  6:  0.559 0.559 0.559 0.559 0.559
  7:  0.559 0.559 0.559 0.559
  8:  0.559 0.559
 
!------------------------------------------------------------------------------!
     1D AXIAL FEEDBACK DISTRIBUTIONS

   N   Br CONC T  COOL RO COOL T  DPLR
   1:   1.0000  0.2860  0.7540  0.5591
   2:   1.0000  0.2860  0.7540  0.5592
   3:   1.0000  0.2860  0.7540  0.5591
   4:   1.0000  0.2860  0.7540  0.5591
   5:   1.0000  0.2860  0.7540  0.5591
   6:   1.0000  0.2860  0.7540  0.5591
   7:   1.0000  0.2860  0.7540  0.5591
   8:   1.0000  0.2860  0.7540  0.5591
   9:   1.0000  0.2860  0.7540  0.5591
  10:   1.0000  0.2860  0.7540  0.5591
  11:   1.0000  0.2860  0.7540  0.5591
  12:   1.0000  0.2860  0.7540  0.5591
  13:   1.0000  0.2860  0.7540  0.5591
  14:   1.0000  0.2860  0.7540  0.5591
  15:   1.0000  0.2860  0.7540  0.5591
  16:   1.0000  0.2860  0.7540  0.5592
!------------------------------------------------------------------------------!
            T/H PARAMETER DISTRIBUTION
            Average Value of Fuel Enthalpy = NaN         [J/Kg]
!------------------------------------------------------------------------------!

     FUEL CENTERLINE TEMPERATURE, [K]                                                
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  NaN    

     3D Maximum Value                      :  0.00000
     Reactor Location of 3D Maximum (X,Y,Z):  0  0  0
     Core    Location of 3D Maximum (X,Y,Z):  1  0 -1

     3D Minimum Value                      : ********
     Reactor Location of 3D minimum (X,Y,Z):  0  0  0
     Core    Location of 3D minimum (X,Y,Z):  1  0 -1

     2D Maximum Value                      :  0.00000
     Reactor Location of 2D Maximum (X,Y)  :  0  0
     Core    Location of 2D Maximum (X,Y)  :  1  0

     2D Minimum Value                      :  0.00000
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.00000
     Reactor Location of the 1D Maximum (Z):        0
     Core    Location of the 1D Maximum (Z):       -1

     1D Minimum                            :  0.00000
     Reactor Location of the 1D minimum (Z):        2
     Core    Location of the 1D minimum (Z):        1
!------------------------------------------------------------------------------!

     CLADDING INNER SURFACE TEMPERATURE, [K]                                         
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  NaN    

     3D Maximum Value                      :  0.00000
     Reactor Location of 3D Maximum (X,Y,Z):  0  0  0
     Core    Location of 3D Maximum (X,Y,Z):  1  0 -1

     3D Minimum Value                      : ********
     Reactor Location of 3D minimum (X,Y,Z):  0  0  0
     Core    Location of 3D minimum (X,Y,Z):  1  0 -1

     2D Maximum Value                      :  0.00000
     Reactor Location of 2D Maximum (X,Y)  :  0  0
     Core    Location of 2D Maximum (X,Y)  :  1  0

     2D Minimum Value                      :  0.00000
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.00000
     Reactor Location of the 1D Maximum (Z):        0
     Core    Location of the 1D Maximum (Z):       -1

     1D Minimum                            :  0.00000
     Reactor Location of the 1D minimum (Z):        2
     Core    Location of the 1D minimum (Z):        1
!------------------------------------------------------------------------------!

     FUEL ENTHALPY, [J/Kg]                                                           
     Scaling Factor                        :  NaN        
     Average Value                         :  NaN    

     3D Maximum Value                      :  NaN    
     Reactor Location of 3D Maximum (X,Y,Z):  0  0  0
     Core    Location of 3D Maximum (X,Y,Z):  1  0 -1

     3D Minimum Value                      :  NaN    
     Reactor Location of 3D minimum (X,Y,Z):  0  0  0
     Core    Location of 3D minimum (X,Y,Z):  1  0 -1

     2D Maximum Value                      :  NaN    
     Reactor Location of 2D Maximum (X,Y)  :  0  0
     Core    Location of 2D Maximum (X,Y)  :  1  0

     2D Minimum Value                      :  NaN    
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  NaN    
     Reactor Location of the 1D Maximum (Z):        0
     Core    Location of the 1D Maximum (Z):       -1

     1D Minimum                            :  NaN    
     Reactor Location of the 1D minimum (Z):        2
     Core    Location of the 1D minimum (Z):        1
!------------------------------------------------------------------------------!
            2D T/H PARAMETER DISTRIBUTION
!------------------------------------------------------------------------------!

     FUEL CENTERLINE TEMPERATURE, [K]                                                

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
  2:  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
  3:  0.000 0.000 0.000 0.000 0.000 0.000 0.000
  4:  0.000 0.000 0.000 0.000 0.000 0.000 0.000
  5:  0.000 0.000 0.000 0.000 0.000 0.000
  6:  0.000 0.000 0.000 0.000 0.000
  7:  0.000 0.000 0.000 0.000
  8:  0.000 0.000
 

     CLADDING INNER SURFACE TEMPERATURE, [K]                                         

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
  2:  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
  3:  0.000 0.000 0.000 0.000 0.000 0.000 0.000
  4:  0.000 0.000 0.000 0.000 0.000 0.000 0.000
  5:  0.000 0.000 0.000 0.000 0.000 0.000
  6:  0.000 0.000 0.000 0.000 0.000
  7:  0.000 0.000 0.000 0.000
  8:  0.000 0.000
 

     FUEL ENTHALPY, [J/Kg]                                                           

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  
  2:  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  
  3:  NaN   NaN   NaN   NaN   NaN   NaN   NaN  
  4:  NaN   NaN   NaN   NaN   NaN   NaN   NaN  
  5:  NaN   NaN   NaN   NaN   NaN   NaN  
  6:  NaN   NaN   NaN   NaN   NaN  
  7:  NaN   NaN   NaN   NaN  
  8:  NaN   NaN  
 
!------------------------------------------------------------------------------!
     1D AXIAL T/H PARAMETER DISTRIBUTIONS 

   N   TF CENT TCL SRF FL ENTH
   1:   0.0000  0.0000 NaN    
   2:   0.0000  0.0000 NaN    
   3:   0.0000  0.0000 NaN    
   4:   0.0000  0.0000 NaN    
   5:   0.0000  0.0000 NaN    
   6:   0.0000  0.0000 NaN    
   7:   0.0000  0.0000 NaN    
   8:   0.0000  0.0000 NaN    
   9:   0.0000  0.0000 NaN    
  10:   0.0000  0.0000 NaN    
  11:   0.0000  0.0000 NaN    
  12:   0.0000  0.0000 NaN    
  13:   0.0000  0.0000 NaN    
  14:   0.0000  0.0000 NaN    
  15:   0.0000  0.0000 NaN    
  16:   0.0000  0.0000 NaN    
!------------------------------------------------------------------------------!
     Temperature Distribution in the Average Fuel Rods, [K]
     Total Power Generated in the Fuel Rods,  [MWt] : 6.80569E-04
!------------------------------------------------------------------------------!
     Fuel Assembly Type:    1
 Fuel Rod Number    1 of    1
 Power of the FRD [MWt] =  6.80569E-04
!------------------------------------------------------------------------------!
          1      2      3      4      5      6      7      8      9     10     11     12     13     14
      POWER   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   CLAD   CLAD   CLAD   COOL
  2: 1.0146    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  3: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  4: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  5: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  6: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  7: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  8: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  9: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 10: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 11: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 12: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 13: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 14: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 15: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 16: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
 17: 0.9997    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
            NEUTRON FLUX
!------------------------------------------------------------------------------!
     Average Neutron Flux (1/cm^3)         : 
     0.10000E+01 0.10000E+01

     Neutron Flux [1/(s cm^2)], group   1                                       
     Scaling Factor                        :  0.10000E+01
     Average Value                         :  1.00000

     3D Maximum Value                      :  1.00000
     Reactor Location of 3D Maximum (X,Y,Z):  1  1  1
     Core    Location of 3D Maximum (X,Y,Z):  1  1  0

     3D Minimum Value                      :  1.00000
     Reactor Location of 3D minimum (X,Y,Z):  3  9  2
     Core    Location of 3D minimum (X,Y,Z):  4  9  1

     2D Maximum Value                      :  1.00000
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  1.00000
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  1.00000
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  1.00000
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!

     Neutron Flux [1/(s cm^2)], group   2                                       
     Scaling Factor                        :  0.10000E+01
     Average Value                         :  1.00000

     3D Maximum Value                      :  1.00000
     Reactor Location of 3D Maximum (X,Y,Z):  1  1  1
     Core    Location of 3D Maximum (X,Y,Z):  1  1  0

     3D Minimum Value                      :  1.00000
     Reactor Location of 3D minimum (X,Y,Z):  3  9  2
     Core    Location of 3D minimum (X,Y,Z):  4  9  1

     2D Maximum Value                      :  1.00000
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  1.00000
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  1.00000
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  1.00000
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!
            2D NEUTRON FLUX DISTRIBUTION
!------------------------------------------------------------------------------!

     Neutron Flux [1/(s cm^2)], group   1                                       

          1     2     3     4     5     6     7     8     9
      ----- ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  2:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  3:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  4:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  5:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  6:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  7:  1.000 1.000 1.000 1.000 1.000 1.000
  8:  1.000 1.000 1.000 1.000 1.000
  9:  1.000 1.000 1.000
 

     Neutron Flux [1/(s cm^2)], group   2                                       

          1     2     3     4     5     6     7     8     9
      ----- ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  2:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  3:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  4:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  5:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  6:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  7:  1.000 1.000 1.000 1.000 1.000 1.000
  8:  1.000 1.000 1.000 1.000 1.000
  9:  1.000 1.000 1.000
 
!------------------------------------------------------------------------------!
     1D AXIAL NEUTRON FLUX DISTRIBUTIONS

  NZ       1     2
   1:  1.000 1.000
   2:  1.000 1.000
   3:  1.000 1.000
   4:  1.000 1.000
   5:  1.000 1.000
   6:  1.000 1.000
   7:  1.000 1.000
   8:  1.000 1.000
   9:  1.000 1.000
  10:  1.000 1.000
  11:  1.000 1.000
  12:  1.000 1.000
  13:  1.000 1.000
  14:  1.000 1.000
  15:  1.000 1.000
  16:  1.000 1.000
  17:  1.000 1.000
  18:  1.000 1.000
!------------------------------------------------------------------------------!
!==============================================================================!
!                                   SKETCH-N                                   !
!  Nodal Neutron Diffusion Code for Solving Steady-State and Kinetics Problems !
!              Version 1.0 (c) Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !
!      Date & Time of Calculation :  1 January  2002   5:04:49 pm              !
!==============================================================================!
                 NUMERICAL RESULTS
                 STEADY-STATE CALCULATIONS
!------------------------------------------------------------------------------!
     k_ef                                  :  1.000002
     Boron Concentration                   :  0.56260E+03 ppm
!------------------------------------------------------------------------------!

    Control Rod Positions (cm, PWR from the bottom, BWR from the top):          

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:   38.   -   38.   -    -    -   38.   - 
  2:    -  401.   -    -    -   38.   -    - 
  3:   38.   -   38.   -   38.   -   38.
  4:    -    -    -  401.   -    -    - 
  5:    -    -   38.   -   38.   - 
  6:    -   38.   -    -    - 
  7:   38.   -   38.   - 
  8:    -    - 
 

    Control Rod Positions (% of the Covered Reactor Core)           :           

         1    2    3    4    5    6    7    8
      ---- ---- ---- ---- ---- ---- ---- ----
  1:   98.   -   98.   -    -    -   98.   - 
  2:    -    0.   -    -    -   98.   -    - 
  3:   98.   -   98.   -   98.   -   98.
  4:    -    -    -    0.   -    -    - 
  5:    -    -   98.   -   98.   - 
  6:    -   98.   -    -    - 
  7:   98.   -   98.   - 
  8:    -    - 
 
  Number of the CR in the Axial Layer:
          0.000  0.000 13.000 13.000 13.000 13.000 13.000 13.000 13.000 13.000
         13.000 13.000 13.000 13.000 13.000 13.000 13.000 14.741
!------------------------------------------------------------------------------!
     POWER DISTRIBUTION            
!------------------------------------------------------------------------------!
     Total Reactor Power                   :  0.69375E-03 MWt
     Average Power Density                 :  0.10308E-03 Wt/cm^3

     NODAL POWER PEAKING FACTORS                                                
     Scaling Factor                        :  0.97008E+04
     Average Value                         :  1.00000

     3D Maximum Value                      :  2.87559
     Reactor Location of 3D Maximum (X,Y,Z):  4  2  9
     Core    Location of 3D Maximum (X,Y,Z):  4  2  8

     3D Minimum Value                      :  0.05746
     Reactor Location of 3D minimum (X,Y,Z):  7  3 17
     Core    Location of 3D minimum (X,Y,Z):  7  3 16

     2D Maximum Value                      :  1.91097
     Reactor Location of 2D Maximum (X,Y)  :  2  4
     Core    Location of 2D Maximum (X,Y)  :  2  4

     2D Minimum Value                      :  0.38234
     Reactor Location of 2D minimum (X,Y)  :  7  3
     Core    Location of 2D minimum (X,Y)  :  7  3

     1D Maximum                            :  1.50284
     Reactor Location of the 1D Maximum (Z):        9
     Core    Location of the 1D Maximum (Z):        8

     1D Minimum                            :  0.14032
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!

     2D RADIAL ASSEMBLY-AVERAGED POWER DENSITY                                  

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.993 1.521 1.071 1.806 1.815 1.007 0.396 0.540
  2:  1.521 1.845 1.661 1.911 1.424 0.559 0.560 0.428
  3:  1.071 1.661 1.015 1.448 0.730 0.726 0.382
  4:  1.806 1.911 1.448 1.435 1.016 0.945 0.541
  5:  1.815 1.424 0.730 1.016 0.542 0.671
  6:  1.007 0.559 0.726 0.945 0.671
  7:  0.396 0.560 0.382 0.541
  8:  0.540 0.428
 
!------------------------------------------------------------------------------!
     1D AXIAL AVERAGE POWER DENSITY

   1:   0.1979
   2:   0.2667
   3:   0.4216
   4:   0.6818
   5:   0.9920
   6:   1.2417
   7:   1.4155
   8:   1.5028
   9:   1.4984
  10:   1.4025
  11:   1.2210
  12:   0.9650
  13:   0.6501
  14:   0.4011
  15:   0.2453
  16:   0.1403
!------------------------------------------------------------------------------!
            FEEDBACK DISTRIBUTION
!------------------------------------------------------------------------------!

     BORON CONCENTRATION, [ppm]                                                      
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.56260

     3D Maximum Value                      :  0.56260
     Reactor Location of 3D Maximum (X,Y,Z):  2  8  3
     Core    Location of 3D Maximum (X,Y,Z):  2  8  2

     3D Minimum Value                      :  0.56260
     Reactor Location of 3D minimum (X,Y,Z):  2  2  2
     Core    Location of 3D minimum (X,Y,Z):  2  2  1

     2D Maximum Value                      :  0.56260
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  0.56260
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.56260
     Reactor Location of the 1D Maximum (Z):        4
     Core    Location of the 1D Maximum (Z):        3

     1D Minimum                            :  0.56260
     Reactor Location of the 1D minimum (Z):        3
     Core    Location of the 1D minimum (Z):        2
!------------------------------------------------------------------------------!

     COOLANT TEMPERATURE, [C]                                                        
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.28600

     3D Maximum Value                      :  0.28601
     Reactor Location of 3D Maximum (X,Y,Z):  1  8  2
     Core    Location of 3D Maximum (X,Y,Z):  1  8  1

     3D Minimum Value                      :  0.28600
     Reactor Location of 3D minimum (X,Y,Z):  1  1 17
     Core    Location of 3D minimum (X,Y,Z):  1  1 16

     2D Maximum Value                      :  0.28600
     Reactor Location of 2D Maximum (X,Y)  :  2  2
     Core    Location of 2D Maximum (X,Y)  :  2  2

     2D Minimum Value                      :  0.28600
     Reactor Location of 2D minimum (X,Y)  :  7  2
     Core    Location of 2D minimum (X,Y)  :  7  2

     1D Maximum                            :  0.28601
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.28600
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!

     COOLANT DENSITY, [g/cm^3]                                                       
     Scaling Factor                        :  0.10000E+01
     Average Value                         :  0.75440

     3D Maximum Value                      :  0.75443
     Reactor Location of 3D Maximum (X,Y,Z):  1  8  2
     Core    Location of 3D Maximum (X,Y,Z):  1  8  1

     3D Minimum Value                      :  0.75436
     Reactor Location of 3D minimum (X,Y,Z):  2  2 17
     Core    Location of 3D minimum (X,Y,Z):  2  2 16

     2D Maximum Value                      :  0.75440
     Reactor Location of 2D Maximum (X,Y)  :  1  1
     Core    Location of 2D Maximum (X,Y)  :  1  1

     2D Minimum Value                      :  0.75440
     Reactor Location of 2D minimum (X,Y)  :  2  2
     Core    Location of 2D minimum (X,Y)  :  2  2

     1D Maximum                            :  0.75443
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.75436
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!

     DOPPLER FUEL TEMPERATURE, [K]                                                   
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.55915

     3D Maximum Value                      :  0.55916
     Reactor Location of 3D Maximum (X,Y,Z):  7  1  3
     Core    Location of 3D Maximum (X,Y,Z):  7  1  2

     3D Minimum Value                      :  0.55915
     Reactor Location of 3D minimum (X,Y,Z):  2  2 17
     Core    Location of 3D minimum (X,Y,Z):  2  2 16

     2D Maximum Value                      :  0.55915
     Reactor Location of 2D Maximum (X,Y)  :  5  1
     Core    Location of 2D Maximum (X,Y)  :  5  1

     2D Minimum Value                      :  0.55915
     Reactor Location of 2D minimum (X,Y)  :  6  2
     Core    Location of 2D minimum (X,Y)  :  6  2

     1D Maximum                            :  0.55916
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.55915
     Reactor Location of the 1D minimum (Z):       17
     Core    Location of the 1D minimum (Z):       16
!------------------------------------------------------------------------------!
            2D FEEDBACK DISTRIBUTION
!------------------------------------------------------------------------------!

     BORON CONCENTRATION, [ppm]                                                      

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.563 0.563 0.563 0.563 0.563 0.563 0.563 0.563
  2:  0.563 0.563 0.563 0.563 0.563 0.563 0.563 0.563
  3:  0.563 0.563 0.563 0.563 0.563 0.563 0.563
  4:  0.563 0.563 0.563 0.563 0.563 0.563 0.563
  5:  0.563 0.563 0.563 0.563 0.563 0.563
  6:  0.563 0.563 0.563 0.563 0.563
  7:  0.563 0.563 0.563 0.563
  8:  0.563 0.563
 

     COOLANT TEMPERATURE, [C]                                                        

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.286 0.286 0.286 0.286 0.286 0.286 0.286 0.286
  2:  0.286 0.286 0.286 0.286 0.286 0.286 0.286 0.286
  3:  0.286 0.286 0.286 0.286 0.286 0.286 0.286
  4:  0.286 0.286 0.286 0.286 0.286 0.286 0.286
  5:  0.286 0.286 0.286 0.286 0.286 0.286
  6:  0.286 0.286 0.286 0.286 0.286
  7:  0.286 0.286 0.286 0.286
  8:  0.286 0.286
 

     COOLANT DENSITY, [g/cm^3]                                                       

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.754 0.754 0.754 0.754 0.754 0.754 0.754 0.754
  2:  0.754 0.754 0.754 0.754 0.754 0.754 0.754 0.754
  3:  0.754 0.754 0.754 0.754 0.754 0.754 0.754
  4:  0.754 0.754 0.754 0.754 0.754 0.754 0.754
  5:  0.754 0.754 0.754 0.754 0.754 0.754
  6:  0.754 0.754 0.754 0.754 0.754
  7:  0.754 0.754 0.754 0.754
  8:  0.754 0.754
 

     DOPPLER FUEL TEMPERATURE, [K]                                                   

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  2:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  3:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  4:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  5:  0.559 0.559 0.559 0.559 0.559 0.559
  6:  0.559 0.559 0.559 0.559 0.559
  7:  0.559 0.559 0.559 0.559
  8:  0.559 0.559
 
!------------------------------------------------------------------------------!
     1D AXIAL FEEDBACK DISTRIBUTIONS

   N   Br CONC T  COOL RO COOL T  DPLR
   1:   0.5626  0.2860  0.7544  0.5592
   2:   0.5626  0.2860  0.7544  0.5592
   3:   0.5626  0.2860  0.7544  0.5592
   4:   0.5626  0.2860  0.7544  0.5592
   5:   0.5626  0.2860  0.7544  0.5592
   6:   0.5626  0.2860  0.7544  0.5592
   7:   0.5626  0.2860  0.7544  0.5592
   8:   0.5626  0.2860  0.7544  0.5592
   9:   0.5626  0.2860  0.7544  0.5592
  10:   0.5626  0.2860  0.7544  0.5592
  11:   0.5626  0.2860  0.7544  0.5592
  12:   0.5626  0.2860  0.7544  0.5592
  13:   0.5626  0.2860  0.7544  0.5592
  14:   0.5626  0.2860  0.7544  0.5592
  15:   0.5626  0.2860  0.7544  0.5592
  16:   0.5626  0.2860  0.7544  0.5592
!------------------------------------------------------------------------------!
            T/H PARAMETER DISTRIBUTION
            Average Value of Fuel Enthalpy =0.125874E+06 [J/Kg]
!------------------------------------------------------------------------------!

     FUEL CENTERLINE TEMPERATURE, [K]                                                
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.55915

     3D Maximum Value                      :  0.55916
     Reactor Location of 3D Maximum (X,Y,Z):  7  1  3
     Core    Location of 3D Maximum (X,Y,Z):  7  1  2

     3D Minimum Value                      :  0.55915
     Reactor Location of 3D minimum (X,Y,Z):  5  4 16
     Core    Location of 3D minimum (X,Y,Z):  5  4 15

     2D Maximum Value                      :  0.55915
     Reactor Location of 2D Maximum (X,Y)  :  5  1
     Core    Location of 2D Maximum (X,Y)  :  5  1

     2D Minimum Value                      :  0.55915
     Reactor Location of 2D minimum (X,Y)  :  6  2
     Core    Location of 2D minimum (X,Y)  :  6  2

     1D Maximum                            :  0.55916
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.55915
     Reactor Location of the 1D minimum (Z):       16
     Core    Location of the 1D minimum (Z):       15
!------------------------------------------------------------------------------!

     CLADDING INNER SURFACE TEMPERATURE, [K]                                         
     Scaling Factor                        :  0.10000E-02
     Average Value                         :  0.55915

     3D Maximum Value                      :  0.55916
     Reactor Location of 3D Maximum (X,Y,Z):  1  1  2
     Core    Location of 3D Maximum (X,Y,Z):  1  1  1

     3D Minimum Value                      :  0.55915
     Reactor Location of 3D minimum (X,Y,Z):  1  1 17
     Core    Location of 3D minimum (X,Y,Z):  1  1 16

     2D Maximum Value                      :  0.55915
     Reactor Location of 2D Maximum (X,Y)  :  2  2
     Core    Location of 2D Maximum (X,Y)  :  2  2

     2D Minimum Value                      :  0.55915
     Reactor Location of 2D minimum (X,Y)  :  1  1
     Core    Location of 2D minimum (X,Y)  :  1  1

     1D Maximum                            :  0.55916
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.55915
     Reactor Location of the 1D minimum (Z):       16
     Core    Location of the 1D minimum (Z):       15
!------------------------------------------------------------------------------!

     FUEL ENTHALPY, [J/Kg]                                                           
     Scaling Factor                        :  0.79444E-05
     Average Value                         :  1.00000

     3D Maximum Value                      :  1.00001
     Reactor Location of 3D Maximum (X,Y,Z):  7  1  3
     Core    Location of 3D Maximum (X,Y,Z):  7  1  2

     3D Minimum Value                      :  0.99999
     Reactor Location of 3D minimum (X,Y,Z):  5  4 16
     Core    Location of 3D minimum (X,Y,Z):  5  4 15

     2D Maximum Value                      :  1.00000
     Reactor Location of 2D Maximum (X,Y)  :  5  1
     Core    Location of 2D Maximum (X,Y)  :  5  1

     2D Minimum Value                      :  1.00000
     Reactor Location of 2D minimum (X,Y)  :  6  2
     Core    Location of 2D minimum (X,Y)  :  6  2

     1D Maximum                            :  1.00001
     Reactor Location of the 1D Maximum (Z):        2
     Core    Location of the 1D Maximum (Z):        1

     1D Minimum                            :  0.99999
     Reactor Location of the 1D minimum (Z):       16
     Core    Location of the 1D minimum (Z):       15
!------------------------------------------------------------------------------!
            2D T/H PARAMETER DISTRIBUTION
!------------------------------------------------------------------------------!

     FUEL CENTERLINE TEMPERATURE, [K]                                                

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  2:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  3:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  4:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  5:  0.559 0.559 0.559 0.559 0.559 0.559
  6:  0.559 0.559 0.559 0.559 0.559
  7:  0.559 0.559 0.559 0.559
  8:  0.559 0.559
 

     CLADDING INNER SURFACE TEMPERATURE, [K]                                         

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  2:  0.559 0.559 0.559 0.559 0.559 0.559 0.559 0.559
  3:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  4:  0.559 0.559 0.559 0.559 0.559 0.559 0.559
  5:  0.559 0.559 0.559 0.559 0.559 0.559
  6:  0.559 0.559 0.559 0.559 0.559
  7:  0.559 0.559 0.559 0.559
  8:  0.559 0.559
 

     FUEL ENTHALPY, [J/Kg]                                                           

          1     2     3     4     5     6     7     8
      ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  2:  1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
  3:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  4:  1.000 1.000 1.000 1.000 1.000 1.000 1.000
  5:  1.000 1.000 1.000 1.000 1.000 1.000
  6:  1.000 1.000 1.000 1.000 1.000
  7:  1.000 1.000 1.000 1.000
  8:  1.000 1.000
 
!------------------------------------------------------------------------------!
     1D AXIAL T/H PARAMETER DISTRIBUTIONS 

   N   TF CENT TCL SRF FL ENTH
   1:   0.5592  0.5592  1.0000
   2:   0.5592  0.5592  1.0000
   3:   0.5592  0.5592  1.0000
   4:   0.5592  0.5592  1.0000
   5:   0.5592  0.5592  1.0000
   6:   0.5592  0.5592  1.0000
   7:   0.5592  0.5592  1.0000
   8:   0.5592  0.5592  1.0000
   9:   0.5592  0.5592  1.0000
  10:   0.5592  0.5592  1.0000
  11:   0.5592  0.5592  1.0000
  12:   0.5592  0.5592  1.0000
  13:   0.5592  0.5592  1.0000
  14:   0.5592  0.5592  1.0000
  15:   0.5592  0.5592  1.0000
  16:   0.5592  0.5592  1.0000
!------------------------------------------------------------------------------!
     Temperature Distribution in the Average Fuel Rods, [K]
     Total Power Generated in the Fuel Rods,  [MWt] : 6.80569E-04
!------------------------------------------------------------------------------!
     Fuel Assembly Type:    1
 Fuel Rod Number    1 of    1
 Power of the FRD [MWt] =  6.80569E-04
!------------------------------------------------------------------------------!
          1      2      3      4      5      6      7      8      9     10     11     12     13     14
      POWER   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   FUEL   CLAD   CLAD   CLAD   COOL
  2: 0.1979  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  3: 0.2667  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  4: 0.4216  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  5: 0.6818  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  6: 0.9920  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  7: 1.2417  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  8: 1.4155  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
  9: 1.5028  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 10: 1.4984  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 11: 1.4025  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 12: 1.2210  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 13: 0.9650  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 14: 0.6501  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 15: 0.4011  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 16: 0.2453  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2    0.0
 17: 0.1403  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.2  559.1  559.1  559.2    0.0
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
            NEUTRON FLUX
!------------------------------------------------------------------------------!
     Average Neutron Flux (1/cm^3)         : 
     0.17569E+09 0.41601E+08

     Neutron Flux [1/(s cm^2)], group   1                                       
     Scaling Factor                        :  0.56918E-08
     Average Value                         :  1.00000

     3D Maximum Value                      :  4.29002
     Reactor Location of 3D Maximum (X,Y,Z):  1  4  9
     Core    Location of 3D Maximum (X,Y,Z):  1  4  8

     3D Minimum Value                      :  0.00072
     Reactor Location of 3D minimum (X,Y,Z):  3  9 18
     Core    Location of 3D minimum (X,Y,Z):  4  9 17

     2D Maximum Value                      :  2.45688
     Reactor Location of 2D Maximum (X,Y)  :  1  4
     Core    Location of 2D Maximum (X,Y)  :  1  4

     2D Minimum Value                      :  0.01507
     Reactor Location of 2D minimum (X,Y)  :  9  3
     Core    Location of 2D minimum (X,Y)  :  9  3

     1D Maximum                            :  1.74462
     Reactor Location of the 1D Maximum (Z):        9
     Core    Location of the 1D Maximum (Z):        8

     1D Minimum                            :  0.02502
     Reactor Location of the 1D minimum (Z):       18
     Core    Location of the 1D minimum (Z):       17
!------------------------------------------------------------------------------!

     Neutron Flux [1/(s cm^2)], group   2                                       
     Scaling Factor                        :  0.24038E-07
     Average Value                         :  1.00000

     3D Maximum Value                      :  5.23094
     Reactor Location of 3D Maximum (X,Y,Z):  2  4  9
     Core    Location of 3D Maximum (X,Y,Z):  2  4  8

     3D Minimum Value                      :  0.00037
     Reactor Location of 3D minimum (X,Y,Z):  3  9 18
     Core    Location of 3D minimum (X,Y,Z):  4  9 17

     2D Maximum Value                      :  3.04796
     Reactor Location of 2D Maximum (X,Y)  :  2  4
     Core    Location of 2D Maximum (X,Y)  :  2  4

     2D Minimum Value                      :  0.00754
     Reactor Location of 2D minimum (X,Y)  :  9  3
     Core    Location of 2D minimum (X,Y)  :  9  3

     1D Maximum                            :  1.70212
     Reactor Location of the 1D Maximum (Z):        9
     Core    Location of the 1D Maximum (Z):        8

     1D Minimum                            :  0.15985
     Reactor Location of the 1D minimum (Z):       18
     Core    Location of the 1D minimum (Z):       17
!------------------------------------------------------------------------------!
            2D NEUTRON FLUX DISTRIBUTION
!------------------------------------------------------------------------------!

     Neutron Flux [1/(s cm^2)], group   1                                       

          1     2     3     4     5     6     7     8     9
      ----- ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.792 2.104 1.925 2.457 2.304 1.314 0.715 0.605 0.092
  2:  2.104 2.354 2.220 2.432 1.891 1.005 0.726 0.484 0.068
  3:  1.925 2.220 1.818 1.934 1.310 1.008 0.584 0.145 0.015
  4:  2.457 2.432 1.934 1.826 1.375 1.159 0.620 0.077
  5:  2.304 1.891 1.310 1.375 0.970 0.766 0.198 0.019
  6:  1.314 1.005 1.008 1.159 0.766 0.193 0.026
  7:  0.715 0.726 0.584 0.620 0.198 0.026
  8:  0.605 0.484 0.145 0.077 0.019
  9:  0.092 0.068 0.015
 

     Neutron Flux [1/(s cm^2)], group   2                                       

          1     2     3     4     5     6     7     8     9
      ----- ----- ----- ----- ----- ----- ----- ----- -----
  1:  1.549 2.029 1.669 2.410 2.893 1.342 0.617 0.622 0.051
  2:  2.029 2.946 2.215 3.048 1.898 0.872 0.643 0.491 0.038
  3:  1.669 2.215 1.582 1.930 1.139 0.963 0.429 0.070 0.008
  4:  2.410 3.048 1.930 2.291 1.351 1.085 0.620 0.043
  5:  2.893 1.898 1.139 1.351 0.844 0.770 0.100 0.010
  6:  1.342 0.872 0.963 1.085 0.770 0.097 0.013
  7:  0.617 0.643 0.429 0.620 0.100 0.013
  8:  0.622 0.491 0.070 0.043 0.010
  9:  0.051 0.038 0.008
 
!------------------------------------------------------------------------------!
     1D AXIAL NEUTRON FLUX DISTRIBUTIONS

  NZ       1     2
   1:  0.032 0.212
   2:  0.162 0.235
   3:  0.306 0.303
   4:  0.490 0.477
   5:  0.791 0.772
   6:  1.152 1.124
   7:  1.441 1.406
   8:  1.643 1.603
   9:  1.745 1.702
  10:  1.739 1.697
  11:  1.628 1.588
  12:  1.417 1.383
  13:  1.120 1.093
  14:  0.755 0.736
  15:  0.466 0.454
  16:  0.284 0.278
  17:  0.132 0.166
  18:  0.025 0.160
!------------------------------------------------------------------------------!
 Number of Thermal-Hydraulics Calculations =      18
 Number of Nonlinear Iterations            =      18
 Number of Source Iteration per Nonlinear  =       8
 Number of SSOR iterations  per Outer      =       2
 n_inter =  2 e_inter =.10000E-07
 Convergence  of k_eff (actual, tolerance) =  0.89109E-05 0.10000E-04
 Convergence  of flux  (actual, tolerance) =  0.16542E-05 0.10000E-04
!------------------------------------------------------------------------------!
 SKETCH has used   7.661016    seconds of CPU time
    This Include: 
       Internal  T/H Module        2.743946     s    35.81700     %
       Nodal Coefficient Update:   2.433499     s    31.76471     %
       CMFD Iterations             2.103025     s    27.45099     %
       XS Updates                 0.1502146     s    1.960767     %
       Input/Output Module        9.0129793E-02 s    1.176473     %
 CHECKING SUM   7.681045     seconds of CPU time
 This include user time, system time =    7.570887      0.1101584    
!------------------------------------------------------------------------------!

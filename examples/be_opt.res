
----------------------------------------------------------------------------------------------------------
                                                  PROGRAM CHAMP
                                                 version 3.08.00
----------------------------------------------------------------------------------------------------------

<<<<<<< HEAD
GIT commit 54040b0797ba5c74d120eeeff62ccd12b5d8f74a, Date:   Mon Jul 4 11:39:02 2016 -0400
Compiled by mussard on Tue Jul 12 12:39:21 CEST 2016 on host antares
Executed by mussard on 2016-07-12 12:40:03 (+0200) on master host 
=======
GIT commit 112fa4cd20d6c214eb733113b751c1d415c9a449, Date:   Tue Jul 5 16:55:14 2016 -0400
Compiled by mussard on Wed Jul  6 09:19:49 EDT 2016 on host antares
Executed by mussard on 2016-07-06 09:23:31 (-0400) on master host 
>>>>>>> master
Executable: /home/mussard/softwares/champ/qmc/champ.exe
Command line arguments: -p -m vmc_mov1
Run in mode >vmc_mov1<


Beginning of control menu --------------------------------------------------------------------------------
 type of system: finite system
 trial energy =   -14.620000

 random number seed = 1837465927472523
 number of steps per block              =         10
 number of blocks after equilibration   =         10
 number of blocks before equilibration  =          1
 number of configurations saved         =     0
 number of walkers initialized to        1

 version of Metropolis  =          6
 step size              =    1.00000
 radial step multiplier =    5.00000
 cos(theta) step size   =    1.00000
 force bias             =    1.00000
End of control menu --------------------------------------------------------------------------------------

Beginning of nuclei menu ---------------------------------------------------------------------------------
 system of dimension 3
 number of atomic center types =   1
 number of atomic centers =     1
 geometry:
                  type  charge            cartesian coordinates
 nucleus #    1:    1     4.0        0.000000    0.000000    0.000000

 type of external potential: nloc=    0
 Nuclear potential energy =     0.0000000
End of nuclei menu ---------------------------------------------------------------------------------------

Beginning of wavefunction menu ---------------------------------------------------------------------------
 numbers of total     electrons =     4
 numbers of spin-up   electrons =     2
 numbers of spin-down electrons =     2
End of wavefunction menu ---------------------------------------------------------------------------------

Beginning of basis menu ----------------------------------------------------------------------------------
 use localized basis functions (ibasis=1)
 numr=   -3
 analytical radial basis functions used
 the analytical basis functions are Slater functions

 number of basis functions =     8
 number of analytical basis functions =     8

 center type #     1
 function #     1:     1s       exponent =     6.285179
 function #     2:     1s       exponent =     3.455497
 function #     3:     2s       exponent =     2.774117
 function #     4:     2s       exponent =     1.192734
 function #     5:     2s       exponent =     0.824535
 function #     6:    2px       exponent =     0.986656
 function #     7:    2py       exponent =     0.986656
 function #     8:    2pz       exponent =     0.986656
 not using recursion for spherical harmonics
End of basis menu ----------------------------------------------------------------------------------------

Beginning of orbitals menu -------------------------------------------------------------------------------
 number of orbitals =     8
 orbital coefficients:
  0.093369  0.915836  0.005504 -0.006041 -0.006733  0.000000  0.000000  0.000000
 -0.014480 -0.163225 -0.083169  0.624494  0.459069  0.000000  0.000000  0.000000
  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  1.000000
  0.000000  0.000000  0.000000  0.000000  0.000000  1.000000  0.000000  0.000000
  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  1.000000  0.000000
 -0.143342  0.128868 -0.647442  3.321945 -2.986689  0.000000  0.000000  0.000000
  0.916555 -2.721666  3.370485 -3.005751  1.736059  0.000000  0.000000  0.000000
 -4.045617  5.658976 -3.236422  1.782572 -0.988422  0.000000  0.000000  0.000000
End of orbitals menu -------------------------------------------------------------------------------------

Beginning of csfs menu -----------------------------------------------------------------------------------
Normalized CSF coefs=  0.982296 -0.187337
Initial CSF rotation coefs=  0.188450
 number of determinants =     4
 determinants have orbitals (spin-up | spin-down):
 det #     1:    1   2  |   1   2
 det #     2:    1   4  |   1   4
 det #     3:    1   5  |   1   5
 det #     4:    1   3  |   1   3
 number of unique spin-up   determinants =     4
 number of unique spin-down determinants =     4
 unique spin-up determinants have orbitals:
 det #     1:    1   2
 det #     2:    1   4
 det #     3:    1   5
 det #     4:    1   3
 unique spin-down determinants have orbitals:
 det #     1:    1   2
 det #     2:    1   4
 det #     3:    1   5
 det #     4:    1   3

 number of CSFs =     2
 CSF coefficients =   0.982296 -0.187337
 sum of square of CSF coefficients =   1.000000
 CSF #     1
 determinants in CSF:   1
 coefficients: 1.00000
 CSF #     2
 determinants in CSF:   2   3   4
 coefficients: 1.00000 1.00000 1.00000

 Orbital occupation information:
 There are   8 total    orbitals
 There are   5 occupied orbitals of indexes:   1   2   3   4   5
 There are   1 closed   orbitals of indexes:   1
 There are   4 active   orbitals of indexes:   2   3   4   5
 There are   0 open     orbitals of indexes:
 There are   3 virtual  orbitals of indexes:   6   7   8
 Last occupied orbital has index   5
 Number of computed orbitals initialized to        5
End of csfs menu -----------------------------------------------------------------------------------------

Beginning of jastrow menu --------------------------------------------------------------------------------
 type of Jastrow factor: ijas=   4
 type of scaling function: isc=   4
 starting spin index: nspin1=   1
 ending   spin index: nspin2=   1
 order of polynomial: nord=   5
 Fock terms: ifock=   0
 analytic Laplacian: ianalyt_lap=   1
 order of e-n polynomial: norda=    5
 order of e-e polynomial: nordb=    5
 order of e-e-n polynomial: nordc=    5
 scale factor: scalek=   0.80000
 a21=   0.00000
 e-n terms: a=  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
 e-e terms: b=  0.500000  1.000000  0.000000  0.000000  0.000000  0.000000
 e-e-n terms: c=  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
End of jastrow menu --------------------------------------------------------------------------------------

Beginning of optimization menu ---------------------------------------------------------------------------
 update of nonlinear parameters in linear optimization method will be done using semiorthogonal derivatives:
 the derivatives will be orthogonalized to [xi Psi_0 +(1-xi) Psi_lin], with xi= 1.00000000E+00
 maximum number of iterations =    1
 stabilization constant: add_diag=  1.0E-08
 fraction of variance: p_var=      0.0000
 energy threshold for convergence =  1.0000E-03
 Requested parameter types: jastrow   csfs      orbitals  exponents 

 Orbital optimization information:
 There are   8 orbitals in the optimization space of indexes:   1   2   3   4   5   6   7   8
 Constructing list of single orbital excitations based on orbital occupancies alone...
 Number of single orbital excitations =          7
 Number of orbital derivatives        =          7
 Number of orbitals involved          =          5
 Constructing list of single orbital excitations checking for additional redundancies...
 Number of single orbital excitations =          7
 Number of orbital derivatives        =          7
 Number of orbitals involved          =          5
 Number of unique spin-up   excited determinants =         18
 Number of unique spin-down excited determinants =         18
 Number of computed orbitals will be        8

 Exponent optimization information:
 Number of exponent parameters =       6
 Exponent parameter #   1 corresponds to basis functions:   1_1s
 Exponent parameter #   2 corresponds to basis functions:   2_1s
 Exponent parameter #   3 corresponds to basis functions:   3_2s
 Exponent parameter #   4 corresponds to basis functions:   4_2s
 Exponent parameter #   5 corresponds to basis functions:   5_2s
 Exponent parameter #   6 corresponds to basis functions:   6_2px   7_2py   8_2pz

 Number of Jastrow parameters:      24
 Number of periodic Jastrow parameters:     0
<<<<<<< HEAD
 Number of CSF parameters:           3
 Number of orbital parameters:       7
 Number of exponent parameters:      6
 Number of geometry parameters:      0
 Total number of parameters:        40
=======
 Number of CSF parameters:           1
 Number of orbital parameters:       7
 Number of exponent parameters:      6
 Number of geometry parameters:      0
 Total number of parameters:        38
>>>>>>> master

End of optimization menu ---------------------------------------------------------------------------------

*************************************** WAVE FUNCTION OPTIMIZATION ***************************************
1 configuration for  4 electrons has been generated by routine sites.
sites:   -0.313319    0.628839   -0.208193   -0.881718   -0.193128    0.281660    0.605008   -0.230834   -0.536871   -0.024506   -0.056022   -0.310910

Optimization will be done with the linear method.

stabilization turned off in optimization.
<<<<<<< HEAD
OPT: optimization of   24 Jastrow,    3 CSF,      7 orbital,    6 exponent and     0 geometry parameters with linear method:
=======
OPT: optimization of   24 Jastrow,    1 CSF,      7 orbital,    6 exponent and     0 geometry parameters with linear method:
>>>>>>> master

Beginning optimization iteration #   1

***************************************** START VMC CALCULATION ******************************************

The following averages will be calculated:
- eloc_av                                            = average of eloc
- dpsi_av                                            = average of dpsi
- dpsi_eloc_av                                       = average of dpsi_eloc
- dpsi_dpsi_av                                       = average of dpsi_dpsi
- deloc_av                                           = average of deloc
- dpsi_deloc_av                                      = average of dpsi_deloc
- dpsi_dpsi_eloc_av                                  = average of dpsi_dpsi_eloc
- eloc_sq_av                                         = average of eloc_sq

The following variances will be calculated:
- object 1328                                        = variance of average gradient_norm
- eloc_av_var                                        = variance of average eloc_av
- object   42                                        = variance of average sigma

The following statistical errors will be calculated:
- gradient_norm_err                                  = statistical error of average gradient_norm
- error_sigma                                        = statistical error of average sigma

<<<<<<< HEAD
Beginning of equilibration (total CPU time is       0.02 s, CPU time since last check is       0.02 s)
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -15.10477 -15.10477(    0) -45.26187(    0)  30.15710(    0)  21.97741(    0)                    0.87500        10   0.00000(    0)
End       of equilibration (total CPU time is       0.02 s, CPU time since last check is       0.00 s)
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -14.46468 -14.46468(    0) -23.09322(    0)   8.62853(    0)  11.98111(    0)                    0.85000        10   0.57041(    0)
 -14.70932 -14.58700( 8649) -23.73700(45522)   9.15000(36873)  11.76154(15526)                    0.87500        20   0.47179( 9862)
 -14.61213 -14.59538( 5806) -27.39881(*****)  12.80343(*****)  13.48930(*****)                    0.85833        30   0.44047( 6498)
 -14.22246 -14.50215( 9173) -24.88218(*****)  10.38003(*****)  12.18645(*****)                    0.85625        40   0.57721(14425)
 -14.20866 -14.44345( 9023) -24.28594(*****)   9.84249(*****)  11.89270(*****)                    0.83500        50   0.63646(12648)
 -14.82661 -14.50731( 9514) -26.91444(*****)  12.40713(*****)  13.16759(*****)                    0.84167        60   0.61442(10559)
 -14.55099 -14.51355( 8176) -25.82345(*****)  11.30990(*****)  12.63426(*****)                    0.84643        70   0.62165( 8954)
 -14.85203 -14.55586( 8176) -26.30283(*****)  11.74697(*****)  13.06038(*****)                    0.84688        80   0.63261( 7831)
 -14.53880 -14.55396( 7269) -25.89373(*****)  11.33976(*****)  12.80172(*****)                    0.84167        90   0.61069( 7246)
 -14.98397 -14.59696( 7710) -26.98448(*****)  12.38751(*****)  13.33591(*****)                    0.84250       100   0.61951( 6541)
=======
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0037404997    0.6225029055  -69.1052682346    6.6519320524  -21.0806767058  -14.4287446535

 /BM/test -0.31331888095519095       0.62883871943766578      -0.20819293062787911      -0.88171829020826586      -0.19312797169391074       0.28166040979961288       0.60500838392687106      -0.23083377482807357      -0.53687051617475234       -2.4506199481476049E-002  -5.6022001822949287E-002 -0.31090953347892053     
 /BM/test  -3.7404997287269108E-003
 /BM/test -0.47400698318092149     
 /BM/test   2.6770310548532965       -4.8424183769256111        1.3951548923744943       -2.1797685010362722      -0.18568410050122783       0.95277160943390038       0.77979812931779624       0.40548149968565445      -0.45606630399215869       0.19925413508194267       0.53183248392991489        4.0302292817303256     
 /BM/test  -28.840281731031205       -10.573838616545174       -3.1444096001630744       -27.050702841605222     
 /BM/test  -69.105268234645777     
 /BM/test  -21.080676705830538     
 /BM/test   5.9874161183893087     
 /BM/test  -14.428744653469259     
 /BM/test -0.69928214252029575       -7.5581908293623545E-002 -0.69557243620268649      -0.77404479608264054      -0.65771454202313973      -0.12272332236179634       -1.1841755951139830       -1.9350953495186682       -2.1580437762207847       -2.0598136950097388      -0.18718373546726075       0.23032663679619336      -0.45976841985682348       0.50359552971986066      -0.44939123839230710      -0.60388188191619285       0.35610820845354496      -0.91802800622495262       0.25694740180485631       -1.9915264007476985       0.11987738041812257       -2.6673034376054865       -2.3140378356136697       -2.8507949943655042       -5.6360291941901508E-002
Beginning of equilibration (total CPU time is       0.02 s, CPU time since last check is       0.02 s)
 /BM/checkBEGIN  0.31664587296845426     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0016475001    0.7187116306  -64.6332187370   18.3180450651  -33.6737768938  -15.3557318288

 /BM/checkEND  0.24863814909575055     
 /BM/checkBEGIN  0.39406869444187365     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0005529207    0.7554206572 -144.2440564563   58.1687070439  -73.8483833245  -15.6796762806

 /BM/checkEND  0.79808546059215857     
 /BM/checkBEGIN  0.67815459097849740     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0011210039    0.7580375764  -55.2611433808   13.8699083220  -28.9562063556  -15.0862980335

 /BM/checkEND  0.82729766733121934     
 /BM/checkBEGIN  0.27458108026582195     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0063892333    0.7113988131  -80.9870697325   26.0619433934  -41.0662607588  -15.0043173654

 /BM/checkEND  0.91174659347303688     
 /BM/checkBEGIN  0.33097075404743137     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0151768412    0.6282167972  -52.6798360994   11.7703649595  -26.2072261287  -14.4368611692

 /BM/checkEND  0.46088521649458514     
 /BM/checkBEGIN  0.88234184120071646     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0256614926    0.6689217880  -61.0835705108   17.1966695226  -32.2635608291  -15.0668913066

 /BM/checkEND  0.78235673269010064     
 /BM/checkBEGIN  0.71662340325085339     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0252618838    0.6930485740 -224.2490274647   98.0438813937 -113.3165347415  -15.2726533479

 /BM/checkEND  0.71904446145888912     
 /BM/checkBEGIN  0.11720876150579684     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0114602153    0.7191465735 -205.6394116609   88.2453981930 -104.2184004202  -15.9730022272

 /BM/checkEND  0.33260434910711112     
 /BM/checkBEGIN  0.33088838211921612     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0712521813    0.5977118350 -107.8763256141   39.2379085508  -54.1965615407  -14.9586529899

 /BM/checkEND  0.40871118735564593     
 /BM/checkBEGIN  0.21602412012035543     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0080220473    0.7097277106 -181.4747346855   76.9051530976  -91.4894168988  -14.5842638012

 /BM/checkEND  0.75285002715464699     
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -15.14183 -15.14183(    0) -59.92363(    0)  44.78180(    0)  29.45321(    0)                    0.92500        10   0.00000(    0)
End       of equilibration (total CPU time is       0.02 s, CPU time since last check is       0.00 s)
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0080220473    0.7097277106 -181.4747346855   76.9051530976  -91.4894168988  -14.5842638012

 /BM/test -0.20105521793083422        2.8335040232265039E-002  -5.1483258803494296E-002  -1.8011550031260477        1.3371592279244644       0.82528231662237950       -2.2724069695289910       -3.8443636892378903      -0.47812498189717823        3.3163442511909613E-002  -4.0722345994028589E-002   7.9777874394237874E-003
 /BM/test  -8.0220473279413486E-003
 /BM/test -0.34287388864560353     
 /BM/test   3.4304404107380835      -0.44835711251147881       0.87312128053349747       0.55957095370198939      -0.14383354692682446      -0.18151389660026090       0.46824734698774123       0.57695879608424394        3.9546696026730901E-002  -2.2309338262020444        2.9564793214312699      -0.54300136755391781     
 /BM/test  -35.367806200600860      -0.75450821363061726      -0.37353824237630623       -147.58637688690840     
 /BM/test  -181.47473468550209     
 /BM/test  -91.489416898837931     
 /BM/test   5.4714692156564197     
 /BM/test  -14.584263801190769     
 /BM/test  -6.9671517611278447E-002  -2.7353602284042466      -0.40566167449547552        2.9544789919886450E-002  0.12894389147513682       0.10352042495585261      -0.91288471988231501      -0.61311692023766196      -0.36553239031151924      -0.32609395402579655        1.7178162064902862        7.5536489249895169        5.3521732844056995        7.9655383761394498        19.075680786557527        17.877897414947519        6.7107187522697602        10.847479468207830        19.765632686873293        34.280462169522956        17.161856402839344        41.391113984180002        30.549444655828356        40.227126703085375        15.885227561479260     
 /BM/checkBEGIN  0.51491572598451185     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.1049091112    0.5991748088 -154.4725435536   62.8442841966  -77.8080557262  -14.9637715296

 /BM/checkEND  0.83109581771175201     
 /BM/checkBEGIN  0.76562837065501199     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0561314509    0.6638498396 -128.2982915928   50.0125488377  -64.9260839825  -14.9135351448

 /BM/checkEND  0.96133178265143826     
 /BM/checkBEGIN  0.64425303761646902     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0083033426    0.6576733349  -41.6968304959    6.7968260732  -21.2068763904  -14.4100503172

 /BM/checkEND  0.63765808152127690     
 /BM/checkBEGIN  0.32956125851055518     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0012209839    0.6229617848  -86.6023055687    6.8375959512  -20.6880221262  -13.8504261750

 /BM/checkEND  0.22627164129955801     
 /BM/checkBEGIN  0.45060917482317464     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0009462316    0.7003247789  -25.2595813273   -0.3920823893  -13.1853965727  -13.5774789620

 /BM/checkEND  0.42352590573214499     
 /BM/checkBEGIN  0.49782480071960222     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0027781708    0.7334182123  -75.6948228945   24.0356295239  -39.4352359957  -15.3996064719

 /BM/checkEND  0.90898119721397208     
 /BM/checkBEGIN  0.73321829910145553     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0010670895    0.7476469131  -56.5195926172   14.5696699667  -29.9010730837  -15.3314031170

 /BM/checkEND  0.36654891237868625     
 /BM/checkBEGIN  0.37197664086824389     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002309650    0.7711085372  -20.3168854079   -2.5264227456  -11.5288929352  -14.0553156809

 /BM/checkEND  0.41986289648285080     
 /BM/checkBEGIN  0.55767275609989397     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0003457216    0.7557309926  -22.6424432675   -1.4106727733  -12.5582254858  -13.9688982592

 /BM/checkEND  0.45658614254316987     
 /BM/checkBEGIN  0.70109141953500753     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0001431541    0.8027110757  -43.2231712592    8.1528402224  -23.1555645112  -15.0027242888

 /BM/checkEND  0.47102553954026050     
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -14.54732 -14.54732(    0) -31.43934(    0)  16.89202(    0)  16.36816(    0)                    0.77500        10   0.62272(    0)
 /BM/checkBEGIN  0.20244328168232428     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002524415    0.7772602547  -33.5148800942    3.4949365193  -18.2333022598  -14.7383657405

 /BM/checkEND  0.16625850977557022     
 /BM/checkBEGIN  0.94883973146898271     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0058363797    0.7389607411  -34.3335288662    4.1759936319  -18.8907716546  -14.7147780227

 /BM/checkEND  0.96089557515365343     
 /BM/checkBEGIN  0.59397013651051012     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0187889163    0.6783274621  -52.3976790194   12.6157869023  -26.7161321264  -14.1003452241

 /BM/checkEND  0.25891961472940395     
 /BM/checkBEGIN  0.72943687653009093     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0165781791    0.6999569757  -53.3008812683   13.3001974194  -28.2004180077  -14.9002205882

 /BM/checkEND  0.94303675411601873     
 /BM/checkBEGIN  0.35216396937835626     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0077087004    0.6915091301  -59.2683776123   15.6940314325  -30.1647237138  -14.4706922814

 /BM/checkEND  0.73992151849126842     
 /BM/checkBEGIN  0.37808780301805101     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0069454512    0.7096795822  -55.2698355881   14.0487264836  -29.0268779182  -14.9781514346

 /BM/checkEND  0.81381624137246078     
 /BM/checkBEGIN  0.33368354591297944     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0062424160    0.7064129354  -60.3605924648   16.6516109319  -31.7570997000  -15.1054887680

 /BM/checkEND  0.93038294379816122     
 /BM/checkBEGIN   5.3466455299744808E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0076570300    0.6979470105  -43.5973236757    8.3490736223  -22.4589695477  -14.1098959254

 /BM/checkEND  0.66283289721790695     
 /BM/checkBEGIN  0.60449345100736096     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0335318793    0.6498669105  -73.9785747786   23.4413503090  -38.4596038011  -15.0182534921

 /BM/checkEND  0.69221329885074212     
 /BM/checkBEGIN   5.3522063349081606E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0037529454    0.6777422886  -39.6129547302    4.7825965780  -19.2611365422  -14.4785399643

 /BM/checkEND  0.50229250943420922     
 -14.66147 -14.60440( 4036) -28.87812(*****)  14.27373(*****)  14.50451(*****)                    0.78750        20   0.50609(11664)
 /BM/checkBEGIN  0.95110331280216442     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0021743349    0.7141492865  -39.5998613711    6.5589281542  -21.2217630660  -14.6628349118

 /BM/checkEND  0.38479618383153635     
 /BM/checkBEGIN  0.28558918656719712     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0062258518    0.6746017643  -50.6464236608   10.7062022567  -24.4197019714  -13.7134997147

 /BM/checkEND   2.8124647683949178E-002
 /BM/checkBEGIN  0.45314982754345934     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0000889424    0.7857966771  -19.5279047334   -2.9713538504  -10.7603401328  -13.7316939831

 /BM/checkEND  0.13378466516784826     
 /BM/checkBEGIN  0.81778628860766744     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002762853    0.7884276996  -23.4035853665   -1.0793777184  -13.0721854518  -14.1515631702

 /BM/checkEND  0.66929406173205663     
 /BM/checkBEGIN  0.17591543941091814     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002060848    0.7996363974  -48.8595108308   10.8458862187  -26.1121117303  -15.2662255117

 /BM/checkEND  0.22875379273612495     
 /BM/checkBEGIN  0.69579146832129268     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0024329955    0.7481825647  -37.7960798512    5.8128311747  -20.6795350818  -14.8667039072

 /BM/checkEND  0.49099089990358991     
 /BM/checkBEGIN  0.27329674161933681     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0041475943    0.7499937679  -48.3579801829   10.7453216936  -25.6737403393  -14.9284186457

 /BM/checkEND  0.40179957198362715     
 /BM/checkBEGIN  0.44287916457475163     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002223540    0.7784112824  -21.0558068247   -2.1601906170  -11.6737698260  -13.8339604430

 /BM/checkEND   8.2260688712590735E-002
 /BM/checkBEGIN  0.28390172471060637     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0005870608    0.7034765254  -24.9374932207   -2.0976194460  -11.5130972227  -13.6107166687

 /BM/checkEND  0.28138410575648365     
 /BM/checkBEGIN  0.52287189595310579     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0000304175    0.7899557107  -21.3257852318   -3.9161799966   -9.4634163441  -13.3795963407

 /BM/checkEND  0.25791132381327131     
 -14.21452 -14.47444(10947) -25.07174(*****)  10.59730(*****)  12.46560(*****)                    0.81667        30   0.57944( 9958)
 /BM/checkBEGIN  0.74193628199262562     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0017729819    0.7586270510  -36.6024678364    5.1529840231  -19.9106670335  -14.7576830104

 /BM/checkEND  0.31613746726626957     
 /BM/checkBEGIN  0.83338701122935177     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0023752726    0.7470226635  -36.2987989176    4.9294843211  -19.2846623667  -14.3551780456

 /BM/checkEND  0.93973881305269558     
 /BM/checkBEGIN  0.59730538316168591     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0060752932    0.7053409701  -76.9680549264   24.7687912607  -39.8153815554  -15.0465902947

 /BM/checkEND  0.24563716190665374     
 /BM/checkBEGIN  0.58046708642854483     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0003558184    0.7266782890  -19.1507348255   -3.0121897332  -10.0995328664  -13.1117225996

 /BM/checkEND   4.6371053274292251E-003
 /BM/checkBEGIN  0.92804535704008018     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002419980    0.7086656081  -18.8452482172   -3.3872630846   -9.5942622487  -12.9815253333

 /BM/checkEND   8.3654985687314110E-002
 /BM/checkBEGIN  0.81280265729725798     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0009108438    0.7228108887  -25.0347817940   -0.4095926428  -13.5262130552  -13.9358056980

 /BM/checkEND  0.91857718627398910     
 /BM/checkBEGIN   1.1759158088072041E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0182924431    0.6700905511  -54.3466920136   13.2975277487  -28.1984792173  -14.9009514685

 /BM/checkEND  0.46609823435217024     
 /BM/checkBEGIN  0.46465517897047803     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0417986860    0.6535007409  -71.2755573441   21.8754846403  -36.4025820125  -14.5270973722

 /BM/checkEND  0.29500520801611785     
 /BM/checkBEGIN  0.92242114154192123     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0370937555    0.6706637013  -95.2814832373   33.9548536008  -48.8082141745  -14.8533605737

 /BM/checkEND  0.98829028263021357     
 /BM/checkBEGIN   7.4173287520356013E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0113019298    0.6831315398  -63.7747333987   17.8040859094  -32.0998873758  -14.2958014664

 /BM/checkEND  0.28135463861123000     
 -14.27657 -14.42497( 9261) -25.24730(*****)  10.82233(*****)  12.45906(*****)                    0.84375        40   0.61507( 7891)
 /BM/checkBEGIN  0.22779018776049398     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0115035247    0.6851001896  -51.3643059466   11.7270629767  -26.5541158558  -14.8270528791

 /BM/checkEND  0.17368154112461909     
 /BM/checkBEGIN  0.53876887885951774     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0160554205    0.6973964661  -54.5590722929   13.8471581068  -28.7278213542  -14.8806632474

 /BM/checkEND  0.20509671331698698     
 /BM/checkBEGIN  0.68120301020227103     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0095150480    0.6623384999  -53.0262906657   12.4025103424  -27.5158403916  -15.1133300492

 /BM/checkEND  0.71166134047173202     
 /BM/checkBEGIN  0.67699166735849658     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0041310467    0.7412552046  -52.3050592781   12.7277034665  -27.7863221938  -15.0586187273

 /BM/checkEND  0.17456093517933979     
 /BM/checkBEGIN  0.74599973110440843     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0293685319    0.6486280539  -70.9547062562   21.6322351154  -36.4767056191  -14.8444705038

 /BM/checkEND  0.45816130310738501     
 /BM/checkBEGIN  0.57227451123879192     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0005047197    0.7337034737  -24.6453679082   -0.5937378616  -12.9081873245  -13.5019251861

 /BM/checkEND  0.20469012822688271     
 /BM/checkBEGIN  0.99716390076051553     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0010020273    0.7342163476  -35.5535176614    4.4071294349  -18.6524348221  -14.2453053871

 /BM/checkEND  0.37014854410521636     
 /BM/checkBEGIN  0.53217528462349506     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0001130965    0.7962433194  -22.3299419453   -1.6978029906  -12.3612425102  -14.0590455009

 /BM/checkEND  0.53214233070162820     
 /BM/checkBEGIN  0.84596668674334197     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0016682089    0.7483804473  -43.7534806231    8.4611522393  -23.4626283837  -15.0014761444

 /BM/checkEND  0.75525225931850670     
 /BM/checkBEGIN  0.33609150893497386     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0036366264    0.6920423494  -71.3712493706   22.0081547806  -37.3957972641  -15.3876424835

 /BM/checkEND  0.74312916580036514     
 -14.69195 -14.47837( 8814) -25.23466(*****)  10.75629(*****)  12.36657(*****)                    0.81500        50   0.61136( 6124)
 /BM/checkBEGIN  0.57930560708505041     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0202063594    0.6841033060  -45.7916607480    9.7139553474  -24.5952600840  -14.8813047367

 /BM/checkEND  0.23470809677662530     
 /BM/checkBEGIN  0.33453041614609802     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0079592838    0.6583372299  -40.3870749761    6.5901491516  -21.5784143871  -14.9882652355

 /BM/checkEND  0.82531873727647209     
 /BM/checkBEGIN  0.26158084667651238     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0058376413    0.6848174631  -34.2042136898    3.2412513915  -17.8084440418  -14.5671926503

 /BM/checkEND  0.43950749376319820     
 /BM/checkBEGIN  0.44067905946588226     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0007488975    0.7829854554  -25.5526001185   -0.0563694559  -14.3115000304  -14.3678694864

 /BM/checkEND  0.78161851088223600     
 /BM/checkBEGIN  0.43306262155860864     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0023276390    0.7588848512  -42.8242859840    8.0619325809  -22.7586550924  -14.6967225115

 /BM/checkEND  0.80502100981207647     
 /BM/checkBEGIN  0.96910167836265160     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0104857413    0.7015706143  -73.9164312039   22.7543200043  -38.0632052097  -15.3088852054

 /BM/checkEND  0.59931562810394823     
 /BM/checkBEGIN  0.80840833130233491     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0055673536    0.7135699886  -89.3563383092   30.6055039685  -46.2641036891  -15.6585997206

 /BM/checkEND  0.35983482526033228     
 /BM/checkBEGIN  0.16258027043690859     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0029494629    0.7260874222  -31.8182776441    2.8146251675  -17.4885608641  -14.6739356966

 /BM/checkEND   2.9311530460756074E-002
 /BM/checkBEGIN  0.74715261025375312     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0113387425    0.6760224100  -44.8053440992    8.7531764955  -23.6039854752  -14.8508089797

 /BM/checkEND  0.84304146604811692     
 /BM/checkBEGIN  0.71969884417976360     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0041244984    0.6902584649  -31.0675471887    2.3456206620  -16.6509766812  -14.3053560192

 /BM/checkEND   3.0998226247117344E-002
 -14.82989 -14.53696( 9086) -25.08094(*****)  10.54398(*****)  12.22099(95302)                    0.80000        60   0.59528( 5252)
 /BM/checkBEGIN  0.67860602595925812     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0091212043    0.6735589617  -41.4612647783    7.2207674814  -21.6947075220  -14.4739400405

 /BM/checkEND  0.62936863275522725     
 /BM/checkBEGIN  0.13363194152723068     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0002867384    0.7394493013  -78.2542740747   20.6289376985  -36.2885713045  -15.6596336060

 /BM/checkEND  0.42863998499002420     
 /BM/checkBEGIN  0.51011869533016707     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0006813496    0.7460630894  -38.5328775278    5.9309478770  -20.6308596717  -14.6999117947

 /BM/checkEND  0.96975793891509099     
 /BM/checkBEGIN  0.50857251081160726     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0024527143    0.7382598010  -49.9580506446   11.3158219960  -26.4969469249  -15.1811249290

 /BM/checkEND  0.51802741296408428     
 /BM/checkBEGIN  0.35142399540118063     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0145186592    0.6959149362  -88.1955226122   30.1475687615  -45.1150166979  -14.9674479364

 /BM/checkEND  0.94453688619908505     
 /BM/checkBEGIN  0.92628101422687692     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0137949531    0.6593564458  -57.8295563140   14.9785984524  -29.4033863971  -14.4247879447

 /BM/checkEND  0.34123194965958348     
 /BM/checkBEGIN  0.15419037466910268     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0026584417    0.7376901433  -47.2963761701   10.1771599070  -24.9617000985  -14.7845401915

 /BM/checkEND  0.20648739676245143     
 /BM/checkBEGIN  0.78046431744695965     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0026279791    0.7378910823  -64.4618661725   18.6110772472  -33.8062969171  -15.1952196699

 /BM/checkEND  0.43037219572878982     
 /BM/checkBEGIN  0.54154067979314746     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0009352003    0.7278686392  -58.8847493364   15.8636177386  -30.8727265132  -15.0091087746

 /BM/checkEND  0.55739016086341664     
 /BM/checkBEGIN  0.12721981541693950     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0196475284    0.6932594494  -77.2923780945   24.9176822694  -39.6189133459  -14.7012310765

 /BM/checkEND  0.99513568810523267     
 -14.90969 -14.59020( 9217) -25.91065(*****)  11.32044(*****)  12.62573(89872)                    0.81429        70   0.58201( 4633)
 /BM/checkBEGIN  0.64802772483733762     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0024922278    0.7235134705  -44.1325764864    8.6188682824  -23.4663307027  -14.8474624203

 /BM/checkEND   3.8125185994100974E-002
 /BM/checkBEGIN  0.17777257583162864     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0063716183    0.6817436309  -42.5214622100    6.3681804281  -21.1148950712  -14.7467146431

 /BM/checkEND  0.69138084363266827     
 /BM/checkBEGIN  0.77681326549850738     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0036018402    0.7192764903  -63.4565403261   17.5759314809  -32.8717210906  -15.2957896097

 /BM/checkEND  0.86134403467259446     
 /BM/checkBEGIN  0.10111552926054301     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0238797419    0.6857326456  -60.8636124975   17.0036396777  -31.9609676897  -14.9573280120

 /BM/checkEND  0.90448306235925813     
 /BM/checkBEGIN  0.86052024374339808     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0592268054    0.6575849258 -115.3437256913   44.0688048775  -58.9470323058  -14.8782274283

 /BM/checkEND  0.73552955303740575     
 /BM/checkBEGIN  0.58795573657807054     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0074887269    0.7163822623  -43.1271969986    8.2480811118  -22.7324043106  -14.4843231989

 /BM/checkEND  0.25261337623741653     
 /BM/checkBEGIN  0.87039509508648294     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0091529964    0.6330214375  -93.4207688093   16.5530943412  -31.5328213475  -14.9797270063

 /BM/checkEND   8.1162661632543376E-002
 /BM/checkBEGIN  0.56784323006249693     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0053705655    0.7003273675  -37.2162112663    5.2918634044  -19.5030531082  -14.2111897038

 /BM/checkEND  0.56302226456501714     
 /BM/checkBEGIN  0.24951528007540347     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0040162842    0.7231425577  -49.1799838604   11.1478817184  -26.1372659726  -14.9893842542

 /BM/checkEND  0.44459706795835174     
 /BM/checkBEGIN  0.90583266430633458     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0060475086    0.7303339099  -67.9170518970   20.3305825588  -35.3153831830  -14.9848006241

 /BM/checkEND  0.16866685336263387     
 -14.83749 -14.62112( 8568) -26.46659(*****)  11.84547(*****)  12.97620(85198)                    0.81875        80   0.55972( 4590)
 /BM/checkBEGIN  0.82962511904745284     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0020068854    0.7154643979  -60.5555319246   15.9903716974  -30.3602437906  -14.3698720932

 /BM/checkEND  0.23909558468711722     
 /BM/checkBEGIN  0.22463816369583256     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0010421986    0.7497126361  -53.9301049253   13.2073575325  -28.2212723590  -15.0139148265

 /BM/checkEND   4.3474199696685645E-002
 /BM/checkBEGIN  0.45090537111076756     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0006187216    0.7544347837  -33.0808508446    3.3210893749  -17.4124558229  -14.0913664480

 /BM/checkEND  0.37815533690546133     
 /BM/checkBEGIN  0.48836744829279510     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0014829441    0.7051707208  -30.3548355272    1.7000450603  -15.0112846721  -13.3112396119

 /BM/checkEND  0.17485541563094031     
 /BM/checkBEGIN  0.36557766726366836     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0332051208    0.6580196110  -84.7436556996   28.4648596041  -43.6626238349  -15.1977642308

 /BM/checkEND  0.42737503600694637     
 /BM/checkBEGIN  0.95597122812853641     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0155853477    0.6735407973  -58.4957205260   15.6658976235  -30.6579492745  -14.9920516510

 /BM/checkEND  0.15136592075482369     
 /BM/checkBEGIN  0.84922449076293915     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0000662004    0.7667883963  -22.8852055477   -1.5571157469  -12.3334055442  -13.8905212911

 /BM/checkEND  0.26410310084865429     
 /BM/checkBEGIN  0.67115413287391235     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0049395031    0.7021110002  -30.6659080601    2.1178226666  -16.4064314256  -14.2886087590

 /BM/checkEND  0.66067973473014874     
 /BM/checkBEGIN  0.61749007444990767     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0072117375    0.6682884998  -59.2305745975    9.9061404128  -24.8690168301  -14.9628764173

 /BM/checkEND  0.33767896556922139     
 /BM/checkBEGIN  0.61729909572850161     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0007248829    0.7630954096  -46.6904099794    9.8790807486  -24.9427805375  -15.0636997890

 /BM/checkEND  0.72642061981134631     
 -14.51819 -14.60968( 7692) -26.23561(*****)  11.62593(*****)  12.86949(76397)                    0.81389        90   0.56463( 4078)
 /BM/checkBEGIN   4.9329460007928816E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0003791748    0.7717377589  -49.1676971247   10.7965475311  -26.0537668642  -15.2572193331

 /BM/checkEND   9.5105469530306408E-002
 /BM/checkBEGIN  0.16872777767858693     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0032234815    0.7398775650  -72.6476957709   22.5977428752  -37.5559514817  -14.9582086065

 /BM/checkEND  0.30742067376676729     
 /BM/checkBEGIN  0.42745615676202320     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0023025512    0.6801117991  -24.4397944372   -0.9538625885  -11.6319152410  -12.5857778295

 /BM/checkEND  0.11380764110237962     
 /BM/checkBEGIN  0.26171791332638250     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0017072006    0.7091301875  -26.9718573918    0.5630857081  -14.1386739547  -13.5755882466

 /BM/checkEND  0.76866441187727830     
 /BM/checkBEGIN  0.29547175212374910     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0091867829    0.6955857740  -55.1536911046   13.7254846020  -28.4951232961  -14.7696386941

 /BM/checkEND  0.58522079241837943     
 /BM/checkBEGIN  0.91050657242275079     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0032010397    0.7437980239  -43.6279761970    8.4116478555  -22.9648807647  -14.5532329092

 /BM/checkEND  0.35020823100425602     
 /BM/checkBEGIN  0.22228129892622306     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0008110311    0.7660439052  -38.2380036702    5.6542979245  -20.3965416358  -14.7422437113

 /BM/checkEND  0.95744826462039256     
 /BM/checkBEGIN   8.7988433461195115E-002
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0009603192    0.7423585209 -104.4247041155   37.9684088824  -53.7632153706  -15.7948064881

 /BM/checkEND  0.11529288164127749     
 /BM/checkBEGIN  0.12982199595973398     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0026210698    0.7433772539  -59.1930913930   15.9497707748  -31.0331414618  -15.0833706870

 /BM/checkEND  0.26827641844486649     
 /BM/checkBEGIN  0.90701582205256059     
hpsi: ifr,psid,exp(psij),d2psi,energy-pe,pe,energy= 1   -0.0027552619    0.7477423706  -34.7976679830    4.3423002536  -19.0708690110  -14.7285687575

 /BM/checkEND   8.1707874532423119E-002
 -14.60487 -14.60920( 6923) -26.26309(*****)  11.65389(*****)  12.85419(68772)                    0.81250       100   0.60073( 5132)
>>>>>>> master
End       of accumulation (total CPU time is       0.04 s, CPU time since last check is       0.02 s)

vmc_mov1      no title
Final results after         100. passes (nstep =     10, nblk =     10)
physical variable  average       rms error    rms er*rt(pass)   sigma                Tcor
total E =        -14.6091980 +-  0.0692275  0.69228  0.60073  0.60073 +-  0.05132    1.33
potential E =    -26.2630875 +-  1.2298950 12.29895
interaction E =    4.2929199 +-  0.0898924  0.89892
jf kinetic E =    12.8541947 +-  0.6877239  6.87724
pb kinetic E =    11.6538895 +-  1.1880825 11.88083
<r> =              0.0000000 +-  0.0000000  0.00000
<r2> =             4.3463938 +-  0.1805169  1.80517
<r3> =                0.0000 +-    0.00000     0.00
<r4> =                  0.00 +-      0.000      0.0
acceptance         0.8125000

Gradient with respect to the parameters:
<<<<<<< HEAD
gradient component #     1 :  2.00899978E-01 (CSF)
gradient component #     2 : -2.22860910E-01 (CSF)
gradient component #     3 :  5.22291295E-02 (CSF)
gradient component #     4 :  9.92291145E-02 (Jastrow)
gradient component #     5 :  6.52394541E-02 (Jastrow)
gradient component #     6 :  4.80058536E-02 (Jastrow)
gradient component #     7 :  4.02186754E-02 (Jastrow)
gradient component #     8 : -3.28415396E-03 (Jastrow)
gradient component #     9 :  7.39206999E-02 (Jastrow)
gradient component #    10 :  9.80947317E-02 (Jastrow)
gradient component #    11 :  1.12914215E-01 (Jastrow)
gradient component #    12 :  1.22743859E-01 (Jastrow)
gradient component #    13 : -6.94523221E-02 (Jastrow)
gradient component #    14 : -2.95485186E-01 (Jastrow)
gradient component #    15 : -1.89561898E-01 (Jastrow)
gradient component #    16 : -3.32154324E-01 (Jastrow)
gradient component #    17 : -6.77453050E-01 (Jastrow)
gradient component #    18 : -6.31415392E-01 (Jastrow)
gradient component #    19 : -2.86334982E-01 (Jastrow)
gradient component #    20 : -3.40887809E-01 (Jastrow)
gradient component #    21 : -7.39376209E-01 (Jastrow)
gradient component #    22 : -1.11204770E+00 (Jastrow)
gradient component #    23 : -6.48158710E-01 (Jastrow)
gradient component #    24 : -1.35885511E+00 (Jastrow)
gradient component #    25 : -9.89786957E-01 (Jastrow)
gradient component #    26 : -1.32785498E+00 (Jastrow)
gradient component #    27 : -6.17349674E-01 (Jastrow)
gradient component #    28 : -1.11942035E-02 (exponent)
gradient component #    29 : -1.70543796E-01 (exponent)
gradient component #    30 :  7.88222108E-03 (exponent)
gradient component #    31 : -5.86896228E-02 (exponent)
gradient component #    32 : -5.76192028E-02 (exponent)
gradient component #    33 : -9.23208354E-02 (exponent)
gradient component #    34 :  9.17072300E-03 (orbitals)
gradient component #    35 :  7.61815408E-02 (orbitals)
gradient component #    36 :  4.55621838E-01 (orbitals)
gradient component #    37 :  8.38868146E-01 (orbitals)
gradient component #    38 : -2.99808698E-02 (orbitals)
gradient component #    39 : -2.72630486E-01 (orbitals)
gradient component #    40 : -5.70319096E-01 (orbitals)
gradient norm :               3.15975635E+00 +-  8.01331636E-01
=======
gradient component #     1 : -3.96062506E-01 (CSF)
gradient component #     2 :  9.49443377E-02 (Jastrow)
gradient component #     3 :  6.00709323E-02 (Jastrow)
gradient component #     4 :  4.02648398E-02 (Jastrow)
gradient component #     5 :  2.86457948E-02 (Jastrow)
gradient component #     6 : -5.27396693E-03 (Jastrow)
gradient component #     7 :  7.89525411E-02 (Jastrow)
gradient component #     8 :  1.03712501E-01 (Jastrow)
gradient component #     9 :  1.18525140E-01 (Jastrow)
gradient component #    10 :  1.27839688E-01 (Jastrow)
gradient component #    11 : -5.50921835E-02 (Jastrow)
gradient component #    12 : -2.64012199E-01 (Jastrow)
gradient component #    13 : -1.57406772E-01 (Jastrow)
gradient component #    14 : -2.96621536E-01 (Jastrow)
gradient component #    15 : -6.05649168E-01 (Jastrow)
gradient component #    16 : -5.67893386E-01 (Jastrow)
gradient component #    17 : -2.58711103E-01 (Jastrow)
gradient component #    18 : -2.89029516E-01 (Jastrow)
gradient component #    19 : -6.61872549E-01 (Jastrow)
gradient component #    20 : -9.95691243E-01 (Jastrow)
gradient component #    21 : -5.85572187E-01 (Jastrow)
gradient component #    22 : -1.22292541E+00 (Jastrow)
gradient component #    23 : -8.93477515E-01 (Jastrow)
gradient component #    24 : -1.19691776E+00 (Jastrow)
gradient component #    25 : -5.59334344E-01 (Jastrow)
gradient component #    26 : -1.06857849E-02 (exponent)
gradient component #    27 : -1.93920577E-01 (exponent)
gradient component #    28 : -4.59822187E-03 (exponent)
gradient component #    29 :  8.72096052E-03 (exponent)
gradient component #    30 : -2.28469799E-02 (exponent)
gradient component #    31 : -6.11122246E-02 (exponent)
gradient component #    32 :  2.64957217E-02 (orbitals)
gradient component #    33 :  1.22208682E-01 (orbitals)
gradient component #    34 :  6.79971649E-01 (orbitals)
gradient component #    35 :  7.67239331E-01 (orbitals)
gradient component #    36 : -2.32529771E-02 (orbitals)
gradient component #    37 : -6.93635274E-03 (orbitals)
gradient component #    38 :  1.06058515E-01 (orbitals)
gradient norm :               2.85952323E+00 +-  9.33805891E-01
>>>>>>> master


Maximum absolute value of local energy derivatives : 9.38736623E+00 (must be zero in VMC within statistical noise except for geometry optimization)

<<<<<<< HEAD
Adding to ovlp_lin  0.0000E+00
 /print_too_much/ovlp           1   1.4194465027519478       -8.1895668112101083E-002  -5.5126377745726929E-002   5.7493847174551194E-003   7.7754310341688893E-003   6.8361763692057487E-003   4.4775711588758149E-003   7.1725367981321647E-003  -7.8514034785576836E-002 -0.10250354491498359      -0.12391860299317603      -0.14299091528403141       -4.8565417787065204E-002  -3.0446514002645975E-002 -0.10015357281313597       -5.1026536871649153E-002  -8.1053825549386893E-002  -5.5530448787290965E-002  -2.5181001835136163E-002 -0.14932252278980440      -0.10752125292389847      -0.13207026035123715       -6.6978249280633673E-002 -0.12220636227804516       -7.6427202029247354E-002 -0.10700347735891658       -5.1447963240161876E-002   6.0671841518274497E-004   3.8337728711019432E-002   1.3040010295121800E-002 -0.14152387036412484        6.4190813266017249E-002   1.2614036747299020E-002   2.7248886763595253E-002  -6.2254302026288823E-002 -0.22287650078691174       -3.3215122778823374E-003   3.3738557117501372E-002 -0.33212081928524195      -0.17707465723313651     
 /print_too_much/ovlp           2  -8.1895668112101083E-002   1.1826400907731551      -0.10259196856203041       -3.1815785425953882E-002  -3.5547115972703852E-002  -3.8986105674053501E-002  -4.1295132946060242E-002   9.9315372639566229E-003 -0.14514674233946681      -0.18653139157335258      -0.21825600163309702      -0.24385309524157162       -2.8073547173892037E-002  0.10271223346199115       -4.9299793784502222E-002   9.7584042559055478E-002  0.19778345875682013       0.20243955544679113       0.10225328617280786       -6.3812719486004518E-002  0.20014725613382406       0.29504562115872202       0.20617272354384822       0.40414092213510100       0.30024516870112938       0.40328382015927511       0.20525130631312505        2.7056704104145256E-003  -1.0119551311385511E-003   5.3417384830346092E-003   4.3249940182508769E-003   9.7712596448555594E-002   8.7271671012205926E-002   5.5704703516729577E-002   1.5475099563862921E-002  0.14542483779441284      -0.22520548249897662        7.7801357346078681E-002  0.17711517448492764       0.29009498403780259     
 /print_too_much/ovlp           3  -5.5126377745726929E-002 -0.10259196856203041       0.61930308782418908        3.7904970706252872E-003  -1.1871437278918595E-003  -6.6668197980344068E-003  -1.1759242352183286E-002   2.2082371814314133E-003  -3.1340026704912649E-002  -4.9437767046790704E-002  -6.5976757289800680E-002  -8.0510851695541952E-002  -1.5905468176031001E-002  -2.1672298998585582E-002  -3.3058175361053932E-002  -2.4458801242283723E-002  -5.4729475196003907E-002  -4.6411901292759072E-002  -1.5939877670357649E-002  -4.9724338161196258E-002  -5.8292539542865729E-002  -8.9719156913200315E-002  -4.5132817141972037E-002 -0.10559198435288852       -7.1155494264412589E-002 -0.10011620694488954       -3.9416453951572450E-002   4.1043245241710775E-005  -2.4811780981881110E-003   6.3920078389896415E-003  -2.1507676485202229E-002   8.8887030393122748E-002  -4.1635222169429706E-002   3.2610558307536577E-002   2.0096431805187519E-002  -2.1730009860346965E-002  -2.9061149580947541E-002  0.13560001277519257      -0.25887105996563853       -5.2646541885758513E-002
 /print_too_much/ovlp           4   5.7493847174551194E-003  -3.1815785425953882E-002   3.7904970706252872E-003   8.4495201448245183E-002   9.2130434877699230E-002   9.5304083835422659E-002   9.6074792700555633E-002  -1.3430149692577276E-002  0.18632337402632260       0.22461628197564920       0.24857090031571971       0.26414786054996853       -2.3533244342708315E-003 -0.26166335418345099       -1.5588360975538507E-002 -0.27239996320580673      -0.52718297863022912      -0.51624567912143959      -0.26146491484353973       -3.5226012301867726E-002 -0.55472451561846015      -0.79197161396390925      -0.53418325216784268       -1.0355344666593425      -0.76318209235168410       -1.0272827373896547      -0.52592579712745646       -1.1152771591630486E-003  -3.1573360694369362E-002  -3.2118260887129824E-003 -0.17348601344488257      -0.14887275790609034       -6.3281003399448726E-002   3.4705061278683808E-003   1.3491247614283663E-002   7.5483999117368916E-002   5.2968861960482450E-002 -0.26099523952882020      -0.23166522550744861       -7.4264497879753311E-002
 /print_too_much/ovlp           5   7.7754310341688893E-003  -3.5547115972703852E-002  -1.1871437278918595E-003   9.2130434877699230E-002  0.10613912242672541       0.11310025311480132       0.11604839295053182       -1.4786913559436066E-002  0.21608314099007231       0.26596136485244415       0.29828448402865604       0.31995868940785499       -9.4448404985936207E-004 -0.28599593687056313       -1.2967542954818612E-002 -0.30029193949448540      -0.57059152427636661      -0.55826200640802881      -0.28796748847791065       -3.1133087141940408E-002 -0.60533165223181129      -0.85072206528734284      -0.58271136381719657       -1.1110693570037142      -0.81877487042311259       -1.1017411382328817      -0.57338259514490630       -3.5491998613168541E-004  -1.6145449123808137E-002  -4.5116788306936995E-003 -0.20548262912888293      -0.18563404490111535       -7.2679148067418653E-002  -1.1789071106905746E-003  -4.3044256073734855E-003   1.2801512318931341E-002   1.1028734253426986E-002 -0.33005152517008585      -0.22874243093287694       -5.7020048378205401E-002
 /print_too_much/ovlp           6   6.8361763692057487E-003  -3.8986105674053501E-002  -6.6668197980344068E-003   9.5304083835422659E-002  0.11310025311480132       0.12261793966409584       0.12727693749958746       -1.5311966306409630E-002  0.23078903386186767       0.28795334633133507       0.32602393825517595       0.35223311496315546       -5.9831511355845635E-004 -0.29582089985113669       -1.2517778934189039E-002 -0.31261696328938626      -0.58735455214024057      -0.57414074582558783      -0.29940644460944554       -3.0469900264872507E-002 -0.62659510507367600      -0.87254352687403980      -0.60257693977969495       -1.1377517418875982      -0.83861413969856358       -1.1278392097813139      -0.59266218773821322        9.6763951887984812E-006  -7.3441362887625705E-003  -5.3179202921987412E-003 -0.21991706356067819      -0.20907875853053581       -7.8079181025322242E-002  -4.7101162121518936E-003  -1.5272967183337816E-002  -2.3021300852381410E-002  -3.6597307951360403E-003 -0.37284520021401069      -0.20647716174499164       -4.6775417847801926E-002
 /print_too_much/ovlp           7   4.4775711588758149E-003  -4.1295132946060242E-002  -1.1759242352183286E-002   9.6074792700555633E-002  0.11604839295053182       0.12727693749958746       0.13326281181122113       -1.5387315680638114E-002  0.23682829042823528       0.29850650052090089       0.34052952184447349       0.37012796982395457       -8.4968694713705162E-004 -0.29785846373214042       -1.3266866486079110E-002 -0.31642256287571513      -0.58999066504850362      -0.57620046987904061      -0.30263315143251646       -3.1681928377246038E-002 -0.63201971394602197      -0.87496623984546318      -0.60703184639586993       -1.1389466525261298      -0.83976214395372040       -1.1287285825352740      -0.59680833999257743        1.7371728377799101E-004  -2.7247628933966261E-003  -5.7789028938300957E-003 -0.22436214646036978      -0.22408075458389498       -8.1346767090746674E-002  -7.0712629602471266E-003  -2.1108020864689658E-002  -4.1020190522080502E-002  -6.9608853588642994E-003 -0.39883805856658516      -0.17707443378348464       -4.2171347558710881E-002
 /print_too_much/ovlp           8   7.1725367981321647E-003   9.9315372639566229E-003   2.2082371814314133E-003  -1.3430149692577276E-002  -1.4786913559436066E-002  -1.5311966306409630E-002  -1.5387315680638114E-002   2.4936374353976137E-003  -3.2894908654412358E-002  -3.8956764857710269E-002  -4.2598221801636083E-002  -4.4870387612230722E-002  -1.7328447067243014E-003   4.1789342035696286E-002  -2.3526020103709078E-003   4.2306932298574296E-002   8.2877555217450993E-002   8.2282830782683725E-002   4.1722922848132615E-002  -2.1527156983270901E-003   8.5717338280622890E-002  0.12330621606169867        8.4447012560908252E-002  0.16415962077183632       0.12141350221644931       0.16354077065717920        8.3839573198295714E-002   1.4206818216520607E-004   4.0400544343260972E-003   6.3700787229128975E-004   2.7985234431663814E-002   2.4364661441758070E-002   9.6085766886140969E-003   7.7558797563950788E-004  -1.2628711653515541E-003  -7.3549627243986587E-003  -5.8953758515878918E-003   4.1936980664338402E-002   4.1728574698687043E-002   1.3084595536881558E-002
 /print_too_much/ovlp           9  -7.8514034785576836E-002 -0.14514674233946681       -3.1340026704912649E-002  0.18632337402632260       0.21608314099007231       0.23078903386186767       0.23682829042823528       -3.2894908654412358E-002  0.47130313246514532       0.57642343159992038       0.64363973781966877       0.68814759882697274        1.7678339028977419E-002 -0.57988507418954782        1.6657114707578913E-002 -0.59909057332300364       -1.1457303125995253       -1.1306447145536396      -0.58418373328274953        3.5766580774527768E-003  -1.2043998121652919       -1.6990904007654564       -1.1751562697943712       -2.2428867713502854       -1.6566034902465816       -2.2297015990685338       -1.1621600631865476       -6.3824645598206664E-004  -2.2863226535048076E-002  -1.0433347992156783E-002 -0.41385947447057414      -0.38744387048906237      -0.14904298681878730       -1.9305647767845435E-002  -2.2160804490884367E-002  -3.2204965104825123E-002   3.1478205857807251E-002 -0.67514877720356503      -0.49541536008044729      -0.14684828027973440     
 /print_too_much/ovlp          10 -0.10250354491498359      -0.18653139157335258       -4.9437767046790704E-002  0.22461628197564920       0.26596136485244415       0.28795334633133507       0.29850650052090089       -3.8956764857710269E-002  0.57642343159992038       0.71548545933165997       0.80728026376257844       0.87002771111571064        1.7696037980885038E-002 -0.69808822509205726        1.0868974946861343E-002 -0.72788293311430152       -1.3775487462826703       -1.3558669405772008      -0.70641334984203752       -1.0940987880090347E-002  -1.4583961889312320       -2.0408388230903540       -1.4171088651854404       -2.6832879548694279       -1.9813702777990727       -2.6651731595969466       -1.3992126663625015       -2.9730402793845534E-004  -1.4771573586118802E-002  -1.3813579582175661E-002 -0.50337399740039768      -0.49944359348523792      -0.19067883857931928       -2.9002776288858584E-002  -4.2508893449996810E-002  -8.5625415206808775E-002   2.6748532542341552E-002 -0.87060516176189617      -0.51145528006789287      -0.15402112662103706     
 /print_too_much/ovlp          11 -0.12391860299317603      -0.21825600163309702       -6.5976757289800680E-002  0.24857090031571971       0.29828448402865604       0.32602393825517595       0.34052952184447349       -4.2598221801636083E-002  0.64363973781966877       0.80728026376257844       0.91793501531108745       0.99544648230042299        1.6475325365069793E-002 -0.77120883533855533        4.2281058639179747E-003 -0.80937070913242337       -1.5208949518323607       -1.4942783099071448      -0.78297874843396187       -2.4919390379139372E-002  -1.6179343535027897       -2.2521586166778889       -1.5676830331374845       -2.9526344669129116       -2.1802319604097136       -2.9310264541377933       -1.5462995801352974       -9.1622143486680407E-005  -8.7931670744127444E-003  -1.6235560907834790E-002 -0.55616635425131822      -0.58019943245112726      -0.22207568333606226       -3.6474776628395766E-002  -5.6411287811639121E-002 -0.12115915802481581        2.9921234993583568E-002  -1.0086352123484543      -0.48475191112827387      -0.15578509379802452     
 /print_too_much/ovlp          12 -0.14299091528403141      -0.24385309524157162       -8.0510851695541952E-002  0.26414786054996853       0.31995868940785499       0.35223311496315546       0.37012796982395457       -4.4870387612230722E-002  0.68814759882697274       0.87002771111571064       0.99544648230042299        1.0850834138190635        1.4832165158873423E-002 -0.81801844008384705       -2.1689829839885988E-003 -0.86274597880515103       -1.6127132143581946       -1.5824081194868995      -0.83266912606210042       -3.7379518037617743E-002  -1.7216975154319414       -2.3876508266975520       -1.6647485530682502       -3.1233589491478710       -2.3065404567940959       -3.0992627709838416       -1.6408732385177700        1.4626258551309057E-005  -4.9403252144176735E-003  -1.7966934990198302E-002 -0.58677358540313751      -0.64027262774597515      -0.24766440725771333       -4.2153830775726042E-002  -6.4732240948559783E-002 -0.14286804011048160        3.6701448030887818E-002  -1.1080945214624578      -0.44057923116713571      -0.15725725015297431     
 /print_too_much/ovlp          13  -4.8565417787065204E-002  -2.8073547173892037E-002  -1.5905468176031001E-002  -2.3533244342708315E-003  -9.4448404985936207E-004  -5.9831511355845635E-004  -8.4968694713705162E-004  -1.7328447067243014E-003   1.7678339028977419E-002   1.7696037980885038E-002   1.6475325365069793E-002   1.4832165158873423E-002   1.4393599507105170E-002   6.6034366172935677E-003   3.2921817005692633E-002   1.4411127153530323E-002   2.2745199888857925E-002   1.4835734410082679E-002   6.3886241918531539E-003   5.2686797600824731E-002   3.2739119788986137E-002   4.2718856092500701E-002   1.9092780616965399E-002   3.6768372205649769E-002   2.4321929041855483E-002   3.1973271865387787E-002   1.4163910645992317E-002   3.9925771314459034E-004   1.0644761849358766E-002  -7.4510969831170434E-004  -2.2747240866450580E-003  -1.3614570830440403E-003   4.6657329798289737E-003  -9.5055366216938650E-003  -8.7144587706557280E-003  -4.9951916191877235E-002  -2.2508713851550587E-002   1.8202111600579135E-003  -3.6308607232273957E-002  -5.1415701835082622E-003
 /print_too_much/ovlp          14  -3.0446514002645975E-002  0.10271223346199115       -2.1672298998585582E-002 -0.26166335418345099      -0.28599593687056313      -0.29582089985113669      -0.29785846373214042        4.1789342035696286E-002 -0.57988507418954782      -0.69808822509205726      -0.77120883533855533      -0.81801844008384705        6.6034366172935677E-003  0.81349380378967595        4.6530405354276994E-002  0.84636642691306463        1.6380742305442624        1.6039593346887955       0.81224759126959611       0.10605990932046261        1.7233733555089970        2.4595387345012085        1.6596017937963978        3.2167101060728669        2.3698496996356369        3.1907768979216371        1.6336383931221548        3.2273424072958856E-003   9.1830804033043473E-002   9.5919545912457327E-003  0.54572908500133233       0.45855841524679630       0.19414416758025266       -1.0986270492913480E-002  -3.1439782416144646E-002 -0.20593450451731043      -0.14768962640392036       0.80624117866810008       0.75231952063015795       0.24651174850348256     
 /print_too_much/ovlp          15 -0.10015357281313597       -4.9299793784502222E-002  -3.3058175361053932E-002  -1.5588360975538507E-002  -1.2967542954818612E-002  -1.2517778934189039E-002  -1.3266866486079110E-002  -2.3526020103709078E-003   1.6657114707578913E-002   1.0868974946861343E-002   4.2281058639179747E-003  -2.1689829839885988E-003   3.2921817005692633E-002   4.6530405354276994E-002   7.7317941943676471E-002   6.5811995102723841E-002  0.11588843656249459        9.6420233010352518E-002   4.6110586805845344E-002  0.12587731931291479       0.14247387443424131       0.19469444340829511       0.10840759926374233       0.21023603957180370       0.14862128049702861       0.19813408494301044        9.6024088754518289E-002   1.1937499562675955E-003   2.9534116920905996E-002  -1.2019289321237095E-003   1.2005717923119497E-002   1.6205696619338267E-002   1.9558306208422666E-002  -2.0424111841876237E-002  -2.1577368300703420E-002 -0.12639582120903148       -7.0973756830844420E-002   3.5587781281042383E-002  -6.2958120985825899E-002   3.1284270731885555E-003
 /print_too_much/ovlp          16  -5.1026536871649153E-002   9.7584042559055478E-002  -2.4458801242283723E-002 -0.27239996320580673      -0.30029193949448540      -0.31261696328938626      -0.31642256287571513        4.2306932298574296E-002 -0.59909057332300364      -0.72788293311430152      -0.80937070913242337      -0.86274597880515103        1.4411127153530323E-002  0.84636642691306463        6.5811995102723841E-002  0.88709496022978129        1.7066879821033609        1.6662484214852782       0.84659654452157440       0.13807730448544042        1.8046039546829604        2.5642342813052181        1.7298734507513558        3.3405525420965319        2.4594937465128197        3.3105042452420150        1.6997281310447363        3.1391809279054245E-003   9.1026857956787111E-002   1.0342473294984944E-002  0.57081556382295240       0.49108848617854339       0.20746655252488733       -1.2466345388219391E-002  -2.6216829146765430E-002 -0.20304082660181549      -0.14714095887221235       0.86738170864703523       0.72871872038298247       0.23550467054455682     
 /print_too_much/ovlp          17  -8.1053825549386893E-002  0.19778345875682013       -5.4729475196003907E-002 -0.52718297863022912      -0.57059152427636661      -0.58735455214024057      -0.58999066504850362        8.2877555217450993E-002  -1.1457303125995253       -1.3775487462826703       -1.5208949518323607       -1.6127132143581946        2.2745199888857925E-002   1.6380742305442624       0.11588843656249459        1.7066879821033609        3.3087617244098055        3.2358545926793454        1.6337164516518641       0.25012346761392479        3.4824718518080431        4.9787782590172469        3.3463019245355383        6.5009768289502290        4.7878944756819806        6.4462060376993975        3.2913984601575805        7.4859564966308356E-003  0.20583952509617515        1.8094595975308714E-002   1.0800744179517308       0.90342289991372615       0.39089796612771643       -3.0382504916672826E-002  -8.2886641774731928E-002 -0.49890111398845471      -0.35431307072474971        1.5853332873446111        1.5053846126485468       0.51671617311079143     
 /print_too_much/ovlp          18  -5.5530448787290965E-002  0.20243955544679113       -4.6411901292759072E-002 -0.51624567912143959      -0.55826200640802881      -0.57414074582558783      -0.57620046987904061        8.2282830782683725E-002  -1.1306447145536396       -1.3558669405772008       -1.4942783099071448       -1.5824081194868995        1.4835734410082679E-002   1.6039593346887955        9.6420233010352518E-002   1.6662484214852782        3.2358545926793454        3.1689992050719695        1.5993882314262748       0.21768758553663403        3.3993419708362183        4.8657285477956691        3.2736998822276746        6.3646835426579287        4.6893674288338616        6.3139328236995595        3.2228893200399966        7.3401768879871876E-003  0.20136274914270097        1.7889974703836964E-002   1.0579870675952945       0.88278875209505725       0.37958008237672158       -2.6616943348986066E-002  -8.2561640167662209E-002 -0.48333579399023374      -0.34611398359932788        1.5461917960106439        1.4988972608497262       0.51458811702516871     
 /print_too_much/ovlp          19  -2.5181001835136163E-002  0.10225328617280786       -1.5939877670357649E-002 -0.26146491484353973      -0.28796748847791065      -0.29940644460944554      -0.30263315143251646        4.1722922848132615E-002 -0.58418373328274953      -0.70641334984203752      -0.78297874843396187      -0.83266912606210042        6.3886241918531539E-003  0.81224759126959611        4.6110586805845344E-002  0.84659654452157440        1.6337164516518641        1.5993882314262748       0.81227318225634804       0.10529723391062618        1.7213513285935704        2.4510508797502837        1.6572388389988788        3.2042156009721054        2.3609570258913521        3.1782226223599821        1.6312263088958758        2.9900112485516245E-003   8.6514473793965374E-002   1.0135631252902289E-002  0.54871235096134874       0.47047389841004272       0.19621066388854724       -8.6333364387967959E-003  -2.5948888996808195E-002 -0.18731649288768920      -0.13862417752396061       0.82817322835585683       0.72255136484856786       0.23371592863525614     
 /print_too_much/ovlp          20 -0.14932252278980440       -6.3812719486004518E-002  -4.9724338161196258E-002  -3.5226012301867726E-002  -3.1133087141940408E-002  -3.0469900264872507E-002  -3.1681928377246038E-002  -2.1527156983270901E-003   3.5766580774527768E-003  -1.0940987880090347E-002  -2.4919390379139372E-002  -3.7379518037617743E-002   5.2686797600824731E-002  0.10605990932046261       0.12587731931291479       0.13807730448544042       0.25012346761392479       0.21768758553663403       0.10529723391062618       0.20722344543996485       0.29599353478123192       0.41038249167108631       0.23869645045301269       0.46516441726421931       0.33278921954148899       0.44472428152584342       0.21783431556657717        2.2891650067535588E-003   5.4093996140862174E-002  -1.5062435043621258E-003   3.5398730093479269E-002   4.3662562595134458E-002   3.9456024693127492E-002  -3.1508204487131472E-002  -3.7481936424940035E-002 -0.21982224735568370      -0.13966693608683203        8.5249155066743729E-002  -7.8649052796604974E-002   2.3324768203307339E-002
 /print_too_much/ovlp          21 -0.10752125292389847       0.20014725613382406       -5.8292539542865729E-002 -0.55472451561846015      -0.60533165223181129      -0.62659510507367600      -0.63201971394602197        8.5717338280622890E-002  -1.2043998121652919       -1.4583961889312320       -1.6179343535027897       -1.7216975154319414        3.2739119788986137E-002   1.7233733555089970       0.14247387443424131        1.8046039546829604        3.4824718518080431        3.3993419708362183        1.7213513285935704       0.29599353478123192        3.6784670814476499        5.2404957496284510        3.5243705261135858        6.8255959206932175        5.0247605327376164        6.7638746032719723        3.4624428208364861        7.3329415456072078E-003  0.20483393585269610        1.9829238898359303E-002   1.1449224179985777       0.97412213689606242       0.41921400441472723       -3.1139690397933240E-002  -7.2412585467685142E-002 -0.48675386649626162      -0.35068488511005524        1.7158438872259834        1.5122714003669551       0.51052834301927974     
 /print_too_much/ovlp          22 -0.13207026035123715       0.29504562115872202       -8.9719156913200315E-002 -0.79197161396390925      -0.85072206528734284      -0.87254352687403980      -0.87496623984546318       0.12330621606169867       -1.6990904007654564       -2.0408388230903540       -2.2521586166778889       -2.3876508266975520        4.2718856092500701E-002   2.4595387345012085       0.19469444340829511        2.5642342813052181        4.9787782590172469        4.8657285477956691        2.4510508797502837       0.41038249167108631        5.2404957496284510        7.5034002855641688        5.0290941056719021        9.7881118212530964        7.2077102867756935        9.7037192871115394        4.9444510029678952        1.2471120748224956E-002  0.33326541576992563        2.5953992086136024E-002   1.5990293123544479        1.3364363833188886       0.58647991725561821       -5.2670487915193243E-002 -0.14611410394111068      -0.84305710229466457      -0.60799081857376258        2.3396482885537204        2.2458588057478055       0.80428933275760550     
 /print_too_much/ovlp          23  -6.6978249280633673E-002  0.20617272354384822       -4.5132817141972037E-002 -0.53418325216784268      -0.58271136381719657      -0.60257693977969495      -0.60703184639586993        8.4447012560908252E-002  -1.1751562697943712       -1.4171088651854404       -1.5676830331374845       -1.6647485530682502        1.9092780616965399E-002   1.6596017937963978       0.10840759926374233        1.7298734507513558        3.3463019245355383        3.2736998822276746        1.6572388389988788       0.23869645045301269        3.5243705261135858        5.0290941056719021        3.3887888666801871        6.5690402653811475        4.8388354813162948        6.5143198823017201        3.3339742276093318        6.9413917135472403E-003  0.19452780710271855        1.9442934565911396E-002   1.1056208326607462       0.93607170958082397       0.39862821926338321       -2.4723255588430526E-002  -7.0044970077891344E-002 -0.45147267005496516      -0.32688922061005621        1.6449731684504587        1.5003308538214246       0.50158909342890290     
 /print_too_much/ovlp          24 -0.12220636227804516       0.40414092213510100      -0.10559198435288852       -1.0355344666593425       -1.1110693570037142       -1.1377517418875982       -1.1389466525261298       0.16415962077183632       -2.2428867713502854       -2.6832879548694279       -2.9526344669129116       -3.1233589491478710        3.6768372205649769E-002   3.2167101060728669       0.21023603957180370        3.3405525420965319        6.5009768289502290        6.3646835426579287        3.2042156009721054       0.46516441726421931        6.8255959206932175        9.7881118212530964        6.5690402653811475        12.798340846487747        9.4285102832534449        12.695213845912122        6.4657461654576309        1.6139294142552052E-002  0.43296599181300144        3.4133827018363538E-002   2.0958167722396990        1.7363237177222999       0.75759921381612605       -6.2512376927876367E-002 -0.19336568133999132       -1.0808318184843113      -0.77396169182601637        3.0343616922650849        3.0287928502192112        1.0723991828953374     
 /print_too_much/ovlp          25  -7.6427202029247354E-002  0.30024516870112938       -7.1155494264412589E-002 -0.76318209235168410      -0.81877487042311259      -0.83861413969856358      -0.83976214395372040       0.12141350221644931       -1.6566034902465816       -1.9813702777990727       -2.1802319604097136       -2.3065404567940959        2.4321929041855483E-002   2.3698496996356369       0.14862128049702861        2.4594937465128197        4.7878944756819806        4.6893674288338616        2.3609570258913521       0.33278921954148899        5.0247605327376164        7.2077102867756935        4.8388354813162948        9.4285102832534449        6.9471916308174286        9.3538670263893096        4.7640963381463166        1.2004740976432088E-002  0.32039665370998371        2.5378493627123300E-002   1.5426080617614977        1.2832370771147610       0.55755274228662188       -4.4224209648327406E-002 -0.14388529503774716      -0.80026663832812350      -0.58101033431423832        2.2407059185722726        2.2213886217035506       0.79168374556772148     
 /print_too_much/ovlp          26 -0.10700347735891658       0.40328382015927511      -0.10011620694488954       -1.0272827373896547       -1.1017411382328817       -1.1278392097813139       -1.1287285825352740       0.16354077065717920       -2.2297015990685338       -2.6651731595969466       -2.9310264541377933       -3.0992627709838416        3.1973271865387787E-002   3.1907768979216371       0.19813408494301044        3.3105042452420150        6.4462060376993975        6.3139328236995595        3.1782226223599821       0.44472428152584342        6.7638746032719723        9.7037192871115394        6.5143198823017201        12.695213845912122        9.3538670263893096        12.594851000733001        6.4138341214396064        1.6062659342438973E-002  0.43038481583963351        3.3943919809709455E-002   2.0786944157202845        1.7211803553418692       0.74926660467566708       -6.0461399682674610E-002 -0.19352095147614212       -1.0732499676581786      -0.77062894613043342        3.0062715327643996        3.0163554340725534        1.0688466965434853     
 /print_too_much/ovlp          27  -5.1447963240161876E-002  0.20525130631312505       -3.9416453951572450E-002 -0.52592579712745646      -0.57338259514490630      -0.59266218773821322      -0.59680833999257743        8.3839573198295714E-002  -1.1621600631865476       -1.3992126663625015       -1.5462995801352974       -1.6408732385177700        1.4163910645992317E-002   1.6336383931221548        9.6024088754518289E-002   1.6997281310447363        3.2913984601575805        3.2228893200399966        1.6312263088958758       0.21783431556657717        3.4624428208364861        4.9444510029678952        3.3339742276093318        6.4657461654576309        4.7640963381463166        6.4138341214396064        3.2820174318517843        6.8595567746338182E-003  0.19189028270683423        1.9249934090153409E-002   1.0884871572893848       0.92093682908330798       0.39034578494853633       -2.2599203496345643E-002  -7.0248471434898208E-002 -0.44366283406014040      -0.32310748414182089        1.6167986343115912        1.4882560618385221       0.49842228389084853     
 /print_too_much/ovlp          28   6.0671841518274497E-004   2.7056704104145256E-003   4.1043245241710775E-005  -1.1152771591630486E-003  -3.5491998613168541E-004   9.6763951887984812E-006   1.7371728377799101E-004   1.4206818216520607E-004  -6.3824645598206664E-004  -2.9730402793845534E-004  -9.1622143486680407E-005   1.4626258551309057E-005   3.9925771314459034E-004   3.2273424072958856E-003   1.1937499562675955E-003   3.1391809279054245E-003   7.4859564966308356E-003   7.3401768879871876E-003   2.9900112485516245E-003   2.2891650067535588E-003   7.3329415456072078E-003   1.2471120748224956E-002   6.9413917135472403E-003   1.6139294142552052E-002   1.2004740976432088E-002   1.6062659342438973E-002   6.8595567746338182E-003   2.4906052618373396E-004   3.3842167790956350E-003  -5.6115585784363481E-005  -7.3745286401682986E-004  -3.7165202355745556E-004   5.5295785797880267E-004   3.8891182140953894E-005  -1.8116650838664356E-003  -9.7318247596942450E-003  -1.8996036786841006E-002  -1.5600223464371839E-003   1.5670484336338659E-003   4.0323493365716778E-003
 /print_too_much/ovlp          29   3.8337728711019432E-002  -1.0119551311385511E-003  -2.4811780981881110E-003  -3.1573360694369362E-002  -1.6145449123808137E-002  -7.3441362887625705E-003  -2.7247628933966261E-003   4.0400544343260972E-003  -2.2863226535048076E-002  -1.4771573586118802E-002  -8.7931670744127444E-003  -4.9403252144176735E-003   1.0644761849358766E-002   9.1830804033043473E-002   2.9534116920905996E-002   9.1026857956787111E-002  0.20583952509617515       0.20136274914270097        8.6514473793965374E-002   5.4093996140862174E-002  0.20483393585269610       0.33326541576992563       0.19452780710271855       0.43296599181300144       0.32039665370998371       0.43038481583963351       0.19189028270683423        3.3842167790956350E-003   8.5555624501225627E-002  -8.9163594102783302E-004  -5.5914246360516812E-003  -2.0318054379265815E-003   1.1195177801538871E-002  -2.5719031094762042E-002  -9.5501080160694057E-002 -0.34485594221561883      -0.15649192620130845       -1.9847238100553295E-002   3.5738591261121624E-002   7.2844397012225476E-002
 /print_too_much/ovlp          30   1.3040010295121800E-002   5.3417384830346092E-003   6.3920078389896415E-003  -3.2118260887129824E-003  -4.5116788306936995E-003  -5.3179202921987412E-003  -5.7789028938300957E-003   6.3700787229128975E-004  -1.0433347992156783E-002  -1.3813579582175661E-002  -1.6235560907834790E-002  -1.7966934990198302E-002  -7.4510969831170434E-004   9.5919545912457327E-003  -1.2019289321237095E-003   1.0342473294984944E-002   1.8094595975308714E-002   1.7889974703836964E-002   1.0135631252902289E-002  -1.5062435043621258E-003   1.9829238898359303E-002   2.5953992086136024E-002   1.9442934565911396E-002   3.4133827018363538E-002   2.5378493627123300E-002   3.3943919809709455E-002   1.9249934090153409E-002  -5.6115585784363481E-005  -8.9163594102783302E-004   7.6017709012113376E-004   7.9341882892578902E-003   1.1389024419652142E-002   8.7011759175881437E-004   1.4301276263485609E-003   1.8784698830940972E-003   2.9003981565240912E-003   1.5075519901047926E-003   2.1574990892015332E-002  -6.8010490549877175E-003  -1.5956648388376612E-002
 /print_too_much/ovlp          31 -0.14152387036412484        4.3249940182508769E-003  -2.1507676485202229E-002 -0.17348601344488257      -0.20548262912888293      -0.21991706356067819      -0.22436214646036978        2.7985234431663814E-002 -0.41385947447057414      -0.50337399740039768      -0.55616635425131822      -0.58677358540313751       -2.2747240866450580E-003  0.54572908500133233        1.2005717923119497E-002  0.57081556382295240        1.0800744179517308        1.0579870675952945       0.54871235096134874        3.5398730093479269E-002   1.1449224179985777        1.5990293123544479        1.1056208326607462        2.0958167722396990        1.5426080617614977        2.0786944157202845        1.0884871572893848       -7.3745286401682986E-004  -5.5914246360516812E-003   7.9341882892578902E-003  0.48332612615888870       0.36254594195827866        5.5244547775279915E-002   5.0874435101612922E-003   6.7581450176763841E-002  0.15006829177798406        5.4419789075567389E-002  0.68449471441286647       0.62143962893340388       0.14224555615132461     
 /print_too_much/ovlp          32   6.4190813266017249E-002   9.7712596448555594E-002   8.8887030393122748E-002 -0.14887275790609034      -0.18563404490111535      -0.20907875853053581      -0.22408075458389498        2.4364661441758070E-002 -0.38744387048906237      -0.49944359348523792      -0.58019943245112726      -0.64027262774597515       -1.3614570830440403E-003  0.45855841524679630        1.6205696619338267E-002  0.49108848617854339       0.90342289991372615       0.88278875209505725       0.47047389841004272        4.3662562595134458E-002  0.97412213689606242        1.3364363833188886       0.93607170958082397        1.7363237177222999        1.2832370771147610        1.7211803553418692       0.92093682908330798       -3.7165202355745556E-004  -2.0318054379265815E-003   1.1389024419652142E-002  0.36254594195827866       0.43898486198020548        8.6213389922147887E-002   2.1293774456850274E-002   4.1911483325970347E-002   8.9948363424144584E-002   1.5882224225553451E-003  0.75842985165523358       0.12143870855565117       0.12722102845875335     
 /print_too_much/ovlp          33   1.2614036747299020E-002   8.7271671012205926E-002  -4.1635222169429706E-002  -6.3281003399448726E-002  -7.2679148067418653E-002  -7.8079181025322242E-002  -8.1346767090746674E-002   9.6085766886140969E-003 -0.14904298681878730      -0.19067883857931928      -0.22207568333606226      -0.24766440725771333        4.6657329798289737E-003  0.19414416758025266        1.9558306208422666E-002  0.20746655252488733       0.39089796612771643       0.37958008237672158       0.19621066388854724        3.9456024693127492E-002  0.41921400441472723       0.58647991725561821       0.39862821926338321       0.75759921381612605       0.55755274228662188       0.74926660467566708       0.39034578494853633        5.5295785797880267E-004   1.1195177801538871E-002   8.7011759175881437E-004   5.5244547775279915E-002   8.6213389922147887E-002  0.20689034790626823        7.5288219436662395E-003  -1.7376533578749673E-002  -9.0965994300284698E-004  -5.2429361326836094E-002  0.10907257277790755       0.10679637735939224        1.7004311435452363E-002
 /print_too_much/ovlp          34   2.7248886763595253E-002   5.5704703516729577E-002   3.2610558307536577E-002   3.4705061278683808E-003  -1.1789071106905746E-003  -4.7101162121518936E-003  -7.0712629602471266E-003   7.7558797563950788E-004  -1.9305647767845435E-002  -2.9002776288858584E-002  -3.6474776628395766E-002  -4.2153830775726042E-002  -9.5055366216938650E-003  -1.0986270492913480E-002  -2.0424111841876237E-002  -1.2466345388219391E-002  -3.0382504916672826E-002  -2.6616943348986066E-002  -8.6333364387967959E-003  -3.1508204487131472E-002  -3.1139690397933240E-002  -5.2670487915193243E-002  -2.4723255588430526E-002  -6.2512376927876367E-002  -4.4224209648327406E-002  -6.0461399682674610E-002  -2.2599203496345643E-002   3.8891182140953894E-005  -2.5719031094762042E-002   1.4301276263485609E-003   5.0874435101612922E-003   2.1293774456850274E-002   7.5288219436662395E-003   4.0548114021958209E-002   5.7375631990177443E-002  0.14398474045002665       -9.6658496121299550E-002   3.3468950109978211E-002  -1.2844341225481634E-002  -2.3281304890665357E-002
 /print_too_much/ovlp          35  -6.2254302026288823E-002   1.5475099563862921E-002   2.0096431805187519E-002   1.3491247614283663E-002  -4.3044256073734855E-003  -1.5272967183337816E-002  -2.1108020864689658E-002  -1.2628711653515541E-003  -2.2160804490884367E-002  -4.2508893449996810E-002  -5.6411287811639121E-002  -6.4732240948559783E-002  -8.7144587706557280E-003  -3.1439782416144646E-002  -2.1577368300703420E-002  -2.6216829146765430E-002  -8.2886641774731928E-002  -8.2561640167662209E-002  -2.5948888996808195E-002  -3.7481936424940035E-002  -7.2412585467685142E-002 -0.14611410394111068       -7.0044970077891344E-002 -0.19336568133999132      -0.14388529503774716      -0.19352095147614212       -7.0248471434898208E-002  -1.8116650838664356E-003  -9.5501080160694057E-002   1.8784698830940972E-003   6.7581450176763841E-002   4.1911483325970347E-002  -1.7376533578749673E-002   5.7375631990177443E-002  0.21730771727820561       0.46025048236048899       -9.0751751020772750E-002  0.11387906964953087        1.5946325152817534E-002  -4.5960122403513382E-002
 /print_too_much/ovlp          36 -0.22287650078691174       0.14542483779441284       -2.1730009860346965E-002   7.5483999117368916E-002   1.2801512318931341E-002  -2.3021300852381410E-002  -4.1020190522080502E-002  -7.3549627243986587E-003  -3.2204965104825123E-002  -8.5625415206808775E-002 -0.12115915802481581      -0.14286804011048160       -4.9951916191877235E-002 -0.20593450451731043      -0.12639582120903148      -0.20304082660181549      -0.49890111398845471      -0.48333579399023374      -0.18731649288768920      -0.21982224735568370      -0.48675386649626162      -0.84305710229466457      -0.45147267005496516       -1.0808318184843113      -0.80026663832812350       -1.0732499676581786      -0.44366283406014040       -9.7318247596942450E-003 -0.34485594221561883        2.9003981565240912E-003  0.15006829177798406        8.9948363424144584E-002  -9.0965994300284698E-004  0.14398474045002665       0.46025048236048899        1.6178243291148613       0.17690101356909144       0.19731064829580100       0.22155919147040606      -0.10099926465259287     
 /print_too_much/ovlp          37  -3.3215122778823374E-003 -0.22520548249897662       -2.9061149580947541E-002   5.2968861960482450E-002   1.1028734253426986E-002  -3.6597307951360403E-003  -6.9608853588642994E-003  -5.8953758515878918E-003   3.1478205857807251E-002   2.6748532542341552E-002   2.9921234993583568E-002   3.6701448030887818E-002  -2.2508713851550587E-002 -0.14768962640392036       -7.0973756830844420E-002 -0.14714095887221235      -0.35431307072474971      -0.34611398359932788      -0.13862417752396061      -0.13966693608683203      -0.35068488511005524      -0.60799081857376258      -0.32688922061005621      -0.77396169182601637      -0.58101033431423832      -0.77062894613043342      -0.32310748414182089       -1.8996036786841006E-002 -0.15649192620130845        1.5075519901047926E-003   5.4419789075567389E-002   1.5882224225553451E-003  -5.2429361326836094E-002  -9.6658496121299550E-002  -9.0751751020772750E-002  0.17690101356909144        1.8469218928130791        4.6734521115659802E-002   3.4554944922683767E-002 -0.10738753166503297     
 /print_too_much/ovlp          38   3.3738557117501372E-002   7.7801357346078681E-002  0.13560001277519257      -0.26099523952882020      -0.33005152517008585      -0.37284520021401069      -0.39883805856658516        4.1936980664338402E-002 -0.67514877720356503      -0.87060516176189617       -1.0086352123484543       -1.1080945214624578        1.8202111600579135E-003  0.80624117866810008        3.5587781281042383E-002  0.86738170864703523        1.5853332873446111        1.5461917960106439       0.82817322835585683        8.5249155066743729E-002   1.7158438872259834        2.3396482885537204        1.6449731684504587        3.0343616922650849        2.2407059185722726        3.0062715327643996        1.6167986343115912       -1.5600223464371839E-003  -1.9847238100553295E-002   2.1574990892015332E-002  0.68449471441286647       0.75842985165523358       0.10907257277790755        3.3468950109978211E-002  0.11387906964953087       0.19731064829580100        4.6734521115659802E-002   1.3761360259720254       0.19421272073921206        8.8840659667667771E-002
 /print_too_much/ovlp          39 -0.33212081928524195       0.17711517448492764      -0.25887105996563853      -0.23166522550744861      -0.22874243093287694      -0.20647716174499164      -0.17707443378348464        4.1728574698687043E-002 -0.49541536008044729      -0.51145528006789287      -0.48475191112827387      -0.44057923116713571       -3.6308607232273957E-002  0.75231952063015795       -6.2958120985825899E-002  0.72871872038298247        1.5053846126485468        1.4988972608497262       0.72255136484856786       -7.8649052796604974E-002   1.5122714003669551        2.2458588057478055        1.5003308538214246        3.0287928502192112        2.2213886217035506        3.0163554340725534        1.4882560618385221        1.5670484336338659E-003   3.5738591261121624E-002  -6.8010490549877175E-003  0.62143962893340388       0.12143870855565117       0.10679637735939224       -1.2844341225481634E-002   1.5946325152817534E-002  0.22155919147040606        3.4554944922683767E-002  0.19421272073921206        2.2497497683329164       0.62245062883372881     
 /print_too_much/ovlp          40 -0.17707465723313651       0.29009498403780259       -5.2646541885758513E-002  -7.4264497879753311E-002  -5.7020048378205401E-002  -4.6775417847801926E-002  -4.2171347558710881E-002   1.3084595536881558E-002 -0.14684828027973440      -0.15402112662103706      -0.15578509379802452      -0.15725725015297431       -5.1415701835082622E-003  0.24651174850348256        3.1284270731885555E-003  0.23550467054455682       0.51671617311079143       0.51458811702516871       0.23371592863525614        2.3324768203307339E-002  0.51052834301927974       0.80428933275760550       0.50158909342890290        1.0723991828953374       0.79168374556772148        1.0688466965434853       0.49842228389084853        4.0323493365716778E-003   7.2844397012225476E-002  -1.5956648388376612E-002  0.14224555615132461       0.12722102845875335        1.7004311435452363E-002  -2.3281304890665357E-002  -4.5960122403513382E-002 -0.10099926465259287      -0.10738753166503297        8.8840659667667771E-002  0.62245062883372881        1.3765110668213698     

Eigenvalues of overlap matrix of current wave function and its first-order derivatives:
overlap eigenvalue #    1:  4.04128409E-11
overlap eigenvalue #    2:  3.49776960E-10
overlap eigenvalue #    3:  9.67592446E-10
overlap eigenvalue #    4:  3.16620630E-09
overlap eigenvalue #    5:  3.47388034E-09
overlap eigenvalue #    6:  1.92027332E-08
overlap eigenvalue #    7:  2.27953971E-08
overlap eigenvalue #    8:  3.29572379E-08
overlap eigenvalue #    9:  5.32986929E-08
overlap eigenvalue #   10:  7.02065459E-08
overlap eigenvalue #   11:  2.35292234E-07
overlap eigenvalue #   12:  7.62822171E-07
overlap eigenvalue #   13:  1.38861392E-06
overlap eigenvalue #   14:  1.41148313E-06
overlap eigenvalue #   15:  2.34421005E-06
overlap eigenvalue #   16:  4.76256831E-06
overlap eigenvalue #   17:  7.99627356E-06
overlap eigenvalue #   18:  1.01576211E-05
overlap eigenvalue #   19:  2.68314758E-05
overlap eigenvalue #   20:  7.14709005E-05
overlap eigenvalue #   21:  8.33006742E-05
overlap eigenvalue #   22:  1.24934863E-04
overlap eigenvalue #   23:  2.55315594E-04
overlap eigenvalue #   24:  3.42642980E-04
overlap eigenvalue #   25:  6.25893518E-04
overlap eigenvalue #   26:  5.08901075E-03
overlap eigenvalue #   27:  1.04544610E-02
overlap eigenvalue #   28:  1.43188844E-02
overlap eigenvalue #   29:  2.36306654E-02
overlap eigenvalue #   30:  7.58736524E-02
overlap eigenvalue #   31:  1.67783496E-01
overlap eigenvalue #   32:  2.86675681E-01
overlap eigenvalue #   33:  5.89610651E-01
overlap eigenvalue #   34:  9.18763273E-01
overlap eigenvalue #   35:  1.00000000E+00
overlap eigenvalue #   36:  1.26205840E+00
overlap eigenvalue #   37:  1.40815386E+00
overlap eigenvalue #   38:  1.88161314E+00
overlap eigenvalue #   39:  2.53191160E+00
overlap eigenvalue #   40:  2.80605511E+00
overlap eigenvalue #   41:  6.43000340E+01
 /print_too_much/amat           1  -20.314907609957185        1.1416950646754782       0.70679847253866113      -0.14806473888535762      -0.17133360431369923      -0.15833376163962096      -0.12196889552897286       -5.1911003877231535E-002  0.73661066418840582        1.1089435302185053        1.4393113232006189        1.7312476558812808       0.45369918912698748       0.58509916066259249       0.86138045740462421       0.72230313512118438        1.2957471388333646        1.1027440835644926       0.52217241093205313        1.1988441992968581        1.4471072065960970        1.9424260149266594        1.1608585022594680        2.2436514554557272        1.5678251582593743        2.1520113407536017        1.0596204622308809        2.4101533780000184E-002 -0.28170836114446873      -0.17816083177125955        1.9741653262038938      -0.87900235960414697      -0.15982549863595949      -0.25472456344757738       0.69057002531384826        1.9003247628866056       -1.0679108693596877      -0.15065161635894708        2.8683479352799441        5.6351283632198683     
 /print_too_much/amat           2   1.1730630596303537       -16.833930477680283        1.4911194179212097       0.31916768304344106       0.51247857236799998       0.58120635995738201       0.61083870382973271      -0.13715351813433646        2.1392394122547991        2.7773090852198767        3.2362955871160821        3.6008765570273580       0.57970190736378813       -1.0453505802369407        1.1672638551065539      -0.90995886561053918       -1.7069110351429564       -1.8847980675386478       -1.0997800022073070        1.7365176274167688       -1.6732910679714599       -2.1912012370594574       -1.9633078960495332       -3.4006570060478856       -2.5483494203234147       -3.4719408453636924       -2.0504894262016888        7.9754592538740343E-002  0.36750884985271121       -9.3162458587020713E-002  0.10812883017400772       -1.3342942988619815       -1.3291886417474204      -0.50323874242141331      -0.22807871174193495      -0.46780069728570561       -6.9017859339275747       -1.2586973131493908      -0.62369577404498333       -6.3631760267580324     
 /print_too_much/amat           3  0.77639855231324284        1.4921315660226082       -8.8068842973272172       -3.1186266722414549E-002   3.4008329914636475E-002  0.10713886884964266       0.17346608859475637       -1.6916036526414741E-002  0.41047674173145926       0.68275239482143713       0.90361526884946397        1.0975486358481674       0.13299787930577212       0.21433155982254867       0.23424549712735621       0.17107231643218945       0.50476345585353588       0.47196385549755815       0.14034212900160431       0.28917600777814056       0.41415403660437877       0.77263722820140734       0.37250624762434781        1.0406577985416900       0.73025201248066918        1.0405256921375301       0.37482405272694946       -2.2392805087231103E-002   8.4246192707063891E-003  -8.6143440315982839E-002  0.22400537253585501       -1.2668795887403439       0.57198025796407026      -0.33803310561753414      -0.15549388967501890      -0.17286049452589031        1.8317625233738171       -1.8156825867584554        2.7451793854896351        1.5518210005593280     
 /print_too_much/amat           4 -0.13346085916828820       0.46863752090370098       -7.2275608201735564E-002  -1.1036064969437138       -1.2604742528164159       -1.3360841744669059       -1.3630630390125971       0.17996364566950840       -2.5670864167907852       -3.1269055460578570       -3.4811842372698689       -3.7137126930936497       -9.7498865919498634E-004   3.4100713460528800       0.11638804899215227        3.5683232543411183        6.7969024210875943        6.6630212738275585        3.4349518669064629       0.30165705796512299        7.1898103835229534        10.126223060990529        6.9451239920029435        13.260317920420897        9.7810555875675220        13.160315711372490        6.8463033535548945        1.9920701423331966E-002  0.28289132934436423        5.0215629457405472E-002   2.4064942314753481        2.1123808358086031       0.88774465617964970       -4.1531526721231615E-002 -0.21106885177282519      -0.41236729992443788      -0.60407233240269065        3.8327008014456077        2.3400241350851401        2.4300213843055012     
 /print_too_much/amat           5 -0.16427640179943265       0.51223913618519656        9.9955596623167720E-004  -1.2689771411191362       -1.4839331189241798       -1.5989596165074660       -1.6509649371023045       0.20798351838689388       -3.0294638170727426       -3.7329983605002850       -4.1954369107655047       -4.5088603181044391       -5.1077359210593976E-003   3.9161585278047171       0.12187595860453948        4.1202164539225956        7.7756171296975696        7.6188365879331954        3.9627774375259541       0.32067037620131789        8.2596863245541190        11.551234880756384        7.9768511553922581        15.110543732449491        11.152732687194923        14.994880360696172        7.8606962538556120        1.9219902858256921E-002  0.21072024394810884        6.6643539751183228E-002   2.8701166169398347        2.6387727356314925        1.0221373667149733        5.9654318532620376E-003  -1.7429142612788384E-002   4.7846661214342667E-002 -0.82637871233711380        4.8075256833533322        2.3660401203273693        2.2559938399860986     
 /print_too_much/amat           6 -0.15017576974223096       0.55571914576117321        8.1919227856138921E-002  -1.3445425433356293       -1.5973900418699727       -1.7409839575089610       -1.8136003933335270       0.22138589317839019       -3.2694257944599281       -4.0630554984241476       -4.5990805076954553       -4.9720956256635986       -9.0814216217725008E-003   4.1430419735465414       0.11939254388639320        4.3751193246110569        8.2037971469088067        8.0376233317418926        4.2077500424067482       0.32009921691793153        8.7391133040615347        12.163949953904112        8.4405964000789488        15.901829758042103        11.744901106498773        15.781081393491412        8.3185640624943513        1.9028217436625416E-002  0.15905870004162370        7.7048084677197606E-002   3.0858178082071630        2.9744796440027530        1.0999853501871553        4.2208639896472722E-002   9.7223203006052472E-002  0.33905628764274143      -0.99992073872598652        5.4056439554303211        2.1516576436525607        2.0497106214835643     
 /print_too_much/amat           7 -0.11520404392605904       0.58601566174387787       0.15716988401741899       -1.3711817598846949       -1.6475017183150831       -1.8115128725071201       -1.9008180959725203       0.22655633796034336       -3.3790203289624974       -4.2279024610359457       -4.8139620231853657       -5.2311703886436627       -1.1095811480416051E-002   4.2181820381747324       0.11699396992731970        4.4674070529105503        8.3367790118688401        8.1680031518576595        4.2971864930685948       0.31654061318974591        8.8994824335666038        12.345672901904098        8.5963820402327418        16.130553815671107        11.922251994264208        16.009912242726770        8.4740669088342440        1.9670168464848148E-002  0.12768013616318363        8.3107710580474073E-002   3.1591088847102347        3.1891080791072639        1.1474000312219732        6.6569168917342214E-002  0.15453644406569289       0.50873993640880144       -1.1690426155544600        5.7652465690941117        1.8356752706107109        1.8783204893777103     
 /print_too_much/amat           8  -9.2496915650186140E-002 -0.14442342762460891       -2.8362534186386869E-002  0.17859261228249579       0.20389954348056888       0.21530371755851474       0.21850661497934759       -2.9445180637466525E-002  0.43077304556442647       0.52781906704055725       0.58715985730872688       0.62447066171969468        5.8693261398607838E-003 -0.55270960723280471       -1.0809599361335770E-002 -0.57795158558811588       -1.1016858228028736       -1.0796989738215559      -0.55605112097690268       -4.2001904085704367E-002  -1.1682404415667491       -1.6438898070679206       -1.1270875122210064       -2.1510201829754294       -1.5841370136382917       -2.1324992147520776       -1.1087553337091858       -3.1090026211103156E-003  -4.1620483783697232E-002  -9.7685935185117083E-003 -0.38608150207553571      -0.34514827577784807      -0.13557526307151221       -6.4971379522566763E-003   2.3324286329503535E-002   4.2740369868558684E-002   8.5877935422601248E-002 -0.61632641629573104      -0.42245427761432419      -0.43255678621154625     
 /print_too_much/amat           9  0.99365565131272993        2.0960542728817764       0.41325195881805094       -2.5810086648497355       -3.0324700682861110       -3.2685067196532622       -3.3711810169581335       0.44025787615961720       -6.4412367313167902       -7.9545338325297950       -8.9437196996899146       -9.6081254884544602      -0.13015059973624543        7.9814705136022734        1.2792416497620529E-002   8.3605379280197756        15.794811663515794        15.511246051333181        8.0761790617587543       0.31523660713802748        16.757283571162809        23.429685502611161        16.239025078356335        30.734136749330769        22.682250501716101        30.505176257340558        16.009721320953371        3.2987893197813478E-002  0.34285856656043401       0.15488033863413514        5.7540396623288013        5.4966167770268850        2.1029772562751163       0.18074543028454049       0.11141436358803603       0.47304728185049605       -1.2675705324309376        9.8378102290371618        5.0272011817347177        5.2768483347949857     
 /print_too_much/amat          10   1.3200309185723547        2.6852973914706162       0.67067518998138709       -3.1528097544614124       -3.7543873932363270       -4.0881686176382432       -4.2525113886176742       0.54521665595773672       -8.0107160085592533       -9.9542152418167191       -11.264686277533023       -12.174065310462176      -0.18268508112806781        9.7352103111148072       -4.9086399040314088E-002   10.216266597663765        19.206699786490091        18.876483805983270        9.8846368872528494       0.26157261064601478        20.400265580232620        28.426362754455344        19.799245085205506        37.309439600143591        27.564783382572649        37.048914968191056        19.537527076835449        3.6488575541667623E-002  0.31402102463698167       0.20231224991410096        7.0284130773745428        7.0898286285361918        2.6924700312526797       0.27792280423894744       0.29838578321219855       0.98151465007127259       -1.6805059258142379        12.629765120669958        5.1636687206271503        5.7301499311395396     
 /print_too_much/amat          11   1.6186503104097452        3.1374272350024914       0.90874037918819950       -3.5131388864171682       -4.2244802766366165       -4.6355752211983718       -4.8536656367178850       0.61145020517236637       -9.0376955827095138       -11.293450253998083       -12.848378998991862       -13.953579916573126      -0.21490444140478004        10.830188411986223       -8.4158572147655875E-002   11.388595781604812        21.326054460272523        20.966270541735405        11.027110384963352       0.23500251195587474        22.682727225517826        31.520957168792314        22.027834014029253        41.368758180146990        30.587736183923418        41.090380573211661        21.747814191784553        3.9297812071217075E-002  0.28213006621455161       0.23601987197433016        7.7901622546599345        8.2390300819712117        3.1377548717913823       0.35307390833977692       0.42174437296710732        1.3532018738358107       -2.0508431028473915        14.589018799174225        4.8370529074234163        5.8484170129630009     
 /print_too_much/amat          12   1.8878782490694899        3.5037490721697484        1.1189478568703235       -3.7473319753692995       -4.5402977887939127       -5.0129765280862948       -5.2772831106474003       0.65460766975746998       -9.7311818047808032       -12.219890539083073       -13.965706017356625       -15.230912289061505      -0.23500861240515825        11.532979843798238      -0.10373672967674175        12.150658191777467        22.678973412053406        22.300435777440349        11.770326147610897       0.22346310908207556        24.152958915722468        33.490556280097401        23.462651562792551        43.943543592714761        32.512385622084764        43.655607950345171        23.172870592021980        4.1609050122680902E-002  0.25662663005728881       0.25996898664121426        8.2390689517155629        9.0935764480033612        3.5010109666977227       0.41017591055350699       0.49018217025118238        1.6069805605505061       -2.3573702927544389        15.993915192810716        4.3115176488403790        5.8544881444904284     
 /print_too_much/amat          13  0.67826931318207451       0.40416201372801086       0.22708857712596295        2.2321320508119739E-003   2.0141735805449118E-005   4.2840656850361668E-003   1.1901720024084822E-002   4.3400987888989537E-003 -0.10594084747997670      -0.13770308960224131      -0.14685415310684186      -0.14235161671960808       -4.9923914759810784E-002  -4.0081700624960137E-003  -8.9449863946528652E-002  -1.6317692130417283E-002  -2.1141373972987543E-002  -1.3353900522535989E-002  -8.8190989449011686E-003 -0.11547761494489928       -2.0088151959930745E-002  -3.2706662590658464E-002  -9.3795717783748678E-003  -2.4732866247790863E-002  -2.9366701340247714E-002  -3.2237710505302442E-002  -1.7341682122023666E-002  -2.9848553841944833E-003  -5.5444545218293909E-002   1.0179803808243531E-002   1.5474402002745968E-002   1.3227903055903718E-002  -6.3383677982673692E-002   9.5183068508158775E-002  0.10480300580708363       0.29749357105119956       0.14359143548677600       -2.5879869961166480E-002  0.38805383293989015       0.28176645648619703     
 /print_too_much/amat          14  0.59568717884452027       -1.5141816522930294       0.36776934382487514        3.4355168307531181        3.9233866352038285        4.1538654426761328        4.2311362034358133      -0.56168208830478594        8.0010791395282990        9.7295624112434940        10.813650074267992        11.515060379164234        7.9631268367220676E-003  -10.642515732689363      -0.34868528273477040       -11.128933662987096       -21.208754301986808       -20.790997831053772       -10.712427323902016      -0.91386967341532266       -22.426752744366826       -31.589550259900896       -21.665938855069903       -41.377831173611838       -30.513852897419984       -41.064325282130440       -21.355642232089231       -6.0263426023113964E-002 -0.84391224664804265      -0.15165165070029774       -7.5647047038717865       -6.5102425864667115       -2.7272569328981380       0.13161267109335817       0.53432060914628865        1.1136348599719510        1.7557836759444285       -11.866340176124528       -7.5871716884860678       -7.9078635054119442     
 /print_too_much/amat          15   1.4023622374775140       0.70626276409818356       0.47482467959735219       0.12913544800773180       0.14495408824851058       0.16526697428415815       0.18759371140569112       -1.5969651639488186E-002   9.7990549542196659E-002   9.6188425580087278E-002  0.12501165648588833       0.17138501161487302       -8.6827736496973729E-002 -0.39410967448579010      -0.15441893018101394      -0.42604369485711635      -0.79728121277932895      -0.78005321652567194      -0.40990901024814669      -0.19331558133478088      -0.81856251674817382       -1.1718060437150242      -0.79443411502331385       -1.5325713501312919       -1.1641941666188522       -1.5496091425112488      -0.81316374990991391       -7.5916023734198118E-003 -0.15427880012438072        1.5045540237938091E-002 -0.20090466423010678      -0.24211522601763025      -0.26982197530675062       0.20501809632527779       0.26788072985502831       0.74395430511399296       0.27593262572610788      -0.51922147952524789       0.72072808082057005       0.28439811123673708     
 /print_too_much/amat          16  0.88309567406169409       -1.4339499120972659       0.40687573766216295        3.5945551251951970        4.1295795544277105        4.3946360828260165        4.4967166572068153      -0.58899361645264792        8.3922294016309937        10.231181700940478        11.407507725021960        12.185603097848235       -2.3982062666261239E-003  -11.128207083214562      -0.37295980051993283       -11.652207932894129       -22.153039775836241       -21.719026638925254       -11.219156181660734      -0.94680599480433347       -23.437669154171161       -32.963703017649465       -22.650884276150396       -43.177892122711164       -31.858760250785149       -42.861014222222934       -22.336828298564249       -6.1546211151986330E-002 -0.84103037994881746      -0.16151923418145034       -7.9355861749372538       -6.9774899108160815       -2.9147525871491311       0.14600871152420580       0.48980446444781522        1.0310590395729884        1.9522688296754129       -12.731394540924708       -7.3562531005946585       -7.7029730159058172     
 /print_too_much/amat          17   1.4744058452107183       -2.9280035838604883       0.90267716876745585        6.8479205194924848        7.7938494242374965        8.2347354783082185        8.3777176467598569       -1.1194961869327880        15.856423702055789        19.241180729549445        21.359057974595508        22.729197269279211        1.2001552477283184E-003  -21.217199710192933      -0.71625154323984486       -22.167540066357883       -42.301996600718056       -41.476899660072498       -21.346328819185473       -1.8435414473613445       -44.691024493136617       -63.022181707237117       -43.186845660150595       -82.576549878255918       -60.901489668375589       -81.963401844921691       -42.582154956006811      -0.12379848542292460       -1.8005300505816551      -0.29055890627638126       -14.963348615360099       -12.825189837895248       -5.4883823360200088       0.34901866846346863        1.2999807853839749        2.8124492920008732        3.3904150841067526       -23.365956985770843       -15.095107750529849       -16.116372171525697     
 /print_too_much/amat          18   1.1139008435270168       -2.9971837807390456       0.78109469427913136        6.7055787519769217        7.6244825097311821        8.0485479658575496        8.1811312846664812       -1.0957426309720963        15.545808761278375        18.865274784225086        20.933077669183895        22.263738592143593        1.3415690721908513E-002  -20.774762764182977      -0.68850047285261784       -21.702054098613264       -41.428668961056367       -40.618270598080080       -20.895509610493402       -1.8036311993056042       -43.772853123264397       -61.738366905592230       -42.292447654486359       -80.886604237769433       -59.645572375379039       -80.277532509084125       -41.691731770531959      -0.12230012074073038       -1.7704393426572143      -0.28735950297054141       -14.646554641358474       -12.528834657872341       -5.3299501009996009       0.31060313567160525        1.2869647376540221        2.7331496184866904        3.3357596317703582       -22.798090071975359       -15.034305984970725       -16.015008043496998     
 /print_too_much/amat          19  0.51828785192890880       -1.5034248799607508       0.28245110436669940        3.4523761236832473        3.9602805507081107        4.2084587802222497        4.3001053765634101      -0.56545173214082556        8.0829944846585509        9.8565058918655577        10.982811019751420        11.721505472874759        1.0514416738758747E-002  -10.686358223714047      -0.34353714065527985       -11.186840675533874       -21.280615226361302       -20.861666951636259       -10.768798356659362      -0.90403264675714112       -22.519647822389445       -31.680930766736207       -21.757381158924229       -41.490449186522213       -30.604841703460252       -41.177849634538070       -21.447470379069394       -6.0033753408712531E-002 -0.81100098840472157      -0.15829000021343428       -7.6183540021301326       -6.6812935229045536       -2.7571377078191830       0.10687433749698960       0.47725406860027531       0.95139594578724185        1.8963482304572032       -12.162774018324239       -7.2971993083951432       -7.6110506306369885     
 /print_too_much/amat          20   2.0956030419454761       0.90834649760824870       0.71777874863728997       0.32496459301775793       0.36823056348947070       0.40899197776770740       0.44798462713975784       -5.1470728507656299E-002  0.48373152643760065       0.54273208413583285       0.63203933119665823       0.73723658745822318      -0.10749364945777540      -0.99884523154937654      -0.18573266045442272       -1.0481326016350780       -1.9898785361096591       -1.9669499979750851       -1.0275191692287817      -0.21846348147121075       -2.0380780589060103       -2.9182315368281371       -2.0071716484491056       -3.8677320198903753       -2.9164254529814007       -3.8979381297558175       -2.0409321937874805       -1.2897697312985207E-002 -0.28032729259789824        1.6668338080418430E-002 -0.54181386379959140      -0.63665622457410487      -0.54622917653555869       0.31681398405809980       0.47238133875228883        1.2898704890481381       0.35869166879876158       -1.2484482678261093       0.97435971385817455        4.8048155585747043E-002
 /print_too_much/amat          21   1.8520876864316711       -2.9533572197779474       0.95558945712901633        7.2481989683721961        8.2908574141125886        8.7949517259835002        8.9778757953230279       -1.1899583707045458        16.852460104115153        20.487176146100207        22.795723611219703        24.312829892723638        2.8935612076264050E-003  -22.450426436041369      -0.73048913594519060       -23.472412193691660       -44.712813775885344       -43.850748746802466       -22.613690449264631       -1.8725925898247855       -47.249827294989217       -66.552272880855739       -45.684593662629311       -87.221634144994539       -64.354456382665830       -86.592261155040561       -45.062998047442569      -0.12647209142393656       -1.8106951118679984      -0.31447622829018884       -15.894892943855709       -13.836212453351017       -5.8869360631877869       0.35657554419786930        1.2088033324202905        2.6362809328875692        3.6753575942579779       -25.233212413834572       -15.185729470379609       -16.207579157670356     
 /print_too_much/amat          22   2.3571076640732995       -4.3808181837757960        1.4670818969363819        10.198461845557626        11.580582239255648        12.218516859697612        12.420746610755756       -1.6703254386889275        23.551078918558272        28.535005576271150        31.650447270654915        33.667709354292057        4.2204077588960942E-003  -31.601439274858677       -1.0498562015094652       -32.987969964568421       -63.016180949405040       -61.804967974001357       -31.784401313153605       -2.7068954294842289       -66.521994324922019       -93.889616674486689       -64.310938578521359       -123.07169422979558       -90.778117030655437       -122.17708515864597       -63.431646126566370      -0.18816250562836956       -2.8104961116229825      -0.42237245136443968       -22.141288738666713       -18.969595574147430       -8.2318138403491741       0.59587049564142958        2.2038456994151225        4.8447812268035166        4.8840270097794942       -34.519237644788490       -22.398345652835769       -24.475783772294029     
 /print_too_much/amat          23   1.2779063734047851       -3.0413625086687652       0.76271274456913585        6.9887588069416360        7.9834531273981666        8.4576631609933752        8.6219617738496730       -1.1451366035927601        16.274978720296577        19.790502953540148        22.006803834134619        23.450661441637177        1.7171291511100195E-002  -21.644729679652222      -0.69737555198600898       -22.629702342647199       -43.125318977387479       -42.285452742379405       -21.792934682863489       -1.8272965071955909       -45.587759903723693       -64.220634824135459       -44.058121211902979       -84.141301169912609       -62.063136040380584       -83.516480782067362       -43.440621662313788      -0.12336397982001479       -1.7470628175443288      -0.30796425051245585       -15.332670736562774       -13.290687234948008       -5.5988766129444816       0.29125966602443176        1.1631003374013149        2.4526555701295849        3.5823494522045678       -24.203143150853904       -15.090136336632121       -15.955812327079080     
 /print_too_much/amat          24   2.3878824064465300       -6.0004791887264890        1.7523066242563943        13.344538471028617        15.124975868001421        15.930141979816039        16.166101679196668       -2.1814325287028300        30.823079806568067        37.333192613872498        41.367209218819596        43.951569400702226        2.4038482084955071E-002  -41.354462504689295       -1.3703902054898744       -43.160653376484966       -82.501786773101273       -80.901770904170419       -41.570426856591155       -3.5872298561656599       -87.104021379346335       -122.97940194562776       -84.176815077536148       -161.16902817530192       -118.84388545966712       -159.96826694336173       -82.995822222244612      -0.24816494762173891       -3.6949615667695519      -0.55577601553063205       -28.987528710653411       -24.637475447311449       -10.634469743296499       0.71683211290765447        2.8947380394196305        6.2691620858330843        6.3618194332093125       -44.806484088345016       -30.261715444347004       -32.715625132970260     
 /print_too_much/amat          25   1.5684000388491015       -4.4586054111903373        1.1939108534450509        9.8296601328377644        11.144380627476503        11.741646797808496        11.919494525362897       -1.6058548944069173        22.721635355875673        27.537746716979441        30.525488469374437        32.442931472074804        1.8332018533780747E-002  -30.454285206892013       -1.0157166564486921       -31.791822729668969       -60.757940063594830       -59.576258011014566       -30.617501550186798       -2.6619143733176855       -64.161320908548987       -90.574615803002743       -61.997269459510299       -118.68539124230263       -87.516514618767346       -117.79696422874922       -61.123610868950287      -0.18472484598524885       -2.7278089000912091      -0.41295538573039808       -21.337129750117086       -18.207746107758968       -7.8270201375599928       0.50941546363858081        2.1546387865487500        4.6145914614239825        4.8146822723236724       -33.079932706785577       -22.182286512125160       -24.065410454204915     
 /print_too_much/amat          26   2.1719893639439318       -5.9898329710867007        1.6709764761530941        13.235214541864170        14.996168487708712        15.790532991758271        16.020727015429273       -2.1618351616406071        30.571418134245285        37.032853118147187        41.031425207354857        43.589226440130844        2.4019268403839167E-002  -41.013144994237919       -1.3690181615549211       -42.807252864017059       -81.831305993452276       -80.238685541921384       -41.224372206891346       -3.5869549562721375       -86.405506141236970       -121.99642646703445       -83.489854240267931       -159.86378060718368       -117.87415910566790       -158.66466754497233       -82.310403223033518      -0.24785305855874329       -3.6759624004413300      -0.55286145590279689       -28.745498223856831       -24.421191642467054       -10.517840002385809       0.69546381450174088        2.8913855498869694        6.2220784987725146        6.3869151959522483       -44.397041601535093       -30.141769727030862       -32.559905195532977     
 /print_too_much/amat          27   1.0576619544245616       -3.0298812611323092       0.67796153948313964        6.8795901102660046        7.8546579045879099        8.3179913084815738        8.4764758327590055       -1.1257489841925619        16.024545038123883        19.491154302434083        21.672006360340184        23.089322812815631        1.7776421024459976E-002  -21.303996877278465      -0.69448949069331078       -22.276440597271804       -42.455834153854958       -41.623670473877425       -21.447289983832128       -1.8244005865842148       -44.889498247362468       -63.238944088061743       -43.372082042780988       -82.838756394264706       -61.095518634210208       -82.215726051899054       -42.756218161810217      -0.12306610917871108       -1.7281571708889412      -0.30500803488469430       -15.090261387438126       -13.074400372693184       -5.4828909301153326       0.26910448879069171        1.1601223819791251        2.4044950083261862        3.6086430442856710       -23.792731265617764       -14.972321315980569       -15.810456363100048     
 /print_too_much/amat          28  -8.8963914350100783E-003  -4.0235498278583276E-002  -6.0133494611186537E-005  -7.4602117249654432E-005  -2.4058211709496700E-004  -1.5869203083390239E-003  -2.7064564416722152E-003  -3.9132898671016998E-004   3.0747184404727457E-003   2.2979637815462449E-003   1.0087426622233629E-003   1.7330582398629880E-004   3.3473141869683858E-003  -9.9299278394505974E-005   7.4889615747009142E-003   2.4818898785674975E-003   2.0374121395546479E-003  -6.4539839961907930E-005  -3.7298101388430482E-005   1.1957356309046159E-002   5.0986521374347948E-003   5.8272487579624221E-003   1.8439432197093764E-003   1.7992271512206337E-003   1.2036069614984879E-003   2.4739773150281308E-004  -3.1396563248012321E-004   8.9513581947126258E-004  -5.3163532962260812E-003   1.7450021910270577E-004   1.2189807833532313E-002   6.0641051467854035E-003  -7.7235648253034909E-003  -2.1196880968141880E-004   2.3319180478981467E-002   5.9091934470779359E-002 -0.11907752636530716        2.0353468558837108E-002  -6.8711438885473639E-003  -4.1734301702968928E-002
 /print_too_much/amat          29 -0.53161441310555058       -2.3331830581927421E-002   4.3425639405154211E-002  0.18569265344203356       0.10773059070460314        5.1673684184077744E-002   1.6596691323467406E-002  -2.7565631408303114E-002  0.18395287193865315       0.13252495931174985        7.9255948081744160E-002   4.2515021101488207E-002  -2.5617213095982999E-002 -0.54923631848738674       -8.6321403501823157E-002 -0.52065820741620372       -1.1898110403622499       -1.1863362561177639      -0.52221708542269130      -0.17851579638719128       -1.1600584850080968       -1.8853545564442800       -1.1338011674515129       -2.5086221381598603       -1.8601521038694995       -2.5119810165392034       -1.1446447833240823       -1.8777230823729726E-002 -0.55170271454920883        4.6146637099984660E-003  0.11014570548991594        4.4242101191397765E-002 -0.15070242761216568       0.27015397565663596        1.0891982839378711        2.5034615999505991      -0.14118808674518091       0.24188341302230690      -0.18490262895163512      -0.97872912449762617     
 /print_too_much/amat          30 -0.18742408815642397       -7.0671376445276685E-002  -9.2354871661575619E-002   5.0744370273429026E-002   6.5943113115080013E-002   7.5463909718415587E-002   8.1111315827664149E-002  -9.4775302966183285E-003  0.15103194543047169       0.19698122389641007       0.22955663283204658       0.25296733357217799        9.0128511439297454E-003 -0.15115398009914957        1.2162329048143738E-002 -0.16150064126671593      -0.29134881340981822      -0.28792350994372212      -0.15820721303764612        1.1677527001932510E-002 -0.31587161824757432      -0.42541481774010759      -0.30868914205018549      -0.55890947995254003      -0.41522347594874498      -0.55596880220571543      -0.30594453883389039        7.5469273808128484E-004   3.2658943317980284E-003  -9.7636427409570469E-003 -0.11702805719187258      -0.16212986371103910       -1.2784196731126501E-002  -1.3644171005764372E-002  -1.7345351842445102E-002  -2.0679382964689470E-002   1.8450345828880943E-003 -0.29330929212998258        1.3787519800779149E-002  0.27723893170657871     
 /print_too_much/amat          31   2.1470968868134932       -1.2146057312152361E-002  0.34090211013488908        2.4973486890882737        2.9258225880758997        3.1390308936700233        3.2121871629822141      -0.39259171405268395        5.7943934105460642        7.0726980885564252        7.8445703747868638        8.2993434969607254       -4.7504266269023016E-002  -7.7545798142253579      -0.35221616983568121       -8.1777205420203636       -15.411295212402669       -15.059925411256186       -7.8222339333759967      -0.79729908649194559       -16.410032487083331       -22.893497557361666       -15.790276629459294       -29.882921911454552       -22.010174693038710       -29.618193517616792       -15.520685760403726       -4.0228052472536427E-002 -0.20664152978681513      -0.12184806423303324       -6.7394081319193981       -5.1644044688545971      -0.77537548328524464       -4.1855199253700459E-002 -0.62124084793647449       -1.4277717730221502        2.1362239131709204       -10.032845024952836       -6.6148868088194508       -5.0003940937477154     
 /print_too_much/amat          32 -0.82386600221540285       -1.3561853992234914       -1.2634848019181835        2.1481510991286479        2.6522463297371179        2.9885086023137850        3.2033888934961063      -0.36019213821693769        5.5682955483287131        7.1142558353329051        8.2270768882330358        9.0592546399950002        5.8308617517119680E-002  -6.5583251918392591      -0.12126953541376874       -7.0031580723629521       -12.902798624257901       -12.647641091977917       -6.7452771945847925      -0.43454484896175688       -13.873042352695277       -19.079207263135437       -13.403036198873538       -24.874766148540218       -18.423352647452333       -24.688277540511393       -13.213118771687311       -4.5323224344144070E-002 -0.16350300544081375      -0.16815426655937241       -5.1334770552166491       -6.2321211121459656       -1.2211153835382595      -0.20099685288980035      -0.28897501891284927       -1.2329156040327494        3.0371215593974124       -10.893359952889449      -0.97937423931843037       -2.8084171895915944     
 /print_too_much/amat          33 -0.17121068633791894       -1.3445413441437710       0.61514638790586940       0.84160576879626492        1.0026083109710713        1.0930587374165095        1.1481115779503548      -0.15188064128660317        2.1420736282402308        2.6796582507890281        3.1063989467999504        3.4651453936198409        3.2422725411499528E-002  -2.6016485111805800        2.2999521594554428E-003  -2.7178684066460219       -5.1327785638470509       -5.0529273509749606       -2.6378697195614560       -4.1558561067392197E-002  -5.4035335358422856       -7.5713216326594841       -5.2670308661214662       -9.9749072583791580       -7.3848536939259297       -9.9248453366073406       -5.2167677677797428        6.1408999371374545E-003  -6.0381908477562002E-002  -6.7959044859806795E-003 -0.75298655073238385       -1.2227941266772466       -2.9117560576025223       -8.2831110575947192E-002  0.19638376678881025       -6.8720459746109078E-003 -0.68296929556223440       -1.6280352683367274      -0.83500981241332262       -2.4043904905146780     
 /print_too_much/amat          34 -0.40186680569427391      -0.77479155972452540      -0.46176173807966336       -2.5726743948814100E-002   4.0116561169705231E-002   8.3796818977303872E-002  0.11114296553135039       -1.0668447948900154E-002  0.27097982456442377       0.40544613715555983       0.50773796905738811       0.58522246099889796       0.10315817673176172        7.8993049404505522E-002  0.21841100519790363        8.1890508756884001E-002  0.25649452897269664       0.22608498154699891        5.1095985097597635E-002  0.33639432839344530       0.23958693508753359       0.47481701616890021       0.18473322690957095       0.56800667948786399       0.40901475432723738       0.55701328005913731       0.17337571739340890        5.7712498999802232E-003  0.29908780600030760       -1.8647173847596715E-002  -7.1875760781646425E-002 -0.30345848570657669      -0.11487885523024562      -0.42542683676654164      -0.57902751253317752       -1.4152179691554823       0.31079238577761803      -0.45670503947351870       0.12546915941064690       0.22066947306744036     
 /print_too_much/amat          35  0.87922012564010399      -0.20002918991062152      -0.29646801799017664        1.0183545740852840E-002  0.18863541456212302       0.30277621327494819       0.36197530294446184       -7.7159436326066401E-003  0.49383025312272644       0.76213532129844264       0.94848894718946297        1.0574052161505656        5.7677731172588731E-002 -0.10824372439922803       0.12855045845060120      -0.19899463067478129       -4.9588651183850141E-002  -1.5226158522709143E-002 -0.16032129501848458       0.21874664772708835      -0.22402079786615325       0.11378116313546482      -0.18696051234494249       0.21978729543411823       0.18289503286314135       0.25362158434639337      -0.14693814870964372        1.7758378174579509E-003  0.83698110923800539       -2.5338706431655613E-002  -1.0020768161433447      -0.66125285584676119       0.21703904956663378      -0.60668661707425420       -2.2789474738446431       -4.2596156801739795        1.9034288814939595       -1.6786713385053569      -0.31973099590534093       0.50934896397142004     
 /print_too_much/amat          36   3.1554606653219603       -1.9217494851736943       0.30460527297903911      -0.22959599544556153       0.26634102619687350       0.55185447840696145       0.69573207546566118        1.6904666630661036E-002  0.85878026498288274        1.4445819374335935        1.8639279653585874        2.1201227640835558       0.23742981153678022       0.52300731051015115       0.58968088684695941       0.34876383741520778        1.5833242933905831        1.6185843246914606       0.40001464447176405        1.0559114039638797        1.2400122385998040        3.0036548296083154        1.2087531315596438        4.0000848336545047        3.0414048744954787        4.0767000821536676        1.3092757368515304        9.8896153444162471E-002   2.6915793913050710       -2.7360784058060528E-002  -2.1226760807575449       -1.3087796042960578       -8.1579001614344900E-003  -1.4995209821756990       -5.1048147450517423       -12.592204773216450       -1.3593440243391366       -2.8333119949745313       -2.7317514895723907      -0.27840202485490484     
 /print_too_much/amat          37  0.11424205290013116        3.2645773039076715       0.39445992919637879       0.35403011900061010       0.15303880425696131        9.4791541947029628E-002   7.5218863626165089E-002  -2.1694719444810620E-002 -0.13648731698660899      -0.35685773044334707      -0.49674910791492349      -0.61653400751974574      -0.34557465324784375       -1.0362660414205185      -0.79570387655179653       -1.1414174508159700       -2.4796384319954621       -2.3108753032503868      -0.94281823723361624       -1.3050204986776772       -2.5764777798169085       -4.1614677391388799       -2.2911234475886246       -5.1252642118222056       -3.7730280863569536       -5.0095381294214612       -2.1325147901572779       -6.4532580927952221E-002 -0.76393275437887320        1.4957644364113697E-002 -0.81788111366548677       -5.4592708883834141E-002  0.73158723671809411        1.0070610415916637       0.50970104359630564        1.5506431975570627        5.5604338743841026      -0.66065012667350764      -0.58120545779590049      -0.22086541574112761     
 /print_too_much/amat          38 -0.32796926854466979       -1.0120212548534697       -1.9350221692473517        3.8427200845891605        4.7486194122181393        5.3437455270167931        5.7087070241176754      -0.63268226928813165        9.7780708683996398        12.456851805281547        14.352092515179354        15.726078903680445        3.6422530014769194E-002  -11.748190030891983      -0.34166952534748207       -12.566784098944970       -23.122810511502617       -22.639127578975863       -12.077581687620302      -0.94355231181783683       -24.884423424994136       -34.179532725078069       -24.003616622129986       -44.518016650402799       -32.960266238469025       -44.178219966503320       -23.656966662922578       -7.8899350881060273E-002 -0.20694125755906528      -0.31111471659344847       -9.7506484751924773       -10.802705363890556       -1.5490185551065938      -0.32166643451111243      -0.93488901607577524       -2.4762585722988528        5.8688175248483105       -19.708354879044244       -2.0830395975571205       -2.5617804373437636     
 /print_too_much/amat          39   4.9807328462037566       -2.5441697416912894        3.8503541350635566        3.0425257436001010        3.1151074362650029        2.8834318222982658        2.5250931509482126      -0.47080051309447751        6.1501142632229495        6.6413644404320893        6.4947862608544975        6.0316132044197017        2.7156604677771257E-002  -9.7675600035500452      -0.24438786517081335       -9.8759519316486983       -19.736323690519153       -19.305981659831467       -9.4425494514121144      -0.71167473412330740       -20.447635155743143       -29.592565187039586       -19.718181469863268       -39.022260144092797       -28.445082213655514       -38.603284954933237       -19.296594784921961        5.0879300061440480E-004 -0.50449518533422633        5.9749063338856337E-002  -8.1393744934844516       -1.6675655824063560       -1.4858849045769240       0.16137327712397764      -0.30272592573371748      -0.86108062439685540       -4.4355763144517777       -3.8021642928295805       -22.956778130673918       -19.900565291831587     
 /print_too_much/amat          40   2.7639306560179402       -4.3959153469275414       0.79543057587974597       0.74545669419639093       0.72594426278250002       0.71049895920665695       0.70048803160217255      -0.13786062023905504        1.7807408382590086        2.0257241482494099        2.1619116050458080        2.2623394825533008        4.9369530015257707E-002  -2.4959223010317531       -1.0163494902015666E-003  -2.4996574936086891       -5.0947483893914827       -5.0074534610241415       -2.4067374237939614      -0.13737445349283375       -5.2209149864615441       -7.7566564243216112       -5.0659990886167661       -10.218195874214150       -7.4962616071555503       -10.110278079026983       -4.9501862508422736       -9.3078132000836630E-002 -0.57576160633031170       0.12545052809213958       -1.6598074840616055       -1.7627110901140837      -0.27886227617995257       0.24952277005908488       0.58176619616524872       0.20803226505022465        2.1769702688356176       -2.2534124592900100       -2.4414014907053936       -21.886946130554421     

Solving generalized eigenvalue equation of linear method with a_diag =  0.0D+00
eigval_srt_ind_to_eigval_ind=     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35    36    37    38    39    40    41
eigval_srt_ind_to_eigval_ind=    31    41    39    38    37    40    28    29    36    35    33    34    32    27    26    30    24    25    23    21    20    22    19    18    13    14    17    16    15    12    11    10     8     9     6     7     4     5     3     1     2
eigval_ind_to_eigval_srt_ind=    40    41    39    37    38    35    36    33    34    32    31    30    25    26    29    28    27    24    23    21    20    22    19    17    18    15    14     7     8    16     1    13    11    12    10     9     5     4     3     6     2
Sorted (complex) eigenvalues:
eigenvalue #    1:   -15.276502 +    0.000000 i
eigenvalue #    2:   -14.596965 +    0.000000 i
eigenvalue #    3:   -14.369736 +    0.000000 i
eigenvalue #    4:   -14.322736 +   -0.129936 i
eigenvalue #    5:   -14.322736 +    0.129936 i
eigenvalue #    6:   -14.144652 +    0.000000 i
eigenvalue #    7:   -14.047581 +    1.768800 i
eigenvalue #    8:   -14.047581 +   -1.768800 i
eigenvalue #    9:   -13.461429 +   -0.023173 i
eigenvalue #   10:   -13.461429 +    0.023173 i
eigenvalue #   11:   -13.258676 +    0.651631 i
eigenvalue #   12:   -13.258676 +   -0.651631 i
eigenvalue #   13:   -12.825652 +    0.000000 i
eigenvalue #   14:   -11.170127 +   -0.884068 i
eigenvalue #   15:   -11.170127 +    0.884068 i
eigenvalue #   16:   -11.048679 +    0.000000 i
eigenvalue #   17:   -10.043894 +    0.852310 i
eigenvalue #   18:   -10.043894 +   -0.852310 i
eigenvalue #   19:    -8.640751 +    0.000000 i
eigenvalue #   20:    -7.942408 +   -3.746713 i
eigenvalue #   21:    -7.942408 +    3.746713 i
eigenvalue #   22:    -6.893737 +    0.000000 i
eigenvalue #   23:    -5.448281 +    0.000000 i
eigenvalue #   24:    -3.358200 +    0.000000 i
eigenvalue #   25:    -3.294472 +    9.071358 i
eigenvalue #   26:    -3.294472 +   -9.071358 i
eigenvalue #   27:    -2.185094 +    0.000000 i
eigenvalue #   28:    -1.271002 +   -5.581547 i
eigenvalue #   29:    -1.271002 +    5.581547 i
eigenvalue #   30:     6.326605 +   -1.202377 i
eigenvalue #   31:     6.326605 +    1.202377 i
eigenvalue #   32:     9.693107 +    0.000000 i
eigenvalue #   33:    12.992228 +    3.807513 i
eigenvalue #   34:    12.992228 +   -3.807513 i
eigenvalue #   35:    20.134110 +   17.370295 i
eigenvalue #   36:    20.134110 +  -17.370295 i
eigenvalue #   37:    39.667929 +    7.945440 i
eigenvalue #   38:    39.667929 +   -7.945440 i
eigenvalue #   39:    44.011257 +    0.000000 i
eigenvalue #   40:    62.755064 +   36.712597 i
eigenvalue #   41:    62.755064 +  -36.712597 i

Reasonable eigenvalue window is:   -15.21648   -14.58696
(sorted) eigenvector with smallest norm of wavefn variation is    2:   -14.5970 +    0.0000 i, norm change=  0.000000
(sorted) eigenvector with largest first coefficient is #          2:   -14.5970 +    0.0000 i, norm change=  0.000000
(sorted) eigenvector with lowest reasonable eigenvalue is #       2:   -14.5970 +    0.0000 i, norm change=  0.000000
selected (sorted) eigenvector is #                                2:   -14.5970 +    0.0000 i, norm change=  0.000000

add_diag, delta_lin1=  0.0E+00  1.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00

add_diag, delta_lin2=  0.0E+00              0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00

CSF parameters variations=   0.0000000   0.0000000   0.0000000
=======
Adding to ovlp_lin  1.0000E-08
 /BM/ovlp           1   1.7807607672868184        6.5249273503136695E-002   7.9022232776588375E-002   8.1359236353403475E-002   7.8002887434610813E-002  -4.5960992401709588E-003   5.5119444377309534E-002   3.5144525351986644E-002   9.1747249277460696E-003  -1.9552521371844378E-002  -3.1584853028315374E-002 -0.21196239220750357       -5.6482843700315755E-002 -0.22493218666117709      -0.42302653608943253      -0.40683719340000124      -0.20873824305372413       -7.2809648803280780E-002 -0.44493963260686087      -0.61638754544271279      -0.42162564105244904      -0.80230122352797650      -0.58603774883203918      -0.79524879908880663      -0.41452344708162059        2.7022547625387468E-003   2.6789182021430318E-002   6.0025210191860766E-003 -0.39420256109380813       -5.2688392992979963E-002  0.13704565617402450        3.7617884510980755E-002  -7.9790296009032979E-002 -0.10319865678939415      -0.17560072553267528      -0.28031107875241290      -0.53139686918457851       0.23920496314741940     
 /BM/ovlp           2   6.5249273503136695E-002   8.6155157224072687E-002   9.2396647383079511E-002   9.5230672163545194E-002   9.6225032814807321E-002  -1.2853612368854517E-002  0.17760930319893120       0.21451107326020491       0.23885279428581896       0.25580419968034107       -6.5022661309934726E-003 -0.26031686880429561       -2.3981147028110605E-002 -0.27316463261002610      -0.52657261221482798      -0.51519850091779062      -0.26175043153891409       -4.7408918274673084E-002 -0.55553436300970560      -0.79272835172740486      -0.53418604584736329       -1.0339228794987605      -0.76374652378179064       -1.0262687087715960      -0.52647274664034427       -1.5087367033199059E-003  -4.6317566913179342E-002  -3.1154688277532927E-003 -0.16442437092405693      -0.14303299563766536       -6.6302188280362745E-002   7.8367336517800867E-003   1.6191952670934595E-002  0.17376610132910453        5.6476098773644234E-002 -0.26207900366862474      -0.16117126407761417       -1.7081715085883742E-002
 /BM/ovlp           3   7.9022232776588375E-002   9.2396647383079511E-002  0.10422026282107227       0.11041747026260396       0.11334364850887368       -1.3905356647120115E-002  0.20085862032093615       0.24671099718814560       0.27759821656579220       0.29937765375086656       -6.2381508762348048E-003 -0.28100605131773193       -2.3249619263879140E-002 -0.29706007840361082      -0.56354310235974481      -0.55062908977519953      -0.28408650122823786       -4.5631683643394538E-002 -0.59860673651463969      -0.84249125586930518      -0.57508876975086309       -1.0974031399545652      -0.81032305128361770       -1.0887245563954480      -0.56632816300324862       -9.5780924301930304E-004  -3.2519907644469048E-002  -3.9294971562369718E-003 -0.19299498245862390      -0.17273054526881548       -7.5616986764080185E-002   6.1599432611259917E-003   1.3352052478130538E-003  0.11560387817351880        2.3278551927173319E-002 -0.31875422679916832      -0.16406748533997728       -1.1029446773254570E-002
 /BM/ovlp           4   8.1359236353403475E-002   9.5230672163545194E-002  0.11041747026260396       0.11890805901142359       0.12335940382206445       -1.4414370518561093E-002  0.21310852156381088       0.26456412134211860       0.29996296947732048       0.32538032280686480       -6.1791956127317249E-003 -0.29030453801665601       -2.2958339702867647E-002 -0.30835826726103477      -0.57948667706135382      -0.56582799659220484      -0.29462965966487786       -4.4732555286543629E-002 -0.61804170924081347      -0.86316985803730972      -0.59347949553901458       -1.1232363772386975      -0.82957004228413211       -1.1141664515955938      -0.58431533980167671       -7.1732622098420740E-004  -2.4706188231380155E-002  -4.4440065961954267E-003 -0.20680537130113308      -0.19255443609613376       -8.1055786472355840E-002   5.2725974854381553E-003  -8.2532161208573873E-003   8.2174581623534859E-002   1.2281660409638073E-002 -0.35514450037641032      -0.14964751208018146       -9.2841047224596629E-003
 /BM/ovlp           5   7.8002887434610813E-002   9.6225032814807321E-002  0.11334364850887368       0.12335940382206445       0.12898798608440304       -1.4625361614571730E-002  0.21892740521703757       0.27378831632312028       0.31227789430636221       0.34044404478859747       -6.2200459701768906E-003 -0.29339651152294266       -2.2881770517432898E-002 -0.31264258592851490      -0.58413584967408383      -0.57021276209832195      -0.29864391767074494       -4.4281749784488511E-002 -0.62452619731971026      -0.86840297928731047      -0.59961053733232461       -1.1292660652582640      -0.83441794073482356       -1.1201623956536650      -0.59040608027652297       -6.2174950184193389E-004  -2.0817891727766324E-002  -4.7268955121501643E-003 -0.21192919812000666      -0.20596517036890716       -8.4363425096324496E-002   4.9848807115722114E-003  -1.3821532016301763E-002   6.6210852376570983E-002   9.6066154725999997E-003 -0.37807691052239856      -0.12868996044915049       -1.1356566458279538E-002
 /BM/ovlp           6  -4.5960992401709588E-003  -1.2853612368854517E-002  -1.3905356647120115E-002  -1.4414370518561093E-002  -1.4625361614571730E-002   2.0946047766894679E-003  -2.8142792046807363E-002  -3.3728015701344649E-002  -3.7434500046984454E-002  -4.0053923994471496E-002  -1.9567364603345627E-005   3.8931239007275575E-002   1.2658115190549712E-003   4.0287880243468965E-002   7.8067707377426299E-002   7.6945519499253656E-002   3.9165094703454528E-002   3.3131399179968657E-003   8.1624849714722281E-002  0.11690557683243696        7.9445153514081657E-002  0.15389984382157618       0.11395654789413978       0.15312969136586041        7.8672696056349967E-002   2.1874579092423041E-004   6.3846422737488268E-003   4.5539286224802536E-004   2.4537494943394883E-002   2.2228625446938294E-002   1.0559332348679836E-002  -4.7141202582215208E-004  -2.0455450201165848E-003  -2.3149399920871880E-002  -7.8527393710174677E-003   3.9762001817733060E-002   2.5058595017918302E-002   7.0158013104201430E-003
 /BM/ovlp           7   5.5119444377309534E-002  0.17760930319893120       0.20085862032093615       0.21310852156381088       0.21892740521703757       -2.8142792046807363E-002  0.40303599553624553       0.49372773437644923       0.55463983955560536       0.59760776370859503       -2.6658833493531375E-003 -0.54101803877875909       -2.4477299918364892E-002 -0.56722223594854171       -1.0799297520615738       -1.0597403719626755      -0.54699083249788316       -5.6613095908510047E-002  -1.1417512693067096       -1.6106411606598101       -1.1042573358184029       -2.1090005678515809       -1.5590909317567849       -2.0949176333136847       -1.0901016005458644       -1.9040192368134784E-003  -5.8378038260699455E-002  -7.6942804760326555E-003 -0.36424331410182087      -0.33734866985694412      -0.15609534837137323        4.5259739175784793E-003   1.3985327040499929E-005  0.19442286586049645        5.0728633150975311E-002 -0.60971266829644044      -0.32873938711830064       -6.4750764908357739E-002
 /BM/ovlp           8   3.5144525351986644E-002  0.21451107326020491       0.24671099718814560       0.26456412134211860       0.27378831632312028       -3.3728015701344649E-002  0.49372773437644923       0.61139162251438173       0.69192369022229627       0.74963029078804766       -5.0780818627531588E-003 -0.65353494769136944       -3.4107988530038824E-002 -0.68930207869731674       -1.3025394866903923       -1.2760354509410945      -0.66273222611373939       -7.5613704424529971E-002  -1.3835917847388970       -1.9402637165785563       -1.3348073981064488       -2.5343632014221100       -1.8731819324737558       -2.5160261265590123       -1.3163621834990806       -1.9872209639420491E-003  -6.0457113464116058E-002  -1.0280619242166222E-002 -0.44383469010096599      -0.42868092978634809      -0.19700913775671641        4.1007455379406844E-003  -9.5212495456360458E-003  0.19654838373970696        4.8823610705441500E-002 -0.77427717404632046      -0.34532149880594365       -5.9535778186598165E-002
 /BM/ovlp           9   9.1747249277460696E-003  0.23885279428581896       0.27759821656579220       0.29996296947732048       0.31227789430636221       -3.7434500046984454E-002  0.55463983955560536       0.69192369022229627       0.78757417319641043       0.85728040072621070       -7.0751534374604574E-003 -0.72709560747296109       -4.1312085720754510E-002 -0.76989126884052439       -1.4477632144022152       -1.4170295474130512      -0.73907290665641767       -8.9346136534825149E-002  -1.5422575553730269       -2.1550469656877738       -1.4858902898741349       -2.8106353669808755       -2.0777207077295543       -2.7896266217861125       -1.4647457152421168       -2.0824924024155922E-003  -6.2247690527128863E-002  -1.2122057097495920E-002 -0.49349134066584766      -0.49714543533710609      -0.22699572689508646        4.1579723468857521E-003  -1.6240037852947697E-002  0.20206644546716479        5.0819187488341688E-002 -0.89446961076028186      -0.32823750918459393       -6.0386717834090664E-002
 /BM/ovlp          10  -1.9552521371844378E-002  0.25580419968034107       0.29937765375086656       0.32538032280686480       0.34044404478859747       -4.0053923994471496E-002  0.59760776370859503       0.74963029078804766       0.85728040072621070       0.93701387715815720       -8.6635159609489421E-003 -0.77772417911414493       -4.6703748262711287E-002 -0.82574507639307626       -1.5475727283269123       -1.5139932119148796      -0.79206394673039426       -9.9382425210819747E-002  -1.6516324485292557       -2.3026068173600152       -1.5901258849247029       -3.0001073611664992       -2.2184902556742827       -2.9774391695218583       -1.5672980643020651       -2.1892592073145256E-003  -6.4529410973569234E-002  -1.3415550576538982E-002 -0.52487107690530088      -0.55089053249647968      -0.25048157183513320        4.5488612491703939E-003  -2.0550482480357513E-002  0.21234283977897639        5.4076141231605312E-002 -0.98544291200926193      -0.29598693945728982       -6.8277199018965007E-002
 /BM/ovlp          11  -3.1584853028315374E-002  -6.5022661309934726E-003  -6.2381508762348048E-003  -6.1791956127317249E-003  -6.2200459701768906E-003  -1.9567364603345627E-005  -2.6658833493531375E-003  -5.0780818627531588E-003  -7.0751534374604574E-003  -8.6635159609489421E-003   7.3462729318203909E-003   1.9217398391112361E-002   1.7492396398191445E-002   2.3972319222879790E-002   4.3444109402141606E-002   3.8791884944521016E-002   1.9261103151432479E-002   2.8688462829872030E-002   5.0444381074205680E-002   6.9520683654234361E-002   4.2299687028297228E-002   8.1214131553235092E-002   5.8500329818512498E-002   7.8315057687460055E-002   3.9330868688614018E-002   2.6700815682987981E-004   8.0612865247525756E-003   3.3443621278866908E-004   8.9835511609656360E-003   7.3321587751846717E-003   4.2467612405948235E-003  -4.9848393702874756E-003  -4.0771103141760759E-003  -3.5708190802726525E-002  -1.7436872269522750E-002   1.7854986422757835E-002  -9.5762213970465720E-003  -2.7020816915212839E-002
 /BM/ovlp          12 -0.21196239220750357      -0.26031686880429561      -0.28100605131773193      -0.29030453801665601      -0.29339651152294266        3.8931239007275575E-002 -0.54101803877875909      -0.65353494769136944      -0.72709560747296109      -0.77772417911414493        1.9217398391112361E-002  0.78969875553786217        7.1115604447300029E-002  0.82887421133719386        1.5956574506884778        1.5606474145218954       0.79372306873175802       0.14061013603069128        1.6843297228329561        2.3998281468877138        1.6191171736486467        3.1297737846257405        2.3109089082361862        3.1059983015381079        1.5951436795180740        4.2724457374266611E-003  0.13002360946020503        9.3996144310138086E-003  0.50777640942554814       0.43374242927445417       0.20178981197372917       -2.2775314856275018E-002  -3.7805401205334765E-002 -0.47518773306243312      -0.15942898862721433       0.79696015919543572       0.51848856084054851        5.3663534151669801E-002
 /BM/ovlp          13  -5.6482843700315755E-002  -2.3981147028110605E-002  -2.3249619263879140E-002  -2.2958339702867647E-002  -2.2881770517432898E-002   1.2658115190549712E-003  -2.4477299918364892E-002  -3.4107988530038824E-002  -4.1312085720754510E-002  -4.6703748262711287E-002   1.7492396398191445E-002   7.1115604447300029E-002   4.3280953917634486E-002   8.3327910669305894E-002  0.15533812401360336       0.14335755478278145        7.1222203246792404E-002   7.2677077281210423E-002  0.17454282519675246       0.24474353863942255       0.15316037110106251       0.29688088968469728       0.21576237101528051       0.28923079418709108       0.14535924450485638        9.2847673641961379E-004   2.6101905146230425E-002   1.2378174471074088E-003   3.2830359649488161E-002   2.9057634350880102E-002   1.7693966990592010E-002  -1.1648596439692131E-002  -1.3925497134854456E-002 -0.11054008404140565       -5.7230902811648124E-002   6.2455511470832364E-002  -9.2783940013369537E-003  -6.3380324268947241E-002
 /BM/ovlp          14 -0.22493218666117709      -0.27316463261002610      -0.29706007840361082      -0.30835826726103477      -0.31264258592851490        4.0287880243468965E-002 -0.56722223594854171      -0.68930207869731674      -0.76989126884052439      -0.82574507639307626        2.3972319222879790E-002  0.82887421133719386        8.3327910669305894E-002  0.87381417416683782        1.6753032137590367        1.6357429818755094       0.83406825776305027       0.16126337922521472        1.7739548972411825        2.5195076240024434        1.7008583372455490        3.2786173810198704        2.4197580068479994        3.2518767796818793        1.6738649377151518        4.2841916745400832E-003  0.13117315838143759        1.0418869620770210E-002  0.53625612575881343       0.46380148614297334       0.21546174164574006       -2.4786448361485983E-002  -3.4540598205628492E-002 -0.47996972220284206      -0.15958068167528716       0.85503472645840795       0.51675565892581909        3.2078709419431251E-002
 /BM/ovlp          15 -0.42302653608943253      -0.52657261221482798      -0.56354310235974481      -0.57948667706135382      -0.58413584967408383        7.8067707377426299E-002  -1.0799297520615738       -1.3025394866903923       -1.4477632144022152       -1.5475727283269123        4.3444109402141606E-002   1.5956574506884778       0.15533812401360336        1.6753032137590367        3.2312655059423809        3.1586910895846074        1.6024288260808675       0.30429707649308924        3.4103543353219550        4.8676363485921854        3.2748060067518452        6.3433847359621893        4.6830055429992115        6.2941567518975319        3.2251593751153109        9.3531317760984223E-003  0.28052651448136584        1.8672886538671940E-002   1.0088412400401023       0.85875102659354585       0.40555195602065547       -5.0043095863939457E-002  -9.1756008091230634E-002  -1.0344601346100777      -0.36743100917667060        1.5778118142983870        1.0403372129172057        9.7433143417628099E-002
 /BM/ovlp          16 -0.40683719340000124      -0.51519850091779062      -0.55062908977519953      -0.56582799659220484      -0.57021276209832195        7.6945519499253656E-002  -1.0597403719626755       -1.2760354509410945       -1.4170295474130512       -1.5139932119148796        3.8791884944521016E-002   1.5606474145218954       0.14335755478278145        1.6357429818755094        3.1586910895846074        3.0903495859144163        1.5671456702271058       0.28386716401912437        3.3297664704714975        4.7571035217679594        3.2015964533585191        6.2056718835278843        4.5827133547755636        6.1593244308051283        3.1548833673215313        9.2428236757255283E-003  0.27675593483689059        1.7996413560182500E-002  0.98586407449761815       0.84016456183433874       0.39454646284731903       -4.7715786714333290E-002  -9.1721259867356131E-002  -1.0207532132537156      -0.36180349263372935        1.5404931137555806        1.0247501548900781       0.11635491364704498     
 /BM/ovlp          17 -0.20873824305372413      -0.26175043153891409      -0.28408650122823786      -0.29462965966487786      -0.29864391767074494        3.9165094703454528E-002 -0.54699083249788316      -0.66273222611373939      -0.73907290665641767      -0.79206394673039426        1.9261103151432479E-002  0.79372306873175802        7.1222203246792404E-002  0.83406825776305027        1.6024288260808675        1.5671456702271058       0.79864657716205223       0.14064539177943480        1.6930026749354568        2.4085240143025430        1.6273511507963860        3.1403964112625431        2.3191114236051504        3.1165669665303994        1.6033252030776595        4.1762684739799216E-003  0.12744833962796062        9.7436951589462417E-003  0.51318205040807952       0.44510414539565590       0.20432427775272621       -2.2436201775289211E-002  -3.4554979713959155E-002 -0.46647373014700594      -0.15416640563438255       0.81755820386029221       0.50098769520070729        5.0809174972001436E-002
 /BM/ovlp          18  -7.2809648803280780E-002  -4.7408918274673084E-002  -4.5631683643394538E-002  -4.4732555286543629E-002  -4.4281749784488511E-002   3.3131399179968657E-003  -5.6613095908510047E-002  -7.5613704424529971E-002  -8.9346136534825149E-002  -9.9382425210819747E-002   2.8688462829872030E-002  0.14061013603069128        7.2677077281210423E-002  0.16126337922521472       0.30429707649308924       0.28386716401912437       0.14064539177943480       0.12379464003319379       0.33785870560274489       0.47761645692381194       0.30094519495704120       0.58661350767905418       0.42756294830310537       0.57339455694156527       0.28749600575256196        1.8888722449566589E-003   5.0893066094149564E-002   2.4119009990404838E-003   6.3234087627705216E-002   5.7030288659658701E-002   3.5694153513386007E-002  -1.8957327840292273E-002  -2.7801936336561833E-002 -0.21085556901137542      -0.11523417046210316       0.11819612624419305       -1.8025547709493717E-003 -0.10236889347838807     
 /BM/ovlp          19 -0.44493963260686087      -0.55553436300970560      -0.59860673651463969      -0.61804170924081347      -0.62452619731971026        8.1624849714722281E-002  -1.1417512693067096       -1.3835917847388970       -1.5422575553730269       -1.6516324485292557        5.0444381074205680E-002   1.6843297228329561       0.17454282519675246        1.7739548972411825        3.4103543353219550        3.3297664704714975        1.6930026749354568       0.33785870560274489        3.6079125762918443        5.1359836297524168        3.4583385161832894        6.6829964981916419        4.9319586839391434        6.6283720326821367        3.4032145783790213        9.4196305337180541E-003  0.28419822902338060        2.0654898814149170E-002   1.0729834532811804       0.91993820906149448       0.43391554228494356       -5.2709175252726004E-002  -8.6323926539157880E-002  -1.0460209884900065      -0.36699571259045127        1.6951968626769620        1.0666056424393853        6.5593129304911457E-002
 /BM/ovlp          20 -0.61638754544271279      -0.79272835172740486      -0.84249125586930518      -0.86316985803730972      -0.86840297928731047       0.11690557683243696       -1.6106411606598101       -1.9402637165785563       -2.1550469656877738       -2.3026068173600152        6.9520683654234361E-002   2.3998281468877138       0.24474353863942255        2.5195076240024434        4.8676363485921854        4.7571035217679594        2.4085240143025430       0.47761645692381194        5.1359836297524168        7.3418807488939626        4.9288500407401159        9.5636573041229553        7.0600412845200253        9.4887317386302925        4.8532933240676641        1.5030033709014040E-002  0.44381171140488407        2.7769600434870645E-002   1.4955585709362840        1.2725466340241844       0.60741064515048837       -7.9113222774646164E-002 -0.15613404959856636       -1.6453521385757639      -0.61527850142554663        2.3365322448221022        1.5502069836079282       0.13889793498398628     
 /BM/ovlp          21 -0.42162564105244904      -0.53418604584736329      -0.57508876975086309      -0.59347949553901458      -0.59961053733232461        7.9445153514081657E-002  -1.1042573358184029       -1.3348073981064488       -1.4858902898741349       -1.5901258849247029        4.2299687028297228E-002   1.6191171736486467       0.15316037110106251        1.7008583372455490        3.2748060067518452        3.2015964533585191        1.6273511507963860       0.30094519495704120        3.4583385161832894        4.9288500407401159        3.3218227638818609        6.4238783541046587        4.7428623082360275        6.3742609821261453        3.2717872525316238        9.0746177305965625E-003  0.27427123243846419        1.9416652853443073E-002   1.0327015846556695       0.88632016822173210       0.41434449139697316       -4.8411331978070619E-002  -8.4053358236812481E-002  -1.0087376317686827      -0.34766345963708334        1.6283053679708803        1.0410202099229255        9.8098232409705943E-002
 /BM/ovlp          22 -0.80230122352797650       -1.0339228794987605       -1.0974031399545652       -1.1232363772386975       -1.1292660652582640       0.15389984382157618       -2.1090005678515809       -2.5343632014221100       -2.8106353669808755       -3.0001073611664992        8.1214131553235092E-002   3.1297737846257405       0.29688088968469728        3.2786173810198704        6.3433847359621893        6.2056718835278843        3.1403964112625431       0.58661350767905418        6.6829964981916419        9.5636573041229553        6.4238783541046587        12.474033637041138        9.2113348512255016        12.380575375167950        6.3296916459991053        1.9606493611595249E-002  0.58074944425390163        3.5382873924713065E-002   1.9527353661157170        1.6561879981346115       0.78648826763948143       -9.9868871015971833E-002 -0.20635161365607790       -2.1508769747141563      -0.79145228121553290        3.0348834478770734        2.0647520094598040       0.23357025685018673     
 /BM/ovlp          23 -0.58603774883203918      -0.76374652378179064      -0.81032305128361770      -0.82957004228413211      -0.83441794073482356       0.11395654789413978       -1.5590909317567849       -1.8731819324737558       -2.0777207077295543       -2.2184902556742827        5.8500329818512498E-002   2.3109089082361862       0.21576237101528051        2.4197580068479994        4.6830055429992115        4.5827133547755636        2.3191114236051504       0.42756294830310537        4.9319586839391434        7.0600412845200253        4.7428623082360275        9.2113348512255016        6.8032471836730712        9.1434578714027452        4.6744710645463101        1.4647553280596276E-002  0.43194325124525329        2.6078804222578467E-002   1.4398477143631583        1.2276430579664992       0.57989318833087811       -7.3714736643430590E-002 -0.15376666621846113       -1.6016409993576755      -0.59515883294914818        2.2469652393885733        1.5082451196478404       0.18386384887565521     
 /BM/ovlp          24 -0.79524879908880663       -1.0262687087715960       -1.0887245563954480       -1.1141664515955938       -1.1201623956536650       0.15312969136586041       -2.0949176333136847       -2.5160261265590123       -2.7896266217861125       -2.9774391695218583        7.8315057687460055E-002   3.1059983015381079       0.28923079418709108        3.2518767796818793        6.2941567518975319        6.1593244308051283        3.1165669665303994       0.57339455694156527        6.6283720326821367        9.4887317386302925        6.3742609821261453        12.380575375167950        9.1434578714027452        12.289171249965754        6.2821650303898764        1.9569717512397777E-002  0.57882350526480941        3.4930435903757218E-002   1.9372571726338501        1.6448542957924568       0.77849032330431789       -9.8763807917123525E-002 -0.20627388114725909       -2.1445275261881083      -0.79072262783842473        3.0121321275322748        2.0482957332544949       0.24597567022086242     
 /BM/ovlp          25 -0.41452344708162059      -0.52647274664034427      -0.56632816300324862      -0.58431533980167671      -0.59040608027652297        7.8672696056349967E-002  -1.0901016005458644       -1.3163621834990806       -1.4647457152421168       -1.5672980643020651        3.9330868688614018E-002   1.5951436795180740       0.14535924450485638        1.6738649377151518        3.2251593751153109        3.1548833673215313        1.6033252030776595       0.28749600575256196        3.4032145783790213        4.8532933240676641        3.2717872525316238        6.3296916459991053        4.6744710645463101        6.2821650303898764        3.2238828604817718        9.0401988157210733E-003  0.27238658300062274        1.8964864879301824E-002   1.0170838984175030       0.87484117682952300       0.40619621190832966       -4.7288364138737972E-002  -8.4020481587656071E-002  -1.0025732808247796      -0.34716372708744525        1.6053260625926953        1.0243615066729341       0.11031010175741640     
 /BM/ovlp          26   2.7022547625387468E-003  -1.5087367033199059E-003  -9.5780924301930304E-004  -7.1732622098420740E-004  -6.2174950184193389E-004   2.1874579092423041E-004  -1.9040192368134784E-003  -1.9872209639420491E-003  -2.0824924024155922E-003  -2.1892592073145256E-003   2.6700815682987981E-004   4.2724457374266611E-003   9.2847673641961379E-004   4.2841916745400832E-003   9.3531317760984223E-003   9.2428236757255283E-003   4.1762684739799216E-003   1.8888722449566589E-003   9.4196305337180541E-003   1.5030033709014040E-002   9.0746177305965625E-003   1.9606493611595249E-002   1.4647553280596276E-002   1.9569717512397777E-002   9.0401988157210733E-003   2.2240464002428766E-004   3.2022368225798855E-003   3.3533403687833411E-005   3.8258321041943245E-004   1.0333153629983694E-003   5.4301830339825588E-004  -3.0651872580181728E-004  -1.3323636128960057E-003  -1.0124196783398676E-002  -1.6880717922200653E-002   1.6464263324223296E-003  -1.1753417431334845E-003   9.8762644170492428E-004
 /BM/ovlp          27   2.6789182021430318E-002  -4.6317566913179342E-002  -3.2519907644469048E-002  -2.4706188231380155E-002  -2.0817891727766324E-002   6.3846422737488268E-003  -5.8378038260699455E-002  -6.0457113464116058E-002  -6.2247690527128863E-002  -6.4529410973569234E-002   8.0612865247525756E-003  0.13002360946020503        2.6101905146230425E-002  0.13117315838143759       0.28052651448136584       0.27675593483689059       0.12744833962796062        5.0893066094149564E-002  0.28419822902338060       0.44381171140488407       0.27427123243846419       0.58074944425390163       0.43194325124525329       0.57882350526480941       0.27238658300062274        3.2022368225798855E-003   9.4492595469680182E-002   5.2927309420464877E-004   1.6567741801140588E-002   2.6113175286511595E-002   2.4236190750606647E-002  -1.5760353166933579E-002  -7.6394126296213338E-002 -0.40769904576748489      -0.15163839367420134        3.7562145966436575E-002   9.6502421682097671E-003   7.0886205341507431E-003
 /BM/ovlp          28   6.0025210191860766E-003  -3.1154688277532927E-003  -3.9294971562369718E-003  -4.4440065961954267E-003  -4.7268955121501643E-003   4.5539286224802536E-004  -7.6942804760326555E-003  -1.0280619242166222E-002  -1.2122057097495920E-002  -1.3415550576538982E-002   3.3443621278866908E-004   9.3996144310138086E-003   1.2378174471074088E-003   1.0418869620770210E-002   1.8672886538671940E-002   1.7996413560182500E-002   9.7436951589462417E-003   2.4119009990404838E-003   2.0654898814149170E-002   2.7769600434870645E-002   1.9416652853443073E-002   3.5382873924713065E-002   2.6078804222578467E-002   3.4930435903757218E-002   1.8964864879301824E-002   3.3533403687833411E-005   5.2927309420464877E-004   4.8928397556146448E-004   6.5272559052962252E-003   7.3394060607717033E-003   2.4526230077471763E-003  -1.6799746172577887E-004   4.6505398063080832E-004  -2.6244887170462266E-003  -3.1128368083337497E-003   1.5267958662707044E-002  -4.7792611027784714E-003  -1.4189822824474133E-002
 /BM/ovlp          29 -0.39420256109380813      -0.16442437092405693      -0.19299498245862390      -0.20680537130113308      -0.21192919812000666        2.4537494943394883E-002 -0.36424331410182087      -0.44383469010096599      -0.49349134066584766      -0.52487107690530088        8.9835511609656360E-003  0.50777640942554814        3.2830359649488161E-002  0.53625612575881343        1.0088412400401023       0.98586407449761815       0.51318205040807952        6.3234087627705216E-002   1.0729834532811804        1.4955585709362840        1.0327015846556695        1.9527353661157170        1.4398477143631583        1.9372571726338501        1.0170838984175030        3.8258321041943245E-004   1.6567741801140588E-002   6.5272559052962252E-003  0.44106410354203196       0.32510095488957674        6.6225343602374392E-002  -6.9094874426314972E-003   4.9716299149192053E-002  -1.0595506423839302E-003   1.5206025985721634E-002  0.63834357815765441       0.46112620284232364        2.5238953501714477E-002
 /BM/ovlp          30  -5.2688392992979963E-002 -0.14303299563766536      -0.17273054526881548      -0.19255443609613376      -0.20596517036890716        2.2228625446938294E-002 -0.33734866985694412      -0.42868092978634809      -0.49714543533710609      -0.55089053249647968        7.3321587751846717E-003  0.43374242927445417        2.9057634350880102E-002  0.46380148614297334       0.85875102659354585       0.84016456183433874       0.44510414539565590        5.7030288659658701E-002  0.91993820906149448        1.2725466340241844       0.88632016822173210        1.6561879981346115        1.2276430579664992        1.6448542957924568       0.87484117682952300        1.0333153629983694E-003   2.6113175286511595E-002   7.3394060607717033E-003  0.32510095488957674       0.38332562643562140       0.10129917034515863       -7.2910016187624511E-003   4.0843711948430284E-002  -7.2457882729069528E-002  -3.4826361040159716E-002  0.66959013342450369        8.3502324366380415E-002   9.9415934760643379E-002
 /BM/ovlp          31  0.13704565617402450       -6.6302188280362745E-002  -7.5616986764080185E-002  -8.1055786472355840E-002  -8.4363425096324496E-002   1.0559332348679836E-002 -0.15609534837137323      -0.19700913775671641      -0.22699572689508646      -0.25048157183513320        4.2467612405948235E-003  0.20178981197372917        1.7693966990592010E-002  0.21546174164574006       0.40555195602065547       0.39454646284731903       0.20432427775272621        3.5694153513386007E-002  0.43391554228494356       0.60741064515048837       0.41434449139697316       0.78648826763948143       0.57989318833087811       0.77849032330431789       0.40619621190832966        5.4301830339825588E-004   2.4236190750606647E-002   2.4526230077471763E-003   6.6225343602374392E-002  0.10129917034515863       0.17031117063023399       -1.9143612661763768E-003  -1.8615823954875155E-002 -0.10098466855818765        9.3467734940035786E-003  0.15329310249943348        5.1410941521035210E-002   3.3140694670378820E-002
 /BM/ovlp          32   3.7617884510980755E-002   7.8367336517800867E-003   6.1599432611259917E-003   5.2725974854381553E-003   4.9848807115722114E-003  -4.7141202582215208E-004   4.5259739175784793E-003   4.1007455379406844E-003   4.1579723468857521E-003   4.5488612491703939E-003  -4.9848393702874756E-003  -2.2775314856275018E-002  -1.1648596439692131E-002  -2.4786448361485983E-002  -5.0043095863939457E-002  -4.7715786714333290E-002  -2.2436201775289211E-002  -1.8957327840292273E-002  -5.2709175252726004E-002  -7.9113222774646164E-002  -4.8411331978070619E-002  -9.9868871015971833E-002  -7.3714736643430590E-002  -9.8763807917123525E-002  -4.7288364138737972E-002  -3.0651872580181728E-004  -1.5760353166933579E-002  -1.6799746172577887E-004  -6.9094874426314972E-003  -7.2910016187624511E-003  -1.9143612661763768E-003   1.6962912637744059E-002   9.5709649243385410E-003   7.2479643798756557E-002   8.2270382684463475E-003  -1.5344003788190279E-002   1.0809806380182333E-002   2.3669615289052457E-002
 /BM/ovlp          33  -7.9790296009032979E-002   1.6191952670934595E-002   1.3352052478130538E-003  -8.2532161208573873E-003  -1.3821532016301763E-002  -2.0455450201165848E-003   1.3985327040499929E-005  -9.5212495456360458E-003  -1.6240037852947697E-002  -2.0550482480357513E-002  -4.0771103141760759E-003  -3.7805401205334765E-002  -1.3925497134854456E-002  -3.4540598205628492E-002  -9.1756008091230634E-002  -9.1721259867356131E-002  -3.4554979713959155E-002  -2.7801936336561833E-002  -8.6323926539157880E-002 -0.15613404959856636       -8.4053358236812481E-002 -0.20635161365607790      -0.15376666621846113      -0.20627388114725909       -8.4020481587656071E-002  -1.3323636128960057E-003  -7.6394126296213338E-002   4.6505398063080832E-004   4.9716299149192053E-002   4.0843711948430284E-002  -1.8615823954875155E-002   9.5709649243385410E-003  0.11603244259796240       0.38443864626246560       -2.8897801622104601E-002   8.3445718955623652E-002   2.2468860439021507E-002   1.0947365008776157E-002
 /BM/ovlp          34 -0.10319865678939415       0.17376610132910453       0.11560387817351880        8.2174581623534859E-002   6.6210852376570983E-002  -2.3149399920871880E-002  0.19442286586049645       0.19654838373970696       0.20206644546716479       0.21234283977897639       -3.5708190802726525E-002 -0.47518773306243312      -0.11054008404140565      -0.47996972220284206       -1.0344601346100777       -1.0207532132537156      -0.46647373014700594      -0.21085556901137542       -1.0460209884900065       -1.6453521385757639       -1.0087376317686827       -2.1508769747141563       -1.6016409993576755       -2.1445275261881083       -1.0025732808247796       -1.0124196783398676E-002 -0.40769904576748489       -2.6244887170462266E-003  -1.0595506423839302E-003  -7.2457882729069528E-002 -0.10098466855818765        7.2479643798756557E-002  0.38443864626246560        1.9305749562074179       0.26556623922908507       -9.4393705159096120E-002  0.13263458765946706        8.7699280548869191E-002
 /BM/ovlp          35 -0.17560072553267528        5.6476098773644234E-002   2.3278551927173319E-002   1.2281660409638073E-002   9.6066154725999997E-003  -7.8527393710174677E-003   5.0728633150975311E-002   4.8823610705441500E-002   5.0819187488341688E-002   5.4076141231605312E-002  -1.7436872269522750E-002 -0.15942898862721433       -5.7230902811648124E-002 -0.15958068167528716      -0.36743100917667060      -0.36180349263372935      -0.15416640563438255      -0.11523417046210316      -0.36699571259045127      -0.61527850142554663      -0.34766345963708334      -0.79145228121553290      -0.59515883294914818      -0.79072262783842473      -0.34716372708744525       -1.6880717922200653E-002 -0.15163839367420134       -3.1128368083337497E-003   1.5206025985721634E-002  -3.4826361040159716E-002   9.3467734940035786E-003   8.2270382684463475E-003  -2.8897801622104601E-002  0.26556623922908507        1.5242653472432188       -6.9257798718902835E-002  0.20041262812835775        4.7707207708180988E-002
 /BM/ovlp          36 -0.28031107875241290      -0.26207900366862474      -0.31875422679916832      -0.35514450037641032      -0.37807691052239856        3.9762001817733060E-002 -0.60971266829644044      -0.77427717404632046      -0.89446961076028186      -0.98544291200926193        1.7854986422757835E-002  0.79696015919543572        6.2455511470832364E-002  0.85503472645840795        1.5778118142983870        1.5404931137555806       0.81755820386029221       0.11819612624419305        1.6951968626769620        2.3365322448221022        1.6283053679708803        3.0348834478770734        2.2469652393885733        3.0121321275322748        1.6053260625926953        1.6464263324223296E-003   3.7562145966436575E-002   1.5267958662707044E-002  0.63834357815765441       0.66959013342450369       0.15329310249943348       -1.5344003788190279E-002   8.3445718955623652E-002  -9.4393705159096120E-002  -6.9257798718902835E-002   1.2343070237749207       0.17810253372165091        1.8272949804257629E-002
 /BM/ovlp          37 -0.53139686918457851      -0.16117126407761417      -0.16406748533997728      -0.14964751208018146      -0.12868996044915049        2.5058595017918302E-002 -0.32873938711830064      -0.34532149880594365      -0.32823750918459393      -0.29598693945728982       -9.5762213970465720E-003  0.51848856084054851       -9.2783940013369537E-003  0.51675565892581909        1.0403372129172057        1.0247501548900781       0.50098769520070729       -1.8025547709493717E-003   1.0666056424393853        1.5502069836079282        1.0410202099229255        2.0647520094598040        1.5082451196478404        2.0482957332544949        1.0243615066729341       -1.1753417431334845E-003   9.6502421682097671E-003  -4.7792611027784714E-003  0.46112620284232364        8.3502324366380415E-002   5.1410941521035210E-002   1.0809806380182333E-002   2.2468860439021507E-002  0.13263458765946706       0.20041262812835775       0.17810253372165091        1.4861240900527999       0.34876003487206353     
 /BM/ovlp          38  0.23920496314741940       -1.7081715085883742E-002  -1.1029446773254570E-002  -9.2841047224596629E-003  -1.1356566458279538E-002   7.0158013104201430E-003  -6.4750764908357739E-002  -5.9535778186598165E-002  -6.0386717834090664E-002  -6.8277199018965007E-002  -2.7020816915212839E-002   5.3663534151669801E-002  -6.3380324268947241E-002   3.2078709419431251E-002   9.7433143417628099E-002  0.11635491364704498        5.0809174972001436E-002 -0.10236889347838807        6.5593129304911457E-002  0.13889793498398628        9.8098232409705943E-002  0.23357025685018673       0.18386384887565521       0.24597567022086242       0.11031010175741640        9.8762644170492428E-004   7.0886205341507431E-003  -1.4189822824474133E-002   2.5238953501714477E-002   9.9415934760643379E-002   3.3140694670378820E-002   2.3669615289052457E-002   1.0947365008776157E-002   8.7699280548869191E-002   4.7707207708180988E-002   1.8272949804257629E-002  0.34876003487206353        1.1638485562934466     

Eigenvalues of overlap matrix of current wave function and its first-order derivatives:
overlap eigenvalue #    1:  1.00194258E-08
overlap eigenvalue #    2:  1.02630603E-08
overlap eigenvalue #    3:  1.05607883E-08
overlap eigenvalue #    4:  1.10732884E-08
overlap eigenvalue #    5:  1.31373440E-08
overlap eigenvalue #    6:  1.55839257E-08
overlap eigenvalue #    7:  2.04226705E-08
overlap eigenvalue #    8:  4.38105638E-08
overlap eigenvalue #    9:  6.50472765E-08
overlap eigenvalue #   10:  7.54174570E-08
overlap eigenvalue #   11:  1.75848918E-07
overlap eigenvalue #   12:  5.47362439E-07
overlap eigenvalue #   13:  8.34297760E-07
overlap eigenvalue #   14:  1.05856553E-06
overlap eigenvalue #   15:  1.28112622E-06
overlap eigenvalue #   16:  3.79233881E-06
overlap eigenvalue #   17:  5.50460912E-06
overlap eigenvalue #   18:  9.09074201E-06
overlap eigenvalue #   19:  1.89399257E-05
overlap eigenvalue #   20:  3.43459623E-05
overlap eigenvalue #   21:  6.29199010E-05
overlap eigenvalue #   22:  9.56501095E-05
overlap eigenvalue #   23:  1.24743848E-04
overlap eigenvalue #   24:  3.33345741E-04
overlap eigenvalue #   25:  5.89964892E-04
overlap eigenvalue #   26:  3.98616108E-03
overlap eigenvalue #   27:  7.22831178E-03
overlap eigenvalue #   28:  1.47261036E-02
overlap eigenvalue #   29:  1.77730242E-02
overlap eigenvalue #   30:  2.80326208E-02
overlap eigenvalue #   31:  1.41318249E-01
overlap eigenvalue #   32:  1.84866791E-01
overlap eigenvalue #   33:  8.12036571E-01
overlap eigenvalue #   34:  1.00000001E+00
overlap eigenvalue #   35:  1.38593571E+00
overlap eigenvalue #   36:  1.60153157E+00
overlap eigenvalue #   37:  2.03882679E+00
overlap eigenvalue #   38:  2.66468804E+00
overlap eigenvalue #   39:  6.26400781E+01
 /BM/amat           1  -25.999876445178206       -1.0887844748503643       -1.1585081449832875       -1.1609427390885914       -1.1082854649468268       -3.6634894202317414E-002 -0.39649002405113498      -0.35672598436032443       -6.7347832889067413E-002  0.33510975051888758       0.95565235267178994        3.3576290506110862        2.2146686832653022        3.9420585103908530        7.2345807055560112        6.6316463439953282        3.3450177029709001        3.6495836617577595        8.2259467783337819        11.278380306230710        7.1521944539374305        13.626663508600194        9.8098878048170235        13.233067474809280        6.7649015216281398        1.8350551610535143E-002   7.0727303157913479E-002 -0.10018346336658757        5.6272702027023005       0.70520053990546872       -1.9467461177168852      -0.33276452898238817       0.80374487786017057       0.77881603396543764       -2.3266527983257266        3.9489129513390200        7.4614668818133492       -2.6649320878413016     
 /BM/amat           2 -0.95287814249602421       -1.0907334298177844       -1.2421663418360942       -1.3190951343806456       -1.3536738438541029       0.17166420458557929       -2.4177819656361854       -2.9630998806057920       -3.3326485994020563       -3.5880403430338674        6.0528211046890198E-002   3.3276371364858379       0.24035293936071628        3.5201149600670503        6.6612506539772731        6.5073983489429850        3.3656563005476894       0.47346145370494075        7.0811295413402942        9.9393221861112480        6.8087637993436809        12.955935339284846        9.5598797914997320        12.848326168767812        6.6999872168191628        8.8079754659259324E-003  0.35555191087968629        4.6011482383945028E-002   2.2970795746710664        2.0438819944776054       0.92255649774282733       -8.4251672501269359E-002 -0.21037985281515925       -1.2892398945575088       0.18365185393175665        3.8101019446666649        1.7000584372626855        1.1919593607108485     
 /BM/amat           3  -1.1471000682692372       -1.2260183787866634       -1.4350454844453262       -1.5464347708943496       -1.6023709198309579       0.19292329884674542       -2.7868364198886955       -3.4513993158093843       -3.9091715821415169       -4.2297967881726812        6.8726955918860438E-002   3.7497794379442899       0.26651256146387831        3.9876729078157402        7.4737339151725477        7.2939181729697751        3.8073936681061906       0.51720420613966045        7.9807381488352416        11.111505390436269        7.6660351710634282        14.465477262620475        10.673030153842362        14.341263913235604        7.5408864791108448       -1.7067368222438706E-003  0.26520420296716640        5.6738479558634695E-002   2.7075826630135431        2.4727395461258799        1.0562846793755933       -6.9142788696402030E-002  -3.4407348675994154E-002 -0.99507836571716357       0.78330191438159602        4.6214338249783484        1.7427936853233339        1.2169747970205460     
 /BM/amat           4  -1.1744117651851986       -1.2912335393629535       -1.5389087695588479       -1.6757658731780560       -1.7496232404447465       0.20301269056042337       -2.9826253216384453       -3.7251680114277939       -4.2459827472758782       -4.6170922654108155        7.5846739738526764E-002   3.9507692593258188       0.28526236682820538        4.2193166431606013        7.8530533977675052        7.6584236841943127        4.0243259569881120       0.54549819093655216        8.4139462745466709        11.649335562410252        8.0754201308706737        15.147972495436857        11.178480513918839        15.015264273796234        7.9419401467454973       -8.4190883676324990E-003  0.20347187274781131        6.3886468328304719E-002   2.9106959792375258        2.7606263359228649        1.1362468358029147       -6.2124044781723065E-002   7.5972113703822064E-002 -0.81653022709164436        1.1646311158221145        5.1407245236101158        1.5819797548602255        1.1637562310610079     
 /BM/amat           5  -1.1191660306339370       -1.3181526337622564       -1.5898935039073754       -1.7450521372894501       -1.8336460601687212       0.20693245784076472       -3.0762107707939705       -3.8687328128026195       -4.4345139205844308       -4.8448649493143323        8.2624825291140958E-002   4.0304291667254271       0.30113600634404702        4.3193057380041813        7.9986208343571832        7.7960085564947308        4.1164059566570499       0.56854196010392477        8.5914623042944473        11.850160212684841        8.2399534453557770        15.392593488513093        11.362684866090850        15.256293589265365        8.1029999349167507       -1.1929446737760196E-002  0.16877525114892489        6.8030286743530621E-002   2.9903970655316021        2.9564137745541474        1.1864539850674041       -6.1154153860351224E-002  0.13819186708186848      -0.72782585419995627        1.3548246259015251        5.4683225235426676        1.3328377521633321        1.1266314957299528     
 /BM/amat           6   6.5082956914211398E-002  0.16391543644264867       0.18796491562327045       0.20029338065126731       0.20610567922259390       -2.5948616807834439E-002  0.37185463489796344       0.45911121266541327       0.51811177265113295       0.55912928139207974       -7.1016122891031626E-003 -0.49995124329298690       -3.3343698065737976E-002 -0.52966703762288658      -0.99958942645332238      -0.97634719532759040      -0.50630948171033585       -6.8797099264606087E-002  -1.0651082732568067       -1.4910597689389218       -1.0235334582262277       -1.9422925719080875       -1.4327281461648966       -1.9254484195733765       -1.0064771297349799       -6.8643792540907383E-004  -4.8949866757336633E-002  -7.0676126050874498E-003 -0.34084722268411144      -0.31718814346132601      -0.14748336961217409        6.0228966183575733E-003   2.6890099606270507E-002  0.18506195881104584       -7.9897583354046620E-002 -0.58017342009984563      -0.24627200269548677      -0.25330887995295093     
 /BM/amat           7 -0.77856711587285599       -2.3568021128744885       -2.7691593651038042       -2.9867392417942922       -3.0956280104965712       0.38535125801710174       -5.5289527067968747       -6.8474496995586263       -7.7587469157590716       -8.4001044336943984        5.2148465230632809E-002   7.2157039898464426       0.34033166084860389        7.6404043829739345        14.334824728417074        14.023806297684867        7.3282312620489654       0.72673555431671988        15.272656787373704        21.269510411721871        14.727612592270715        27.775978837252499        20.504627666561582        27.555246339843304        14.504699626628497       -8.9451617386950844E-003  0.46466802770173943       0.11404055691996919        5.0842258339595521        4.8242636602112219        2.1811140847170947       -6.2690999342484816E-002  -3.4146718933856679E-002  -1.8020540510400314        2.0342794046199106        8.8638868771950925        3.2891002242780640        3.1884540808998825     
 /BM/amat           8 -0.47910828048663334       -2.8797490425154582       -3.4246868771627246       -3.7208489285033957       -3.8779869567553664       0.47646142179387135       -6.8625417670014341       -8.5405536264557185       -9.7235631053101521       -10.571145806371169        4.7196688827328498E-002   8.8165072157455420       0.36900421722568677        9.3480359840220899        17.468675545612339        17.096316622800888        8.9744833422048860       0.79503315615172376        18.629570434316722        25.863658036121365        17.981728831205952        33.788739231773540        24.956860430871103        33.528841578521785        17.719530886246591       -1.9620153499152027E-002  0.47472099440013649       0.15038742915561135        6.2179009923390964        6.1394050437496039        2.7571190431410972       -6.3635043531129143E-002   7.2905999679271805E-002  -1.9526253509289955        2.8976960687468587        11.229344982354823        3.4605729396348819        3.3684935211427884     
 /BM/amat           9  -9.0319788104537224E-002  -3.2252257976561634       -3.8688102978338867       -4.2286200051576364       -4.4293740864063835       0.53592655072418305       -7.7605680687927903       -9.7066039841850493       -11.102210699577263       -12.118827904103918        5.2157822207222149E-002   9.8663527943737641       0.40732489966824598        10.481791458239831        19.520119280670169        19.103082429016663        10.063586307006032       0.87027823296348394        20.845057627536264        28.867407090051714        20.121697545667487        37.700615376771829        27.858197389113730        37.413822192058895        19.832607189248165       -2.7209958204736484E-002  0.47832205962577312       0.17653780933679489        6.9328013598581171        7.1281079268230982        3.1826674920579037       -6.9410829449453149E-002  0.14730603396795106       -2.0798386148198231        3.4907634254372724        12.955526543491660        3.2356761910080873        3.4468449964893004     
 /BM/amat          10  0.33946559642626406       -3.4643335809933178       -4.1824429662266303       -4.5945161077958279       -4.8344600428962163       0.57673182170187121       -8.3918511620128839       -10.544576807779553       -12.112302797866482       -13.271837516459927        6.1924933457297437E-002   10.586015758895007       0.44871482738607898        11.268427676644901        20.925396171828481        20.475053649380264        10.816958532209947       0.94551706188066231        22.373822763605517        30.924878692681624        21.593635078279064        40.366594829604757        29.840057272845982        40.061034462671209        21.285825365132020       -3.2173632036713418E-002  0.48624185703601874       0.19522506613299651        7.3892305338407489        7.9057383743763765        3.5182517225711596       -7.8573355647569287E-002  0.19404317514512393       -2.2118801578660214        3.8839510975501135        14.264086219722035        2.8112887715082060        3.5599561517382483     
 /BM/amat          11  0.47099392275672125        7.1201905421022796E-002   7.8807323595569065E-002   8.4109664958732483E-002   8.7725222135942249E-002  -8.0580656433327999E-003  0.10171354182401662       0.12192193490536039       0.14047670335209969       0.15592031949277552       -3.1719383793615739E-002 -0.21730936299802778       -6.9521543791058099E-002 -0.23922500651625975      -0.44701919922216571      -0.42759130833239223      -0.21969959891543728      -0.10541703273244063      -0.47985600618249813      -0.67313771790651100      -0.44834249756006950      -0.85718312121403351      -0.63235359049289952      -0.84781591155694214      -0.43885075768454740       -2.9004690389301119E-003  -4.9697911083435353E-002  -3.3441804568775328E-003 -0.14300267165989269      -0.11043545103415574       -5.9416962000780893E-002   4.9064045629202170E-002   5.2206531109762921E-002  0.18558698075306662       0.21554673393452850      -0.24347386279139643       -2.4704960575431978E-002  0.57331883896592073     
 /BM/amat          12   3.0955021639061409        3.3207385405380117        3.7937523432605764        4.0294200751927978        4.1317719326965134      -0.52423053591867630        7.3906562523914143        9.0485098598163756        10.163280380025956        10.925939284603910      -0.17726448060798461       -10.156067039108031      -0.71062607305499448       -10.741324989076674       -20.314664419460211       -19.844056558381737       -10.268502215618135       -1.4005877430485247       -21.597468940295826       -30.288322950253466       -20.768716028176666       -39.490433286440066       -29.130320110479193       -39.159403455986691       -20.433660100773501       -1.8292942224127251E-002 -0.99646412258768358      -0.13932957205354909       -7.0895102347171362       -6.1985924424631449       -2.8083906657741053       0.24554579021210513       0.50599636805039916        3.6123884846995624       -1.0462155948659220       -11.595486214050741       -5.4955518879176726       -3.7045202893487001     
 /BM/amat          13  0.85088726691647598       0.26672894045007334       0.29550275121469000       0.31190559890014963       0.32082296909970898       -3.9802839794051983E-002  0.49325898544669533       0.58702531611144926       0.66385255646994068       0.72368954749449443       -5.5034371758728773E-002 -0.81399648996269969      -0.12723141427826901      -0.86641809218942833       -1.6451881037687632       -1.6003672190838754      -0.82125351252320811      -0.19489594450299741       -1.7360581249280425       -2.4563124382818158       -1.6638774012471851       -3.1938640853936726       -2.3625376335641448       -3.1719676519231563       -1.6414773375141809       -8.2415559798298954E-003 -0.16166624301055099       -1.4559724154004729E-002 -0.49820486432184602      -0.42812869378410601      -0.24513207651574104       0.11471895956080164       0.17529047090638394       0.61429739362729052       0.48990790768943970      -0.86856400798826394      -0.20381757249859156        1.2865096040354480     
 /BM/amat          14   3.2876297393604972        3.5038169922311830        4.0236310354513796        4.2870698000174254        4.4062433051138683      -0.55417500563724820        7.8274523395797928        9.6012781435357422        10.804754364768691        11.634559853779818      -0.19047578374466281       -10.718916355541563      -0.74893700192736423       -11.346727498812809       -21.421269184157737       -20.923583733627311       -10.846864435502521       -1.4629258513845724       -22.786162392403281       -31.911501567089768       -21.913291326259099       -41.603985928590468       -30.694057925144559       -41.257626745168722       -21.562961107493173       -1.4943179104123239E-002  -1.0014481867933482      -0.15243931313753226       -7.5051998746384463       -6.6343715183995480       -3.0006571947722813       0.26625403921410751       0.46997193483647814        3.6642570521972448       -1.2565227471543725       -12.419023463701741       -5.5641254348090730       -3.3726151187192333     
 /BM/amat          15   6.1906055250923355        6.6599907308167374        7.5763696746200004        8.0284312170610370        8.2199344508783803       -1.0527024172647061        14.756943448125345        18.031156784706425        20.231074879837308        21.734713940119047      -0.35678388423244733       -20.362171901177554       -1.4252695508151123       -21.515177233961868       -40.753475173943201       -39.818115518933041       -20.575091060742828       -2.8064383306866811       -43.287861500064061       -60.787227320024414       -41.640225590750191       -79.281574125302328       -58.486521904930328       -78.624735927679950       -40.974870274002811       -4.6940905532293131E-002  -2.1223552342465188      -0.27703360589379483       -14.086008348883793       -12.272631131971940       -5.6406732559088901       0.53456687565186245        1.2076859794980028        7.6879542967849765       -1.5507941874233495       -22.957666895342182       -11.093656920693400       -7.1125969316120319     
 /BM/amat          16   5.9471532275793342        6.5083349316032759        7.3977468537338940        7.8369254027513788        8.0232104595990830       -1.0272752736029904        14.415509976081239        17.618886815744453        19.769052418835091        21.238398220278892      -0.34930096512024811       -19.891439185942453       -1.4023272673364513       -21.019810782673396       -39.821015767107085       -38.905056576042512       -20.099261079388413       -2.7689788783655986       -42.302854224354462       -59.412501825149974       -40.686104120927766       -77.477321348256936       -57.154136771464422       -76.832035880993530       -40.032492180492156       -4.8091035109625091E-002  -2.0972934552773936      -0.26852417661358247       -13.755079141351775       -12.003766492039734       -5.4881938730114213       0.51165128917867242        1.2052006739574126        7.5893583744133153       -1.4503007634380651       -22.427243602550096       -10.835006671150834       -7.3614270363997960     
 /BM/amat          17   3.0440720586628154        3.3514382991725933        3.8440770890337550        4.0945185125025620        4.2083969761607420      -0.52863932946966896        7.4843245098761884        9.1868679030080429        10.340259095653188        11.135483727412332      -0.18289097142700714       -10.245880790685099      -0.72574246476014004       -10.848887836642344       -20.484270721503773       -20.006122128819221       -10.368682799491605       -1.4250563499023847       -21.796278819972571       -30.530114565158282       -20.954525471185253       -39.791112843946081       -29.355324447790782       -39.456382529038805       -20.615994244366878       -1.6150005971086734E-002 -0.97663612094193630      -0.14394370175184201       -7.1729523806429025       -6.3638589930075220       -2.8462861261817571       0.24306275020028922       0.46816140176452992        3.5662667593149582       -1.1517681355417846       -11.886019310506066       -5.3046571307557766       -3.6160418883125018     
 /BM/amat          18   1.1093782340779133       0.52355726982110551       0.57765282130522622       0.60585829971180505       0.61914063717695256       -8.4794220618523750E-002   1.0339431581348326        1.2229902040727341        1.3743570562372329        1.4901659109592211       -6.6494918190088514E-002  -1.5973568177790174      -0.16109956862696351       -1.6783373907246859       -3.2103325886877165       -3.1420994412763101       -1.6094039503328412      -0.24678630304985821       -3.3624280152432506       -4.7791305900743346       -3.2538728557095413       -6.2634852994366952       -4.6383309601641427       -6.2306576090183654       -3.2199724590321588       -1.5428984827065062E-002 -0.31281679705946469       -2.9793789614120186E-002 -0.94518935833105666      -0.83525593000944043      -0.49263273698201937       0.18670883100236499       0.34851100866311446        1.1970008227897073       0.80713404440378733       -1.6535446193173937      -0.47837262348346954        2.0422116717336944     
 /BM/amat          19   6.5160109320540061        7.0659748972051410        8.0730864101163942        8.5750682580173709        8.7935037870752311       -1.1211225830435319        15.728797392997992        19.239965344625425        21.612837855911685        23.241935093088824      -0.36825911606766187       -21.613201833492333       -1.4746518459108422       -22.844257840844765       -43.216284902991944       -42.229028687809503       -21.852244577821892       -2.8952837733787460       -45.913430730071887       -64.407049565443288       -44.180611021665356       -84.019497362904843       -61.989448353897046       -83.330880676740463       -43.483404036580332       -4.0775495532806794E-002  -2.1484699185375251      -0.30315452756351252       -15.006225847535198       -13.155122832502350       -6.0368077956568111       0.56222115936074391        1.1476283068024049        7.8274501445915252       -2.0181087050704001       -24.633089420328410       -11.511383691785927       -6.7949920495770435     
 /BM/amat          20   9.0338011369486182        9.9531915595884577        11.287698277801034        11.941262173296506        12.213173259477459       -1.5770141977002656        22.001206791672647        26.842981675363770        30.097212684211549        32.321894451580079      -0.52337905473854762       -30.420384628070654       -2.1046441500978474       -32.115171185605490       -60.903968505303055       -59.522198998348465       -30.725995135950598       -4.1496083720385339       -64.640917672852780       -90.864830249837553       -62.206413646240797       -118.55463910954290       -87.466962001021841       -117.58540299773752       -61.223894309624839       -8.1986866185111751E-002  -3.3172819402471032      -0.41284067996342844       -20.881386689890554       -18.186591651966079       -8.4442761906275994       0.84092150511985897        2.0372045383774813        12.031743039118089       -1.7236329991758481       -33.998985542653649       -16.585337892187546       -10.302076548275480     
 /BM/amat          21   6.1623746276981493        6.7913137723848305        7.7524783693581698        8.2325470009047130        8.4421283688023490       -1.0741713981190568        15.108934130745462        18.493464137877812        20.776039127682587        22.342089469930325      -0.35876871924007198       -20.763729798181551       -1.4417970973144594       -21.953585449439000       -41.533394251077233       -40.578075424869859       -20.993783987221601       -2.8410549981538549       -44.141097550172312       -61.924926688003893       -42.458929724823996       -80.754838797741371       -59.576546414838717       -80.084988816765673       -41.780895561632960       -4.1475814897784735E-002  -2.0859422755437431      -0.28724003805125137       -14.425954620654283       -12.668932332462679       -5.7660402584630077       0.51983578922176843        1.1159763645707423        7.5792422133386488       -1.8964193344064704       -23.682272984053057       -11.077046166524061       -7.2524151780421278     
 /BM/amat          22   11.743583280872675        12.975082249502815        14.694460715585828        15.533588206696685        15.879061219071094       -2.0507310355858968        28.649499589752015        34.953848101197948        39.175863802170021        42.054336924778269      -0.68731907412391680       -39.647553925850445       -2.7758555469855675       -41.858726370447663       -79.408075569774638       -77.598090462893126       -40.039104109268273       -5.4917097577689713       -84.292353820713700       -118.51827684628194       -81.095306730309233       -154.60345139606704       -114.05254475699226       -153.32712496005536       -79.801606519551626      -0.11109140236148685       -4.3658589667503458      -0.52961523641098407       -27.233726882583927       -23.658768678439287       -10.933687788288463        1.0655031251892746        2.6871614392856671        15.777190928428636       -2.1459214759778931       -44.198180575638624       -21.869581052975764       -14.589596037514532     
 /BM/amat          23   8.5716404146837277        9.5751773150607509        10.847419323816881        11.472560956872172        11.734298499463836       -1.5110950268208616        21.139842085621041        25.812074450673411        28.947215138571067        31.090208723006491      -0.51776779481920365       -29.248116321250080       -2.0771640843908665       -30.892064523710868       -58.585387347886673       -57.243543401217309       -29.543187370095573       -4.1041995160360916       -62.207061622391279       -87.449331872486738       -59.834925353131794       -114.04976229197253       -84.136895743549744       -113.10385385523576       -58.876277499682146       -8.3593799296537430E-002  -3.2390120569429239      -0.39134618568649937       -20.081191298453142       -17.537830856618445       -8.0643019266522380       0.78804826932360605        2.0039609487185599        11.729069269431090       -1.5200950306633949       -32.725212025368059       -15.915306533128373       -10.875992955662085     
 /BM/amat          24   11.634796141604044        12.871393235317807        14.574364017262528        15.406936000506585        15.751063052401387       -2.0317101172094469        28.407199618789107        34.668444127520701        38.861524687042660        41.721242068743720      -0.69117081576308426       -39.323691467182684       -2.7811544701756521       -41.525160470472450       -78.770290409341001       -76.968345833868014       -39.713671602692543       -5.5004018144575291       -83.628596826059777       -117.58183154962208       -80.443355757984577       -153.35923309556853       -113.13243144183062       -152.08749588458863       -79.154374057943642      -0.11201398747681957       -4.3502088216685655      -0.52403671618025760       -27.013126695144550       -23.495495954356649       -10.823752568226960        1.0549607423061478        2.6858159703911992        15.723046215649999       -2.0627854429560415       -43.874038274566651       -21.633391358198207       -14.704100597567475     
 /BM/amat          25   6.0528074508256502        6.6866395022627856        7.6311363698237322        8.1045077299197601        8.3126523706305537       -1.0550017069820372        14.864327977670108        18.205149882989140        20.458348824776785        22.005278000442544      -0.36252979265246549       -20.436778489156204       -1.4468398988147757       -21.616729822450033       -40.889531175809864       -39.942401213494804       -20.665187320364311       -2.8493041303781546       -43.470841938117104       -60.979554622103784       -41.800747030565603       -79.499005774091302       -58.647855318268761       -78.833819242506280       -41.127474832019516       -4.2456267123400337E-002  -2.0704988085170171      -0.28166641562384631       -14.203427542888649       -12.503523221885668       -5.6539606829125502       0.50907495772941280        1.1152724669316321        7.5254567393176153       -1.8087486269316742       -23.354536200571850       -10.839823640270701       -7.3612152695901436     
 /BM/amat          26  -3.9623236501765038E-002   7.2514706587408287E-003   8.7502483611448072E-003   8.8437958943346894E-003   8.6124735680862749E-003  -1.3623839156416242E-003   1.7600753181513194E-002   2.2183083974548937E-002   2.5686016062062784E-002   2.8337062822483251E-002   2.1438579123302104E-004  -2.1006950655898561E-002  -1.7100024775763834E-005  -2.2418352125314595E-002  -4.1013225110264256E-002  -4.0007255356102278E-002  -2.1362386428489928E-002   2.3481060866853073E-004  -4.4009581189448019E-002  -5.9083669441008069E-002  -4.2686142141470595E-002  -7.7819929795650500E-002  -5.6852523244989017E-002  -7.6869218437019299E-002  -4.1646158448698012E-002  -3.9123991259432397E-005  -7.6048404815430917E-003  -9.7919017947154310E-004  -5.0911537684695104E-003  -1.4912671860524671E-002  -7.2641506539118401E-003   3.1751449853096797E-003   1.8426901215235603E-002   5.3033524020880782E-002  -4.1779032799871385E-002  -2.3405008506789636E-002   9.6109513530924541E-003   4.5157427708628795E-002
 /BM/amat          27 -0.37690318228238034       0.38958108627112009       0.32166290501529238       0.28809938721716249       0.27121720325844290       -5.6844716838516614E-002  0.61758778087018706       0.70009835715978530       0.76555910669914196       0.82283493166108723       -5.2985822618154668E-002  -1.1164678044968361      -0.17946944810453377       -1.1461938847298612       -2.3576793988335116       -2.3167769376554439       -1.1045004987700100      -0.35250779013688655       -2.4249740868174614       -3.6735324941308996       -2.3357050365290961       -4.7926058019910132       -3.5586030657824068       -4.7662821133925721       -2.3078659488187654       -4.8399530906802532E-002 -0.68542303717279163       -1.0612645344122465E-002 -0.23264458500813631      -0.37478146687752872      -0.33764884878441254       0.16127149352757511       0.91359192685847923        2.5221905035551462        2.3629069481414686      -0.55338729611814164      -0.13657918099331590       0.23905572027768901     
 /BM/amat          28  -8.7279368423741760E-002   4.4467455141740254E-002   5.7233504567687506E-002   6.4054761079760519E-002   6.7611646825635391E-002  -7.8275694055838542E-003  0.11795273419468670       0.15069041931133845       0.17472991352761316       0.19196267726246644        2.0540038939800285E-003 -0.13588809475860245        1.0971101524835415E-003 -0.14516739543888330      -0.26362850097509161      -0.25855662015364178      -0.14004921767899189        2.0952302420934599E-004 -0.28348323239280337      -0.38321026282498877      -0.27529906187916720      -0.50259322880511537      -0.37167204976528173      -0.49920979336900945      -0.27184380366256961        1.2341571346143713E-003  -3.8160184925000503E-004  -6.3396064733994683E-003  -9.9334243910705652E-002 -0.10797059724313106       -3.4658327786665971E-002   2.0159174984381795E-003  -4.1442444590329974E-003   2.5009972422890953E-002  -6.9453749383217467E-002 -0.21004630149117401       -1.6852359437139836E-002  0.29492182763144154     
 /BM/amat          29   5.7496207716269492        2.2869781418590471        2.7118217867561705        2.9200209433974269        3.0039571368553353      -0.33754410308418675        5.0473946974061086        6.2206637031145027        6.9588424327278151        7.4182353785897863      -0.19809667283192839       -7.0502652163976043      -0.64525374177020178       -7.5371154377725293       -14.054016295878098       -13.665223913846910       -7.1488158336987784       -1.2007560739746546       -15.090262121408543       -20.882258854610456       -14.412851307763313       -27.085992439598364       -19.932068033873726       -26.812824380868495       -14.139624741998118        1.9918034700790135E-002 -0.18572880394997740       -9.4336551913037359E-002  -6.1670280801709900       -4.6315784102616586      -0.90873150198207431        7.1918622374887331E-002 -0.52533228069274174       0.43558795505248460       -2.3281481690655927       -9.2680129045501172       -4.9565966673413095       -3.2812594476316792     
 /BM/amat          30  0.70565086982445724        1.9687673490678146        2.4415220983204011        2.7441610973796298        2.9471659347801045      -0.30975149414206105        4.7290955280345415        6.0756077241594575        7.0981091600177422        7.8900080216138520      -0.14070973830938466       -5.9914443913385549      -0.48967491280471009       -6.4787750660165102       -11.850915674062447       -11.542300771559638       -6.1698753337955923      -0.90731280910428724       -12.805514301360205       -17.523686759613078       -12.267922173141306       -22.696275428182265       -16.787095515099409       -22.496913980949486       -12.067885381744663        1.9426397959978681E-002 -0.18215348757444577      -0.10906064897476452       -4.6002746534220895       -5.5091615784243215       -1.4409508619841802        9.9545854609510043E-002 -0.41595951476643361       0.84337768087790943       -2.1346305041704037       -9.7147554594670602      -0.30054382971987603       -3.0544461133940519     
 /BM/amat          31  -2.0212639194329509       0.88284609486492327        1.0415630202583435        1.1331504462128934        1.1912911629171634      -0.16433853166453913        2.2517402073138788        2.8053372453449583        3.2354031113929196        3.5767505358272147        3.4041831479860965E-002  -2.7014629714320320        2.9875310723421200E-003  -2.8368034700019442       -5.3362037601041470       -5.2495195571364466       -2.7484477346967870       -4.2032915984336172E-002  -5.6357898428343258       -7.8786786225639966       -5.4925049520215206       -10.365478002249997       -7.6767784410620497       -10.306312110179796       -5.4312603588748578        1.1619751218022761E-002 -0.20013615145956903       -3.6020033249127989E-002 -0.91311124109552089       -1.4668233776124608       -2.3997866550184268        3.6568105317080581E-002  0.18910798684476726        1.0291099499352165       -1.2572611853123878       -2.2498341340041028      -0.37613520133465150      -0.93720913864983757     
 /BM/amat          32 -0.55093860151523011       -8.1135239233164430E-002  -7.0948610251519162E-002  -7.0057129665318213E-002  -7.3496662330286275E-002  -1.5748295731515399E-003  -3.8669474215694617E-002  -5.8271174019391553E-002  -7.4756658166608064E-002  -8.9589482533584813E-002   7.9029677937243054E-002  0.24126653852888308       0.19152984761024877       0.28898363958851864       0.54696723373094580       0.49800699327893894       0.23901261035442828       0.32140696592737833       0.61784729592766285       0.88497800049926356       0.53028300865386180        1.0541204617400650       0.76494164565756662        1.0211542359170258       0.49600885640715753        7.4366117933316982E-003  0.13020472747802186        1.0461749536704501E-003  0.11793530797794120       0.12184944402779743        3.4614569133766947E-002 -0.17494933064883803      -0.10755174366148182      -0.43300258151908999      -0.52547387875873008       0.22694718537914918        2.1093504775282490E-004 -0.55647647759696839     
 /BM/amat          33   1.1294884947289727       -6.6637856399912515E-002   9.3847960383826104E-002  0.18316106990028252       0.23370987805825652        1.0170075938616491E-002  0.12279514661116624       0.24947567421995395       0.33903143940834068       0.39505787187666130        1.2154644721790939E-003  0.11495178884417923        4.5428304207012715E-002   4.5924976219233138E-002  0.37555826839864748       0.40410383657174775        7.4646221882703201E-002  0.13288757183828634       0.25956618193149428       0.75301639525989961       0.28432884531380953        1.0417021242218292       0.79738681884898455        1.0611879149886410       0.30366274843920982        4.4611684785131639E-002  0.70114900891836429       -8.2821005953617136E-003 -0.70303119352791421      -0.59842182549954293       0.25312962900224922       -9.1153653477021326E-002  -1.2712552222395281       -2.6852260893408646       -2.2621678509417347       -1.2154111218367796      -0.17165923112296211      -0.41478527800814213     
 /BM/amat          34   1.4305394377857339       -1.5094151646440888       -1.0735308958238192      -0.89282586360009331      -0.82788857173930774       0.20270492254077888       -1.9787857281648513       -2.1784140084920782       -2.3951962651151195       -2.6319383200759745       0.29489579178738218        4.1821680013993374       0.93426436905603405        4.2703397804644432        9.0486181315079151        8.9108867574119053        4.1307784945965995        1.8105753980588539        9.2140424336325637        14.382451786819924        8.8642927741157749        18.709502875629017        13.968527790570786        18.643466690355559        8.7939579172384050       0.25307441183152868        3.4328789905649648        3.5153109756055942E-002  0.11124909368910725        1.0705837286538185        1.4325708722338806      -0.73941703874357667       -4.4857695016306325       -12.421584601275853       -13.905324577267248        1.3409778807720516      -0.87216949720705550       -3.2261699657868910     
 /BM/amat          35   2.6101253121228298       0.14736112201854767       -8.8294523118394530E-002 -0.15235988446758522      -0.15995651812409298       -5.1513147825144046E-006 -0.18927656806389592      -0.38563754978505749      -0.54173341561803889      -0.64768104056696996       -8.3619528618223582E-002 -0.40675323824158438      -0.26733185308109642      -0.33946007028726166       -1.1044554790166368       -1.1221516660183939      -0.36113326748537578      -0.58457784107337174      -0.98781694772567352       -2.0944903868027733      -0.94538561365447560       -2.6994743273592210       -2.0825015228839749       -2.7331357860167298      -0.98562844515339698       -5.3037464915727639E-002 -0.48629675006101403        7.9825850400194520E-002 -0.18682216380175290       0.52069455368405360      -0.16495141862213872       -8.9268595407390156E-002  -4.0713828343431790E-002   2.2287669273546307E-002   5.5599647523897877       0.84803733281825211       -1.3305617658431648       -6.5797564058439626     
 /BM/amat          36   4.0076836892096770        3.6472529227349408        4.5302192516900677        5.0693434199275496        5.4067522391869485      -0.55474806987880632        8.5844652632044980        10.994487475695523        12.769029219281283        14.092666348826599      -0.31690931388793275       -11.114624838794903       -1.0242947263194087       -12.042053564813498       -21.994646423443392       -21.385368818290093       -11.433325159735981       -1.8586135116323610       -23.814984646269505       -32.512622337141387       -22.757377114732016       -42.049408711169669       -31.060043710118276       -41.653925457289418       -22.361735830133089        4.4375856743992061E-002 -0.23783293598036431      -0.21971380287295195       -9.0886669846956387       -9.6284407252575654       -2.1738459609922525       0.19021748209688633      -0.85069702547276316        1.2784939555596644       -4.1397111112203406       -17.788878931631487       -1.7499720064186062       -2.1842258277091040     
 /BM/amat          37   7.8496866454289664        2.1653601211553739        2.1621275370505644        1.9892749388723845        1.7328441758342872      -0.35051179776788460        4.3181299900825838        4.5633681285236918        4.3437368422040867        3.9099541245342273        5.2009119311475505E-002  -6.9249092153563483      -0.10235765943311473       -6.9753584360571290       -14.021950820670039       -13.764126698077828       -6.7156899376049779      -0.42906333709327965       -14.487594570322003       -21.099794587202638       -14.029754018614074       -27.891879195417346       -20.390325466443763       -27.639166311892140       -13.774199032868662       -1.3207726321222090E-002 -0.42304161081292757        5.7672571947053275E-002  -6.0046256069703032       -1.0194688313166611      -0.61046522450837637      -0.14363955853254082      -0.31631333642820852      -0.34159406094205846       -1.6104696135807364       -3.0547527433983421       -14.610340627049613       -15.288737508115990     
 /BM/amat          38  -3.5532154320243152        6.4746661419626927E-002   5.0633035608366642E-002   9.0501402786939036E-002  0.15754500741422556       -4.6931320229080219E-002  0.40944556454988462       0.56287456895239940       0.74187417131560973       0.95325997354680403        9.2913780725612405E-002 -0.17661373248916823       0.16232430554750787      -0.14935243686595079      -0.36251580584711784      -0.41104981634847637      -0.19902156837377505       0.19046494334079012      -0.32732715105495569      -0.61892576536179833      -0.38866419189489509      -0.87060544962553399      -0.70390691693493812      -0.89519306924699782      -0.41513521487347482       -1.4394985471016197E-002  -6.3156209782962214E-002  0.13796436387563343       0.14669967473901652       -1.2674246429873268      -0.47298814942735168      -0.21000936364346745      -0.19136445880224090      -0.17531813713108613       -2.1658165013786359       -1.0088731090789054        1.9963598226066583       -24.182373010346900     

Solving generalized eigenvalue equation of linear method with a_diag =  1.0D-08
eigval_srt_ind_to_eigval_ind=     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35    36    37    38    39
eigval_srt_ind_to_eigval_ind=    32    39    35    36    38    37    34    33    31    30    29    28    27    26     7     8    22    23    25    24    20    21    19    18    17    16    14    15    12    11    13    10     9     6     5     3     2     4     1
eigval_ind_to_eigval_srt_ind=    39    37    36    38    35    34    15    16    33    32    30    29    31    27    28    26    25    24    23    21    22    17    18    20    19    14    13    12    11    10     9     1     8     7     3     4     6     5     2
Sorted (complex) eigenvalues:
eigenvalue #    1:   -14.679489 +    0.000000 i
eigenvalue #    2:   -14.609198 +    0.000000 i
eigenvalue #    3:   -14.410678 +    0.170708 i
eigenvalue #    4:   -14.410678 +   -0.170708 i
eigenvalue #    5:   -13.897604 +   -0.102730 i
eigenvalue #    6:   -13.897604 +    0.102730 i
eigenvalue #    7:   -13.417048 +   -0.204235 i
eigenvalue #    8:   -13.417048 +    0.204235 i
eigenvalue #    9:   -12.774874 +   -0.275914 i
eigenvalue #   10:   -12.774874 +    0.275914 i
eigenvalue #   11:   -11.358591 +    0.000000 i
eigenvalue #   12:   -10.852153 +    0.000000 i
eigenvalue #   13:   -10.108378 +    0.000000 i
eigenvalue #   14:    -9.653332 +    0.000000 i
eigenvalue #   15:    -9.240154 +   12.116860 i
eigenvalue #   16:    -9.240154 +  -12.116860 i
eigenvalue #   17:    -7.076850 +    0.786420 i
eigenvalue #   18:    -7.076850 +   -0.786420 i
eigenvalue #   19:    -6.359880 +   -0.170831 i
eigenvalue #   20:    -6.359880 +    0.170831 i
eigenvalue #   21:    -2.209525 +    1.039794 i
eigenvalue #   22:    -2.209525 +   -1.039794 i
eigenvalue #   23:    -1.486356 +   -2.214335 i
eigenvalue #   24:    -1.486356 +    2.214335 i
eigenvalue #   25:     1.106132 +    0.000000 i
eigenvalue #   26:     1.790527 +    0.000000 i
eigenvalue #   27:     1.834064 +    1.387192 i
eigenvalue #   28:     1.834064 +   -1.387192 i
eigenvalue #   29:     1.903859 +   -6.886607 i
eigenvalue #   30:     1.903859 +    6.886607 i
eigenvalue #   31:     4.874874 +    0.000000 i
eigenvalue #   32:     6.040722 +   -2.975772 i
eigenvalue #   33:     6.040722 +    2.975772 i
eigenvalue #   34:     9.346986 +    0.000000 i
eigenvalue #   35:    21.853763 +    0.000000 i
eigenvalue #   36:    27.033324 +  -34.462280 i
eigenvalue #   37:    27.033324 +   34.462280 i
eigenvalue #   38:    42.080369 +    0.000000 i
eigenvalue #   39:    59.534737 +    0.000000 i

Reasonable eigenvalue window is:   -15.20993   -14.59920
(sorted) eigenvector with smallest norm of wavefn variation is    2:   -14.6092 +    0.0000 i, norm change=  0.000000
(sorted) eigenvector with largest first coefficient is #          2:   -14.6092 +    0.0000 i, norm change=  0.000000
(sorted) eigenvector with lowest reasonable eigenvalue is #       1:   -14.6795 +    0.0000 i, norm change=  Infinity
Warning: lowest_eigval, largest_1st_coef, smallest_norm eigenvecs are:   1:  -14.679 +   0.000 i   2:  -14.609 +   0.000 i   2:  -14.609 +   0.000 i
selected (sorted) eigenvector is #                                2:   -14.6092 +    0.0000 i, norm change=  0.000000

add_diag, delta_lin1=  1.0E-08  1.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00

add_diag, delta_lin2=  1.0E-08              0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00

CSF parameters variations=   0.0000000
>>>>>>> master
Jastrow parameters variations=   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Exponent parameters variations=   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Orbital rotations parameters variations=   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
CSFs coefficients:
csfs
 csf_coef
 CSF norm:    1.00000000    
<<<<<<< HEAD
     0.94953691    -0.18108898    -0.18108898    -0.18108898
=======
     0.98229579    -0.18733653
>>>>>>> master
 end
end
Jastrow parameters:
jastrow
 parameters
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (a(iparmj),iparmj=1,nparma)
  0.500000000000       1.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (b(iparmj),iparmj=1,nparmb)
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (c(iparmj),iparmj=1,nparmc)
 end
end
Orbital coefficients:
orbitals
 coefficients
  9.33688390E-02  9.15836262E-01  5.50412700E-03 -6.04136300E-03 -6.73312300E-03  0.00000000E+00  0.00000000E+00  0.00000000E+00 (coef(i,j),j=1,nbasis)
 -1.44796690E-02 -1.63225017E-01 -8.31691630E-02  6.24493932E-01  4.59069428E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00
 -1.43342395E-01  1.28867516E-01 -6.47441737E-01  3.32194479E+00 -2.98668907E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  9.16554740E-01 -2.72166599E+00  3.37048517E+00 -3.00575054E+00  1.73605913E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
 -4.04561705E+00  5.65897624E+00 -3.23642173E+00  1.78257227E+00 -9.88421902E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
 end
end
basis
 basis_functions
  1
1s         6.28517900
1s         3.45549700
2s         2.77411700
2s         1.19273400
2s         0.82453500
2px        0.98665600
2py        0.98665600
2pz        0.98665600
 end
end


OPT: iter    energy         error      diff          sigma                grad norm        p_var  nxt err  nxt stab
<<<<<<< HEAD
OPT:  1   -14.5969650 +-  0.0770999   0.00000  0.61951 +-  0.06541     3.15976 +-  0.80133 0.000  0.03855  0.0E+00
=======
OPT:  1   -14.6091980 +-  0.0692275   0.00000  0.60073 +-  0.05132     2.85952 +-  0.93381 0.000  0.03461  1.0E-08
>>>>>>> master

For next iteration #   2 new wave function:
CSFs coefficients:
csfs
 csf_coef
<<<<<<< HEAD
     0.94953691    -0.18108898    -0.18108898    -0.18108898
=======
     0.98229579    -0.18733653
>>>>>>> master
 end
end
Jastrow parameters:
jastrow
 parameters
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (a_new(iparmj),iparmj=1,nparma)
  0.500000000000       1.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (b_new(iparmj),iparmj=1,nparmb)
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (c_new(iparmj),iparmj=1,nparmc)
 end
end
Orbital coefficients:
orbitals
 coefficients
  9.33688390E-02  9.15836262E-01  5.50412700E-03 -6.04136300E-03 -6.73312300E-03  0.00000000E+00  0.00000000E+00  0.00000000E+00 (coef_new(i,j),j=1,nbasis)
 -1.44796690E-02 -1.63225017E-01 -8.31691630E-02  6.24493932E-01  4.59069428E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00
 -1.43342395E-01  1.28867516E-01 -6.47441737E-01  3.32194479E+00 -2.98668907E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  9.16554740E-01 -2.72166599E+00  3.37048517E+00 -3.00575054E+00  1.73605913E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
 -4.04561705E+00  5.65897624E+00 -3.23642173E+00  1.78257227E+00 -9.88421902E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
 end
end
basis
 basis_functions
  1
1s         6.28517900
1s         3.45549700
2s         2.77411700
2s         1.19273400
2s         0.82453500
2px        0.98665600
2py        0.98665600
2pz        0.98665600
 end
end

 
Optimization ended after   2 iterations.
Warning: Convergence not reached.
optimization: Maximun number of iterations   1 reached.

Wave function at final iteration:
CSFs coefficients:
csfs
 csf_coef
<<<<<<< HEAD
     0.94953691    -0.18108898    -0.18108898    -0.18108898
=======
     0.98229579    -0.18733653
>>>>>>> master
 end
end
Jastrow parameters:
jastrow
 parameters
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (a_new(iparmj),iparmj=1,nparma)
  0.500000000000       1.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (b_new(iparmj),iparmj=1,nparmb)
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (c_new(iparmj),iparmj=1,nparmc)
 end
end
Orbital coefficients:
orbitals
 coefficients
  9.33688390E-02  9.15836262E-01  5.50412700E-03 -6.04136300E-03 -6.73312300E-03  0.00000000E+00  0.00000000E+00  0.00000000E+00 (coef_new(i,j),j=1,nbasis)
 -1.44796690E-02 -1.63225017E-01 -8.31691630E-02  6.24493932E-01  4.59069428E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00
 -1.43342395E-01  1.28867516E-01 -6.47441737E-01  3.32194479E+00 -2.98668907E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  9.16554740E-01 -2.72166599E+00  3.37048517E+00 -3.00575054E+00  1.73605913E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
 -4.04561705E+00  5.65897624E+00 -3.23642173E+00  1.78257227E+00 -9.88421902E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
 end
end
basis
 basis_functions
  1
1s         6.28517900
1s         3.45549700
2s         2.77411700
2s         1.19273400
2s         0.82453500
2px        0.98665600
2py        0.98665600
2pz        0.98665600
 end
end


OPT: the best wave function was found at iteration #   1  -14.6091980 +-  0.0692275
Best wave function:
CSFs coefficients:
csfs
 csf_coef
<<<<<<< HEAD
     0.94953691    -0.18108898    -0.18108898    -0.18108898 (csf_coef_best(icsf),icsf=1,ncsf)
=======
     0.98229579    -0.18733653 (csf_coef_best(icsf),icsf=1,ncsf)
>>>>>>> master
 end
end
Jastrow parameters:
jastrow
 parameters
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (a_best(iparmj),iparmj=1,nparma)
  0.500000000000       1.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (b_best(iparmj),iparmj=1,nparmb)
   0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000       0.00000000000     (c_best(iparmj),iparmj=1,nparmc)
 end
end
Orbital coefficients:
orbitals
 coefficients
  9.33688390E-02  9.15836262E-01  5.50412700E-03 -6.04136300E-03 -6.73312300E-03  0.00000000E+00  0.00000000E+00  0.00000000E+00 (coef_best(i,j),j=1,nbasis)
 -1.44796690E-02 -1.63225017E-01 -8.31691630E-02  6.24493932E-01  4.59069428E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  1.00000000E+00  0.00000000E+00
 -1.43342395E-01  1.28867516E-01 -6.47441737E-01  3.32194479E+00 -2.98668907E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  9.16554740E-01 -2.72166599E+00  3.37048517E+00 -3.00575054E+00  1.73605913E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
 -4.04561705E+00  5.65897624E+00 -3.23642173E+00  1.78257227E+00 -9.88421902E-01  0.00000000E+00  0.00000000E+00  0.00000000E+00
 end
end
basis
 basis_functions
  1
1s         6.28517900
1s         3.45549700
2s         2.77411700
2s         1.19273400
2s         0.82453500
2px        0.98665600
2py        0.98665600
2pz        0.98665600
 end
end


Restoring best wavefunction.

Exit of menu.


Some warnings were encountered.

Total CPU time is 0000h00m00s04c
The program ended normally.

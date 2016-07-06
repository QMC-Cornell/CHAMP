
----------------------------------------------------------------------------------------------------------
                                                  PROGRAM CHAMP
                                                 version 3.08.00
----------------------------------------------------------------------------------------------------------

GIT commit 112fa4cd20d6c214eb733113b751c1d415c9a449, Date:   Tue Jul 5 16:55:14 2016 -0400
Compiled by mussard on Wed Jul  6 09:19:49 EDT 2016 on host antares
Executed by mussard on 2016-07-06 09:23:31 (-0400) on master host 
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

Beginning of linearresponse menu -------------------------------------------------------------------------
 Requested parameter types: jastrow   csfs      orbitals  exponents 

 Orbital parameter information:
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

 Exponent parameter information:
 Number of exponent parameters =       6
 Exponent parameter #   1 corresponds to basis functions:   1_1s
 Exponent parameter #   2 corresponds to basis functions:   2_1s
 Exponent parameter #   3 corresponds to basis functions:   3_2s
 Exponent parameter #   4 corresponds to basis functions:   4_2s
 Exponent parameter #   5 corresponds to basis functions:   5_2s
 Exponent parameter #   6 corresponds to basis functions:   6_2px   7_2py   8_2pz

 Number of Jastrow parameters:      24
 Number of periodic Jastrow parameters:     0
 Number of CSF parameters:           1
 Number of orbital parameters:       7
 Number of exponent parameters:      6
 Number of geometry parameters:      0
 Total number of parameters:        38

End of linearresponse menu -------------------------------------------------------------------------------

************************************** LINEAR-RESPONSE CALCULATION ***************************************
1 configuration for  4 electrons has been generated by routine sites.
sites:   -0.313319    0.628839   -0.208193   -0.881718   -0.193128    0.281660    0.605008   -0.230834   -0.536871   -0.024506   -0.056022   -0.310910

***************************************** START VMC CALCULATION ******************************************

The following averages will be calculated:
- dpsi_av                                            = average of dpsi
- dpsi_eloc_av                                       = average of dpsi_eloc
- dpsi_dpsi_av                                       = average of dpsi_dpsi
- dpsi_dpsi_eloc_av                                  = average of dpsi_dpsi_eloc
- deloc_av                                           = average of deloc
- dpsi_deloc_av                                      = average of dpsi_deloc
- eloc_av                                            = average of eloc
- eloc_sq_av                                         = average of eloc_sq

The following variances will be calculated:
- eloc_av_var                                        = variance of average eloc_av
- object   42                                        = variance of average sigma

The following statistical errors will be calculated:
- error_sigma                                        = statistical error of average sigma

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
Beginning of equilibration (total CPU time is       0.05 s, CPU time since last check is       0.05 s)
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
End       of equilibration (total CPU time is       0.05 s, CPU time since last check is       0.00 s)
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
End       of accumulation (total CPU time is       0.08 s, CPU time since last check is       0.03 s)

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
overlap eigenvalue #   34:  1.38593571E+00
overlap eigenvalue #   35:  1.60153157E+00
overlap eigenvalue #   36:  2.03882679E+00
overlap eigenvalue #   37:  2.66468804E+00
overlap eigenvalue #   38:  6.26400781E+01
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

Solving generalized eigenvalue equation with a_diag =  1.0D-08
eigval_srt_ind_to_eigval_ind=     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35    36    37    38    39    40    41    42    43    44    45    46    47    48    49    50    51    52    53    54    55    56    57    58    59    60    61    62    63    64    65    66    67    68    69    70    71    72    73    74    75    76
eigval_srt_ind_to_eigval_ind=    68    30    32    31    69    70    35    36    73    74    75    76    38    37    72    71    34    33    67    29    28    66    65    27    26    64    46    45     7     8    24    25    62    63    61    60    22    23    20    21    58    59    57    56    18    19    55    17    16    54    53    52    15    14    11    12    50    49    13    51    47    48    10     9    44     6    43     5    41    42     3     4     2    40    39     1
eigval_ind_to_eigval_srt_ind=    76    73    71    72    68    66    29    30    64    63    55    56    59    54    53    49    48    45    46    39    40    37    38    31    32    25    24    21    20     2     4     3    18    17     7     8    14    13    75    74    69    70    67    65    28    27    61    62    58    57    60    52    51    50    47    44    43    41    42    36    35    33    34    26    23    22    19     1     5     6    16    15     9    10    11    12
Sorted (complex) (unique) eigenvalues:   38
eigenvalue #    1:   -14.679489 +    0.000000 i (    2)
eigenvalue #    2:   -14.410678 +   -0.170708 i (    2)
eigenvalue #    3:   -14.410678 +    0.170708 i (    2)
eigenvalue #    4:   -13.897604 +    0.102730 i (    2)
eigenvalue #    5:   -13.897604 +   -0.102730 i (    2)
eigenvalue #    6:   -13.417048 +    0.204235 i (    2)
eigenvalue #    7:   -13.417048 +   -0.204235 i (    2)
eigenvalue #    8:   -12.774874 +   -0.275914 i (    2)
eigenvalue #    9:   -12.774874 +    0.275914 i (    2)
eigenvalue #   10:   -11.358591 +    0.000000 i (    2)
eigenvalue #   11:   -10.852153 +    0.000000 i (    2)
eigenvalue #   12:   -10.108378 +    0.000000 i (    2)
eigenvalue #   13:    -9.653332 +    0.000000 i (    2)
eigenvalue #   14:    -9.240154 +  -12.116860 i (    2)
eigenvalue #   15:    -9.240154 +   12.116860 i (    2)
eigenvalue #   16:    -7.076850 +    0.786420 i (    2)
eigenvalue #   17:    -7.076850 +   -0.786420 i (    2)
eigenvalue #   18:    -6.359880 +   -0.170832 i (    2)
eigenvalue #   19:    -6.359880 +    0.170832 i (    2)
eigenvalue #   20:    -2.209526 +    1.039794 i (    2)
eigenvalue #   21:    -2.209526 +   -1.039794 i (    2)
eigenvalue #   22:    -1.486357 +   -2.214335 i (    2)
eigenvalue #   23:    -1.486357 +    2.214335 i (    2)
eigenvalue #   24:     1.106132 +    0.000000 i (    2)
eigenvalue #   25:     1.790535 +    0.000000 i (    2)
eigenvalue #   26:     1.834060 +   -1.387192 i (    2)
eigenvalue #   27:     1.834060 +    1.387192 i (    2)
eigenvalue #   28:     1.903859 +    6.886606 i (    2)
eigenvalue #   29:     1.903859 +   -6.886606 i (    2)
eigenvalue #   30:     4.874876 +    0.000000 i (    2)
eigenvalue #   31:     6.040722 +    2.975774 i (    2)
eigenvalue #   32:     6.040722 +   -2.975774 i (    2)
eigenvalue #   33:     9.346986 +    0.000000 i (    2)
eigenvalue #   34:    21.853763 +    0.000000 i (    2)
eigenvalue #   35:    27.033324 +   34.462280 i (    2)
eigenvalue #   36:    27.033324 +  -34.462280 i (    2)
eigenvalue #   37:    42.080369 +    0.000000 i (    2)
eigenvalue #   38:    59.534738 +    0.000000 i (    2)

Exit of menu.


Total CPU time is 0000h00m00s08c
The program ended normally.

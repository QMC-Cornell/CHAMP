
----------------------------------------------------------------------------------------------------------
                                                  PROGRAM CHAMP
                                                 version 3.08.00
----------------------------------------------------------------------------------------------------------

GIT commit 0fe703577f13497b4775ed40e68b29f1426c6df9, Date:   Tue Jun 21 17:34:40 2016 -0400
Compiled by mussard on Thu Jun 23 14:12:17 EDT 2016 on host antares
Executed by mussard on 2016-06-23 14:12:25 (-0400) on master host 
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
Normalized CSF coefs=  0.949537 -0.181089 -0.181089 -0.181089
Initial CSF rotation coefs=  0.184198  0.184198  0.184198
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

 number of CSFs =     4
 CSF coefficients =   0.949537 -0.181089 -0.181089 -0.181089
 sum of square of CSF coefficients =   1.000000
 CSF #     1
 determinants in CSF:   1
 coefficients: 1.00000
 CSF #     2
 determinants in CSF:   4
 coefficients: 1.00000
 CSF #     3
 determinants in CSF:   3
 coefficients: 1.00000
 CSF #     4
 determinants in CSF:   2
 coefficients: 1.00000

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
 Requested parameter types: orbitals  

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

 Number of Jastrow parameters:       0
 Number of periodic Jastrow parameters:     0
 Number of CSF parameters:           0
 Number of orbital parameters:       7
 Number of exponent parameters:      0
 Number of geometry parameters:      0
 Total number of parameters:         7

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
- d2psi_av                                           = average of d2psi
- d2psi_eloc_av                                      = average of d2psi_eloc
- eloc_av                                            = average of eloc
- eloc_sq_av                                         = average of eloc_sq

The following variances will be calculated:
- eloc_av_var                                        = variance of average eloc_av
- object   42                                        = variance of average sigma

The following statistical errors will be calculated:
- error_sigma                                        = statistical error of average sigma

Beginning of equilibration (total CPU time is       0.06 s, CPU time since last check is       0.06 s)
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -15.10477 -15.10477(    0) -45.26187(    0)  30.15710(    0)  21.97741(    0)                    0.87500        10   0.00000(    0)
End       of equilibration (total CPU time is       0.06 s, CPU time since last check is       0.00 s)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.8890089401339843E-002 -0.42456732033242556       -1.2298504012193168      -0.82146560470276309      -0.87353647292629277      -0.20454596794362798       -2.9037645445181390E-003
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.38339161350641143      -0.71211375278351630       -15.108507887040176        23.643427581338152      -0.21059592652824244      -0.61732808107092585       0.29983179052467507     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.58842627633642763       0.55289420888811969        3.1629727757966313       0.62903228694009228       0.72997418669064085        1.4790601592084371       -1.9449177270559537     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    2.7963714374949102        1.2101929775821680        23.853852739205873       -7.6409148216079652        1.1339867367662897      -0.63531150467720199       -9.8692269715564258     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.26189135812890624      -0.10467077032168362        1.5587652420253624        1.3750033351569946       -1.8241380457589016        6.1143594561143759      -0.61858156710991574     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.97621062303881057       -2.1112690177886937        14.176281217062796        11.313127503024999       -2.8284763838136389        27.903697197453127       -66.333275449228722     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   5.1137466671781921E-002   1.8433773914650082E-002   9.7264911993985209E-002   1.6256902632709835        1.2878058167993309       0.83670483335315249       -1.1779486314080854     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    7.2030613771144961E-003 -0.24813883676579898        2.4145350938638357        29.665459768245093       0.85261246774622834       -1.4807059669592912       -5.2833691088507093     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.1757217261607986E-002 -0.15963555797651330       0.67382950694165411        1.8150787481743351        5.9458603497492858E-002 -0.30472533351339071      -0.12371923099699154     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.27390714852284287       -1.4528584604281853        7.7523532005990434        21.579101376710074       0.11352820539683284       -1.1671993536281706       0.21174270461074393     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   9.3534885752314820E-002 -0.41409101528536896      -0.32334939121790857        1.1142874671503438       -1.7595407538954415      -0.48639645557278199       0.30437174130327665     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.32013488632852966       -1.9803836731399551       -2.7775607494680519        31.302516747904797      -0.84873537185691861        7.8429945489855102E-002  0.46779389419769712     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.16248738723373962        7.5959032356893788E-002 -0.56040825793596727       0.81287033081922966        2.1329866695757262       0.59350788234542717      -0.97328874849618086     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.73776867258435885       0.84450751698263238       -5.8083684888327660        34.582525293722448        1.6181069013457932      -0.58300667421928742       -9.6815240579678914     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -4.3178010047982031E-002 -0.36537857917769162      -0.60483420608035965       0.96958215403612413      -0.11811979337849166       0.44586072629785084      -0.70239911751512507     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.10967186457152980       -1.0484579124501539       -5.6159037829663356        36.469793451200715       0.22817387208921994      -0.75660694277853802       -3.5257460270490499     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.8386896453366004E-002 -0.39052930870101554      -0.54299679252653055      -0.77207186702371100      -0.66642788127327945        3.5228894851689301        3.2484644955531032     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.45468496979883605       -1.1572330087675884       -7.3608693420765050        5.1511342578539718       -1.6138074610081259        18.853731096639468        20.225150912874298     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.7536267321573704E-002 -0.41443632260517982       -1.1266989954609492      -0.95421600358237990       -1.0071411910475181      -0.59750708918166007       0.35375662132166741     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.30792130398871792       -1.2705299855426693       -11.336816143081299       0.16353523742559631      -0.54879516222158342        5.6435468599049896E-002   3.2159562031342744E-002
 warning: d2psi not yet implemented for orbitals/orbitals parameters
    enow      eave  (eerr )    peave (peerr)    tpbave(tpberr    tjfave(tjferr    fave  (ferr)     accave     iter     sigma
 -14.46468 -14.46468(    0) -23.09322(    0)   8.62853(    0)  11.98111(    0)                    0.85000        10   0.57041(    0)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.1881060929892419E-002 -0.25967056140061662      -0.86079185075744946       0.31109177271020994        2.2966105895098524E-002 -0.26228285181169247       -9.3367894946785449E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.38687443074649358      -0.56785783557007585       -9.4386959513610087        30.806104266729204        6.7838115205496388E-002 -0.92010514295875123      -0.44859605890290660     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.12921256042447560       0.97050475377541490        2.2158293178373802       0.23210123850687348        3.0879619106792485E-002 -0.24367200282293444       -8.4401354648039967E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.44941055443529399        3.5611401114412606        13.386647841344313        8.0648235095697558       0.10024595975712665      -0.74279448606079745      -0.60416120473941959     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.30921150094837063       0.82089127932016936        1.7763021194501756       0.62691976869212895       -1.1930407301095969      -0.83831440359983256       0.52265586647671636     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   -1.1095220905168617        2.7188228709960649        11.398747629101905        12.222823530295187      -0.71362152223192821       0.23868692606186784       0.24050865863693585     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.28590090034072957      -0.36956080174189515       -1.1056978447789991       0.35876168089115373      -0.50785806879500439       -2.4553951772329099        1.2102295032304520     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    1.3212049294234289       -1.0700315389800763       -10.496056965634111        28.436522925722436      -0.88226545132174938       -2.3237473627821772        2.0802065378162125     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   8.0642346245825572E-002 -0.24275688561628214      -0.84312961289900334        1.2926497889073263      -0.11055146492067652       -1.2804555540328875       0.63582123584344330     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.15787915858805088      -0.92308453967795945       -6.5994049643119901        39.396799517461702      -0.31853913169221099      -0.83265507785272008        1.0154526958016867     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -1.6379500075541864E-003  -1.9206672137144436E-002 -0.92497881903006796        1.8498890461674092E-002   1.4693039278463769       0.30564128989372108      -0.74055433348264033     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.12480791046156725       0.68606847505526247       -11.223884387319478        32.809701274814692       0.88931360028522743       -1.5741897345898308       -4.1313022144455038     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -9.7620039063314873E-002 -0.40507349816149435      -0.95978530349146607       0.16616575206895057       -1.0353608962353296       -4.1324792998518267E-002  -8.6576146431985718E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.39635010113453006       -1.1006463689794963       -11.077410413244889        34.256348114249569      -0.19128427766634865      -0.54804433449227141       0.36972119858926705     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.8096862434716280E-002  0.12774406363430679        1.8754515444332904       -1.2082735112087515       -1.3459787547018474       -7.1217823053887255E-002   3.3788319317764798E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.40236645521453585       0.73663472159100896        7.5701407747790936       -4.9959100492076312      -0.38005868417909827        7.4404727974519705E-002   7.1855744939108263E-002
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.13253309518674927      -0.15995706778912919      -0.76848981540642469      -0.71294471688387173        1.0148811510873332      -0.34093830643612999      -0.18441524132178563     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.47430199883278856      -0.37379893057431374       -6.4292884260763197       -12.909949205535060       0.44219476475660763       -1.8458702265946429       0.81768079245704284     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.3893395098569392E-002  0.16473276206648227        1.8664541386576916E-002  0.42463041925485301      -0.77282159171563314       0.18707128014427174      -0.28767066210905501     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.25813823197891012       0.15685018478348273       -1.3890568384701487        19.332897085757313       -5.3204316721664952E-002 -0.65446940490802297      -0.67711913726272810     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.70932 -14.58700( 8649) -23.73700(45522)   9.15000(36873)  11.76154(15526)                    0.87500        20   0.47179( 9862)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -6.9324567226803760E-002 -0.14673799415976424       -1.2703329661592784       -1.3262778839808425       0.90459929294695218      -0.82412050039957141       0.15263736998555802     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.26314006043012411       0.15722157790875382       -13.016410623541462       -6.5854947041861740       0.31058260422961409       -1.7275294943557036       0.63633951762584939     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   1.9072746813994915E-002 -0.41549391534739927       -1.3248782085189954       -2.3071512340537943      -0.68803587304025060      -0.38117580039627652        2.7594217417108158E-004
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.16545566013642668      -0.70465044699682078       -15.171445431089159       -17.433417857913373      -0.20718287811531819       -1.1707530276459741       0.43977044448407593     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -6.7341511425200007E-002 -0.57617750535137524       -1.2516963900941684       -1.5391661190041728       -1.7190908999782415      -0.16398500856397097       0.12326708981202916     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.30906577235290533       -1.4117924402427047       -14.787260057061836        3.0022184480646992      -0.68230764286934664       0.24802752416012785      -0.13331545079797008     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -2.2364628180173720E-003  0.14308234303829973        1.0983748852843767      -0.33918376579669285        1.0789454325774028      -0.53283236164484482       -2.4668828065921133E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    3.8033285121813497E-002  0.64860398181126266        4.0210938619854444        7.5131179898547371       0.48446976122563939       -1.4589940153487970      -0.79715525480865390     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   1.8740690594789985E-002   4.4872857151650167E-002  0.36471667112467987        1.2590161833324403      -0.71610545349868138       0.34544357751447963      -0.54743815812684582     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.18110257454233153      -0.49132429607845096        3.1753592656121619        25.055259538564421        1.7463329744777396E-002 -0.89900082742053167       -2.2872052415622290     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.7105459558523062E-002  0.34611242122451785       0.20228011252340097       0.42660875323713510       -1.1755838397218787       0.74685147292175624      -0.57868201215187320     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.29667470443995159       0.78628389960359413       0.27485848108590188        17.217584065845081        5.8837985367432101E-002 -0.38115979499119096       -1.2634747898808782     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -9.0939021886091015E-002  -8.2019446024277395E-002  -1.0144136961133987       -2.0744149553632711        1.4209730801566580        1.5488965690315133      -0.70757217284429108     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.26572393344144429       0.72699069646067327       -12.628171162865033       -16.418745501729894       0.48525917861189555        2.1632796197850603       -6.7850910948918166     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -2.8198602160851228E-002 -0.35259790061289997       -1.4265068549196109       -1.4879424044983371      -0.38518286077691888       -1.1223670458516091       0.58961062907688078     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.11364358375366954      -0.45168012468131158       -16.117198944483839        6.0241372724951905      -0.43769685896953753      -0.46982270667363990       0.45812156184517216     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.3562529036897112E-002 -0.37047260719535552       -1.4108151202874892       -2.9938753908356763      -0.34512868187434936      -0.99771604009479964       0.51107354959480900     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.16504492518402517       -1.0439397659969156       -12.340266543700887       -60.791167344304498      -0.38602848120358596      -0.37218405922811476      -0.11325312953357668     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -4.8067545953935878E-002 -0.21286637816463466       -1.1818564261205005        3.5444831895225980E-002  0.36337198025743661       -1.0362015316974715       0.40607856252658997     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.25315363990924566      -0.13323712788022726       -12.823327650070315        32.505700810074586       -1.9707774937925460E-002  -1.1381773094553471       0.75713614131516449     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.61213 -14.59538( 5806) -27.39881(*****)  12.80343(*****)  13.48930(*****)                    0.85833        30   0.44047( 6498)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.87341336162116201        1.9576732452594252        4.6753270329338754       -3.9761509644684612        1.0073626903337327        1.4338804800990916       -2.4537061167747658E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    3.8903470113224738        10.386310770047265        15.672515427117725       -16.542114772447732        8.2013662627430685E-003   5.6459953117737580       -12.300582620373241     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.7201367330552190E-002 -0.22529001322773598       0.12551899960799134       0.38782844937219246      -0.85508475321119137      -0.62173038456389740       0.36226443181366291     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.34869382902636914       -1.0192289175640150      -0.99510481866839762        20.100410725657795      -0.47127166730210757        1.1313454190290079E-002  0.29500081887252699     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.10118546997093433       0.30193541360034037       0.44239579956175989        1.9048214992601211       -1.8036497481450597       0.32026070830403797      -0.32523531818481044     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.35212764189917933      -0.42281804496738912        6.8483636986619505        22.690239342768152      -0.28669729004066496      -0.63632555200372853       0.45911593748122681     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   6.6854225591786637E-002  -2.6648586134695923E-002  -8.8754368719643845E-002   1.8134748503809681        2.2054893676111833       -1.7616331718596172       0.10790992071121733     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.58203864778872538      -0.21363938024400375        1.1378971630482466        30.582611130569266        1.7437415771093780       -5.1279148683616862        1.8166727726177641     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.1031244178444530E-002 -0.17771264925718172      -0.41506711624300902        1.5025631839886495        1.0648484032965526       0.52583470472236593      -0.86965877578111872     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.20306808513184893      -0.62422217202747154       -2.0742795736814861        35.211663837981732       0.56197033419064790      -0.41022156656007330       -6.9315982892219257     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -3.0266205279999955E-002 -0.18472522607736758        1.0218476944225783        1.7546708795995820       0.21031530725572734       -1.0416720096771870       0.35304523628545470     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   -3.0500314237861259E-002  -1.4903172908648594        10.434340271982155        18.133093116328229       -9.8990899200401671E-003  -1.2823118795835304        1.0531174928654092     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.10366957245380555      -0.25661388735687141      -0.55261393154048355        1.4394395944063778      -0.61042738622501647      -0.81693027748177460       0.46139956515138597     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.33810445310061704       -1.2527034446898087       -4.1024963043644229        37.555499460687273      -0.44048644810375381      -0.16093283503540096       0.54776994562251047     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11261018929190472      -0.13721501340941408       0.34847253836597308        1.6879106306264382      -0.46478694235343981      -0.65450443189530727       0.36509629723811737     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.45450531351933660       -1.2631481772947293        4.4075973400078849        26.284951061093231      -0.31865493184819493       -9.7217273881819860E-002  0.40095531781886629     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11357143822616299      -0.12825138061978050       0.34647014314156810        1.6875289783800438      -0.58190704220079315      -0.64325901724081003       0.37564156614472738     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.44798691444019156       -1.2517341627242879        4.4063593380538233        26.292156877946947      -0.37172957095965586       -2.9291989394401507E-003  0.34912252849960040     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11303335281726860        8.6904538567917716E-002   1.8577603107759897       0.57827831966263543        1.1343782591905196       0.77552474866479715      -0.77198338028673641     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.52292043419427303        2.4130247808164845E-002   11.795124234098054        12.328019381092007       0.49088054174193879       0.71053887531640847       -9.4477672541458055     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.22246 -14.50215( 9173) -24.88218(*****)  10.38003(*****)  12.18645(*****)                    0.85625        40   0.57721(14425)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   8.2682103317820346E-002 -0.18704928018705053       -3.2892659357205194E-002  0.34398877746521478      -0.57349343103044581      -0.94028675015855590       0.53201902134863810     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.26250252102707283      -0.88502496831021449       -1.5181046594597505        15.685999435092887      -0.47800968463436871      -0.13805787359544927       0.39667006589119930     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.9804275310885328E-002 -0.29920618820799455       -1.2844772612872315      -0.35950929147439681      -0.16721865679105566       -1.0983080973640988       0.55745975038449891     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.18706048942271933      -0.37932718287207673       -14.500322042355531        28.870357597558371      -0.31881734442773019      -0.58982057045471215       0.68850660882556558     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -9.5905591094174361E-002 -0.29263547550148716      -0.81775614660640106       0.92794348430157580      -0.54424840903095073      -0.67045703536728130       0.38837685784931314     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.44490389956699322       -1.0101663917594146       -7.8138836243119059        38.058812455038520      -0.37195141816713601       -5.4142283135323745E-002  0.39032936458840861     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.13906345407808601      -0.50614146709966512      -0.65972304908084312        7.8312157742772368E-003  -1.4821474103317624        2.1001011793756947       -1.3261099577646849     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.56858441867707454       -1.4079390207730411       -9.5114716910321313        34.076656739689952      -0.19019520403813653        1.8400384077834502       -9.3539311023086427     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.16328657124161788      -0.19509054149120275       0.94147686051016022       0.26695723665553994      -0.11103976560578471       0.35024643282574042      -0.82414724967165831     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.54987705619659344      -0.80069156997718627        5.1758630418313771        12.971601941385327       0.37011594997792990       -1.8461282838436079      -0.74695900942659021     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.20696650600751104       0.47692302945953968        2.8180458175179313       0.83681252660122019       0.65122113608115495        9.6437183823624900E-002 -0.63664243262762021     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.83824631859442500       0.86207227259711372        21.301570713526573       -1.9854889724251639       0.52906365046656889       -1.5864812941342401       -2.1777576432907471     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -4.0822482700329640E-003 -0.14281733472598465      -0.69431088356306070       0.95696244734969638       0.64253081247201749       0.30035827058736581      -0.52879248818305391     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   -5.5498214723662488E-002 -0.25695532205250515       -6.8219348186575450        38.700402153448877       0.35506397552595631      -0.19167496575375914       -5.5796601551662617     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -3.5443395806609122E-002 -0.32554960339961242      -0.92011312816808311      -0.21669702806756291      -0.50257459057177112      -0.86378339160684514       0.47742657689788570     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.23221290066188088       -1.0015603652180745       -9.1459144245689306        13.453867661489681      -0.41744681111900578      -0.18558337562171812       0.33037019840395199     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -6.8121733523721820E-002 -0.30013833423187569      -0.87042538853757878       0.68878183554023920      -0.50227964713092776      -0.86322420200393835       0.48470810615673676     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.31113665819200070      -0.93465137126684972       -8.7730383828248613        37.295113562958356      -0.41829228703651522      -0.18013533725587935       0.52642079148597076     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.5915108437363081E-002  -8.3522699439812984E-002  0.63028409720804202        1.9026707713352786      -0.89497641881707057      -0.50626277076144399       0.32746005257976429     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.25856639450246005       -1.4118698191644867        7.7433851979649475        21.952952437131060      -0.46632064300496667       0.20165513751398906       0.19291049386205514     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.20866 -14.44345( 9023) -24.28594(*****)   9.84249(*****)  11.89270(*****)                    0.83500        50   0.63646(12648)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   5.0112143679467890E-002 -0.23079364569401320       -1.2693353602708386       -1.1794592480524684       0.44822126177579746      -0.20127023604541919      -0.12930366214986458     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.24375237186034351       0.29492456242948323       -15.702258500813983        17.131928331785254       0.26828525161761430      -0.90415108137906719       -1.4020157834505509     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.27420216530556613       0.54204765781073638        1.4384138566276268      -0.29552912255337760       0.77705041879023196      -0.59266086798802731       -4.8416453926315335E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.88559040384080889        2.2522800436236832        5.7680723399239691        9.9896263334505520       0.47721758711553885       -1.4910592292326608       -2.1852408664127245     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   7.0183621714438787E-002 -0.28053138051034543      -0.72632230988651647       -1.9702251938339237       0.37430072521400248        1.1276688257690493       -1.1243514381394599     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.40627291300654178       -2.7528819389182901        7.3828599277041587       -140.74820163843214       0.34537493139099934        1.3360142931709222       -13.777817442311196     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -2.6107222929630150E-002 -0.68278677836682189      -0.76602505730651238       -1.9522338446387317       -3.0885479695038338       0.68738140282522442      -0.38186614815770159     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.10994518399067527       -5.1245102997908187        9.1231083008241285       -155.67328465655308      -0.79171087774515814       0.39046162291674824       -1.4851145471326876     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   8.1975586523063235E-002 -0.52424387864331712      -0.83769623436596830       -1.5132403038359434       -2.8198880234566261      -0.49752437171677250       0.31487341126647117     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.53798533220045186       -4.0365642273474656        3.8085419413437211       -108.81416456445137       -1.5590152958112300       0.42765601995881086       -1.6551906059731802     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.33551686607588344       -4.4414730991737782E-002  -1.2691753389812817       -2.7157124174447631        1.6858470259796037      -0.27638230416421788      -0.33610228071961568     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    1.7262241460435601       0.30552133115928853       -9.3639349207457663       -67.896930757538826        1.4792119901147203       -1.5293431909496671       -5.7880411624254648     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   5.8544784647261434E-002 -0.41033207392881327       -1.3139054985971879      -0.78203821577380195       -1.0134563059405148       -1.3289027802103359       0.73284997753102366     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.39518699153486536       -1.0929497349422912       -13.925842442627451        12.739943791597533      -0.85607625277136579      -0.47823831546877177       0.61838017374848309     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.4848917748404709E-002  -7.0877961674390172E-002   1.3217280220822389        1.6860502750566684        1.3760973549242250      -0.61134051103440712      -0.31455264533675126     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.23253089218640560       -1.0026583074056716        12.547093870575019        14.061730907761040       0.80875883385573411       -2.7614688921173358      -0.11808210686490044     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -6.1741731387295176E-002 -0.27435880511405558      -0.34708363372339590        1.4957221451664573       0.23256845972680956      -0.28942755535027564      -0.23965262617968466     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.20891473438627842       -1.1973196543828164       -1.9122973026442551        34.289700710414692       0.21396510311184394       -1.4288729763390158      -0.81909725757453555     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.5158378681549615E-002 -0.10850265058368200      -0.28165546561165372        1.4168907132450070        1.4066340899951144        1.0933836522520481       -1.3046765693842570     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.31829780188061774      -0.36916468986049605       -1.5712650349044881        33.508871577320839       0.81953636483351811      -0.68649163638800281       -8.4170716987246497     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.82661 -14.50731( 9514) -26.91444(*****)  12.40713(*****)  13.16759(*****)                    0.84167        60   0.61442(10559)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.15898584722692694       0.15911793903926416        1.7150130212791601        1.4903358713051427      -0.30864385443578701       0.71388216429366291      -0.83419863714780662     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.44426704102550035      -0.53104251835949334        15.425448288405947        8.3426139843989464       0.25879990540453263      -0.84164906230938585       -1.9498330720274362     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   7.7513360464955966E-002  0.60593384880929513        2.0325346157605897        1.4359870837587705      -0.38572429843301531        1.4061381876339580       -1.1028917811154433     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    8.7875858215510444E-002   1.0145064915451358        18.099680853510073        3.2021372046379346       0.14949726357018250        1.1249682749243430       -8.0466018298440130     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -2.7631721315405228E-002 -0.36250543176795236      -0.65557071643887455        1.0018460888977441      -0.90653054314937653       -4.2913466967466035E-002 -0.20876932926910488     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.12796459370432609       -1.3780334035434849       -6.4037938817281308        37.402023267919503      -0.13281434786730584      -0.99995281500546485       0.17995391485939527     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   7.9848752924032550E-003   3.0448107594041600E-003  0.75121341823436039        1.8317875004322364        1.1696078718264287      -0.53158488719211838      -0.43298913012389983     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.20779276375921085      -0.67187363143370860        8.3573481529382718        20.371179406222598       0.87137895468967075       -2.7218314911968822       -1.6827710495721433     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.12370599426663831       -4.9887323611534168E-002  0.40641789183805443       0.61668450362608929      -0.75801444692578301       0.64427528720378180      -0.55200654841029340     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.36224552011526889      -0.66566371663073920        1.9059680641236849        22.262707689038717        6.2822180796435964E-002 -0.13046837846422110       -2.7917028662230066     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.23437671775944630       -1.2540742818079948E-002 -0.71348112684492027       0.77378837104328868        1.4049244626821875       -2.0493899662437429       0.65321424111880411     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.90445337164825279       0.20399591460618355       -6.3116821201815805        29.414932479541772       0.75554619897339825       -3.5793666800886568        2.1643510076670238     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -6.8383890130066710E-002  -9.5658585661670584E-002   1.2991356905893305        1.1552478775823007        1.8105288053545108        6.0324081945898937        6.9935121790216703     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.12289995685416720      -0.92647941037875936        12.034647752121412        8.5697769100169605       -2.6316986035374432        46.968660949040384       -2.7030154292045547     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   3.0325758184297853E-002  -3.3026093151813075E-002  0.13384320897455551        1.7056568048733578        1.8369187434677476      -0.76809329745525612      -0.40064878744044291     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.23288181230568303      -0.34792248598913339        2.7486762864506047        27.935844803605296        1.2843192358117872       -3.8516530210016655       0.38041791709485773     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -4.1657657626140383E-002 -0.36173811616151719        1.5698189167813357        1.2385640018008777      -0.53447503408206432      -0.63649211743214962        1.6617612115566099E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.10414049788978341       -2.1500965675961501        12.603146206303173        14.174250458949338      -0.11357487510776489       -1.5994145577165693        1.2308645569491206     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   1.2282926740246181        2.3909811952489570        4.3673593418350176       -4.6158285497433180       0.67682329891367166       -1.7172256509876462       0.70668679320259908     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    5.6261610198207563        12.980351172512650        12.296587820499495       -25.772180479001502       0.26738988111886508       -1.7355754485501964        1.3306691245830711     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.55099 -14.51355( 8176) -25.82345(*****)  11.30990(*****)  12.63426(*****)                    0.84643        70   0.62165( 8954)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.17969665066374957       -1.7697708338936902E-002  0.34707731338286540      -0.54947736705491279       0.69704727957674639       -1.6913223688622230       0.67312353360326171     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.66349030459665503      -0.16483444102805936        1.5835045429574783       -12.420449494216102       0.12391750549239297       -1.9856812848090633        1.4379189189022130     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -1.5910056484398239E-002  -5.9878435579394346E-002 -0.13116293820902980        1.3333016866756195        1.1079341536185561        2.1710166894334786      -0.70456448415528117     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.21188237984831085      -0.33072647097657953      -0.68551013704578434        32.074927475371062      -0.14600908015332353        7.8808055945163851       -26.055406984246790     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.19309515509883196       0.14271625635143240        1.8948536836115621        1.5256141544197448       0.30553172852960281       -1.0323423330235386       0.46824460356041597     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.58618484215289524      -0.54786417267203147        16.917985958365829        5.9658355262626710       -6.9607626302249981E-002 -0.73163660426534516       0.79889852640040970     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -1.3793567771165114E-002 -0.21581532384427812      -0.65712055894772148       0.83869425955682586       -5.9485303942451294E-002 -0.85219731620238492       0.42799491009429758     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.24879661756211147      -0.76231448365617860       -6.0227329919639798        31.963262248301906      -0.20571975275199611      -0.43637987014025043       0.57038653705474351     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.46265451050436879       0.74437603528165230        3.2995210777996928       -1.2435829324350349        1.3817933850444017       -7.7754317156077529E-002 -0.52966279958204898     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    1.8986398929415860        3.7413032286050547        15.382643177253176        1.8383356521292360       0.86518310877482985       -1.8904628417270604       -2.3143106247452288     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.4906919001113986E-002 -0.15687210890042924      -0.35999596766953612       0.75831655070063830       0.66474768012263985      -0.76747825222258503       0.15498648361445261     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.37252879387188065      -0.50401491995348779       -3.9336892193395370        27.606296973174643       0.18727769432993804       -1.4304803656229341       0.72700891929939360     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.10545235927870798       -4.9122207643402553E-003 -0.13547017245062232       -1.0940871175871405      -0.79640731085584626      -0.57026438005156310       0.28137404632527685     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.34070118583863540       -1.2314735651075008        4.1397501408291202       -67.235338748135078      -0.37802125387893609      -0.28767905056954973       0.29548613125010265     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11363825563338609      -0.40976107970600884       -1.1131050592427103      -0.90585927341229644      -0.74207105841934962      -0.14115596502074701       -5.4306256366185940E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.43016239379702220      -0.94441072296671003       -12.690178135216591        9.5177052798229926      -0.13457661463933429      -0.61385983787218390        4.4027408574892429E-002
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -3.6335666422471297E-002 -0.13579623948250000      -0.28924144795705087       0.21197773510179621       0.83003572605301623        4.8643968331861362        5.6327213198504511     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.21917753111284916      -0.69209840754992558       -6.6823411369419577E-002   1.4972957732692844       -4.5744868819134803        52.303397491269791       -74.358105018961538     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   1.8827670115813985E-002  -4.2840820276933153E-003 -0.44754977206640079       -1.4326752188889165        1.9128605976175561        4.7141914019557056        3.4842515708476820     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.23899008544632144       0.78265685670661211       -6.6156943693363877       -11.658937021833053      -0.47549552827225761        19.821214745095897        45.347612198495398     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.85203 -14.55586( 8176) -26.30283(*****)  11.74697(*****)  13.06038(*****)                    0.84688        80   0.63261( 7831)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -7.2022375026738639E-002 -0.19893301538713443       -1.2658048014526737       -1.4817326499293988       0.57353829195061234       -8.8246752538842524E-002  -7.4759905422048398E-002
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.37903238375465592       0.31691003896402281       -15.663892509464841        9.0097049142928025       0.21116163045386172      -0.43425827335100453      -0.76532965594336166     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   2.5959896780178588E-002   1.5358529034601061       0.78718880295183158        1.4035606375421379       -2.9501758718284443        1.0079315634433428      -0.55971191058441028     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   -2.7014563524330161E-002   3.5861085031041289        9.3637830571604663        15.498802546174664      -0.35690536386575589       0.13699435177412142       0.15259182571618446     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.12398080550637007      -0.54973456995427383      -0.84113030132033484      -0.48746213168991953       -1.5310622839451480       0.35153758876528313      -0.33085726279998995     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.44093546283890150       -1.6072101742786946       -10.100635238577505        14.924793310026525      -0.14766487719439442      -0.59447348466754990       0.18451787700707759     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.13461427240660784      -0.11785604835233239        9.1592162129316854E-002   1.5744118091428410       0.43787589319512255      -0.89397783841488909       0.32904965370822598     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.50100170265794719      -0.83817110369442926        1.7633075358641981        29.372875184118040        2.6425667685121442E-002  -1.0148272172882757       0.69191170897658283     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.20431293558293934       0.16704158102830263        1.0114625504136161        1.1338741197251969       0.15242567934108381      -0.87362108979133957       0.31626745489686375     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.72014735615968861       0.11442257737252627        7.9673205349621359        17.880030735990992       -3.0821626006236801E-002 -0.94250784886317041       0.44465630582351040     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11137604829689229      -0.54627444930394164      -0.96439756946680810      -0.13136827160043746       -2.2940552224783590       0.42605747854248388      -0.20344643952147509     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.40887252586069928       -1.8454108544046117       -11.463545130167709        28.854523716579234      -0.61380348135456730       0.30172492182454880        4.7681955178055679E-002
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.11911332315892952       0.40786483040148813        1.4172589374075282       0.22843528913578892      -0.27454197446815604      -0.81795901768155588       0.43724014105024667     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.36613108825047558        1.3719726771635918        7.5236586055690911        13.770468410547572      -0.27925408456749046      -0.19604314585322585       0.37735405687358364     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.7717812789312108E-002 -0.16106136535850069      -0.82424550213576020       0.69047340187566941       0.75022848115247942       0.54769187877508807      -0.52527266451007792     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.37383056969469808      -0.11512186772727623       -8.7399270705006362        38.683111675373226       0.24352442633841365       0.76637359978756003       -7.4753312057504839     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.1683273703666667E-002 -0.14056404434483549      -0.53299715495196387        1.1655671225031246       0.75033363905737649       0.54888915113218484      -0.52340417274760176     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.35982560749282888      -0.39482145672017011       -4.8389351469573443        36.967907892165528       0.24410575746583160       0.75858119046200467       -7.5039121326438680     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -8.5680286215120804E-002   1.3390597692262834E-002  0.45704153982836188        1.2394862945208933       -1.6265563375691272      -0.44468989703545986       0.32085168051594776     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.26254123898215781      -0.76439048084588768        3.7384458802065290        24.324307107678045      -0.78904695648531187       0.38749155179455452       0.16192186103163100     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.53880 -14.55396( 7269) -25.89373(*****)  11.33976(*****)  12.80172(*****)                    0.84167        90   0.61069( 7246)
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.18267465939855984      -0.63945807134022137      -0.65549293109649420      -0.87143484749697742       -3.2293022967029841      -0.68047433000684387       0.46456534423929696     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.90965763084664408       -3.4029920369351037       -2.6447262233170532       -37.303989010573204       -2.0211204059359322       0.55146751908455605      -0.92677139635888839     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -2.6642461586305521E-002 -0.30920238876301553      -0.41707404266243192       -1.8478159062820911       -3.7482885508718582E-002   5.2657926387351204E-003 -0.37844795197703512     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.22805111061388542       -4.5769143596233119        21.006733803400451       -221.22348943947031       0.12534318258902150      -0.85770139958746616       -3.3099295675756242     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:   8.6752694596217168E-003  -3.6447811902396328E-002  -1.1219630018187532       -1.4702730428129740        1.5363848836995100      -0.13345110149764339      -0.44984465764730103     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:    3.2164744412649604E-002  0.47789998119610760       -10.928093266432045       -18.467530806013070       0.92580383482475348       -2.1700614935283911       -1.0889008385539420     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -5.5446597535623594E-002  -9.2157639587902088E-002 -0.47341761081469497      -0.85668663063652739        1.1322276314854487        1.5754000975133806       -1.1200059245305338     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.32933836998395438      -0.34667773449176614       -3.3501080182228380       -22.814253199607855       0.43958886115077866        2.6301934706919972       -14.154267981989276     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.13166385772621972       -9.6918783150165902E-002  0.68480377605181519        1.2328455598229055E-002  -1.0851347070216211      -0.32231660907870602       0.21037427565738503     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.33838152618166156      -0.51641294560324347        2.6445266897039712        10.230416605397732      -0.46045962672216173       0.19897132774647353        5.9907187595752281E-002
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -3.7795251049248607E-002 -0.19957677167217011      -0.20386259540050033        4.0381698982556551E-002  0.35616420514481856       -6.5996823331246054E-002 -0.35151096636702500     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.24447600441692838      -0.68249954878490982       -2.9813792424204046        10.757379191777430       0.25729691328461202       -1.1603850386049890       -1.1571395283109847     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.15713850920510325      -0.13548034223206923       0.40860982505338606        1.5470535819116453       0.58747122415888864        2.6690865995750706        1.3659273260093752     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.64690002399516011       -1.0400112020744168        4.3503596920174745        25.692428532210197       -1.7943586223285894        19.585699802828231       -43.769280803846577     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb: -0.11268715947206885      -0.12744870915697173      -0.53463487553473510      -0.84706375361285047        1.0196857102235843        3.6480263755273072        1.2601242674406521     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.40525310748428306      -0.11191003403410681       -5.3177725951575301       -7.2430566266301293       -1.7512282311337035        24.095304393392546       -46.227132226914357     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  0.19483173747915103       0.14069791295054351        1.8869495925884414        1.5357599253292482       0.53850265345321735       -5.0419327660344584E-002 -0.56115468279235181     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:   0.70897995653978130      -0.55533159659457176        17.040704995485875        5.7779504887949917       0.49280408665266640       -1.7343529291262407       -1.4923722395584869     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 /BM/logicals:  F F F F T F F
 /BM/dpsi_orb:  -3.3573526257620750E-002 -0.24745974240919447      -0.90665356555551224      -0.14366305881796501       -6.8133530005979687E-002 -0.93626651473824041       0.45630441343956368     
 / BM/logicals:  F F F F T F F
 /BM/deloc_orb:  -0.22398424335149539      -0.68477278518303109       -8.8441434139665329        13.607971223997003      -0.21845129642451200      -0.58396716319699893       0.63832918843799524     
 warning: d2psi not yet implemented for orbitals/orbitals parameters
 -14.98397 -14.59696( 7710) -26.98448(*****)  12.38751(*****)  13.33591(*****)                    0.84250       100   0.61951( 6541)
End       of accumulation (total CPU time is       0.08 s, CPU time since last check is       0.02 s)

vmc_mov1      no title
Final results after         100. passes (nstep =     10, nblk =     10)
physical variable  average       rms error    rms er*rt(pass)   sigma                Tcor
total E =        -14.5969650 +-  0.0770999  0.77100  0.61951  0.61951 +-  0.06541    1.55
potential E =    -26.9844781 +-  2.3529745 23.52974
interaction E =    4.5036037 +-  0.2268205  2.26820
jf kinetic E =    13.3359136 +-  1.1894450 11.89445
pb kinetic E =    12.3875132 +-  2.2939148 22.93915
<r> =              0.0000000 +-  0.0000000  0.00000
<r2> =             4.3709380 +-  0.1698931  1.69893
<r3> =                0.0000 +-    0.00000     0.00
<r4> =                  0.00 +-      0.000      0.0
acceptance         0.8425000

Adding to ovlp_lin  1.0000E-08
 /BM/ovlp           1   4.0548124021958211E-002   5.7375631990177443E-002  0.14398474045002665       -9.6658496121299550E-002   3.3468950109978211E-002  -1.2844341225481634E-002  -2.3281304890665357E-002
 /BM/ovlp           2   5.7375631990177443E-002  0.21730772727820560       0.46025048236048899       -9.0751751020772750E-002  0.11387906964953087        1.5946325152817534E-002  -4.5960122403513382E-002
 /BM/ovlp           3  0.14398474045002665       0.46025048236048899        1.6178243391148612       0.17690101356909144       0.19731064829580100       0.22155919147040606      -0.10099926465259287     
 /BM/ovlp           4  -9.6658496121299550E-002  -9.0751751020772750E-002  0.17690101356909144        1.8469219028130790        4.6734521115659802E-002   3.4554944922683767E-002 -0.10738753166503297     
 /BM/ovlp           5   3.3468950109978211E-002  0.11387906964953087       0.19731064829580100        4.6734521115659802E-002   1.3761360359720254       0.19421272073921206        8.8840659667667771E-002
 /BM/ovlp           6  -1.2844341225481634E-002   1.5946325152817534E-002  0.22155919147040606        3.4554944922683767E-002  0.19421272073921206        2.2497497783329163       0.62245062883372881     
 /BM/ovlp           7  -2.3281304890665357E-002  -4.5960122403513382E-002 -0.10099926465259287      -0.10738753166503297        8.8840659667667771E-002  0.62245062883372881        1.3765110768213698     

Eigenvalues of overlap matrix of current wave function and its first-order derivatives:
overlap eigenvalue #    1:  1.91217022E-02
overlap eigenvalue #    2:  6.53386454E-02
overlap eigenvalue #    3:  9.88819081E-01
overlap eigenvalue #    4:  1.27578486E+00
overlap eigenvalue #    5:  1.70436860E+00
overlap eigenvalue #    6:  2.02214300E+00
overlap eigenvalue #    7:  2.64942308E+00
 dens_mat_wfdet_bld: dens_mat_wfdet=   2.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        1.8032406822777609        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        6.5586439240746267E-002   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        6.5586439240746267E-002   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        6.5586439240746267E-002   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
 dens_mat_wfdet_bld: dens_mat_wfdet_up_trace=   1.9999999999999998     
 dens_mat_wfdet_bld: dens_mat_wfdet_dn_trace=   1.9999999999999998     
 dens_mat_wfdet_bld: dens_mat_wfdet_trace=   3.9999999999999996     
 /BM/amat           1  0.16645270922901539       0.25848257782671591       0.68652224357789704       -1.1001282962718597        3.1840052903980121E-002  -6.2019239535582443E-002 -0.11916691890494666     
 /BM/amat           2  0.23082347328563924       0.89308580920190384        2.4586444888104451       0.57872875077055630       -1.6382547880669396E-002  -8.6963046230628338E-002 -0.16152933282377679     
 /BM/amat           3  0.60221923055768034        1.6134454239326823        11.023120434364445        1.2228738740950016        4.6824626639950573E-002  0.50234026736392945       -1.7526847531036684     
 /BM/amat           4 -0.40385964045781408      -0.81499908712709757        4.1328610959912009        32.519888194250399        2.1532040996616875E-002  -7.6808137178456248E-002  -1.7883974537970195     
 /BM/amat           5  0.16687865786638639       0.72739977454891225       0.40387804931562910        6.5509996925184346       0.37905463361692071       0.75187668403058128       -1.2649764401583574     
 /BM/amat           6  -2.6115121822251702E-002  -6.9957976059004889E-002   2.3730111325394647       -3.9311789938343336      -0.96724801124187865        9.8827405775700115       -10.814675266457133     
 /BM/amat           7  -9.0313621913302139E-002  -8.9112100629948121E-002  -1.2662504631985390       0.60943823077972570      -0.95660846210460382        6.6444885346690601       -1.7940621597526167     
 /BM/bmat           1  -1.6057050485233840E+253  -1.8883480343298559E+253  -1.0242101328837837E+254  -8.8479228183629966E+253  -6.7394674523250949E+254  -1.8517886422631951E+255  -5.1083243111364564E+253
 /BM/bmat           2  -1.8883480343298559E+253  -2.3730428142479352E+254  -1.2864342221788863E+255  -4.9613784257235972E+254   1.0233320446856355E-003   2.6199820991876379E+254  -1.5164160847938008E+255
 /BM/bmat           3  -1.0242101328837837E+254  -1.2864342221788863E+255   2.1654156090729448E+254   3.4895760555665111E+255   1.1474005484302364E-002  -4.2707262782437682E+254  -2.3279796464890380E+254
 /BM/bmat           4  -8.8479228183629966E+253  -4.9613784257235972E+254   3.4895760555665111E+255  -6.4638378887935509E+254  -2.1963556771236840E+255   6.8785975475868169E+254  -5.2221969109911281E-003
 /BM/bmat           5  -6.7394674523250949E+254   1.0233320446856355E-003   1.1474005484302358E-002  -2.1963556771236840E+255   1.0086198231588820E+254  -3.0075054944264199E+253   1.4342978991614510E+254
 /BM/bmat           6  -1.8517886422631951E+255   2.6199820991876379E+254  -4.2707262782437682E+254   6.8785975475868169E+254  -3.0075054944264199E+253   6.2852576249349066E+254  -2.5340955314344414E+254
 /BM/bmat           7  -5.1083243111364564E+253  -1.5164160847938008E+255  -2.3279796464890380E+254  -5.2221969109911559E-003   1.4342978991614510E+254  -2.5340955314344414E+254   8.4098027035561773E+253

Solving generalized eigenvalue equation with a_diag =  1.0D-08
eigval_srt_ind_to_eigval_ind=     1     2     3     4     5     6     7     8     9    10    11    12    13    14
eigval_srt_ind_to_eigval_ind=     1     2     5     7     9    14    11    12    13    10     8     6     4     3
eigval_ind_to_eigval_srt_ind=     1     2    14    13     3    12     4    11     5    10     7     8     9     6
Sorted (complex) (unique) eigenvalues:   14
eigenvalue #    1: ************ +    0.000000 i (    1)
eigenvecEQU     1:    -1.000000      0.000000
eigenvecEQU     1:     0.163031     -0.000000
eigenvecEQU     1:     0.081298      0.000000
eigenvecEQU     1:    -0.062242      0.000000
eigenvecEQU     1:    -0.035941     -0.000000
eigenvecEQU     1:    -0.083072     -0.000000
eigenvecEQU     1:     0.040643      0.000000
eigenvalue #    2: ************ +    0.000000 i (    1)
eigenvecOPP     2:     1.000000     -0.000000
eigenvecOPP     2:    -0.890619      0.000000
eigenvecOPP     2:     0.255610      0.000000
eigenvecOPP     2:     0.048158     -0.000000
eigenvecOPP     2:    -0.024441     -0.000000
eigenvecOPP     2:    -0.136540     -0.000000
eigenvecOPP     2:     0.149585      0.000000
eigenvalue #    3: ************ +    0.000000 i (    1)
eigenvecOPP     3:     1.000000      0.000000
eigenvecOPP     3:     0.399912      0.000000
eigenvecOPP     3:    -0.233430     -0.000000
eigenvecOPP     3:     0.024250      0.000000
eigenvecOPP     3:    -0.089148      0.000000
eigenvecOPP     3:    -0.067601      0.000000
eigenvecOPP     3:    -0.009387      0.000000
eigenvalue #    4: ************ +    0.000000 i (    1)
eigenvecEQU     4:     0.376539      0.000000
eigenvecEQU     4:     1.000000     -0.000000
eigenvecEQU     4:    -0.288653      0.000000
eigenvecEQU     4:     0.291973      0.000000
eigenvecEQU     4:     0.041971     -0.000000
eigenvecEQU     4:    -0.036324     -0.000000
eigenvecEQU     4:     0.244568      0.000000
eigenvalue #    5: ************ +    0.000000 i (    1)
eigenvecEQU     5:    -1.000000     -0.000000
eigenvecEQU     5:    -0.707449      0.000000
eigenvecEQU     5:    -0.430034     -0.000000
eigenvecEQU     5:     0.409786     -0.000000
eigenvecEQU     5:     0.263671      0.000000
eigenvecEQU     5:    -0.275828      0.000000
eigenvecEQU     5:    -0.330276     -0.000000
eigenvalue #    6: ************ +    0.000000 i (    1)
eigenvecOPP     6:    -0.375119      0.000000
eigenvecOPP     6:    -0.846914     -0.000000
eigenvecOPP     6:    -0.204297      0.000000
eigenvecOPP     6:    -0.720359     -0.000000
eigenvecOPP     6:     1.000000     -0.000000
eigenvecOPP     6:    -0.323350      0.000000
eigenvecOPP     6:     0.667027     -0.000000
eigenvalue #    7: ************ +    0.000000 i (    1)
eigenvecOPP     7:     0.138985      0.000000
eigenvecOPP     7:    -0.128589      0.000000
eigenvecOPP     7:    -0.720549     -0.000000
eigenvecOPP     7:     0.042568      0.000000
eigenvecOPP     7:    -1.000000     -0.000000
eigenvecOPP     7:     0.394299      0.000000
eigenvecOPP     7:     0.730014     -0.000000
eigenvalue #    8: ************ +    0.000000 i (    1)
eigenvecEQU     8:     0.138985      0.000000
eigenvecEQU     8:    -0.128589      0.000000
eigenvecEQU     8:    -0.720549     -0.000000
eigenvecEQU     8:     0.042568      0.000000
eigenvecEQU     8:    -1.000000     -0.000000
eigenvecEQU     8:     0.394299      0.000000
eigenvecEQU     8:     0.730014      0.000000
eigenvalue #    9: ************ +    0.000000 i (    1)
eigenvecEQU     9:     0.375119     -0.000000
eigenvecEQU     9:     0.846914      0.000000
eigenvecEQU     9:     0.204297     -0.000000
eigenvecEQU     9:     0.720359      0.000000
eigenvecEQU     9:    -1.000000     -0.000000
eigenvecEQU     9:     0.323350     -0.000000
eigenvecEQU     9:    -0.667027      0.000000
eigenvalue #   10: ************ +    0.000000 i (    1)
eigenvecOPP    10:    -1.000000      0.000000
eigenvecOPP    10:    -0.707449      0.000000
eigenvecOPP    10:    -0.430034     -0.000000
eigenvecOPP    10:     0.409786     -0.000000
eigenvecOPP    10:     0.263671     -0.000000
eigenvecOPP    10:    -0.275828     -0.000000
eigenvecOPP    10:    -0.330276      0.000000
eigenvalue #   11: ************ +    0.000000 i (    1)
eigenvecOPP    11:    -0.376539      0.000000
eigenvecOPP    11:    -1.000000      0.000000
eigenvecOPP    11:     0.288653     -0.000000
eigenvecOPP    11:    -0.291973      0.000000
eigenvecOPP    11:    -0.041971     -0.000000
eigenvecOPP    11:     0.036324     -0.000000
eigenvecOPP    11:    -0.244568     -0.000000
eigenvalue #   12: ************ +    0.000000 i (    1)
eigenvecEQU    12:     1.000000     -0.000000
eigenvecEQU    12:     0.399912      0.000000
eigenvecEQU    12:    -0.233430      0.000000
eigenvecEQU    12:     0.024250     -0.000000
eigenvecEQU    12:    -0.089148     -0.000000
eigenvecEQU    12:    -0.067601     -0.000000
eigenvecEQU    12:    -0.009387      0.000000
eigenvalue #   13: ************ +    0.000000 i (    1)
eigenvecEQU    13:    -1.000000      0.000000
eigenvecEQU    13:     0.890619     -0.000000
eigenvecEQU    13:    -0.255610     -0.000000
eigenvecEQU    13:    -0.048158      0.000000
eigenvecEQU    13:     0.024441      0.000000
eigenvecEQU    13:     0.136540      0.000000
eigenvecEQU    13:    -0.149585     -0.000000
eigenvalue #   14: ************ +    0.000000 i (    1)
eigenvecOPP    14:     1.000000     -0.000000
eigenvecOPP    14:    -0.163031      0.000000
eigenvecOPP    14:    -0.081298     -0.000000
eigenvecOPP    14:     0.062242     -0.000000
eigenvecOPP    14:     0.035941      0.000000
eigenvecOPP    14:     0.083072      0.000000
eigenvecOPP    14:    -0.040643     -0.000000

Exit of menu.


Total CPU time is 0000h00m00s08c
The program ended normally.

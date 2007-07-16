(*compute cuthi and cutlo for include file parameters.h*)
dblmax=1.7976931348623157 10^308
dblmin=2.2250738585072014 10^(-308)
dblepsilon=2.2204460492503131 10^(-16)
cuthi=Sqrt[dblmax/dblepsilon]
cutlo=Sqrt[dblmin]
Print["dblmax ",dblmax]
Print["dblmin ",dblmin]
Print["dblepsilon ",dblepsilon]
Print["cuthi ",cuthi]
Print["cutlo ",cutlo]

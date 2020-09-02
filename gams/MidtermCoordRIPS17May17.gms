$TITLE RIPS HYDRO THERMAL COORD MODEL, 17 MAY 2017, MICHAEL CRAIG

*Turn off output in .lst file
$offlisting
$offsymxref offsymlist

Options
         optcr = 1E-3
         reslim = 1000
         limcol = 0
         limrow = 0
         threads = 0
         solprint = silent
*         solvelink = 5
         ;

Sets
         egu                                     "electricity generators"
         windegu(egu)                            "wind electricity generators"
         solaregu(egu)                           "solar electricity generators"
         hydroegu(egu)                           "hydroelectric generators"
         h                                       "hours"
         z                                       "demand zones"
         l                                       "lines between pairs of zones"
                 ;

alias(h,hh);

Parameters
*Unit-specific parameters
         pCapac(egu,h)                           "capacity of egu (GW)"
         pOpcost(egu,h)                          "plant operating cost (thousands$/GWh)"
*Max hourly generation for renewables
         pMaxgenwind(z,h)                          "maximum hourly generation by all wind generators"
         pMaxgensolar(z,h)                         "maximum hourly generation by all solar generators"
*Monthly max generation for hydro
         pMaxgenhydro(hydroegu)                  "maximum monthly hydro generation (MWh)"
*Zonal parameters
         pDemand(z,h)                              "electricity demand (GWh)"
         pEguzones(egu)                          "zone EGU is in
         pLinesources(l)                         "which zone line carries power from, i.e. source zone for line"
         pLinesinks(l)                           "which zone line carries power to, i.e. sink zone for line"
         pLinecapacs(l,h)                        "transmission limit per line (GW/hr)"
*Diagnostic parameters
         pModelstat
         pSolvestat
                 ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load egu, windegu, solaregu, hydroegu, h, z, l, pCapac, pOpcost
$load pMaxgenwind, pMaxgensolar, pMaxgenhydro, pDemand
$load pEguzones, pLinesources, pLinesinks, pLinecapacs
$gdxin

Variables
         vTotalopcost                            "total cost of power generation (thousands $)"
                 ;

Positive Variables
         vGen(egu,h)                             "power generation at plant egu at end of hour h (GW)"
         vLineflow(l,h)                          "flow over lines per hour (GW)"
                 ;

Equations
         objfunc                                 "define objective function to be minimized"
         meetdemand(z,h)                               "must meet electric demand"
         maxwindgen(z,h)                           "restrict wind generation to maximum aggregate output"
         maxsolargen(z,h)                          "restrict solar generation to maxmimum aggregate output"
         maxmonthlyhydrogen(hydroegu)
         genlimit(egu,h)                 "limit generation + spin reserves to max capac"
                  ;

******************OBJECTIVE FUNCTION******************
*Minimize total operational cost
objfunc .. vTotalopcost =e= sum((egu,h), vGen(egu,h)*pOpcost(egu,h));
******************************************************

******************ZONAL DEMAND AND RESERVE CONSTRAINTS******************
*Demand requirement
meetdemand(z,h).. sum(egu$[pEguzones(egu)=z],vGen(egu,h)) - sum(l$[pLinesources(l)=z],vLineflow(l,h))
         + sum(l$[pLinesinks(l)=z],vLineflow(l,h)) =e= pDemand(z,h);
***********************************************************

******************GENERATION CONSTRAINTS******************
*Enforce max generation on all wind generators
maxwindgen(z,h).. pMaxgenwind(z,h) =g= sum(windegu$[pEguzones(windegu)=z],vGen(windegu,h));

*Enforce max generation on all solar generators
maxsolargen(z,h).. pMaxgensolar(z,h) =g= sum(solaregu$[pEguzones(solaregu)=z],vGen(solaregu,h));

*Enforce max generation on monthly hydro generation
maxmonthlyhydrogen(hydroegu).. sum(h,vGen(hydroegu,h)) =l= pMaxhydrogen(hydroegu);

*Enforce maximum capacity on generation
genlimit(egu,h).. vGen(egu,h) =l= pCapac(egu,h);
**********************************************************

*****************LINE FLOW LIMITS*************************
*Limit max value of line flow to line capacity (already bounded at 0 since declared as positive variable)
vLineflow.up(l,h)=pLinecapacs(l,h);
**********************************************************

model ripsCoord /all/;
solve ripsCoord using mip minimizing vTotalopcost;

pModelstat = unitcommitment.Modelstat;
pSolvestat = unitcommitment.solvestat;

$TITLE RIPS UNIT COMMITMENT, 17 MAY 2017, MICHAEL CRAIG

*Turn off output in .lst file
$offlisting
$offsymxref offsymlist

Options
         optcr = 1E-3
         reslim = 1000
         limcol = 0
         limrow = 0
         solprint = off
         ;

Sets
         z                                       "demand zones"
         l                                       "lines between pairs of zones"
         h                                       "hours"
         egu                                     "electricity generators"
         windegu(egu)                            "wind electricity generators"
         solaregu(egu)                           "solar electricity generators"
         hydroegu(egu)                           "hydroelectric generators"
         pumphydroegu(egu)                       "Pumped Storage generators"
         ;

alias(h,hh);

Parameters
*Time-varying unit-specific parameters
         pCapac(egu,h)                           "capacity of egu (GW)"
*Unit-specific parameters
         pOpcost(egu)                            "plant operating cost (thousands$/GWh)"
         pMinload(egu)                           "minimum load of EGU (GW)"
         pRamprate(egu)                          "ramp rate of EGU, assumed to be the same up & down (GW/hr)"
         pStartupfixedcost(egu)                  "start up cost of EGU (thousands$)"
         pMindowntime(egu)                       "MDT (hours)"
         pRegeligible(egu)                       "whether eligible to provide reg reserves (1) or not (0)"
         pFlexeligible(egu)                      "whether eligible to provide flex reserves (1) or not (0)"
         pConteligible(egu)                      "whether eligible to provide contigency reserves (1) or not (0)"
*Max hourly generation for renewables
         pMaxgenwind(z,h)                          "maximum hourly generation by all wind generators"
         pMaxgensolar(z,h)                         "maximum hourly generation by all solar generators"
*Daily max generation for hydro
         pMaxgenhydro(hydroegu)                  "maximum daily hydro generation (MWh)"
*Carried unit-specific parameters from last optimization
         pOnoroffinitial(egu)                    "whether plant is set to on or off in hour 1 from prior optimization"
         pMdtcarriedhours(egu)                   "MDT carried over from last optimization (hours)"
         pGenabovemininitial(egu)                "gen above min load from last period of prior optimization (GW)"
*Zonal parameters
         pDemand(z,h)                            "electricity demand (GWh)"
         pRegupreserves(z,h)                     "required hourly up regulation reserves (GW)"
         pRegdownreserves(z,h)                   "required hourly down regulation reserves (GW)"
         pFlexreserves(z,h)
         pContreserves(z,h)
         pEguzones(egu)                          "zone EGU is in"
         pLinesources(l)                         "which zone line carries power from, i.e. source zone for line"
         pLinesinks(l)                           "which zone line carries power to, i.e. sink zone for line"
         pLinecapacs(l)                          "transmission limit per line (GW/hr)"
*Scalars
         pRampratetoregreservescalar             "converts timeframe that ramp rate is given in to reg reserve provision timeframe"
         pRampratetoflexreservescalar
         pRampratetocontreservescalar
         pMaxregupoffer(egu)
         pMaxflexoffer(egu)
         pMaxcontoffer(egu)
         pCnse                                   "cost of non-served energy (thousands$/GWh)"
         pCO2price                               "cost of CO2"
*Pumped Hydro Parameters
         pInitsoc(pumphydroegu)                  "initial state of charge of each pumped hydro unit [GWh]"
         pMaxsoc(pumphydroegu)                   "max state of charge of each pumped hydro unit [GWh]"
         pEfficiency(pumphydroegu)               "efficiency of each pumped hydro unit"
*Diagnostic parameters
         pModelstat
         pSolvestat
                 ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load z, l, h, egu, windegu, solaregu, hydroegu, pumphydroegu
$load pDemand, pEguzones, pCapac, pOpcost, pMinload, pRamprate, pStartupfixedcost, pMindowntime
$load pOnoroffinitial, pMdtcarriedhours, pGenabovemininitial, pRegeligible, pFlexeligible, pConteligible
$load pMaxgensolar, pMaxgenwind, pMaxgenhydro, pRampratetoregreservescalar, pRegupreserves, pRegdownreserves
$load pRampratetoflexreservescalar, pFlexreserves, pRampratetocontreservescalar, pContreserves
$load pCnse, pCO2price, pLinecapacs, pLinesources, pLinesinks
$load pInitsoc, pMaxsoc, pEfficiency
$gdxin

*DEFINE PARAMETERS
pMaxregupoffer(egu) = pRegeligible(egu)*pRamprate(egu)*pRampratetoregreservescalar;
pMaxflexoffer(egu) = pFlexeligible(egu)*pRamprate(egu)*pRampratetoflexreservescalar;
pMaxcontoffer(egu) = pConteligible(egu)*pRamprate(egu)*pRampratetocontreservescalar;

Variables
         vTotalopcost                            "total cost of power generation (thousands $)"
                 ;

Positive Variables
         vGen(egu,h)                             "power generation at plant egu at end of hour h (GW)"
         vGenabovemin(egu,h)                     "power generation above minimum stable load (GW)"
         vRegup(egu,h)                           "regulation up reserves provided (GW)"
         vFlex(egu,h)
         vCont(egu,h)
         vNse(z,h)                              "nonserved energy (GW)"
         vSoc(pumphydroegu,h)                    "state of charge for pumped hydro (GW)"
         vCharge(pumphydroegu,h)                 "charged energy for pumped hydro (GW)"
         vLineflow(l,h)                         "Line Flow in GW"
         ;

Binary Variables
         vTurnon(egu,h)                          "indicates whether plant decides to turn on (1) or not (0) in hour h"
         vTurnoff(egu,h)                         "indicates whether plant decides to turn off (1) or not (0) in hour h"
         vOnoroff(egu,h)                         "indicates whether plant is up (1) or down (0) in hour h"
                 ;

Equations
         objfunc                                 "define objective function to be minimized"
         meetdemand(z,h)                         "must meet electric demand"
         meetflexreserves(z,h)                   "meet hourly spinning reserve requirements"
         meetcontreserves(z,h)
         meetregupreserves(z,h)                  "meet hourly regulation reserve requirements"
         maxwindgen(z,h)                         "restrict wind generation to maximum aggregate output"
         maxsolargen(z,h)                        "restrict solar generation to maxmimum aggregate output"
         maxdailyhydrogen(hydroegu)
         definegenabovemin(egu,h)                "establish relationship between Gen (total gen) and Genabovemin (gen just above min stable load)"
         rampconstraintup(egu,h)                 "ramping up constraint for t>1"
         rampconstraintdown(egu,h)               "ramping down constraint for t>1"
         rampconstraintupinitial(egu,h)          "ramping up constraint at t=1"
         rampconstraintdowninitial(egu,h)        "ramping down constraint at t=1"
         statusofplant(egu,h)                    "balance whether thermal plant is on or off with whether shutting down or starting up"
         determineloadabovemin(egu,h)            "determine what each thermal unit's generation is above its minimum load. Constraints Genabovemin to be between max and min capacity"
         enforcemindowntime(egu,h)               "make sure plant, once it turns off, doesn't turn back on before MDT passes"
         enforcemindowntimecarryover(egu,h)      "enforce MDT from turn off decisions in prior optimization"
         flexreservelimit(egu,h)                 "limit spin reserves provided by generator to multiple of ramp rate"
         contreservelimit(egu,h)
         regupreservelimit(egu,h)                "limit reg reserves provided by generator to multiple of ramp rate"
         genplusresuplimit(egu,h)                "limit generation + spin reserves to max capac"
         genandsoc(pumphydroegu,h)               "restrict gen by pump hydro to state of charge"
         maxsto(pumphydroegu,h)                  "set max soc limit"
         limitcharging(pumphydroegu,h)           "limit charging to capacity times efficiency"
         defsoc(pumphydroegu,h)                  "link state of charge and charging and discharging"
         ;

* PRINT ZONES SET IN LST FILE IN ORDER TO CHECK IF ORDER IS CORRECT!
display z;

* PRINT h SET IN LST FILE IN ORDER TO CHECK IF ORDER IS CORRECT!
display h;

******************OBJECTIVE FUNCTION******************
*Minimize total operational cost
objfunc .. vTotalopcost =e= sum((z,h),pCnse*vNse(z,h))
                 + sum((egu,h), vGen(egu,h)*pOpcost(egu)+pStartupfixedcost(egu)*vTurnon(egu,h));
******************************************************

******************ZONAL DEMAND AND RESERVE CONSTRAINTS******************
*Demand requirement
meetdemand(z,h).. sum(egu$[pEguzones(egu)=ORD(z)],vGen(egu,h))
                  - sum(pumphydroegu$[pEguzones(pumphydroegu)=ORD(z)],vCharge(pumphydroegu,h))
                  - sum(l$[pLinesources(l)=ORD(z)],vLineflow(l,h))
                  + sum(l$[pLinesinks(l)=ORD(z)],vLineflow(l,h)) + vNse(z,h) =e= pDemand(z,h);

*Reserve requirements
meetflexreserves(z,h)  .. sum(egu$[pEguzones(egu)=ORD(z)],vFlex(egu,h)) =g= pFlexreserves(z,h);
meetcontreserves(z,h) .. sum(egu$[pEguzones(egu)=ORD(z)],vCont(egu,h)) =g= pContreserves(z,h);
meetregupreserves(z,h) .. sum(egu$[pEguzones(egu)=ORD(z)],vRegup(egu,h)) =g= pRegupreserves(z,h);
***********************************************************

******************GENERATION CONSTRAINTS******************
*Enforce max generation on all wind generators
maxwindgen(z,h).. pMaxgenwind(z,h) =g= sum(windegu$[pEguzones(windegu)=ORD(z)],vGen(windegu,h));

*Enforce max generation on all solar generators
maxsolargen(z,h).. pMaxgensolar(z,h) =g= sum(solaregu$[pEguzones(solaregu)=ORD(z)],vGen(solaregu,h));

*Enforce max generation on daily hydro generation
maxdailyhydrogen(hydroegu).. sum(h,vGen(hydroegu,h)) =l= pMaxgenhydro(hydroegu);

*Constrain plants to generate below their max capacity
definegenabovemin(egu,h).. vGen(egu,h) =e= vOnoroff(egu,h)*pMinload(egu)+vGenabovemin(egu,h);

*Establish relationship between gen above min load, gen output, and min load
determineloadabovemin(egu,h) .. vGenabovemin(egu,h) =l= (pCapac(egu,h)-pMinload(egu))*vOnoroff(egu,h);
**********************************************************

*****************LINE FLOW LIMITS*************************
*Limit max value of line flow to line capacity (already bounded at 0 since declared as positive variable)
vLineflow.up(l,h)=pLinecapacs(l);
**********************************************************

******************RESERVE PROVISION CONSTRAINTS******************
*Spin and reg reserves limited by ramp rate
regupreservelimit(egu,h)$[pMaxregupoffer(egu)>0] .. vRegup(egu,h) =l= pMaxregupoffer(egu)*vOnoroff(egu,h);
flexreservelimit(egu,h)$[pMaxflexoffer(egu)>0] .. vFlex(egu,h) =l= pMaxflexoffer(egu)*vOnoroff(egu,h);
contreservelimit(egu,h)$[pMaxcontoffer(egu)>0] .. vCont(egu,h) =l= pMaxcontoffer(egu)*vOnoroff(egu,h);

vRegup.fx(egu,h)$[pMaxregupoffer(egu)=0] = 0;
vFlex.fx(egu,h)$[pMaxflexoffer(egu)=0] = 0;
vCont.fx(egu,h)$[pMaxcontoffer(egu)=0] = 0;

*Limit generation + up (reg + spin) reserves to max capacity
genplusresuplimit(egu,h) .. vGen(egu,h) + vFlex(egu,h) + vCont(egu,h) + vRegup(egu,h) =l= pCapac(egu,h);
*****************************************************************

******************RAMPING CONSTRAINTS******************
*Ensure plants are limited to their ramping speed
rampconstraintup(egu,h)$[ORD(h)>1] .. (vGenabovemin(egu,h) + vFlex(egu,h) + vCont(egu,h) + vRegup(egu,h)) - vGenabovemin(egu,h-1) =l= pRamprate(egu);
rampconstraintdown(egu,h)$[ORD(h)>1] .. (vGenabovemin(egu,h-1) - vGenabovemin(egu,h)) =l= pRamprate(egu);

*Enforce ramp rates in hour 1 of optimization window
rampconstraintupinitial(egu,h)$[ORD(h)=1] .. (vGenabovemin(egu,h) + vFlex(egu,h) + vCont(egu,h) + vRegup(egu,h)) - pGenabovemininitial(egu) =l= pRamprate(egu);
rampconstraintdowninitial(egu,h)$[ORD(h)=1] .. (pGenabovemininitial(egu) - vGenabovemin(egu,h)) =l= pRamprate(egu);
*******************************************************

******************ON/OFF CONSTRAINTS******************
*Constrains status of plant per whether it's on/off, turning on, or shutting down
statusofplant(egu,h) .. vOnoroff(egu,h) =e= pOnoroffinitial(egu)$[ORD(h)=1]+vOnoroff(egu,h-1)$[ORD(h)>1]+vTurnon(egu,h)-vTurnoff(egu,h);

*Limit plant to not start up until it reaches its min down time
enforcemindowntime(egu,h)$[ORD(h)>pMdtcarriedhours(egu)] .. 1-vOnoroff(egu,h) =g= sum(hh$[ORD(hh)<=ORD(h) and ORD(hh)>(ORD(h)-pMindowntime(egu))],vTurnoff(egu,hh));

*Enforce MDT hours carried over from last optimization
enforcemindowntimecarryover(egu,h)$[ORD(h)<=pMdtcarriedhours(egu)] .. vOnoroff(egu,h) =l= 0;
******************************************************

******************STORAGE CONSTRAINTS******************
*Limit generation to state of charge
genandsoc(pumphydroegu,h).. vGen(pumphydroegu,h) =l= pInitsoc(pumphydroegu)$[ord(h)=1] +
                                                     vSoc(pumphydroegu,h-1);
*
*Link state of charge, charging, and discharging
defsoc(pumphydroegu,h).. vSoc(pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=1]
                                                  + vSoc(pumphydroegu,h-1)$[ord(h)>1]
                                                  - vGen(pumphydroegu,h)
                                                  + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,h);
*
*Limit state of charge to maximum storage capacity
maxsto(pumphydroegu,h).. vSoc(pumphydroegu,h) =l= pMaxsoc(pumphydroegu);
*
*Limit rate of charging to capacity times efficiency
limitcharging(pumphydroegu,h).. vCharge(pumphydroegu,h) =l= pCapac(pumphydroegu,h);
******************************************************

model ripsUC /all/;
solve ripsUC using mip minimizing vTotalopcost;

pModelstat = ripsUC.Modelstat;
pSolvestat = ripsUC.solvestat;

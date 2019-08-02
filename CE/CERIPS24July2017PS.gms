$TITLE CAPACITY EXPANSION WITH CHRONOLOGICAL DEMAND THROUGH PYTHON API, 27 MAY 2016, MICHAEL CRAIG

*$offsymxref offsymlist

Options
         optcr = 5E-2
         reslim = 80000
         limcol = 0
         limrow = 0
         solprint = off
         ;

Sets
         egu                             existing generators
         windegu(egu)                    existing wind generators
         solaregu(egu)                   existing solar generators
         hydroegu(egu)                   existing hydro generators
         pumphydroegu(egu)               existing pumped hydro generators
*         type                            plant types (no cooling info)
         tech                            candidate plant types for new construction (with cooling info)
*         tech2d(type, tech)              2d set mapping types to techs
         techcurtailed(tech)             plant types that can be curtailed for new construction
         techrenew(tech)                 renewable plant types for new construction
         technotcurtailed(tech)          plant types that are not curtailed and not renewables for new construction
         c                               cells that new techs can be placed in
         g                               gcm scenarios
         h                               hours (1d set. over all gcm scenarios)
         h2(g,h)                         hours in each gcm scenario
         dispatchh(g,h)                  hours that are not peak hour
         springh(g,h)                    hours representing spring
         summerh(g,h)                    hours representing summer
         winterh(g,h)                    hours representing winter
         fallh(g,h)                      hours representing fall
         specialh(g,h)                   hours representing special periods
         z                               zones
         l                               lines
         peakh(g,z,h)                    hours for each zone with peak net demand and respective gcm
         ;

Parameters
*SIZE PARAMETERS [GW]
         pCapac(g,egu,h)                                hourly capacity of existing generators accounting for curtailments [GW]
         pCapactech(tech)                               nameplate capacity of new builds for cost calculations [GW]
         pCapactechcurtailed(g,c,techcurtailed,h)       hourly capacity of new builds accounting for curtailments [GW]
*TIME- AND SPACE-VARYING OP COSTS AND HEAT RATES [MMBtu/GWh]
         pOpcost(egu)                                   total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
         pOpcosttechcurt(c, techcurtailed)              total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
         pOpcosttechnotcurt(z, technotcurtailed)        total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
         pOpcosttechrenew(z, techrenew)                 total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
*NEW TECH FIXED COST PARAMETERS
         pFom(tech)                      fixed O&M cost [thousand USD per GW per yr]
         pOcc(tech)                      overnight capital cost [thousand USD per GW]
*EMISSIONS RATES [short ton/GWh]
         pCO2emrate(egu)                                CO2 emissions rate of existing generators [short tons per GWh]
         pCO2emratetechcurt(c, techcurtailed)           CO2 emissions rate of potential new generators [short tons per GWh]
         pCO2emratetechnotcurt(z, technotcurtailed)     CO2 emissions rate of potential new generators [short tons per GWh]
         pCO2emratetechrenew(z, techrenew)              CO2 emissions rate of potential new generators [short tons per GWh]
*EMISSIONS CAP AND COST
         pCO2emcap                       CO2 annual emissions cap [10^3 short tons]
*HOURLY CAPACITY FACTORS FOR RENEWABLES
         pCf(g,z,techrenew,h)              hourly capacity factors for potential new renewables in each zone
         pMaxgenwind(g,z,h)                max hourly generation for existing wind [GWh]
         pMaxgensolar(g,z,h)               max hourly generation for existing solar [GWh]
*MAX GENERATION VALUES FOR HYDRO POWER PER SEASON
         pMaxhydrogenspr(g,hydroegu)       max generation by each hydro unit in spring [GWh]
         pMaxhydrogensum(g,hydroegu)       max generation by each hydro unit in summer [GWh]
         pMaxhydrogenwin(g,hydroegu)       max generation by each hydro unit in winter [GWh]
         pMaxhydrogenfal(g,hydroegu)       max generation by each hydro unit in fall [GWh]
         pMaxhydrogenspe(g,hydroegu)       max generation by each hydro unit in special hours [GWh]
*INITIAL STATE OF CHARGE FOR PUMPED HYDRO
         pInitsoc(pumphydroegu)          initial state of charge of each pumped hydro unit [GWh]
         pMaxsoc(pumphydroegu)           max state of charge of each pumped hydro unit [GWh]
         pEfficiency(pumphydroegu)       efficiency of each pumped hydro unit
*FINANCIAL PARAMETERS
         pR                              discount rate
         pLife(tech)                     lifetime of tech [years]
         pCrf(tech)                      capital recovery factor
*BUILD LIMITS ON NEW TECHS
*         pNmax(z,type)                   max number techs built per zone and plant type (no cooling info)
         pNmaxReBlock(z, techrenew)      max number RE techs built per aggregated block and zone
*ZONAL PARAMETERS
         pDemand(g,z,h)                  hourly electricity demand [GWh]
         pPlanningreserve(z)             planning margin reserve capacity [GW]
*         pPeakhtozone(peakh)             peakh for each zone
         pEguzones(egu)                  zone EGU is in
         pCellzones(c)                   zone that each cell eligible for new tech builds is in
         pLinesources(l)                 which zone line carries power from
         pLinesinks(l)                   which zone line carries power to
         pLinecapacs(l)                  transmission limit per line (GW)
*WEIGHT TO SCALE UP VAR COSTS AND EMISSIONS FROM REPRESENTATIVE SEASONAL HOURS TO ENTIRE SEASON
         pWeightspring                   weight for spring hours
         pWeightsummer                   weight for summer hours
         pWeightfall                     weight for fall hours
         pWeightwinter                   weight for winter hours
*DIAGNOSTIC PARAMETERS
         pModelstat                      model status whether optimal solution achieved
         pSolvestat                      solver status whether terminated normally
*OTHER PARAMETERS
         pNgcm                           number of gcms considered
         pHourIniSpring(g)                           number of gcms considered
         pHourIniSummer(g)                           number of gcms considered
         pHourIniFall(g)                           number of gcms considered
         pHourIniWinter(g)                           number of gcms considered
         pHourIniSpecial(g)                           number of gcms considered
         pPreviousHourSpecial(g, h)
         ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
*******$load egu, windegu, solaregu, hydroegu, pumphydroegu, type, tech, tech2d, techcurtailed, techrenew, technotcurtailed, c, g
$load egu, windegu, solaregu, hydroegu, pumphydroegu, tech, techcurtailed, techrenew, technotcurtailed, c, g
$load h, h2, springh, summerh, winterh, fallh, specialh, z, l, peakh
$load pCapac, pCapactech, pCapactechcurtailed, pOpcost, pOpcosttechcurt, pOpcosttechnotcurt, pOpcosttechrenew
$load pFom, pOcc, pCO2emrate, pCO2emratetechcurt, pCO2emratetechnotcurt, pCO2emratetechrenew, pCO2emcap, pCf
$load pMaxgenwind, pMaxgensolar, pMaxhydrogensum, pMaxhydrogenspr, pMaxhydrogenwin, pMaxhydrogenfal, pMaxhydrogenspe
$load pDemand, pEguzones, pCellzones, pLinesources, pLinesinks, pLinecapacs
$load pR, pLife, pNmaxReBlock, pPlanningreserve, pWeightspring, pWeightsummer, pWeightfall, pWeightwinter
$load pInitsoc, pMaxsoc, pEfficiency
$gdxin

*pPeakhtozone,

*DEFINE NON-PEAK HOURS THAT WILL BE INCLUDED IN DEMAND=SUPPLY EQUATION
*dispatchh(g,h) = winterh(h) + summerh(h) + specialh(h) + fallh(h) + springh(h);

*CALCULATE CAPITAL RECOVERY FACTOR
pCrf(tech) = pR / (1 - (1 / ( (1 + pR)**pLife(tech))));

*GET NUMBER OF GCMS USED IN THE SIMULATION
pNgcm = card(g);

*GET POSITION OF INITIAL HOUR OF EACH SEASON IN EACH GCM (NECESSARY FOR PUMPED HYDRO SIMULATION)
pHourIniSpring(g) = smin(h$springh(g,h), ord(h));
pHourIniSummer(g) = smin(h$summerh(g,h), ord(h));
pHourIniWinter(g) = smin(h$winterh(g,h), ord(h));
pHourIniFall(g) = smin(h$fallh(g,h), ord(h));
pHourIniSpecial(g) = smin(h$specialh(g,h), ord(h));

* CREATE ANOTHER NAME ('Alias') FOR SET h
Alias(h, hh);

* FOR EACH OF THE SPECIAL HOURS, GET THE PREVIOUS HOUR (NECESSARY FOR PUMPED HYDRO SIMULATION)
pPreviousHourSpecial(g, h) = smax(hh$[specialh(g,hh) and ord(hh) < ord(h)], ord(hh));

Variable
         vZ                              obj func [billion $ per yr]
         vIc                             total investment costs for new plants = fixed O&M + overnight capital costs [thousand USD per yr]
         vVc                             total variable costs for new and existing plants = variable O&M + fuel + emission costs [thousand USD per yr]
         vVcspring                       variable costs for spring hours
         vVcsummer                       variable costs for summer hours
         vVcfall                         variable costs for fall hours
         vVcwinter                       variable costs for winter hours
         vVcspecial                      variable costs for special hours
         ;

Positive variables
         vPtechcurtailed(g,c,techcurtailed,h)                      hourly electricity generation by new techs that can be curtailed by cell c and hour h [GWh]
         vPtechrenew(g,z,techrenew,h)                              hourly electricity generation by new renewable techs [GWh]
         vPtechnotcurtailed(g,z,technotcurtailed,h)                hourly electricity generation by new plants [GWh]
         vPegu(g,egu,h)                                            hourly electricity generation by existing plants [GWh]
         vLineflow(g,l,h)                                          flow over lines per hour (GW)
         vCO2emsannual                                           co2 emissions in entire year from new and existing plants [10^3 short ton]
         vCO2emssummer                                           co2 emissions in summer from new and existing plants [10^3 short ton]
         vCO2emsspring                                           co2 emissions in spring from new and existing plants [10^3 short ton]
         vCO2emswinter                                           co2 emissions in winter from new and existing plants [10^3 short ton]
         vCO2emsfall                                             co2 emissions in fall from new and existing plants [10^3 short ton]
         vCO2emsspecial                                          co2 emissions in special hours from new and existing plants [10^3 short ton]
         vSoc(g,pumphydroegu,h)                                    state of charge of pumped hydro units [GWh]
         vCharge(g,pumphydroegu,h)                                 charging by pumped hydro units [GWh]
         ;

Integer variable
         vNrenew(z,techrenew)                                    number of newly constructed renewable plants
         vNnotcurtailed(z,technotcurtailed)                      number of newly constructed plants of types that cant be curtailed in zone z
         vNcurtailed(c,techcurtailed)                            number of newly constructed plants of types that can be curtailed in cell c
         ;

Equations
         objfunc                                                 objective function = sum investment and variable costs
         investmentcost                                          calculate investment costs = fixed O&M + annualized capital costs
         varcosttotal                                            sum annual variable costs
         varcostspring                                           calculate spring variable costs
         varcostsummer                                           calculate summer variable costs
         varcostwinter                                           calculate winter variable costs
         varcostfall                                             calculate fall variable costs
         varcostspecial                                          calculate special hour variable costs
         meetdemand(g,z,h)                                       meet supply with demand
         meetreservemargin(z)                                    meet planning reserve requirement with installed capacity
         curtailedtechgen(g,c,techcurtailed,h)                   restrict electricity generation by new curtailed plants to number built and hourly capacities
         notcurtailedtechgen(g,z,technotcurtailed,h)             restrict elec gen by not curtailed plants to number built and constant capacity
         renewtechgen(g,z,techrenew,h)                           restrict electricity generation by new renewables to number built and capacity and capacity factor
*         maxzonalbuild(z,type)                                   max number plant types built per zone
         egugen(g,egu,h)                                         restrict electricity generation by existing generators to hourly capacities
         eguwindgen(g,z,h)                                       restrict electricity generation by existing wind generation to maximum aggregate output per zone
         egusolargen(g,z,h)                                      restrict electricity generation by existing solar generation to maximum aggregate output per zone
         hydrogenspr(g,hydroegu)                                   restrict total max elec gen by hydro units in spring
         hydrogensum(g,hydroegu)                                   restrict total max elec gen by hydro units in summer
         hydrogenwin(g,hydroegu)                                   restrict total max elec gen by hydro units in winter
         hydrogenfal(g,hydroegu)                                   restrict total max elec gen by hydro units in fall
         hydrogenspe(g,hydroegu)                                   restrict total max elec gen by hydro units in special hours
         co2emsannual                                            sum annual co2 emissions
         co2emssummer                                            calculate summer co2 emissions
         co2emsspring                                            calculate spring co2 emissions
         co2emswinter                                            calculate winter co2 emissions
         co2emsfall                                              calculate fall co2 emissions
         co2emsspecial                                           calculate special hour co2 emissions
         enforceco2emissionscap                                  restrict total co2 emissions to cap
         genandsoc(g,pumphydroegu,h)                             restrict gen by pump hydro to state of charge
         maxsto(g,pumphydroegu,h)                                set max soc limit
         limitcharging(g,pumphydroegu,h)                         limit charging to capacity times efficiency
         defsocspr(g,pumphydroegu,h)                         link state of charge and charging and discharging
         defsocsum(g,pumphydroegu,h)                         link state of charge and charging and discharging
         defsocfal(g,pumphydroegu,h)                           link state of charge and charging and discharging
         defsocwin(g,pumphydroegu,h)                         link state of charge and charging and discharging
         defsocspe(g,pumphydroegu,h)                        link state of charge and charging and discharging
         ;

* PRINT ZONES SET IN LST FILE IN ORDER TO CHECK IF ORDER IS CORRECT!
display z;

* PRINT NUMBER OF GCMS SCENARIOS CONSIDERED IN LST FILE
display pNgcm;


******************PUT FILE******************
* defines a put file to write some of the model parameters
* File file_name "this defines a specific external file"  / /Users/kiko/report.txt /;
* put file_name;
******************************************************

******************OBJECTIVE FUNCTION******************
*Objective: minimize fixed + variable costs (vIc and vVc are in thousands $ per year. vZ is in billion $ per year.)
objfunc..                vZ =e= vIc/1e6 + vVc/1e6;
******************************************************

******************CALCULATE COSTS******************
*Fixed costs = annual fixed O&M + fixed annualized capital costs
investmentcost..         vIc =e= sum((c,techcurtailed),vNcurtailed(c,techcurtailed)*pCapactech(techcurtailed)*(pFom(techcurtailed)+pOcc(techcurtailed)*pCrf(techcurtailed)))
                               + sum((z,techrenew),vNrenew(z,techrenew)*pCapactech(techrenew)*(pFom(techrenew)+pOcc(techrenew)*pCrf(techrenew)))
                               + sum((z,technotcurtailed),vNnotcurtailed(z,technotcurtailed)*pCapactech(technotcurtailed)*(pFom(technotcurtailed)+pOcc(technotcurtailed)*pCrf(technotcurtailed)));

*Variable costs = electricity generation costs by new and existing plants
varcosttotal..      vVc =e= vVcspring + vVcsummer + vVcwinter + vVcfall + vVcspecial;

varcostspring..     vVcspring =e= (pWeightspring/pNgcm)*sum((g,h)$[springh(g,h)],sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pOpcosttechcurt(c, techcurtailed))
                                                                                 + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pOpcosttechnotcurt(z, technotcurtailed))
                                                                                 + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pOpcosttechrenew(z, techrenew))
                                                                                 + sum(egu,vPegu(g,egu,h)*pOpcost(egu)));

varcostsummer..     vVcsummer =e= (pWeightsummer/pNgcm)*sum((g,h)$[summerh(g,h)],sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pOpcosttechcurt(c, techcurtailed))
                                                                                 + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pOpcosttechnotcurt(z, technotcurtailed))
                                                                                 + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pOpcosttechrenew(z, techrenew))
                                                                                 + sum(egu,vPegu(g,egu,h)*pOpcost(egu)));

varcostwinter..     vVcwinter =e= (pWeightwinter/pNgcm)*sum((g,h)$[winterh(g,h)],sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pOpcosttechcurt(c, techcurtailed))
                                                                                 + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pOpcosttechnotcurt(z, technotcurtailed))
                                                                                 + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pOpcosttechrenew(z, techrenew))
                                                                                 + sum(egu,vPegu(g,egu,h)*pOpcost(egu)));

varcostfall..       vVcfall =e= (pWeightfall/pNgcm)*sum((g,h)$[fallh(g,h)],sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pOpcosttechcurt(c, techcurtailed))
                                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pOpcosttechnotcurt(z, technotcurtailed))
                                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pOpcosttechrenew(z, techrenew))
                                                                           + sum((egu),vPegu(g,egu,h)*pOpcost(egu)));

varcostspecial..    vVcspecial =e= (1/pNgcm)*sum((g,h)$[specialh(g,h)],sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pOpcosttechcurt(c, techcurtailed))
                                                                        + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pOpcosttechnotcurt(z, technotcurtailed))
                                                                        + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pOpcosttechrenew(z, techrenew))
                                                                        + sum(egu,vPegu(g,egu,h)*pOpcost(egu)));
***************************************************

******************SYSTEM-WIDE GENERATION AND CAPACITY CONSTRAINTS******************
*Demand = generation by new and existing plants in each zone
meetdemand(g,z,h) $ h2(g,h)..   sum((techcurtailed,c)$[pCellzones(c)=ORD(z)],vPtechcurtailed(g,c,techcurtailed,h))
                                + sum(technotcurtailed,vPtechnotcurtailed(g,z,technotcurtailed,h))
                                + sum(techrenew,vPtechrenew(g,z,techrenew,h))
                                + sum(egu$[pEguzones(egu)=ORD(z)],vPegu(g,egu,h))
                                - sum(l$[pLinesources(l)=ORD(z)],vLineflow(g,l,h))
                                + sum(l$[pLinesinks(l)=ORD(z)],vLineflow(g,l,h))
                                - sum(pumphydroegu$[pEguzones(pumphydroegu)=ORD(z)],vCharge(g,pumphydroegu,h)) =e= pDemand(g,z,h);

*Total installed capacity must exceed peak demand + planning reserve margin in each zone
meetreservemargin(z)..  pPlanningreserve(z) =l= sum((g,h)$[peakh(g,z,h)], sum(egu$[pEguzones(egu)=ORD(z)],pCapac(g,egu,h))
                                                                         + sum((techcurtailed,c)$[pCellzones(c)=ORD(z)],pCapactechcurtailed(g,c,techcurtailed,h)*vNcurtailed(c,techcurtailed))
                                                                         + sum(technotcurtailed,pCapactech(technotcurtailed)*vNnotcurtailed(z,technotcurtailed))
                                                                         + sum(techrenew,pCapactech(techrenew)*vNrenew(z,techrenew)*pCf(g,z,techrenew,h)));
***********************************************************************************

*****************LINE FLOW LIMITS*************************
*Limit max value of line flow to line capacity (already bounded below at 0 since declared as positive variable)
vLineflow.up(g,l,h) $ h2(g,h) = pLinecapacs(l);
**********************************************************

******************GENERATION CONSTRAINTS ON NEW UNITS******************
*Generation by techs that can be curtailed is limited by curtailed capacity
curtailedtechgen(g,c,techcurtailed,h)$h2(g,h)..          vPtechcurtailed(g,c,techcurtailed,h) =l= vNcurtailed(c,techcurtailed)*pCapactechcurtailed(g,c,techcurtailed,h);

*Generation by techs that cant be curtailed is limited by their capacity
notcurtailedtechgen(g,z,technotcurtailed,h)$h2(g,h)..    vPtechnotcurtailed(g,z,technotcurtailed,h) =l= vNnotcurtailed(z,technotcurtailed)*pCapactech(technotcurtailed);

*Renewable generation limited by capacity factor and nameplate capacity
renewtechgen(g,z,techrenew,h)$h2(g,h)..                  vPtechrenew(g,z,techrenew,h) =l= vNrenew(z,techrenew)*pCapactech(techrenew)*pCf(g,z,techrenew,h);
***********************************************************************

******************BUILD DECISIONS******************
* Limit total number builds to input value
* build limits are imposed by zone and plant type. Plant types do not take into account cooling type, just fuel source.
* So a upper bound of 10 on type 'Coal Steam', mean that the model can build [10 Coal/OT] or [5 Coal/OT + 5 Coal/RC] or ...
* Syntax of the right hand side conditional assignment:
* eq.. (expression) $ (logical condition) <======> if(logical condition) then {expression}
*
*maxzonalbuild(z,type)..  sum(c$[pCellzones(c)=ORD(z)], sum(techcurtailed$[tech2d(type,techcurtailed)], vNcurtailed(c,techcurtailed))) $ (sum(techcurtailed$[tech2d(type,techcurtailed)],1) > 0) +
*                         sum(techrenew$[tech2d(type,techrenew)], vNrenew(z,techrenew)) $ (sum(techrenew$[tech2d(type,techrenew)],1) > 0) +
*                         sum(technotcurtailed$[tech2d(type,technotcurtailed)], vNnotcurtailed(z,technotcurtailed)) $ (sum(technotcurtailed$[tech2d(type,technotcurtailed)],1) > 0) =l=
*                         pNmax(z,type);


* Fix an upper bound on integer build variables (otherwise GAMS automatically sets it to 100)
* (note that the constraint above will still be applied)
vNnotcurtailed.up(z,technotcurtailed) = 5000;
vNcurtailed.up(c,techcurtailed) = 10;

*Define upper bound for each block of renewable (wind or solar)
vNrenew.up(z,techrenew) = pNmaxReBlock(z, techrenew);
*vNrenew.up(z,techrenew) = 5000;
***************************************************

******************GENERATION CONSTRAINTS ON EXISTING UNITS******************
*Enforce max hourly generation limit based on hourly capacity
egugen(g,egu,h)$h2(g,h)..          vPegu(g,egu,h) =l= pCapac(g,egu,h);

*Enforce max generation on existing wind and solar plants
eguwindgen(g,z,h)$h2(g,h)..        pMaxgenwind(g,z,h) =g= sum(windegu$[pEguzones(windegu)=ORD(z)],vPegu(g,windegu,h));
egusolargen(g,z,h)$h2(g,h)..       pMaxgensolar(g,z,h) =g= sum(solaregu$[pEguzones(solaregu)=ORD(z)],vPegu(g,solaregu,h));

*Enforce max generation on total hydro generation
hydrogenspr(g,hydroegu)..   pMaxhydrogenspr(g,hydroegu) =g= sum(h$springh(g,h),vPegu(g,hydroegu,h));
hydrogensum(g,hydroegu)..   pMaxhydrogensum(g,hydroegu) =g= sum(h$summerh(g,h),vPegu(g,hydroegu,h));
hydrogenwin(g,hydroegu)..   pMaxhydrogenwin(g,hydroegu) =g= sum(h$winterh(g,h),vPegu(g,hydroegu,h));
hydrogenfal(g,hydroegu)..   pMaxhydrogenfal(g,hydroegu) =g= sum(h$fallh(g,h),vPegu(g,hydroegu,h));
hydrogenspe(g,hydroegu)..   pMaxhydrogenspe(g,hydroegu) =g= sum(h$specialh(g,h),vPegu(g,hydroegu,h));
****************************************************************************

******************CO2 EMISSIONS CONSTRAINT******************
*Co2 emissions = electricity generation * co2 emissions rate
*(each seasonal component is in thousands tons. Convert vCO2emsannual to million short tons)
co2emsannual..   vCO2emsannual =e= 1e-3 * (vCO2emsspring + vCO2emssummer + vCO2emswinter + vCO2emsfall + vCO2emsspecial);

co2emsspring..   vCO2emsspring =e= (pWeightspring/pNgcm)*sum((g,h)$[springh(g,h)],sum(egu,vPegu(g,egu,h)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pCO2emratetechcurt(c, techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pCO2emratetechnotcurt(z, technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pCO2emratetechrenew(z, techrenew)));

co2emssummer..   vCO2emssummer =e= (pWeightsummer/pNgcm)*sum((g,h)$[summerh(g,h)],sum(egu,vPegu(g,egu,h)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pCO2emratetechcurt(c, techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pCO2emratetechnotcurt(z, technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pCO2emratetechrenew(z, techrenew)));

co2emswinter..   vCO2emswinter =e= (pWeightwinter/pNgcm)*sum((g,h)$[winterh(g,h)],sum(egu,vPegu(g,egu,h)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pCO2emratetechcurt(c, techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pCO2emratetechnotcurt(z, technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pCO2emratetechrenew(z, techrenew)));

co2emsfall..     vCO2emsfall =e= (pWeightfall/pNgcm)*sum((g,h)$[fallh(g,h)],sum(egu,vPegu(g,egu,h)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pCO2emratetechcurt(c, techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pCO2emratetechnotcurt(z, technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pCO2emratetechrenew(z, techrenew)));

co2emsspecial..  vCO2emsspecial =e= (1/pNgcm)*sum((g,h)$[specialh(g,h)],sum(egu,vPegu(g,egu,h)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(g,c,techcurtailed,h)*pCO2emratetechcurt(c, techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(g,z,technotcurtailed,h)*pCO2emratetechnotcurt(z, technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(g,z,techrenew,h)*pCO2emratetechrenew(z, techrenew)));

*Meet emissions cap (change RHS to million short tons)
enforceco2emissionscap.. vCO2emsannual =l= pCO2emcap/1e3;
************************************************************

******************STORAGE CONSTRAINTS******************
*Limit generation to state of charge
genandsoc(g,pumphydroegu,h)$h2(g,h).. vPegu(g,pumphydroegu,h) =l= vSoc(g,pumphydroegu,h);

*Link state of charge, charging, and discharging
defsocspr(g,pumphydroegu,h)$springh(g,h).. vSoc(g,pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=pHourIniSpring(g)]
                                                                      + vSoc(g,pumphydroegu,h-1)$[ord(h)>pHourIniSpring(g)]
                                                                      - vPegu(g,pumphydroegu,h)
                                                                      + pEfficiency(pumphydroegu) * vCharge(g,pumphydroegu,h);

defsocsum(g,pumphydroegu,h)$summerh(g,h)..   vSoc(g,pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=pHourIniSummer(g)]
                                                                      + vSoc(g,pumphydroegu,h-1)$[ord(h)>pHourIniSummer(g)]
                                                                      - vPegu(g,pumphydroegu,h)
                                                                      + pEfficiency(pumphydroegu) * vCharge(g,pumphydroegu,h);

defsocfal(g,pumphydroegu,h)$fallh(g,h)..     vSoc(g,pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=pHourIniFall(g)]
                                                                      + vSoc(g,pumphydroegu,h-1)$[ord(h)>pHourIniFall(g)]
                                                                      - vPegu(g,pumphydroegu,h)
                                                                      + pEfficiency(pumphydroegu) * vCharge(g,pumphydroegu,h);

defsocwin(g,pumphydroegu,h)$winterh(g,h)..   vSoc(g,pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=pHourIniWinter(g)]
                                                                      + vSoc(g,pumphydroegu,h-1)$[ord(h)>pHourIniWinter(g)]
                                                                      - vPegu(g,pumphydroegu,h)
                                                                      + pEfficiency(pumphydroegu) * vCharge(g,pumphydroegu,h);

defsocspe(g,pumphydroegu,h)$specialh(g,h)..  vSoc(g,pumphydroegu,h) =e= pInitsoc(pumphydroegu)$[ord(h)=pHourIniSpecial(g)]
                                                                      + sum(hh$[ord(hh)=pPreviousHourSpecial(g, h)], vSoc(g,pumphydroegu,hh))
                                                                      - vPegu(g,pumphydroegu,h)
                                                                      + pEfficiency(pumphydroegu) * vCharge(g,pumphydroegu,h);

*Limit state of charge to maximum storage capacity
maxsto(g,pumphydroegu,h)$h2(g,h) .. vSoc(g,pumphydroegu,h) =l= pMaxsoc(pumphydroegu);

*Limit rate of charging to capacity times efficiency
limitcharging(g,pumphydroegu,h)$h2(g,h) .. vCharge(g,pumphydroegu,h) =l= pCapac(g,pumphydroegu,h);
***************************************************

Model expansion includes all equations /all/;
solve expansion using mip minimizing vZ;

pModelstat = expansion.Modelstat;
pSolvestat = expansion.solvestat;

$TITLE CAPACITY EXPANSION WITH CHRONOLOGICAL DEMAND THROUGH PYTHON API, 27 MAY 2016, MICHAEL CRAIG

*$offlisting
*$offsymxref offsymlist

Options
         optcr = 5E-1
         reslim = 20000
         limcol = 0
         limrow = 0
*         solprint = silent
         ;

Sets
         egu                             existing generators
         windegu(egu)                    existing wind generators
         solaregu(egu)                   existing solar generators
         hydroegu(egu)                   existing hydro generators
         pumphydroegu(egu)               existing pumped hydro generators
         tech                            candidate plant types for new construction
         techcurtailed(tech)             plant types that can be curtailed for new construction
         techrenew(tech)                 renewable plant types for new construction
         technotcurtailed(tech)          plant types that are not curtailed and not renewables for new construction
         c                               cells that new techs can be placed in
         h                               hours
         dispatchh(h)                    hours that are not peak hour
         peakh(h)                        hours for each zone with peak net demand
         springh(h)                      hours representing spring
         summerh(h)                      hours representing summer
         winterh(h)                      hours representing winter
         fallh(h)                        hours representing fall
         specialh(h)                     hours representing special periods
         z                               zones
         l                               lines
         ;

Parameters
*SIZE PARAMETERS [GW]
         pCapac(egu,h)                   hourly capacity of existing generators accounting for curtailments [GW]
         pCapactech(tech)                nameplate capacity of new builds for cost calculations [GW]
         pCapactechcurtailed(c,techcurtailed,h)   hourly capacity of new builds accounting for curtailments [GW]
*TIME- AND SPACE-VARYING OP COSTS AND HEAT RATES [MMBtu/GWh]
         pHr(egu)                        heat rate of existing generators [MMBtu per GWh]
         pHrtech(tech)                   heat rate of new builds [MMBtu per GWh]
         pOpcost(egu)                    total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
         pOpcosttech(tech)               total operational cost [thousand USD per GWh] = VOM + FuelCost*HR + EmsCost*EmsRate*HR
*NEW TECH FIXED COST PARAMETERS
         pFom(tech)                      fixed O&M cost [thousand USD per GW per yr]
         pOcc(tech)                      overnight capital cost [thousand USD per GW]
*EMISSIONS RATES [short ton/MMBtu]
         pCO2emrate(egu)                 CO2 emissions rate of existing generators [short ton per MMBtu]
         pCO2emratetech(tech)            CO2 emissions rate of potential new generators [short ton per MMBtu]
*EMISSIONS CAP AND COST
         pCO2emcap                       CO2 annual emissions cap [short tons]
*HOURLY CAPACITY FACTORS FOR RENEWABLES
         pCf(z,techrenew,h)              hourly capacity factors for potential new renewables in each zone
         pMaxgenwind(z,h)                max hourly generation for existing wind [GWh]
         pMaxgensolar(z,h)               max hourly generation for existing solar [GWh]
*MAX GENERATION VALUES FOR HYDRO POWER PER SEASON
         pMaxhydrogenspr(hydroegu)       max generation by each hydro unit in spring [GWh]
         pMaxhydrogensum(hydroegu)       max generation by each hydro unit in summer [GWh]
         pMaxhydrogenwin(hydroegu)       max generation by each hydro unit in winter [GWh]
         pMaxhydrogenfal(hydroegu)       max generation by each hydro unit in fall [GWh]
         pMaxhydrogenspe(hydroegu)       max generation by each hydro unit in special hours [GWh]
*INITIAL STATE OF CHARGE FOR PUMPED HYDRO
         pInitsoc(pumphydroegu)          initial state of charge of each pumped hydro unit [GWh]
         pMaxsoc(pumphydroegu)           max state of charge of each pumped hydro unit [GWh]
         pEfficiency(pumphydroegu)       efficiency of each pumped hydro unit
*FINANCIAL PARAMETERS
         pR                              discount rate
         pLife(tech)                     lifetime of tech [years]
         pCrf(tech)                      capital recovery factor
*BUILD LIMITS ON NEW TECHS
         pNmax(z,tech)                   max number techs built per zone and tech
*ZONAL PARAMETERS
         pDemand(z,h)                    hourly electricity demand [GWh]
         pPlanningreserve(z)             planning margin reserve capacity [GW]
         pPeakhtozone(peakh)             peakh for each zone
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
         ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load egu, windegu, solaregu, hydroegu, pumphydroegu, tech, techcurtailed, techrenew, technotcurtailed, c
$load h, peakh, springh, summerh, winterh, fallh, specialh, z, l
$load pCapac, pCapactech, pCapactechcurtailed, pHr, pHrtech, pOpcost, pOpcosttech
$load pFom, pOcc, pCO2emrate, pCO2emratetech, pCO2emcap, pCf, pMaxgenwind, pMaxgensolar
$load pMaxhydrogensum, pMaxhydrogenspr, pMaxhydrogenwin, pMaxhydrogenfal, pMaxhydrogenspe
$load pDemand, pEguzones, pCellzones, pLinesources, pLinesinks, pLinecapacs
$load pR, pLife, pNmax, pPlanningreserve, pPeakhtozone, pWeightspring, pWeightsummer, pWeightfall, pWeightwinter
$load pInitsoc, pMaxsoc, pEfficiency
$gdxin

*DEFINE NON-PEAK HOURS THAT WILL BE INCLUDED IN DEMAND=SUPPLY EQUATION
dispatchh(h) = winterh(h) + summerh(h) + specialh(h) + fallh(h) + springh(h);

*CALCULATE CAPITAL RECOVERY FACTOR
pCrf(tech) = pR / (1 - (1 / ( (1 + pR)**pLife(tech))));

Variable
         vZ                              obj func [thousand USD per yr]
         vIc                             total investment costs for new plants = fixed O&M + overnight capital costs [thousand USD per yr]
         vVc                             total variable costs for new and existing plants = variable O&M + fuel + emission costs [thousand USD per yr]
         vVcspring                       variable costs for spring hours
         vVcsummer                       variable costs for summer hours
         vVcfall                         variable costs for fall hours
         vVcwinter                       variable costs for winter hours
         vVcspecial                      variable costs for special hours
         ;

Positive variables
         vPtechcurtailed(c,techcurtailed,h)                      hourly electricity generation by new techs that can be curtailed by cell c and hour h [GWh]
         vPtechrenew(z,techrenew,h)                              hourly electricity generation by new renewable techs [GWh]
         vPtechnotcurtailed(z,technotcurtailed,h)                hourly electricity generation by new plants [GWh]
         vPegu(egu,h)                                            hourly electricity generation by existing plants [GWh]
         vLineflow(l,h)                                          flow over lines per hour (GW)
         vCO2emsannual                                           co2 emissions in entire year from new and existing plants [short ton]
         vCO2emssummer                                           co2 emissions in summer from new and existing plants [short ton]
         vCO2emsspring                                           co2 emissions in spring from new and existing plants [short ton]
         vCO2emswinter                                           co2 emissions in winter from new and existing plants [short ton]
         vCO2emsfall                                             co2 emissions in fall from new and existing plants [short ton]
         vCO2emsspecial                                          co2 emissions in special hours from new and existing plants [short ton]
         vSoc(pumphydroegu,h)                                    state of charge of pumped hydro units [GWh]
         vCharge(pumphydroegu,h)                                 charging by pumped hydro units [GWh]
         ;

Integer variable
         vNcurtailed(c,techcurtailed)                            number of newly constructed plants of types that can be curtailed in cell c
         vNrenew(z,techrenew)                                    number of newly constructed renewable plants
         vNnotcurtailed(z,technotcurtailed)                      number of newly constructed plants of types that cant be curtailed in zone z
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
         meetdemand(z,h)                                         meet supply with demand
         meetreservemargin(z)                                    meet planning reserve requirement with installed capacity
         curtailedtechgen(c,techcurtailed,h)                     restrict electricity generation by new curtailed plants to number built and hourly capacities
         notcurtailedtechgen(z,technotcurtailed,h)               restrict elec gen by not curtailed plants to number built and constant capacity
         renewtechgen(z,techrenew,h)                             restrict electricity generation by new renewables to number built and capacity and capacity factor
         maxzonalbuildcurtailed(z,techcurtailed)                 max number curtailed techs built per zone
         maxzonalbuildrenew(z,techrenew)                         max number renew techs built per zone
         maxzonalbuildnotcurtailed(z,technotcurtailed)           max number noncurtailed techs built per zone
         egugen(egu,h)                                           restrict electricity generation by existing generators to hourly capacities
         eguwindgen(z,h)                                         restrict electricity generation by existing wind generation to maximum aggregate output per zone
         egusolargen(z,h)                                        restrict electricity generation by existing solar generation to maximum aggregate output per zone
         hydrogenspr(hydroegu)                                   restrict total max elec gen by hydro units in spring
         hydrogensum(hydroegu)                                   restrict total max elec gen by hydro units in summer
         hydrogenwin(hydroegu)                                   restrict total max elec gen by hydro units in winter
         hydrogenfal(hydroegu)                                   restrict total max elec gen by hydro units in fall
         hydrogenspe(hydroegu)                                   restrict total max elec gen by hydro units in special hours
         co2emsannual                                            sum annual co2 emissions
         co2emssummer                                            calculate summer co2 emissions
         co2emsspring                                            calculate spring co2 emissions
         co2emswinter                                            calculate winter co2 emissions
         co2emsfall                                              calculate fall co2 emissions
         co2emsspecial                                           calculate special hour co2 emissions
         enforceco2emissionscap                                  restrict total co2 emissions to cap
         genandsoc(pumphydroegu,h)                               restrict gen by pump hydro to state of charge
         maxsto(pumphydroegu,h)                                  set max soc limit
         limitcharging(pumphydroegu,h)                           limit charging to capacity times efficiency
         defsocspr(pumphydroegu,springh)                         link state of charge and charging and discharging
         defsocsum(pumphydroegu,summerh)                         link state of charge and charging and discharging
         defsocfal(pumphydroegu,fallh)                           link state of charge and charging and discharging
         defsocwin(pumphydroegu,winterh)                         link state of charge and charging and discharging
         defsocspe(pumphydroegu,specialh)                        link state of charge and charging and discharging
         ;

******************OBJECTIVE FUNCTION******************
*Objective: minimize fixed + variable costs
objfunc..                vZ =e= vIc + vVc;
******************************************************

******************CALCULATE COSTS******************
*Fixed costs = annual fixed O&M + fixed annualized capital costs
investmentcost..         vIc =e= sum((c,techcurtailed),vNcurtailed(c,techcurtailed)*pCapactech(techcurtailed)*(pFom(techcurtailed)+pOcc(techcurtailed)*pCrf(techcurtailed)))
                               + sum((z,techrenew),vNrenew(z,techrenew)*pCapactech(techrenew)*(pFom(techrenew)+pOcc(techrenew)*pCrf(techrenew)))
                               + sum((z,technotcurtailed),vNnotcurtailed(z,technotcurtailed)*pCapactech(technotcurtailed)*(pFom(technotcurtailed)+pOcc(technotcurtailed)*pCrf(technotcurtailed)));

*Variable costs = electricity generation costs by new and existing plants
varcosttotal..           vVc =e= vVcspring + vVcsummer + vVcwinter + vVcfall + vVcspecial;
varcostspring..          vVcspring =e= pWeightspring*sum(springh,sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,springh)*pOpcosttech(techcurtailed))
                                                               + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,springh)*pOpcosttech(technotcurtailed))
                                                               + sum((techrenew,z),vPtechrenew(z,techrenew,springh)*pOpcosttech(techrenew))
                                                               + sum(egu,vPegu(egu,springh)*pOpcost(egu)));
varcostsummer..          vVcsummer =e= pWeightsummer*sum(summerh,sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,summerh)*pOpcosttech(techcurtailed))
                                                               + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,summerh)*pOpcosttech(technotcurtailed))
                                                               + sum((techrenew,z),vPtechrenew(z,techrenew,summerh)*pOpcosttech(techrenew))
                                                               + sum(egu,vPegu(egu,summerh)*pOpcost(egu)));
varcostwinter..          vVcwinter =e= pWeightwinter*sum(winterh,sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,winterh)*pOpcosttech(techcurtailed))
                                                               + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,winterh)*pOpcosttech(technotcurtailed))
                                                               + sum((techrenew,z),vPtechrenew(z,techrenew,winterh)*pOpcosttech(techrenew))
                                                               + sum(egu,vPegu(egu,winterh)*pOpcost(egu)));
varcostfall..            vVcfall =e= pWeightfall*sum(fallh,sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,fallh)*pOpcosttech(techcurtailed))
                                                         + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,fallh)*pOpcosttech(technotcurtailed))
                                                         + sum((techrenew,z),vPtechrenew(z,techrenew,fallh)*pOpcosttech(techrenew))
                                                         + sum((egu),vPegu(egu,fallh)*pOpcost(egu)));
varcostspecial..         vVcspecial =e= sum(specialh,sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,specialh)*pOpcosttech(techcurtailed))
                                                   + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,specialh)*pOpcosttech(technotcurtailed))
                                                   + sum((techrenew,z),vPtechrenew(z,techrenew,specialh)*pOpcosttech(techrenew))
                                                   + sum(egu,vPegu(egu,specialh)*pOpcost(egu)));
***************************************************

******************SYSTEM-WIDE GENERATION AND CAPACITY CONSTRAINTS******************
*Demand = generation by new and existing plants in each zone
meetdemand(z,dispatchh)..       sum((techcurtailed,c)$[pCellzones(c)=ORD(z)],vPtechcurtailed(c,techcurtailed,dispatchh))
                              + sum(technotcurtailed,vPtechnotcurtailed(z,technotcurtailed,dispatchh))
                              + sum(techrenew,vPtechrenew(z,techrenew,dispatchh))
                              + sum(egu$[pEguzones(egu)=ORD(z)],vPegu(egu,dispatchh))
                              - sum(l$[pLinesources(l)=ORD(z)],vLineflow(l,dispatchh))
                              + sum(l$[pLinesinks(l)=ORD(z)],vLineflow(l,dispatchh))
                              - sum(pumphydroegu,vCharge(pumphydroegu,dispatchh)) =e= pDemand(z,dispatchh);

*Total installed capacity must exceed peak demand + planning reserve margin in each zone
meetreservemargin(z)..  pPlanningreserve(z) =l= sum(peakh$[pPeakhtozone(peakh)=ORD(z)], sum(egu$[pEguzones(egu)=ORD(z)],pCapac(egu,peakh))
                                                                                      + sum((techcurtailed,c)$[pCellzones(c)=ORD(z)],pCapactechcurtailed(c,techcurtailed,peakh)*vNcurtailed(c,techcurtailed))
                                                                                      + sum(technotcurtailed,pCapactech(technotcurtailed)*vNnotcurtailed(z,technotcurtailed))
                                                                                      + sum(techrenew,pCapactech(techrenew)*vNrenew(z,techrenew)*pCf(z,techrenew,peakh)));
***********************************************************************************

*****************LINE FLOW LIMITS*************************
*Limit max value of line flow to line capacity (already bounded at 0 since declared as positive variable)
vLineflow.up(l,h)=pLinecapacs(l);
**********************************************************

******************GENERATION CONSTRAINTS ON NEW UNITS******************
*Generation by techs that can be curtailed is limited by curtailed capacity
curtailedtechgen(c,techcurtailed,h)..                    vPtechcurtailed(c,techcurtailed,h) =l= vNcurtailed(c,techcurtailed)*pCapactechcurtailed(c,techcurtailed,h);

*Generation by techs that cant be curtailed is limited by their capacity
notcurtailedtechgen(z,technotcurtailed,h)..              vPtechnotcurtailed(z,technotcurtailed,h) =l= vNnotcurtailed(z,technotcurtailed)*pCapactech(technotcurtailed);

*Renewable generation limited by capacity factor and nameplate capacity
renewtechgen(z,techrenew,h)..                            vPtechrenew(z,techrenew,h) =l= vNrenew(z,techrenew)*pCapactech(techrenew)*pCf(z,techrenew,h);
***********************************************************************

******************BUILD DECISIONS******************
*Limit total number builds to input value
maxzonalbuildcurtailed(z,techcurtailed)..               sum(c$[pCellzones(c)=ORD(z)],vNcurtailed(c,techcurtailed)) =l= pNmax(z,techcurtailed);
maxzonalbuildrenew(z,techrenew)..                       vNrenew(z,techrenew) =l= pNmax(z,techrenew);
maxzonalbuildnotcurtailed(z,technotcurtailed)..         vNnotcurtailed(z,technotcurtailed) =l= pNmax(z,technotcurtailed);

*Change upper bound on build variables
vNcurtailed.up(c,techcurtailed) = 1000;
vNrenew.up(z,techrenew) = 1000;
vNnotcurtailed.up(z,technotcurtailed) = 1000;
***************************************************

******************GENERATION CONSTRAINTS ON EXISTING UNITS******************
*Enforce max hourly generation limit based on hourly capacity
egugen(egu,h)..          vPegu(egu,h) =l= pCapac(egu,h);

*Enforce max generation on existing wind and solar plants
eguwindgen(z,h)..        pMaxgenwind(z,h) =g= sum(windegu$[pEguzones(windegu)=ORD(z)],vPegu(windegu,h));
egusolargen(z,h)..       pMaxgensolar(z,h) =g= sum(solaregu$[pEguzones(solaregu)=ORD(z)],vPegu(solaregu,h));

*Enforce max generation on total hydro generation
hydrogenspr(hydroegu).. pMaxhydrogenspr(hydroegu) =g= sum(springh,vPegu(hydroegu,springh));
hydrogensum(hydroegu).. pMaxhydrogensum(hydroegu) =g= sum(summerh,vPegu(hydroegu,summerh));
hydrogenwin(hydroegu).. pMaxhydrogenwin(hydroegu) =g= sum(winterh,vPegu(hydroegu,winterh));
hydrogenfal(hydroegu).. pMaxhydrogenfal(hydroegu) =g= sum(fallh,vPegu(hydroegu,fallh));
hydrogenspe(hydroegu).. pMaxhydrogenspe(hydroegu) =g= sum(specialh,vPegu(hydroegu,specialh));
****************************************************************************

******************CO2 EMISSIONS CONSTRAINT******************
*Co2 emissions = electricity generation * co2 emissions rate
co2emsannual..   vCO2emsannual =e= vCO2emsspring + vCO2emssummer + vCO2emswinter + vCO2emsfall + vCO2emsspecial;
co2emsspring..   vCO2emsspring =e= pWeightspring*sum(springh,sum(egu,vPegu(egu,springh)*pHr(egu)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,springh)*pHrtech(techcurtailed)*pCO2emratetech(techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,springh)*pHrtech(technotcurtailed)*pCO2emratetech(technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(z,techrenew,springh)*pHrtech(techrenew)*pCO2emratetech(techrenew)));
co2emssummer..   vCO2emssummer =e= pWeightsummer*sum(summerh,sum(egu,vPegu(egu,summerh)*pHr(egu)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,summerh)*pHrtech(techcurtailed)*pCO2emratetech(techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,summerh)*pHrtech(technotcurtailed)*pCO2emratetech(technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(z,techrenew,summerh)*pHrtech(techrenew)*pCO2emratetech(techrenew)));
co2emswinter..   vCO2emswinter =e= pWeightwinter*sum(winterh,sum(egu,vPegu(egu,winterh)*pHr(egu)*pCO2emrate(egu))
                                                           + sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,winterh)*pHrtech(techcurtailed)*pCO2emratetech(techcurtailed))
                                                           + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,winterh)*pHrtech(technotcurtailed)*pCO2emratetech(technotcurtailed))
                                                           + sum((techrenew,z),vPtechrenew(z,techrenew,winterh)*pHrtech(techrenew)*pCO2emratetech(techrenew)));
co2emsfall..     vCO2emsfall =e= pWeightfall*sum(fallh,sum(egu,vPegu(egu,fallh)*pHr(egu)*pCO2emrate(egu))
                                                     + sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,fallh)*pHrtech(techcurtailed)*pCO2emratetech(techcurtailed))
                                                     + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,fallh)*pHrtech(technotcurtailed)*pCO2emratetech(technotcurtailed))
                                                     + sum((techrenew,z),vPtechrenew(z,techrenew,fallh)*pHrtech(techrenew)*pCO2emratetech(techrenew)));
co2emsspecial..  vCO2emsspecial =e= sum(specialh,sum(egu,vPegu(egu,specialh)*pHr(egu)*pCO2emrate(egu))
                                               + sum((techcurtailed,c),vPtechcurtailed(c,techcurtailed,specialh)*pHrtech(techcurtailed)*pCO2emratetech(techcurtailed))
                                               + sum((technotcurtailed,z),vPtechnotcurtailed(z,technotcurtailed,specialh)*pHrtech(technotcurtailed)*pCO2emratetech(technotcurtailed))
                                               + sum((techrenew,z),vPtechrenew(z,techrenew,specialh)*pHrtech(techrenew)*pCO2emratetech(techrenew)));

*Meet emissions cap
enforceco2emissionscap.. vCO2emsannual =l= pCO2emcap;
************************************************************

******************STORAGE CONSTRAINTS******************
*Limit generation to state of charge
genandsoc(pumphydroegu,h) .. vPegu(pumphydroegu,h) =l= vSoc(pumphydroegu,h);

*Link state of charge, charging, and discharging
defsocspr(pumphydroegu,springh) .. vSoc(pumphydroegu,springh) =e= pInitsoc(pumphydroegu)$[ORD(springh)=1] + vSoc(pumphydroegu,springh-1)$[ORD(springh)>1] - vPegu(pumphydroegu,springh)
                                                                                 + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,springh);
defsocsum(pumphydroegu,summerh) .. vSoc(pumphydroegu,summerh) =e= pInitsoc(pumphydroegu)$[ORD(summerh)=1] + vSoc(pumphydroegu,summerh-1)$[ORD(summerh)>1] - vPegu(pumphydroegu,summerh)
                                                                                 + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,summerh);
defsocfal(pumphydroegu,fallh) .. vSoc(pumphydroegu,fallh) =e= pInitsoc(pumphydroegu)$[ORD(fallh)=1] + vSoc(pumphydroegu,fallh-1)$[ORD(fallh)>1] - vPegu(pumphydroegu,fallh)
                                                                                 + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,fallh);
defsocwin(pumphydroegu,winterh) .. vSoc(pumphydroegu,winterh) =e= pInitsoc(pumphydroegu)$[ORD(winterh)=1] + vSoc(pumphydroegu,winterh-1)$[ORD(winterh)>1] - vPegu(pumphydroegu,winterh)
                                                                                 + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,winterh);
defsocspe(pumphydroegu,specialh) .. vSoc(pumphydroegu,specialh) =e= pInitsoc(pumphydroegu)$[ORD(specialh)=1] + vSoc(pumphydroegu,specialh-1)$[ORD(specialh)>1] - vPegu(pumphydroegu,specialh)
                                                                                 + pEfficiency(pumphydroegu) * vCharge(pumphydroegu,specialh);

*Limit state of charge to maximum storage capacity
maxsto(pumphydroegu,h) .. vSoc(pumphydroegu,h) =l= pMaxsoc(pumphydroegu);

*Limit rate of charging to capacity times efficiency
limitcharging(pumphydroegu,h) .. vCharge(pumphydroegu,h) =l= pCapac(pumphydroegu,h);
***************************************************

Model expansion includes all equations /all/;
solve expansion using mip minimizing vZ;

pModelstat = expansion.Modelstat;
pSolvestat = expansion.solvestat;

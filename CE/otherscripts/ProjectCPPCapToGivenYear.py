#Michael Craig
#October 4, 2016
#Project CPP CO2 emissions cap to given year (forward or backward) based on 2022
#& 2030 limits.  

#Inputs: year to project cap to, 2022 and 2030 CPP CO2 emissions cap (short tons)
#Outputs: projected CO2 emissions cap (short tons)
def setCppCapBasedOnCurrYear(currYear,co2Cpp2022Limit,co2Cpp2030Limit):
    (deltaYear,deltaEms) = (2030-2022,co2Cpp2022Limit-co2Cpp2030Limit)
    emsReduxPerYear = deltaEms/deltaYear
    diffCurrYearFrom2022 = currYear-2022
    return co2Cpp2022Limit - diffCurrYearFrom2022*emsReduxPerYear

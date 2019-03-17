def calc_stream_temperature_1_plant(T_stream, T_eq, k_coeff, q_power, Q, depth,stream_vel, impulse_fxn2, distx, sub_distx, dist_impulse, time_stepx):
    
    """
    
    Equation to calculate the stream temperature
    based on the analytical solution of the Edenger et al. (1968), equation #19
    
    """

    timex = (distx + sub_distx/(len(steps_per_km)))/stream_vel # calculate distance traveled in seconds
    time_impulse_sec = dist_impulse/stream_vel # convert location of power plant effluent impulse (km) to seconds based on velocity
    k_coeff = k_coeff * kcal_to_J # convert kcal/sec*m2*C to Joules/sec*m2*C
    alpha = k_coeff/(water_density * c_p_j * depth) # alpha parameter from solution
    beta = (q_power * a_const)/(water_density * c_p_j * Q) # calculate beta parameter
    power_x  = impulse_fxn2[distx]*beta*exp((-1)*alpha[0]*(timex - time_impulse_sec)) # calculate power plant term
    stream_x = (T_stream - T_eq)*exp((-1)*alpha[0]*timex) + T_eq + power_x   #  calculate stream temperature at timex
    
    return stream_x



# --------- input variables ------------
# variable list: Q(m3/sec), depth(m), stream_velocity(km/sec) T_stream - initial(deg C), dry bulb temp (deg C), vapor_pressure (mb), q_ns (kcal/m2*sec),
# q_na (kcal/m2*sec),  wind (m/sec), location of power plant effluent discharge (km), distance between calculations (km)
input_vars = pd.DataFrame([100,  5.27e-4,  28.481531,  30, 900,  0.02, 0.05,  2, 4, 0.1])
input_vars = input_vars.transpose()
input_vars.columns = ['Q','stream_velocity', 'T_stream_i','T_dry_bulb','vap_press','q_ns','q_na','wind', 'power_plant_km', 'dist_between']
a_const = 1
heat_x = 20500  # constant heat rate in Btu/kWh (max observed in SERC = 20,500)

# ----- establish time series with power plant impulse ---------
km1 = 0 # first distance
km2 = 30 # final distance
impulse_distance_series = impseq2(int(input_vars['power_plant_km'][0]),km1,km2)

# ------- flow and power plant capacity scenarios ------------
power_plant_range = [1000,2000,3000,4000] # power plant capacity (MW)
Q_range = [100,1000,10000, 100000] # discharge range (m3/sec)

# ------- data frames and lists to write data -------
T_stream_tot = [] # list of stream temperature time series
scenarios_tot_Q = [] # list of streamflow scenarios
scenarios_tot_P = [] # list of power plant capacity scenarios
T_equil = [] # list of equilibrium temperatures

# ------------------- loop through all power plant capacity scenarios ------------------
for power_x in power_plant_range:
    
    # --------- loop through all flow scenarios ---------
    for flow_x in Q_range:
        
        T_stream_list2 = [] # empty list to write stream T time series

        # ------------- calculate surface energy, T_equil., and k_coeff -------    
        q_surf, T_eq, k_coeff = calc_surf_energy(input_vars['T_stream_i'],input_vars['T_dry_bulb'],
                input_vars['vap_press'],input_vars['q_ns'],input_vars['q_na'],
                input_vars['wind'])
        T_equil.append(T_eq)
        # -------------------- calculate power plant waste heat output ----------------
        #q_power = power_plant_efficiency(power_x, 0.35)
        waste_heat = waste_heat_fxn(power_x, heat_x) # output [J/sec]
        
        # --------------- calculate flow depth --------------------
        depth_x = a_d * ((flow_x*m3sec_to_ft3sec)**b_d) # depth in feet
        depth_x = depth_x /m_to_ft  # depth in meters

        # -------------- print scenarios -----------
        print('power: ',power_x, ' flow: ',flow_x, ' depth: ',depth_x)   
        
        # ------------- loop through each hour of time series -----------
        for i in range(len(impulse_distance_series)):

            steps_per_km = list(range(0,int(1/input_vars['dist_between']))) # list with total time steps in each hour
            
            # update steps_per_hr if more than one time step per hour
            if len(steps_per_km) < 1:
                steps_per_km =  list([1])
            
            # loop through each time step in the hour and calculate stream temperature
            for j in steps_per_km:
                T_stream_new = calc_stream_temperature_1_plant(input_vars['T_stream_i'], T_eq, k_coeff, waste_heat, 
                    flow_x, depth_x, input_vars['stream_velocity'][0],
                    impulse_distance_series, i, j,int(input_vars['power_plant_km'][0]),input_vars['dist_between'] )
                
                T_stream_list2.append(T_stream_new)

        # --------- write each scenario ------
        scenarios_tot_Q.append(flow_x)
        scenarios_tot_P.append(power_x)
        # ------- append new time series of stream temperatures --------
        T_stream_tot.append(T_stream_list2)




def calc_surf_energy(T_i,dbt,ea, q_ns, q_na, wind):
    
    """
    Equation to calculate surface energy, equilibrium temperature, and constant k
    Based on River Basin Model (RBM) 'Energy.f90' subroutine. 
    
    """
    
    # do linear fit of temperature between two T_fit
    T_fit = []
    T_fit.append(T_i - 1.0)
    T_fit.append(T_i + 1.0)
    
    q_fit = [] # q_fit empty list

    # loop between two T_fit temperatures
    for i in range(0,2):
        # calculate energy from evaporation
        e0 = 2.1718e8*exp((-4157.0) /(T_fit[i] + 239.09)) # calc saturation vap pressure
        lvp = 597 - (0.57 * T_fit[i])
                        # calc latent heat of vapor
        q_evap=1000.*lvp*evap_coeff*wind
        if q_evap.values < 0: 
            q_evap=0
        q_evap=q_evap*(e0-ea)

        # calculate energy tranfser from convection
        rb=pf*(dbt-T_fit[i])
        q_conv=rb*q_evap

        # calculate loss of land surface longwave radiation
        q_ws=6.693e-2 + 1.471e-3 * T_fit[i]

        # calculate all energy components (units = kcal/sec*m2) 
        q_fit.append(q_ns + q_na - q_ws - q_evap + q_conv)

#
#     q=AT+B
#
#     Linear fit over the range of 2.0 deg C.
#     These results can be used to estimate the "equilibrium" 
#     temperature and linear rate constant.
#
    A = (q_fit[0] - q_fit[1]) / (T_fit[0] - T_fit[1]) #units: kcal/sec * m2 * deg C
    q_surf=0.5*(q_fit[0]+q_fit[1]) # units: kcal/sec * m2
    # B=(q_surf/A)-(T_fit[0]+T_fit[1])/2 # incorrect - in original Energy.f90 code
    B = q_surf - (A *(T_fit[0]+T_fit[1])/2)

    # calculate equilibrium temperature and k_coefficient 
    T_equil = (-1)*(B/A)  # deg C
    k_coeff = A
    
    return q_surf, T_equil, k_coeff
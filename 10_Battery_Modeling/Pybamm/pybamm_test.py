import  pybamm
#print(list(pybamm.parameter_sets)) #shows list of available parameter sets
#x =  input("Enter battery type: ")
batt_type_dict =  { #dictionary of available types of batteries, #need to go back to supply more descriptive comments 
                 "BBOXX lead-acid":  "Sulzer2019",
                 "Enertech":  "Ai2020",
                 "LG _Chen":  "Chen2020",
                 "graphite/silicon":  "Chen2020_composite",
                 "Kokam SLPB 75106100_Ecker":  "Ecker2015",
                 "Graphite negative electrode":  "Graphite negative electrode parameters",
                 "Graphite halfcell": "Ecker2015_graphite_halfcell",
                 "Graphite electrode":  "Graphite electrode parameters",
                 "NMC622":  "MSMR_Example",
                 "Kokam SLPB78205130H_Marquis":  "Marquis 2019",
                 "NMC532 pouch":  "Mohtat2020",
                 "Lithium plating":  "SEI parameters", 
                 "NCA pouch": "NCA_Kim2011",
                 "LG M50_Okane":  "OKane2022",
                 "LG M50_graphite+SiOx negative electrode":  "OKane2022_graphite_SiOx_halfcell",
                 "LG M50_ORegan":  "ORegan2022",
                 "LFP": "Prada2013",
                 "LCO":  "Ramadass2004", #use with caution
                 "Kokam SLPB78205130H half-cell":  "Xu2019",
                 "1C discharge":  "1C discharge from full",
                 } #general thought is that user can input type of battery with a more user-friendly name

#chemistry =  pybamm.parameter_sets(batt_type_dict(x)) #sets battery chemistry, currently have to set as a value
#parameter_values = pybamm.ParameterValues(y)
current_profile =  "" #inputted current profile for aircraft simulation, preferrable to use a constant time step
#Experiment class?
lithium_model =  pybamm.lithium_ion.DFN()
lithium_simulation =  pybamm.Simulation(lithium_model)
#lithium_simulation.solve([0, 10])
#lithium_simulation.plot()
print()
import  pandas as  pd
import  os


def save_data_stick_fixed(CG_bat_1, 
                          CG_bat_2,
                          vehicle_tag,
                          optimized_vehicle, 
                          output_stick_fixed, 
                          planform_optimization_problem,
                          elapsed_time_stick_fixed
                          ) :
    
    data_stick_fixed = pd.DataFrame({ 
                                     'AoA'                            : [output_stick_fixed[0]],
                                     'mw_span'                        : [output_stick_fixed[1]],
                                     'mw_AR'                          : [output_stick_fixed[2]],
                                     'mw_root_twist'                  : [output_stick_fixed[3]],
                                     'mw_tip_twist'                   : [output_stick_fixed[4]],
                                     'vt_span'                        : [output_stick_fixed[5]],
                                     'vt_AR'                          : [output_stick_fixed[6]],
                                     'ht_span'                        : [output_stick_fixed[7]],
                                     'ht_AR'                          : [output_stick_fixed[8]],
                                     'ht_tip_twist'                   : [output_stick_fixed[9]],     
                                     'CD'                             : [planform_optimization_problem.summary['CD'                        ]],
                                     'CM_residual'                    : [planform_optimization_problem.summary['CM_residual'               ]],
                                     'spiral_criteria'                : [planform_optimization_problem.summary['spiral_criteria'           ]],
                                     'static_margin'                  : [planform_optimization_problem.summary['static_margin'             ]],
                                     'CM_alpha'                       : [planform_optimization_problem.summary['CM_alpha'                  ]],
                                     'phugoid_damping_ratio'          : [planform_optimization_problem.summary['phugoid_damping_ratio'     ]],
                                     'short_period_damping_ratio'     : [planform_optimization_problem.summary['short_period_damping_ratio']],
                                     'dutch_roll_frequency'           : [planform_optimization_problem.summary['dutch_roll_frequency'      ]],
                                     'dutch_roll_damping_ratio'       : [planform_optimization_problem.summary['dutch_roll_damping_ratio'  ]],
                                     'spiral_doubling_time'           : [planform_optimization_problem.summary['spiral_doubling_time'      ]],
                                     'CD'                             : [planform_optimization_problem.summary['CD'                        ]],
                                     'CM_residual'                    : [planform_optimization_problem.summary['CM_residual']],
                                     'spiral_criteria'                : [planform_optimization_problem.summary['spiral_criteria']],
                                     'static_margin'                  : [planform_optimization_problem.summary['static_margin']],
                                     'CM_alpha'                       : [planform_optimization_problem.summary['CM_alpha']],
                                     'F_x_residual'                   : [planform_optimization_problem.summary['F_x_residual']],
                                     'F_y_residual'                   : [planform_optimization_problem.summary['F_y_residual']],
                                     'F_z_residual'                   : [planform_optimization_problem.summary['F_z_residual']],
                                     'M_x_residual'                   : [planform_optimization_problem.summary['M_x_residual']],
                                     'M_y_residual'                   : [planform_optimization_problem.summary['M_y_residual']],
                                     'M_z_residual'                   : [planform_optimization_problem.summary['M_z_residual']],
                                     'CN_beta'                        : [planform_optimization_problem.summary['CN_beta']],
                                     'CL_beta'                        : [planform_optimization_problem.summary['CL_beta']],
                                     'phugoid_damping_ratio'          : [planform_optimization_problem.summary['phugoid_damping_ratio']],
                                     'short_period_damping_ratio'     : [planform_optimization_problem.summary['short_period_damping_ratio']],
                                     'dutch_roll_frequency'           : [planform_optimization_problem.summary['dutch_roll_frequency']],
                                     'dutch_roll_damping_ratio'       : [planform_optimization_problem.summary['dutch_roll_damping_ratio']],
                                     'spiral_doubling_time'           : [planform_optimization_problem.summary['spiral_doubling_time']],
                                     'longmodes1'                     : [planform_optimization_problem.summary['longmodes1']],
                                     'longmodes2'                     : [planform_optimization_problem.summary['longmodes2']],
                                     'longmodes3'                     : [planform_optimization_problem.summary['longmodes3']],
                                     'longmodes4'                     : [planform_optimization_problem.summary['longmodes4']],
                                     'phugoidfreq'                    : [planform_optimization_problem.summary['phugoidfreq']],
                                     'shortPeriodFreqHz'              : [planform_optimization_problem.summary['shortPeriodFreqHz']],
                                     'phugoidtime'                    : [planform_optimization_problem.summary['phugoidtime']],
                                     'shortPeriodTime'                : [planform_optimization_problem.summary['shortPeriodTime']],
                                     'latmodes1'                      : [planform_optimization_problem.summary['latmodes1']],
                                     'latmodes2'                      : [planform_optimization_problem.summary['latmodes2']],
                                     'latmodes3'                      : [planform_optimization_problem.summary['latmodes3']],
                                     'latmodes4'                      : [planform_optimization_problem.summary['latmodes4']],
                                     'lat_dr_time'                    : [planform_optimization_problem.summary['lat_dr_time']],
                                     'lat_rollsubfreq'                : [planform_optimization_problem.summary['lat_rollsubfreq']],
                                     'rollSubsistenceTimeConstant'    : [planform_optimization_problem.summary['rollSubsistenceTimeConstant']],
                                     'rollSubsistenceDamping'         : [planform_optimization_problem.summary['rollSubsistenceDamping']],
                                     'spiralFreqHz'                   : [planform_optimization_problem.summary['spiralFreqHz']],
                                     'spiralDamping'                  : [planform_optimization_problem.summary['spiralDamping']],
                                     'run_time'                       : [elapsed_time_stick_fixed],
                                     'prop module 1 origin x '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin[0][0]],
                                     'prop module 1 origin y '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin[0][1]],
                                     'prop module 1 origin z '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin[0][2]],
                                     'prop module 2 origin x '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin[0][0]],
                                     'prop module 2 origin y '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin[0][1]],
                                     'prop module 2 origin z '        : [optimized_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin[0][2]],
                                     'lift module 1 origin x '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin[0][0]],
                                     'lift module 1 origin y '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin[0][1]],
                                     'lift module 1 origin z '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin[0][2]],
                                     'lift module 2 origin x '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin[0][0]],
                                     'lift module 2 origin y '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin[0][1]],
                                     'lift module 2 origin z '        : [optimized_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin[0][2]],
                                     'aircraft CG x '                 : [optimized_vehicle.mass_properties.center_of_gravity[0,0]],
                                     'aircraft CG y '                 : [optimized_vehicle.mass_properties.center_of_gravity[0,1]],
                                     'aircraft CG z '                 : [optimized_vehicle.mass_properties.center_of_gravity[0,2]],
                                     'aircraft MOI xx '               : [optimized_vehicle.mass_properties.moments_of_inertia.tensor[0,0]],
                                     'aircraft MOI yy '               : [optimized_vehicle.mass_properties.moments_of_inertia.tensor[1,1]],  
                                     'aircraft MOI zz '               : [optimized_vehicle.mass_properties.moments_of_inertia.tensor[2,2]],
                                     }) 
    
    cg_x1   =  str(CG_bat_1[0]).replace('.', "")
    cg_y1   =  str(CG_bat_1[1]).replace('.', "")
    cg_z1   =  str(CG_bat_1[2]).replace('.', "")
    cg_x2   =  str(CG_bat_2[0]).replace('.', "")
    cg_y2   =  str(CG_bat_2[1]).replace('.', "")
    cg_z2   =  str(CG_bat_2[2]).replace('.', "")
    
    def new_filename(filename):
        return filename.replace('[', '').replace(']', '').replace(' ', '_').replace('.', '_')

    raw_file_name = f"{cg_x1}_{cg_y1}_{cg_z1}_{cg_x2}_{cg_y2}_{cg_z2}_Baseline_{vehicle_tag}_Opt_Results"
    excel_file_name = new_filename(raw_file_name) 

    current_dir =  os.path.dirname(os.path.abspath(__file__))
    
    # Go one folder back and into Raw_Data
    load_dir = os.path.join(current_dir)
    
    excel_file = os.path.join(load_dir, excel_file_name + '.xlsx')
    
    with pd.ExcelWriter(excel_file, mode='w', engine='openpyxl') as writer:
        data_stick_fixed.to_excel(writer, sheet_name='Stick_Fixed', index=False)

    return
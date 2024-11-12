import  pandas as  pd


def Optimization_Data_2_CSV(output_stick_fixed, 
                            output_elevator_sizing, 
                            output_aileron_and_rudder_sizing, 
                            output_flap_sizing, 
                            planform_optimization_problem, 
                            elevator_sizing_optimization_problem, 
                            aileron_rudder_sizing_optimization_problem, 
                            flap_sizing_optimization_problem,
                            elapsed_time_stick_fixed,
                            elapsed_time_elevator_sizing,
                            elapsed_time_aileron_and_rudder_sizing,
                            elapsed_time_flap_sizing):
    
    data_stick_fixed = pd.DataFrame({
                                     'mw_span'                   : [output_stick_fixed[0]],   
                                     'mw_AR'                     : [output_stick_fixed[1]],   
                                     'mw_root_twist'             : [output_stick_fixed[2]],   
                                     'mw_tip_twist'              : [output_stick_fixed[3]],   
                                     'vt_span'                   : [output_stick_fixed[4]],   
                                     'vt_AR'                     : [output_stick_fixed[5]],  
                                     'ht_span'                   : [output_stick_fixed[6]],                            
                                     'ht_AR'                     : [output_stick_fixed[7]],                            
                                     'AoA'                       : [output_stick_fixed[8]],           
                                     'CD'                        : [planform_optimization_problem.summary['CD'                        ]],
                                     'CL_residual'               : [planform_optimization_problem.summary['CL_residual'               ]],
                                     'CM_residual'               : [planform_optimization_problem.summary['CM_residual'               ]],
                                     'static_margin'             : [planform_optimization_problem.summary['static_margin'             ]],
                                     'CM_alpha'                  : [planform_optimization_problem.summary['CM_alpha'                  ]],
                                     'phugoid_damping_ratio'     : [planform_optimization_problem.summary['phugoid_damping_ratio'     ]],
                                     'short_period_damping_ratio': [planform_optimization_problem.summary['short_period_damping_ratio']],
                                     'dutch_roll_frequency'      : [planform_optimization_problem.summary['dutch_roll_frequency'      ]],
                                     'dutch_roll_damping_ratio'  : [planform_optimization_problem.summary['dutch_roll_damping_ratio'  ]],
                                     'spiral_doubling_time'      : [planform_optimization_problem.summary['spiral_doubling_time'      ]],
                                     'spiral_criteria'           : [planform_optimization_problem.summary['spiral_criteria'           ]],
                                     'run_time'                  : [elapsed_time_stick_fixed]
                                    }) 
    
    data_elevator_sizing = pd.DataFrame({
                                         'ht_elevator_chord_fraction' : [output_elevator_sizing[0]],   
                                         'ht_elevator_span_frac_start': [output_elevator_sizing[1]],   
                                         'ht_elevator_span_frac_end'  : [output_elevator_sizing[2]],      
                                         'elevator_surface_area'      : [elevator_sizing_optimization_problem.summary['elevator_surface_area'   ]],
                                         'elevator_push_deflection'   : [elevator_sizing_optimization_problem.summary['elevator_push_deflection']],
                                         'elevator_pull_deflection'   : [elevator_sizing_optimization_problem.summary['elevator_pull_deflection']]
                                        })
    
    data_aileron_and_rudder_sizing = pd.DataFrame({
                                                   'mw_aileron_chord_fraction'   : [output_aileron_and_rudder_sizing[0]],   
                                                   'mw_aileron_span_frac_start'  : [output_aileron_and_rudder_sizing[1]],   
                                                   'mw_aileron_span_frac_end'    : [output_aileron_and_rudder_sizing[2]],      
                                                   'vs_rudder_chord_fraction'    : [output_aileron_and_rudder_sizing[3]],
                                                   'vs_rudder_span_frac_start'   : [output_aileron_and_rudder_sizing[4]],
                                                   'vs_rudder_span_frac_end'     : [output_aileron_and_rudder_sizing[5]],
                                                   'aileron_rudder_surface_area' : [aileron_rudder_sizing_optimization_problem.summary['aileron_rudder_surface_area' ]],
                                                   'aileron_roll_deflection'     : [aileron_rudder_sizing_optimization_problem.summary['aileron_roll_deflection'     ]],
                                                   'rudder_roll_deflection'      : [aileron_rudder_sizing_optimization_problem.summary['rudder_roll_deflection'      ]],
                                                   'aileron_crosswind_deflection': [aileron_rudder_sizing_optimization_problem.summary['aileron_crosswind_deflection']],
                                                   'rudder_crosswind_deflection' : [aileron_rudder_sizing_optimization_problem.summary['rudder_crosswind_deflection' ]]
                                                  })
    
    data_flap_sizing = pd.DataFrame({
                                     'mw_flap_chord_fraction'   : [output_flap_sizing[0]],   
                                     'mw_flap_span_frac_start'  : [output_flap_sizing[1]],   
                                     'mw_flap_span_frac_end'    : [output_flap_sizing[2]],      
                                     'flap_surface_area'        : [flap_sizing_optimization_problem['flap_surface_area']],
                                     'flap_criteria'            : [flap_sizing_optimization_problem['flap_criteria'    ]]
                                    })    
    
    with pd.ExcelWriter('Output_Single_Engine_Piston_Baseline.xlsx', mode='a', engine='openpyxl') as writer:
        data_stick_fixed.to_excel(writer, sheet_name='Stick_Fixed', index=False)
        data_elevator_sizing.to_excel(writer, sheet_name='Elevator_Sizing', index=False)
        data_aileron_and_rudder_sizing.to_excel(writer, sheet_name='Aileron_Rudder_Sizing', index=False)
        data_flap_sizing.to_excel(writer, sheet_name='Flap_Sizing', index=False)
    
    return
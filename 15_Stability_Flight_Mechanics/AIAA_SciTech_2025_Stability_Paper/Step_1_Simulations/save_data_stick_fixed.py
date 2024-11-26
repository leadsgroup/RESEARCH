import  pandas as  pd


def save_data_stick_fixed(CG_bat_1, 
                          CG_bat_2,
                          vehicle_tag, 
                          output_stick_fixed, 
                          planform_optimization_problem,
                          elapsed_time_stick_fixed
                          ) :
    
    data_stick_fixed = pd.DataFrame({ 
                                     'mw_root_twist'             : [output_stick_fixed[0]],   
                                     'mw_tip_twist'              : [output_stick_fixed[1]],   
                                     'vt_span'                   : [output_stick_fixed[2]],   
                                     'vt_AR'                     : [output_stick_fixed[3]],  
                                     'ht_span'                   : [output_stick_fixed[4]],                            
                                     'ht_AR'                     : [output_stick_fixed[5]],                            
                                     'AoA'                       : [output_stick_fixed[6]],           
                                     'CD'                        : [planform_optimization_problem.summary['CD'                        ]],
                                     'CM_residual'               : [planform_optimization_problem.summary['CM_residual'               ]],
                                     'spiral_criteria'           : [planform_optimization_problem.summary['spiral_criteria'           ]],
                                     'static_margin'             : [planform_optimization_problem.summary['static_margin'             ]],
                                     'CM_alpha'                  : [planform_optimization_problem.summary['CM_alpha'                  ]],
                                     'CL_stick_fixed_residual'   : [planform_optimization_problem.summary['CL_stick_fixed_residual'   ]],
                                     'phugoid_damping_ratio'     : [planform_optimization_problem.summary['phugoid_damping_ratio'     ]],
                                     'short_period_damping_ratio': [planform_optimization_problem.summary['short_period_damping_ratio']],
                                     'dutch_roll_frequency'      : [planform_optimization_problem.summary['dutch_roll_frequency'      ]],
                                     'dutch_roll_damping_ratio'  : [planform_optimization_problem.summary['dutch_roll_damping_ratio'  ]],
                                     'spiral_doubling_time'      : [planform_optimization_problem.summary['spiral_doubling_time'      ]],
                                     'run_time'                  : [elapsed_time_stick_fixed]
                                    }) 
    
    
    cg_x1   =  str(CG_bat_1[0]).replace('.', "_")
    cg_y1   =  str(CG_bat_1[1]).replace('.', "_")
    cg_z1   =  str(CG_bat_1[2]).replace('.', "_")
    cg_x2   =  str(CG_bat_2[0]).replace('.', "_")
    cg_y2   =  str(CG_bat_2[1]).replace('.', "_")
    cg_z2   =  str(CG_bat_2[2]).replace('.', "_")

    excel_file_name =  cg_x1 + '' +  cg_y1 +'' +  cg_z1 +'' + cg_x2 + '' +  cg_y2 +'' +  cg_z2  +  '_Baseline+' + vehicle_tag+  '_Opt_Results.xlsx'
    with pd.ExcelWriter(excel_file_name, mode='a', engine='openpyxl') as writer:
        data_stick_fixed.to_excel(writer, sheet_name='Stick_Fixed', index=False)
    
    return
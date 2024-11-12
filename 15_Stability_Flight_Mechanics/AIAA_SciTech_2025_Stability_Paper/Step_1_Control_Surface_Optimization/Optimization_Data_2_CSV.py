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
    
    df = pd.DataFrame({
    'mw_root_twist': [output_stick_fixed[0]], 
    'mw_tip_twist': [output_stick_fixed[1]], 
    'mw_dihedral': [output_stick_fixed[2]], 
    'vt_span': [output_stick_fixed[3]], 
    'vt_taper': [output_stick_fixed[4]], 
    'c_g_x': [output_stick_fixed[5]],
    'CD': [planform_optimization_problem.summary['CD']],
    'CM_residual': [planform_optimization_problem.summary['CM_residual']],
    'spiral_criteria': [planform_optimization_problem.summary['spiral_criteria']],
    'static_margin': [planform_optimization_problem.summary['static_margin']],
    'CM_alpha': [planform_optimization_problem.summary['CM_alpha']],
    'phugoid_damping_ratio': [planform_optimization_problem.summary['phugoid_damping_ratio']],
    'short_period_damping_ratio': [planform_optimization_problem.summary['short_period_damping_ratio']],
    'dutch_roll_frequency': [planform_optimization_problem.summary['dutch_roll_frequency']],
    'dutch_roll_damping_ratio': [planform_optimization_problem.summary['dutch_roll_damping_ratio']],
    'spiral_doubling_time': [planform_optimization_problem.summary['spiral_doubling_time']],
    'run_time': [elapsed_time_stick_fixed]
    }) 
    
    df.to_csv('Output_Single_Engine_Piston_Baseline.csv')
    
    return
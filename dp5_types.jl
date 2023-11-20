struct DP5Report 
    x_final::Float64
    idid::Int64

    num_func_evals::Int64
    num_computed_steps::Int64
    num_accepted_steps::Int64
    num_rejected_steps::Int64
end

struct DP5Options
    # originally in work[1] - work[7]
    uround
    safety_factor
    step_size_selection_one
    step_size_selection_two
    beta
    maximal_step_size
    initial_step_size

    # originally in iwork[1] - iwork[4]
    maximum_allowed_steps
    # choice_coefficients - not being used anymore
    # print error messages - could be boolean
    stiffness_test_activation_step
    # number dense components - will never be used
end
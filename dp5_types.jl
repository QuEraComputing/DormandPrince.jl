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
    print_error_messages
    stiffness_test_activation_step
end

DP5Options() = DP5Options(
    2.3e-16, #uround
    0.9,     # safety_factor
    0.2,     # step_size_selection_one
    10.0,    # step_size_selection_two
    0.04,    # beta
    0.0,     # maximal step size - default to 0.0, later set to xend - x
    0.0,     # initial step size - default to 0.0, trigger hinit later
    100000,  # maximum number of allowed steps 
    true,    # whether or not error messages should be printed
    1000 # stiffness test activated after step J * this number
)
function Base.show(io::IO, mime::MIME"text/plain", solver::DP5Solver)
    tab(n) = " "^(n + get(io, :indent, 0))

    printstyled(io, tab(0), "DP5Solver Object", color=:underline)
    println(io, tab(0), "\n")

    print(io, tab(0), "Current x: ")
    printstyled(io, solver.vars.x, "\n", color=:blue)
    print(io, tab(0), "Current y: ")
    printstyled(io, solver.y, "\n", color=:green)

    println(io, tab(0), "") 

    show(io, mime, solver.options)
end


function Base.show(io::IO, mime::MIME"text/plain", solver::DP8Solver)
    tab(n) = " "^(n + get(io, :indent, 0))

    printstyled(io, tab(0), "DP8Solver Object", color=:underline)
    println(io, tab(0), "\n")

    print(io, tab(0), "Current x: ")
    printstyled(io, solver.vars.x, "\n", color=:blue)
    print(io, tab(0), "Current y: ")
    printstyled(io, solver.y, "\n", color=:green)

    println(io, tab(0), "") 

    show(io, mime, solver.options)
end


function Base.show(io::IO, ::MIME"text/plain", options::Options)

    tab(n) = " "^(n + get(io, :indent, 0))
    println(io, tab(0), "Options: ")

    print(io, tab(4), "Absolute Tolerance: ")
    printstyled(io, options.atol, "\n", color=:light_magenta)
    print(io, tab(4), "Relative Tolerance: ")
    printstyled(io, options.rtol, "\n", color=:light_magenta)

    print(io, tab(4), "uround: ")
    printstyled(io, options.uround, "\n", color=:light_magenta)
    print(io, tab(4), "Safety Factor: ")
    printstyled(io, options.safety_factor, "\n", color=:light_magenta)
    print(io, tab(4), "Step Size Selection 1: ")
    printstyled(io, options.step_size_selection_one, "\n", color=:light_magenta)
    print(io, tab(4), "Step Size Selection 2: ")
    printstyled(io, options.step_size_selection_two, "\n", color=:light_magenta)
    print(io, tab(4), "β: ")
    printstyled(io, options.beta, "\n", color=:light_magenta)
    print(io, tab(4), "Maximal step size: ")
    printstyled(io, options.maximal_step_size, "\n", color=:light_magenta)
    print(io, tab(4), "Initial Step Size: ")
    printstyled(io, options.initial_step_size, "\n", color=:light_magenta)
    print(io , tab(4), "Maximum Allowed Steps: ")
    printstyled(io, options.maximum_allowed_steps, "\n", color=:light_magenta)
    print(io, tab(4), "Print Error Messages: ")
    printstyled(io, options.print_error_messages, "\n", color=:light_magenta)
    print(io, tab(4), "Stiffness Test Activation Step: ")
    printstyled(io, options.stiffness_test_activation_step, "\n", color=:light_magenta)

    
end

function Base.show(io::IO, ::MIME"text/plain", report::Report)
    tab(n) = " "^(n + get(io, :indent, 0))
    printstyled(io, tab(0), "Integration Report", color=:underline)
    println(io, tab(0), "\n")

    print(io, tab(4), "x_final: ")
    printstyled(io, report.x_final, "\n", color=:blue)

    print(io, tab(4), "checks: ")
    if report.checks == INPUT_CHECKS_SUCCESSFUL
        printstyled(io, report.checks, "\n", color=:green)
    else
        printstyled(io, report.checks, "\n", color=:red)
    end

    print(io, tab(4), "idid: ")
    if report.idid == COMPUTATION_SUCCESSFUL
        printstyled(io, report.idid, "\n", color=:green)
    else
        printstyled(io, report.idid, "\n", color=:red)
    end

    print(io, tab(4), "Function Evaluations: ")
    printstyled(io, report.num_func_evals, "\n", color=:light_magenta)
    print(io, tab(4), "Computed Steps: ")
    printstyled(io, report.num_computed_steps, "\n", color=:light_magenta)
    print(io, tab(4), "Accepted Steps: ")
    printstyled(io, report.num_accepted_steps, "\n", color=:green)
    print(io, tab(4), "Rejected Steps: ")
    printstyled(io, report.num_rejected_steps, "\n", color=:red)
end
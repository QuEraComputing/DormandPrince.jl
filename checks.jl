function check_max_allowed_steps(max_allowed_steps, print_error_messages)
    if max_allowed_steps < 0
        print_error_messages ? println("Maximum allowed steps cannot be negative") : nothing
        return true
    else
        return false
    end
end

function check_uround(uround, print_error_messages)
    if (uround <= 1e-35) || (uround >= 1.0) 
        print_error_messages ? println("Unsupported uround, your uround was: ",uround) : nothing
        return true
    else
        return false
    end
end

function check_beta(beta, print_error_messages)
    if beta > 0.2
        print_error_messages ? println("Curious input for beta: ", beta) : nothing
        return true
    else
        return false
    end
end

function check_safety_factor(safety_factor, print_error_messages)
    if (safety_factor >= 1.0) || (safety_factor <= 1e-4)
        print_error_messages ? println("Curious input for safety factor", safety_factor) : nothing
        return true
    else
        return false
    end
end
# include("types.jl")

# make enums for every error type and return that instead of printing errors
function check_max_allowed_steps(options::Options)
    if options.maximum_allowed_steps < 0
        return false
    else
        return true
    end
end

function check_uround(options::Options)
    if (options.uround <= 1e-35) || (options.uround >= 1.0) 
        return false
    else
        return true
    end
end

function check_beta(options::Options)
    if options.beta > 0.2
        return false
    else
        return true
    end
end

function check_safety_factor(options::Options)
    if (options.safety_factor >= 1.0) || (options.safety_factor <= 1e-4)
        return false
    else
        return true
    end
end

function variable_chosen_topology(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)
    z_c = _PM.var(pm,nw)[:z_c] = JuMP.@variable(pm.model,
        base_name="$(nw)_z_c", binary = true, start = 1.0)
    
    report 
end


function add_z_c_to_solution!(pm::_PM.AbstractPowerModel, result::Dict{String,Any})
    solution = get!(result, "solution", Dict{String,Any}())
    sol_nw   = get!(solution, "nw", Dict{String,Any}())

    # data network keys are often Strings: "1","2",...
    for nw_key in sort(collect(keys(pm.data["nw"])))
        # Convert to Int for _PM.var indexing (matches your earlier working approach)
        nw_int = parse(Int, nw_key)

        # Ensure per-nw dict exists
        nw_sol = get!(sol_nw, nw_key, Dict{String,Any}())

        # Only add if the variable exists
        if haskey(_PM.var(pm, nw_int), :z_c)
            nw_sol["z_c"] = JuMP.value(_PM.var(pm, nw_int)[:z_c])
        end
    end

    return result
end
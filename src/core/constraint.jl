function sos_binary_constraint(pm::_PM.AbstractPowerModel)
    #nws = keys(_PM.ref(pm, :nw))
    nws = collect(1:length(pm.data["nw"])) # Vector{String}
    println(nws)
    JuMP.@constraint(pm.model, sum(_PM.var(pm, nw)[:z_c] for nw in nws) == 1)
end

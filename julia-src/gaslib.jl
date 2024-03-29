import XMLDict
import ZipFile


"Parses GasLib data appearing in a compressed ZIP directory."
function parse_gaslib(zip_path::Union{IO,String})
    # Read in the compressed directory and obtain file paths.
    zip_reader = ZipFile.Reader(zip_path)
    file_paths = [x.name for x in zip_reader.files]

    # Parse the topology XML file.
    fid = findfirst(x -> occursin(".net", x), file_paths)
    topology_xml = _parse_xml_file(zip_reader, fid)

    # Parse the compressor XML file.
    fid = findfirst(x -> occursin(".cs", x), file_paths)
    compressor_xml = fid !== nothing ? _parse_xml_file(zip_reader, fid) : Dict()

    # Parse the combined decisions XML file.
    fid = findfirst(x -> occursin(".cdf", x), file_paths)
    cd_xml = fid !== nothing ? _parse_xml_file(zip_reader, fid) : Dict()

    # Parse the nomination XML file(s).
    nomination_xml_dict = _parse_all_nomination_data(file_paths, zip_reader)

    # Compute bulk data from averages of network data.
    density = round(_compute_gaslib_density(topology_xml); digits=4)
    @info "mean density of all sources: $density"
    temperature = round(_compute_gaslib_temperature(topology_xml); digits=4)
    @info "mean temperature of all sources: $temperature"
    molar_mass = round(_compute_gaslib_molar_mass(topology_xml); digits=4)
    @info "mean molar mass of all sources: $molar_mass"
    sound_speed = sqrt(8.314 * temperature * inv(molar_mass))

    # Per the approximation in "Approximating Nonlinear Relationships for
    # Optimal Operation of Natural Gas Transport Networks" by Kazda and Li.
    isentropic_exponent = 1.29 - 5.8824e-4 * (temperature - 273.15)

    # Create a dictionary for all components.
    junctions = _read_gaslib_junctions(topology_xml)
    pipes = _read_gaslib_pipes(topology_xml, junctions, density)
    short_pipes = _read_gaslib_short_pipes(topology_xml, density)
    valves = _read_gaslib_valves(topology_xml, density)
    resistors = _read_gaslib_resistors(topology_xml, density)
    loss_resistors = _read_gaslib_loss_resistors(topology_xml, density)

    control_valves = _read_gaslib_control_valves(topology_xml, density)
    compressors = _read_gaslib_compressors(
        topology_xml,
        compressor_xml,
        temperature,
        8.314 * inv(molar_mass),
        isentropic_exponent,
        density,
    )
    exits = _read_gaslib_sinks(topology_xml, nomination_xml_dict, density)
    entries = _read_gaslib_sources(topology_xml, nomination_xml_dict, density)

    decision_groups = !isempty(cd_xml) ? _read_gaslib_cd(cd_xml) : Dict{String,Any}()

    # Get additional metadata.
    name = topology_xml["information"]["title"]
    gas_specific_gravity = 1000.0 * molar_mass * inv(28.9626)
    compressibility_factor = min(1.0, sound_speed^2 * molar_mass * inv(8.314 * temperature))

    # Build the master data dictionary.
    data = Dict{String,Any}(
        "compressors" => compressors,
        "exits" => exits,
        "nodes" => junctions,
        "entries" => entries,
        "pipes" => pipes,
        "control_valves" => control_valves,
        "resistors" => resistors,
        "loss_resistors" => loss_resistors,
        "short_pipes" => short_pipes,
        "valves" => valves,
        "is_si_units" => 1,
        "temperature" => temperature,
        "name" => name,
        "R" => 8.314,
        "gas_specific_gravity" => gas_specific_gravity,
        "specific_heat_capacity_ratio" => isentropic_exponent,
        "gas_molar_mass" => molar_mass,
        "decision_groups" => decision_groups
    )
    
    
    # Assign nodal IDs in place of string IDs.
    data = _correct_ids(data)

    # Return the dictionary.
    return data
end

function _parse_all_nomination_data(file_paths, zip_reader)
    nomination_data = Dict{String,XMLDict.XMLDictElement}()
    fids = findall(x -> occursin(".scn", x), file_paths)
    for fid in fids 
        nomination_xml = _parse_xml_file(zip_reader, fid)
        nomination_path = file_paths[fid]
        nomination_data[nomination_path] = nomination_xml
    end 
    return nomination_data
end 

_get_report(nomination_info::Vector{Any}) = 
    filter(x -> x["type"] == "report", nomination_info)

function _parse_xml_file(zip_reader, path_index)
    xml_str = ZipFile.read(zip_reader.files[path_index], String)
    return XMLDict.parse_xml(xml_str)
end


function _correct_ids(data::Dict{String,<:Any})
    new_data = deepcopy(data)
    junction_names = sort(collect(keys(data["nodes"])))
    junction_mapping = Dict(k => i for (i, k) in enumerate(junction_names))

    for (junction_name, junction) in data["nodes"]
        i = junction_mapping[junction_name]
        new_data["nodes"][string(i)] = junction
        new_data["nodes"][string(i)]["id"] = i
        delete!(new_data["nodes"], junction_name)
    end

    for node_type in ["entries", "exits"]
        for p in keys(data[node_type])
            node_names = sort(collect(keys(data[node_type][p])))
            node_mapping = Dict(k => i for (i, k) in enumerate(node_names))

            for (node_name, _) in data[node_type][p]
                i = node_mapping[node_name]
                new_data[node_type][p][string(i)] = data[node_type][p][node_name]
                new_data[node_type][p][string(i)]["id"] = i
                new_data[node_type][p][string(i)]["node_id"] = junction_mapping[node_name]
                new_data[node_type][p][string(i)]["name"] = node_name
                delete!(new_data[node_type][p], node_name)
            end
        end 
    end

    compressor_mapping = Dict()
    control_valve_mapping = Dict() 
    valve_mapping = Dict()

    for edge_type in [
        "compressors",
        "pipes",
        "resistors",
        "control_valves",
        "short_pipes",
        "valves",
        "loss_resistors",
    ]
        edge_names = sort(collect(keys(data[edge_type])))
        edge_mapping = Dict(k => a for (a, k) in enumerate(edge_names))
        (edge_type == "compressors") && (compressor_mapping = edge_mapping)
        (edge_type == "control_valves") && (control_valve_mapping = edge_mapping)
        (edge_type == "valves") && (valve_mapping = edge_mapping)

        for (edge_name, edge) in data[edge_type]
            edge_id = edge_mapping[edge_name]
            fr_junction, to_junction = edge["fr_node"], edge["to_node"]
            edge["fr_node"] = junction_mapping[fr_junction]
            edge["to_node"] = junction_mapping[to_junction]
            if (edge_type == "compressors") 
                edge["fuel_node"] = junction_mapping[edge["fuel_node"]]
            end 
            edge["id"] = edge_id
            new_data[edge_type][string(edge_id)] = edge
            delete!(new_data[edge_type], edge_name)
        end
    end

    
    for (dg_key, dg_val) in new_data["decision_groups"]
        for (_, d_val) in dg_val 
            for (_, decision) in d_val 
                for component in decision
                    if (component["component_type"] == "compressor") 
                        old_id = component["id"]
                        component["id"] = compressor_mapping[old_id]
                    end 
                    if (component["component_type"] == "valve") 
                        old_id = component["id"]
                        if !haskey(valve_mapping, old_id)
                            # this is to accomodate bug in GasLib 4197 - valve_22
                            delete!(new_data["decision_groups"], dg_key)
                        else 
                            component["id"] = valve_mapping[old_id]
                        end
                    end 
                    if (component["component_type"] == "control_valve") 
                        old_id = component["id"]
                        component["id"] = control_valve_mapping[old_id]
                    end 
                end 
            end 
        end 
    end 

    return new_data
end


function _compute_node_temperature(node::XMLDict.XMLDictElement)
    if uppercase(node["gasTemperature"][:unit]) == "CELSIUS"
        return 273.15 + parse(Float64, node["gasTemperature"][:value])
    elseif uppercase(node["gasTemperature"][:unit]) == "K"
        return parse(Float64, node["gasTemperature"][:value])
    end
end


function _compute_gaslib_density(topology::XMLDict.XMLDictElement)
    # only source has the field density
    node_types = ["source"]
    node_xml = vcat([get(topology["nodes"], x, []) for x in node_types]...)
    nodes = filter(x -> "normDensity" in collect(keys(x)), node_xml)
    sum_density = sum([parse(Float64, node["normDensity"][:value]) for node in nodes])
    return sum_density * inv(length(nodes)) # Return the mean density.
end


function _compute_gaslib_temperature(topology::XMLDict.XMLDictElement)
    # only source has the field temperature
    node_types = ["source"]
    node_xml = vcat([get(topology["nodes"], x, []) for x in node_types]...)
    nodes = filter(x -> "gasTemperature" in collect(keys(x)), node_xml)
    sum_temperature = sum([_compute_node_temperature(node) for node in nodes])
    return sum_temperature * inv(length(nodes)) # Return the mean temperature.
end


function _compute_gaslib_molar_mass(topology::XMLDict.XMLDictElement)
    # only source has the field molar mass
    node_types = ["source"]
    node_xml = vcat([get(topology["nodes"], x, []) for x in node_types]...)
    nodes = filter(x -> "molarMass" in collect(keys(x)), node_xml)
    sum_molar_mass = sum([parse(Float64, node["molarMass"][:value]) for node in nodes])
    return 1.0e-3 * sum_molar_mass * inv(length(nodes)) # Return the mean.
end


function _get_component_dict(data)
    return data isa Array ? Dict{String,Any}(x[:id] => x for x in data) :
           Dict{String,Any}(x[:id] => x for x in [data])
end


function _get_compressor_entry(
    compressor,
    stations,
    T::Float64,
    R::Float64,
    kappa::Float64,
    density::Float64,
)
    fr_junction, to_junction = compressor[:from], compressor[:to]
    fuel_gas_node = compressor[:fuelGasVertex]
    inlet_p_min = _parse_gaslib_pressure(compressor["pressureInMin"]) 
    outlet_p_max = _parse_gaslib_pressure(compressor["pressureOutMax"])
    flow_min = density * _parse_gaslib_flow(compressor["flowMin"])
    flow_max = density * _parse_gaslib_flow(compressor["flowMax"])
    bypass_required = :internalBypassRequired in keys(compressor) ?
        parse(Int, compressor[:internalBypassRequired]) : 1

    drag_in = "dragFactorIn" in keys(compressor) ? 
        parse(Float64, compressor["dragFactorIn"][:value]) : NaN 
    diameter_in = "diameterIn" in keys(compressor) ?
        _parse_gaslib_length(compressor["diameterIn"]) : NaN 
    pressure_loss_in = "pressureLossIn" in keys(compressor) ? 
        _parse_gaslib_pressure(compressor["pressureLossIn"]) : NaN 
    
    
    drag_out = "dragFactorOut" in keys(compressor) ? 
        parse(Float64, compressor["dragFactorOut"][:value]) : NaN 
    diameter_out = "diameterOut" in keys(compressor) ?
        _parse_gaslib_length(compressor["diameterOut"]) : NaN 
    pressure_loss_out = "pressureLossOut" in keys(compressor) ? 
        _parse_gaslib_pressure(compressor["pressureLossOut"]) : NaN 
    

    c_ratio_min, c_ratio_max = 1.0, outlet_p_max * inv(inlet_p_min)
    
    # Calculate the maximum power.
    exp = kappa * inv(kappa - 1.0)
    H_max = R * T * 1.0 * exp * (c_ratio_max^inv(exp) - 1.0)

    # Assume a worst-case efficiency of 0.1 in the computation of power_max.
    power_max = H_max * abs(flow_max) * inv(0.1)

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "fuel_node" => fuel_gas_node,
        "name" => compressor[:id],
        "internal_bypass_required" => bypass_required,
        "min_inlet_pressure" => inlet_p_min,
        "max_outlet_pressure" => outlet_p_max,
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "drag_in" => drag_in, 
        "drag_out" => drag_out, 
        "diameter_in" => diameter_in, 
        "diameter_out" => diameter_out, 
        "pressure_loss_in" => pressure_loss_in, 
        "pressure_loss_out" => pressure_loss_out,
        "c_ratio_min" => c_ratio_min,
        "c_ratio_max" => c_ratio_max,
        "power_max" => power_max
    )
end


function _get_sink_entry(sink, density::Float64)
    if sink["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", sink["flow"])
        withdrawal_min = density * _parse_gaslib_flow(sink["flow"][min_id])
        max_id = findfirst(x -> x[:bound] == "upper", sink["flow"])
        withdrawal_max = density * _parse_gaslib_flow(sink["flow"][max_id])
    else
        withdrawal_min = density * _parse_gaslib_flow(sink["flow"])
        withdrawal_max = density *_parse_gaslib_flow(sink["flow"])
    end

    if haskey(sink, "pressure")
        if sink["pressure"] isa Array 
            min_id = findfirst(x -> x[:bound] == "lower", sink["pressure"])
            pressure_min = _parse_gaslib_pressure(sink["pressure"][min_id])
            max_id = findfirst(x -> x[:bound] == "upper", sink["pressure"])
            pressure_max = _parse_gaslib_pressure(sink["pressure"][max_id])
        else 
            pressure_min = _parse_gaslib_pressure(sink["pressure"])
            pressure_max = _parse_gaslib_pressure(sink["pressure"])
        end 
    else 
        pressure_min, pressure_max = NaN, NaN
    end 

    return Dict{String,Any}(
        "node_id" => sink[:id],
        "min_withdrawal" => withdrawal_min,
        "max_withdrawal" => withdrawal_max,
        "min_pressure" => pressure_min, 
        "max_pressure" => pressure_max,
    )
end


function _get_junction_entry(junction)
    lat_sym = :geoWGS84Lat in keys(junction) ? :geoWGS84Lat : :x
    lat = parse(Float64, junction[lat_sym])
    lon_sym = :geoWGS84Long in keys(junction) ? :geoWGS84Long : :y
    lon = parse(Float64, junction[lon_sym])

    elevation = parse(Float64, junction["height"][:value])
    # default pressure units in barg (gauge pressure); 1 barg = (1 + 1.01325) * 1e5 Pa
    p_min = _parse_gaslib_pressure(junction["pressureMin"]) 
    p_max = _parse_gaslib_pressure(junction["pressureMax"]) 

    return Dict{String,Any}(
        "node_id" => junction[:id], 
        "name" => junction[:id],
        "x_coord" => lat,
        "y_coord" => lon,
        "min_pressure" => p_min,
        "max_pressure" => p_max,
        "elevation" => elevation,
    )
end


function _parse_gaslib_length(entry)
    if isapprox(parse(Float64, entry[:value]), 0.0, atol = 1e-6)
        return 0.0
    elseif entry[:unit] == "m"
        return parse(Float64, entry[:value])
    elseif entry[:unit] == "km"
        return parse(Float64, entry[:value]) * 1000.0
    elseif entry[:unit] == "mm"
        return parse(Float64, entry[:value]) * inv(1000.0)
    end
end

function _parse_gaslib_pressure(entry)
    if isapprox(parse(Float64, entry[:value]), 0.0, atol = 1e-6)
        return 0.0
    elseif entry[:unit] == "bar" 
        return parse(Float64, entry[:value]) * 1e5
    elseif entry[:unit] == "barg"
        return (parse(Float64, entry[:value]) + 1.01325) * 1e5
    elseif entry[:unit] == "Pa"
        return parse(Float64, entry[:value])
    end
end

function _parse_gaslib_flow(entry)
    if isapprox(parse(Float64, entry[:value]), 0.0, atol = 1e-2)
        return 0.0
    elseif entry[:unit] == "m_cube_per_s" 
        return parse(Float64, entry[:value])
    elseif entry[:unit] == "m_cube_per_hour"
        return parse(Float64, entry[:value]) * inv(3600.0)
    elseif entry[:unit] == "1000m_cube_per_hour"
        return parse(Float64, entry[:value]) * inv(3.6)
    end
end


function _get_pipe_entry(pipe, junctions, density::Float64)
    fr_junction, to_junction = pipe[:from], pipe[:to]
    p_min = min(junctions[fr_junction]["min_pressure"], junctions[to_junction]["min_pressure"])
    p_max = max(junctions[fr_junction]["max_pressure"], junctions[to_junction]["max_pressure"])

    if "pressureMin" in keys(pipe)
        p_min = max(p_min, _parse_gaslib_pressure(pipe["pressureMin"])) 
    end

    if "pressureMax" in keys(pipe)
        p_max = min(p_max, _parse_gaslib_pressure(pipe["pressureMax"]))
    end

    diameter = _parse_gaslib_length(pipe["diameter"])
    length = _parse_gaslib_length(pipe["length"])
    roughness = _parse_gaslib_length(pipe["roughness"])

    # Compute the friction factor as per ``Evaluating Gas Network Capacities'' Eq (2.19).
    friction_factor = (2.0 * log10(diameter * inv(roughness)) + 1.138)^(-2)

    # Determine bidirectionality.
    flow_min = density * _parse_gaslib_flow(pipe["flowMin"])
    flow_max = density * _parse_gaslib_flow(pipe["flowMax"])


    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "diameter" => diameter,
        "length" => length,
        "min_pressure" => p_min,
        "max_pressure" => p_max,
        "friction_factor" => friction_factor,
        "min_flow" => flow_min, 
        "max_flow" => flow_max,
        "name" => pipe[:id]
    )
end


function _get_loss_resistor_entry(loss_resistor, density::Float64)
    fr_junction, to_junction = loss_resistor[:from], loss_resistor[:to]
    flow_min = density * _parse_gaslib_flow(loss_resistor["flowMin"])
    flow_max = density * _parse_gaslib_flow(loss_resistor["flowMax"])
    p_loss = _parse_gaslib_pressure(loss_resistor["pressureLoss"])

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "pressure_loss" => p_loss,
        "name" => loss_resistor[:id]
    )
end


function _get_short_pipe_entry(short_pipe, density::Float64)
    fr_junction, to_junction = short_pipe[:from], short_pipe[:to]
    flow_min = density * _parse_gaslib_flow(short_pipe["flowMin"])
    flow_max = density * _parse_gaslib_flow(short_pipe["flowMax"])

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "name" => short_pipe[:id], 
        "min_flow" => flow_min, 
        "max_flow" => flow_max
    )
end


function _get_source_entry(source, density::Float64)
    if source["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", source["flow"])
        injection_min = density * _parse_gaslib_flow(source["flow"][min_id]) 
        max_id = findfirst(x -> x[:bound] == "upper", source["flow"])
        injection_max = density * _parse_gaslib_flow(source["flow"][max_id]) 
    else
        injection_min = density * _parse_gaslib_flow(source["flow"]) 
        injection_max = density * _parse_gaslib_flow(source["flow"]) 
    end

    if haskey(source, "pressure")
        if source["pressure"] isa Array 
            min_id = findfirst(x -> x[:bound] == "lower", source["pressure"])
            pressure_min = _parse_gaslib_pressure(source["pressure"][min_id])
            max_id = findfirst(x -> x[:bound] == "upper", source["pressure"])
            pressure_max = _parse_gaslib_pressure(source["pressure"][max_id])
        else 
            pressure_min = _parse_gaslib_pressure(source["pressure"])
            pressure_max = _parse_gaslib_pressure(source["pressure"])
        end
    else 
        pressure_min, pressure_max = NaN, NaN 
    end  

    return Dict{String,Any}(
        "node_id" => source[:id],
        "min_injection" => injection_min,
        "max_injection" => injection_max,
        "min_pressure" => pressure_min, 
        "max_pressure" => pressure_max,
    )
end


function _get_resistor_entry(resistor, density::Float64)
    fr_junction, to_junction = resistor[:from], resistor[:to]
    flow_min = density * _parse_gaslib_flow(resistor["flowMin"])
    flow_max = density * _parse_gaslib_flow(resistor["flowMax"])

    diameter = _parse_gaslib_length(resistor["diameter"])
    drag = parse(Float64, resistor["dragFactor"][:value])

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "drag" => drag,
        "diameter" => diameter,
        "name" => resistor[:id]
    )
end


function _get_valve_entry(valve, density::Float64)
    fr_junction, to_junction = valve[:from], valve[:to]
    flow_min = density * _parse_gaslib_flow(valve["flowMin"])
    flow_max = density * _parse_gaslib_flow(valve["flowMax"])
    pressure_differential_max = _parse_gaslib_pressure(valve["pressureDifferentialMax"])

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "max_pressure_differential" => pressure_differential_max,
        "name" => valve[:id]
    )
end

function _get_control_valve_entry(control_valve, density::Float64)
    fr_junction, to_junction = control_valve[:from], control_valve[:to]
    flow_min = density * _parse_gaslib_flow(control_valve["flowMin"])
    flow_max = density * _parse_gaslib_flow(control_valve["flowMax"])

    pressure_differential_min = "pressureDifferentialMin" in keys(control_valve) ?
        _parse_gaslib_pressure(control_valve["pressureDifferentialMin"]) : NaN
    pressure_differential_max = "pressureDifferentialMax" in keys(control_valve) ?
        _parse_gaslib_pressure(control_valve["pressureDifferentialMax"]) : NaN
    
    pressure_in_min = _parse_gaslib_pressure(control_valve["pressureInMin"])
    pressure_out_max = _parse_gaslib_pressure(control_valve["pressureOutMax"])

    drag_in = "dragFactorIn" in keys(control_valve) ? 
        parse(Float64, control_valve["dragFactorIn"][:value]) : NaN 
    diameter_in = "diameterIn" in keys(control_valve) ?
        _parse_gaslib_length(control_valve["diameterIn"]) : NaN 
    pressure_loss_in = "pressureLossIn" in keys(control_valve) ? 
        _parse_gaslib_pressure(control_valve["pressureLossIn"]) : NaN 
    
    drag_out = "dragFactorOut" in keys(control_valve) ? 
        parse(Float64, control_valve["dragFactorOut"][:value]) : NaN 
    diameter_out = "diameterOut" in keys(control_valve) ?
        _parse_gaslib_length(control_valve["diameterOut"]) : NaN 
    pressure_loss_out = "pressureLossOut" in keys(control_valve) ? 
        _parse_gaslib_pressure(control_valve["pressureLossOut"]) : NaN 

    bypass_required = :internalBypassRequired in keys(control_valve) ?
        parse(Int, control_valve[:internalBypassRequired]) : 1

    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "name" => control_valve[:id],
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "internal_bypass_required" => bypass_required,
        "min_pressure_differential" => pressure_differential_min, 
        "max_pressure_differential" => pressure_differential_max,
        "min_pressure_in" => pressure_in_min, 
        "max_pressure_out" => pressure_out_max, 
        "drag_in" => drag_in, 
        "drag_out" => drag_out, 
        "diameter_in" => diameter_in, 
        "diameter_out" => diameter_out, 
        "pressure_loss_in" => pressure_loss_in, 
        "pressure_loss_out" => pressure_loss_out
    )
end

function _get_decision_component_entry(comp, id)
    components_in_decision = []
    component_name = Dict{String,String}(
        "valve" => "valve", 
        "controlValve" => "control_valve", 
        "compressorStation" => "compressor"
    )
    for (component_type, decisions) in comp
        if isa(decisions, Vector) 
            for decision in decisions 
                component = Dict{String,Any}(
                    "component_type" => component_name[component_type], 
                    "id" => decision[:id], 
                )
                (haskey(decision, :value)) && (component["value"] = decision[:value])
                (haskey(decision, :flowDirection)) && (component["flow_direction"] = decision[:flowDirection])
                (haskey(decision, :mode)) && (component["mode"] = decision[:mode])
                push!(components_in_decision, component)
            end 
        else     
            component = Dict{String,Any}(
                    "component_type" => component_name[component_type], 
                    "id" => decisions[:id], 
                )
                (haskey(decisions, :value)) && (component["value"] = decisions[:value])
                (haskey(decisions, :flowDirection)) && (component["flow_direction"] = decisions[:flowDirection])
                (haskey(decisions, :mode)) && (component["mode"] = decisions[:mode])
                push!(components_in_decision, component)
        end 
    end 
    return components_in_decision
end 

function _get_decision_group_entry(dg)
    # this can be a vector of decisions or just one decision 
    components = Dict{String,Any}()
    
    if ~isa(dg["decision"], Vector)
        id = dg["decision"][:id]
        comp = [(x, val) for (x, val) in dg["decision"] if x != :id]
        components[id] = _get_decision_component_entry(comp, id)
    else 
        for decision in dg["decision"]
            id = decision[:id]
            comp = [(x, val) for (x, val) in decision if x != :id]
            components[id] = _get_decision_component_entry(comp, id)
        end 
    end 
    components_with_new_ids = Dict{String,Any}(
        string(i) => components[k] for (i, k) in enumerate(keys(components)))
    return Dict{String,Any}("decisions" => components_with_new_ids)
end 

function _read_gaslib_compressors(
    topology,
    compressor_stations,
    T::Float64,
    R::Float64,
    kappa::Float64,
    density::Float64,
)
    if "compressorStation" in keys(topology["connections"])
        compressors = _get_component_dict(topology["connections"]["compressorStation"])
        stations = "compressorStation" in keys(compressor_stations) ?
            compressor_stations["compressorStation"] : nothing
        return Dict{String,Any}(
            i => _get_compressor_entry(x, stations, T, R, kappa, density)
            for (i, x) in compressors
        )
    else
        return Dict{String,Any}()
    end
end


function _read_gaslib_sinks(
    topology::XMLDict.XMLDictElement,
    nomination_dict::Dict{String,XMLDict.XMLDictElement},
    density::Float64,
)   
    sinks = Dict{String,Any}() 
    for (p, nominations) in nomination_dict
        node_xml = _get_component_dict(get(nominations["scenario"], "node", []))

        sink_ids = [x[:id] for x in get(topology["nodes"], "sink", [])]
        sink_xml = filter(x -> x.second[:id] in sink_ids, node_xml)
        sink_data = Dict{String,Any}(i => _get_sink_entry(x, density) for (i, x) in sink_xml)

        sinks[p] = sink_data
    end 
    return sinks
end


function _read_gaslib_junctions(topology::XMLDict.XMLDictElement)
    node_types = ["innode", "sink", "source"]
    node_xml = vcat([get(topology["nodes"], x, []) for x in node_types]...)
    return Dict{String,Any}(x[:id] => _get_junction_entry(x) for x in node_xml)
end


function _read_gaslib_pipes(
    topology::XMLDict.XMLDictElement,
    junctions::Dict{String,<:Any},
    density::Float64,
)
    pipe_xml = _get_component_dict(get(topology["connections"], "pipe", []))
    return Dict{String,Any}(
        i => _get_pipe_entry(x, junctions, density) for (i, x) in pipe_xml
    )
end


function _read_gaslib_loss_resistors(topology::XMLDict.XMLDictElement, density::Float64)
    loss_resistor_xml = _get_component_dict(get(topology["connections"], "resistor", []))
    loss_resistor_xml = filter(x -> "pressureLoss" in collect(keys(x.second)), loss_resistor_xml)
    return Dict{String,Any}(
        i => _get_loss_resistor_entry(x, density) for (i, x) in loss_resistor_xml
    )
end


function _read_gaslib_sources(
    topology::XMLDict.XMLDictElement,
    nominations_dict::Dict{String,XMLDict.XMLDictElement},
    density::Float64,
)   
    sources = Dict{String,Any}()
    for (p, nominations) in nominations_dict
        node_xml = _get_component_dict(get(nominations["scenario"], "node", []))

        # Collect source nodes with positive injections.
        source_ids = [x[:id] for x in get(topology["nodes"], "source", [])]
        source_xml = filter(x -> x.second[:id] in source_ids, node_xml)
        source_data = Dict{String,Any}(i => _get_source_entry(x, density) for (i, x) in source_xml)
    
        sources[p] = source_data
    end 
    return sources
end

function _read_gaslib_control_valves(topology::XMLDict.XMLDictElement, density::Float64)
    control_valve_xml = _get_component_dict(get(topology["connections"], "controlValve", []))
    return Dict{String,Any}(
        i => _get_control_valve_entry(x, density) for (i, x) in control_valve_xml
    )
end

function _read_gaslib_resistors(topology::XMLDict.XMLDictElement, density::Float64)
    resistor_xml = _get_component_dict(get(topology["connections"], "resistor", []))
    resistor_xml = filter(x -> "dragFactor" in collect(keys(x.second)), resistor_xml)
    resistor_xml = filter(x -> parse(Float64, x.second["dragFactor"][:value]) > 0.0, resistor_xml)
    return Dict{String,Any}(i => _get_resistor_entry(x, density) for (i, x) in resistor_xml)
end


function _read_gaslib_short_pipes(topology::XMLDict.XMLDictElement, density::Float64)
    # Parse all short pipes in the network.
    short_pipe_xml = _get_component_dict(get(topology["connections"], "shortPipe", []))
    short_pipes = Dict{String,Any}(
        i => _get_short_pipe_entry(x, density) for (i, x) in short_pipe_xml
    )

    # Resistors with drag equal to zero should also be considered short pipes.
    resistor_xml = _get_component_dict(get(topology["connections"], "resistor", []))
    resistor_xml = filter(x -> "dragFactor" in collect(keys(x.second)), resistor_xml)
    resistor_xml = filter(x -> 
        isapprox(parse(Float64, x.second["dragFactor"][:value]), 0.0, atol=0.001), 
        resistor_xml)
    resistors = Dict{String,Any}(i => _get_short_pipe_entry(x, density) for (i, x) in resistor_xml)

    # Return the merged dictionary of short pipes and resistors.
    return merge(short_pipes, resistors)
end


function _read_gaslib_valves(topology::XMLDict.XMLDictElement, density::Float64)
    valve_xml = _get_component_dict(get(topology["connections"], "valve", []))
    return Dict{String,Any}(i => _get_valve_entry(x, density) for (i, x) in valve_xml)
end

function _read_gaslib_cd(cd::XMLDict.XMLDictElement)
    cd_xml = _get_component_dict(get(cd, "decisionGroup", []))
    decision_groups = Dict{String,Any}(i => _get_decision_group_entry(x) for (i, x) in cd_xml)
    return Dict{String,Any}(string(i) => decision_groups[k] for (i, k) in enumerate(keys(decision_groups)))
end 
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
    fids = findall(x -> occursin(".scn", x), file_paths)

    # TODO: If length(fids) > 1, treat the input as a multinetwork.
    if length(fids) > 1 # If multiple scenarios are defined.
        # Parse and use the last nominations file when sorted by name.
        nomination_path = sort(file_paths[fids])[end]
        fid = findfirst(x -> occursin(nomination_path, x), file_paths)
        nomination_xml = _parse_xml_file(zip_reader, fid)

        # Print a warning message stating that the above file is being used.
        @warn "Multiple nomination file paths found " *
            "in GasLib data. Selecting last nomination file " *
            "(i.e., \"$(file_paths[fid])\") " *
            "after sorting by name."
    
    else # If only one nominations scenario is defined...
        # Parse and use the only nominations file available.
        nomination_xml = _parse_xml_file(zip_reader, fids[end])
    end

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

    manual_control_valves = _read_gaslib_manual_control_valves(topology_xml, density)
    automated_control_valves = _read_gaslib_automated_control_valves(topology_xml, density)
    compressors = _read_gaslib_compressors(
        topology_xml,
        compressor_xml,
        temperature,
        8.314 * inv(molar_mass),
        isentropic_exponent,
        density,
    )
    deliveries = _read_gaslib_deliveries(topology_xml, nomination_xml, density)
    receipts = _read_gaslib_receipts(topology_xml, nomination_xml, density)

    combined_decisions = _read_gaslib_cd(cd_xml)
    
    return (junctions = junctions, pipes = pipes, short_pipes=short_pipes, 
        valves = valves, resistors = resistors, loss_resistors = loss_resistors, 
        manual_control_valves = manual_control_valves, 
        automated_control_valves = automated_control_valves, 
        compressors = compressors, deliveries = deliveries, 
        receipts = receipts, combined_decisions = cd_xml)


    # Add auxiliary nodes for bidirectional control_valves.
    # _add_auxiliary_junctions!(junctions, control_valves)

    # Get additional metadata.
    name = topology_xml["information"]["title"]
    gas_specific_gravity = 1000.0 * molar_mass * inv(28.9626)
    compressibility_factor = min(1.0, sound_speed^2 * molar_mass * inv(8.314 * temperature))

    # Build the master data dictionary.
    data = Dict{String,Any}(
        "compressor" => compressors,
        "delivery" => deliveries,
        "junction" => junctions,
        "receipt" => receipts,
        "pipe" => pipes,
        "control_valve" => control_valves,
        "resistor" => resistors,
        "loss_resistor" => loss_resistors,
        "short_pipe" => short_pipes,
        "valve" => valves,
        "is_si_units" => 1,
        "is_english_units" => 0,
        "is_per_unit" => 0,
        "sound_speed" => sound_speed,
        "temperature" => temperature,
        "name" => name,
        "R" => 8.314,
        "compressibility_factor" => compressibility_factor,
        "gas_specific_gravity" => gas_specific_gravity,
        "standard_density" => density,
        "specific_heat_capacity_ratio" => isentropic_exponent,
        "gas_molar_mass" => molar_mass,
    )

    # Assign nodal IDs in place of string IDs.
    data = _correct_ids(data)

    # Return the dictionary.
    return data
end


function _parse_xml_file(zip_reader, path_index)
    xml_str = ZipFile.read(zip_reader.files[path_index], String)
    return XMLDict.parse_xml(xml_str)
end


function _correct_ids(data::Dict{String,<:Any})
    new_data = deepcopy(data)
    junction_names = sort(collect(keys(data["junction"])))
    junction_mapping = Dict(k => i for (i, k) in enumerate(junction_names))

    for (junction_name, junction) in data["junction"]
        i = junction_mapping[junction_name]
        new_data["junction"][string(i)] = junction
        new_data["junction"][string(i)]["id"] = i
        new_data["junction"][string(i)]["index"] = i
        delete!(new_data["junction"], junction_name)
    end

    for node_type in ["delivery", "receipt"]
        node_names = sort(collect(keys(data[node_type])))
        node_mapping = Dict(k => i for (i, k) in enumerate(node_names))

        for (node_name, _) in data[node_type]
            i = node_mapping[node_name]
            new_data[node_type][string(i)] = data[node_type][node_name]
            new_data[node_type][string(i)]["id"] = i
            new_data[node_type][string(i)]["index"] = i
            new_data[node_type][string(i)]["junction_id"] = junction_mapping[node_name]
            delete!(new_data[node_type], node_name)
        end
    end

    for edge_type in [
        "compressor",
        "pipe",
        "resistor",
        "control_valve",
        "short_pipe",
        "valve",
        "loss_resistor",
    ]
        edge_names = sort(collect(keys(data[edge_type])))
        edge_mapping = Dict(k => a for (a, k) in enumerate(edge_names))

        for (edge_name, edge) in data[edge_type]
            edge_id = edge_mapping[edge_name]
            fr_junction, to_junction = edge["fr_junction"], edge["to_junction"]
            edge["fr_junction"] = junction_mapping[fr_junction]
            edge["to_junction"] = junction_mapping[to_junction]
            edge["id"] = edge["index"] = edge_id
            new_data[edge_type][string(edge_id)] = edge
            delete!(new_data[edge_type], edge_name)
        end
    end

    return new_data
end


function _add_auxiliary_junctions!(junctions, control_valves)
    new_junctions = Dict{String,Any}()
    new_control_valves = Dict{String,Any}()

    for (a, control_valve) in control_valves
        if control_valve["bypass_required"] == 1
            fr_junction = junctions[control_valve["fr_junction"]]
            to_junction = junctions[control_valve["to_junction"]]

            junction_aux_name = a * "_aux_junction"
            junction_aux = deepcopy(fr_junction)

            control_valve_reverse_name = a * "_reverse"
            control_valve_reverse = deepcopy(control_valve)
            control_valve_reverse["fr_junction"] = junction_aux_name
            control_valve_reverse["flow_min"] = -control_valve["flow_max"]
            control_valve_reverse["flow_max"] = -control_valve["flow_min"]
            control_valve["to_junction"] = junction_aux_name

            push!(new_junctions, junction_aux_name => junction_aux)
            push!(new_control_valves, control_valve_reverse_name => control_valve_reverse)
        end
    end

    # junctions = _IM.update_data!(junctions, new_junctions)
    # control_valves = _IM.update_data!(control_valves, new_control_valves)
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
        "name" => compressor[:id],
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


function _get_delivery_entry(delivery, density::Float64)
    if delivery["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", delivery["flow"])
        withdrawal_min = density * _parse_gaslib_flow(delivery["flow"][min_id])
        max_id = findfirst(x -> x[:bound] == "upper", delivery["flow"])
        withdrawal_max = density * _parse_gaslib_flow(delivery["flow"][max_id])
    else
        withdrawal_min = density * _parse_gaslib_flow(delivery["flow"])
        withdrawal_max = density *_parse_gaslib_flow(delivery["flow"])
    end

    if haskey(delivery, "pressure")
        if delivery["pressure"] isa Array 
            min_id = findfirst(x -> x[:bound] == "lower", delivery["pressure"])
            pressure_min = _parse_gaslib_pressure(delivery["pressure"][min_id])
            max_id = findfirst(x -> x[:bound] == "upper", delivery["pressure"])
            pressure_max = _parse_gaslib_pressure(delivery["pressure"][max_id])
        else 
            pressure_min = _parse_gaslib_pressure(delivery["pressure"])
            pressure_max = _parse_gaslib_pressure(delivery["pressure"])
        end 
    else 
        pressure_min, pressure_max = NaN, NaN
    end 

    is_dispatchable = withdrawal_min != withdrawal_max ? 1 : 0

    return Dict{String,Any}(
        "node_id" => delivery[:id],
        "min_withdrawal" => withdrawal_min,
        "max_withdrawal" => withdrawal_max,
        "nominal_withdrawal" => withdrawal_max,
        "min_pressure" => pressure_min, 
        "max_pressure" => pressure_max,
        "is_dispatchable" => is_dispatchable
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
        "node_name" => junction[:id],
        "x_coord" => lat,
        "y_coord" => lon,
        "min_pressure" => p_min,
        "max_pressure" => p_max,
        "elevation" => elevation,
        "slack_bool" => 0
    )
end


function _parse_gaslib_length(entry)
    if isapprox(parse(Float64, entry[:value]), 0.0, atol = 1e-2)
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
    if isapprox(parse(Float64, entry[:value]), 0.0, atol = 1e-2)
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
    flow_max = density * _parse_gaslib_flow(loss_resistor["flowMin"])
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
        "max_flow" => flow_max, 
        "name" => short_pipe[:id]
    )
end


function _get_receipt_entry(receipt, density::Float64)
    if receipt["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", receipt["flow"])
        injection_min = density * _parse_gaslib_flow(receipt["flow"][min_id]) 
        max_id = findfirst(x -> x[:bound] == "upper", receipt["flow"])
        injection_max = density * density * _parse_gaslib_flow(receipt["flow"][max_id]) 
    else
        injection_min = density * _parse_gaslib_flow(receipt["flow"]) 
        injection_max = density * _parse_gaslib_flow(receipt["flow"]) 
    end

    if haskey(receipt, "pressure")
        if receipt["pressure"] isa Array 
            min_id = findfirst(x -> x[:bound] == "lower", receipt["pressure"])
            pressure_min = _parse_gaslib_pressure(receipt["pressure"][min_id])
            max_id = findfirst(x -> x[:bound] == "upper", receipt["pressure"])
            pressure_max = _parse_gaslib_pressure(receipt["pressure"][max_id])
        else 
            pressure_min = _parse_gaslib_pressure(receipt["pressure"])
            pressure_max = _parse_gaslib_pressure(receipt["pressure"])
        end
    else 
        pressure_min, pressure_max = NaN, NaN 
    end  

    is_dispatchable = injection_min != injection_max ? 1 : 0

    return Dict{String,Any}(
        "node_id" => receipt[:id],
        "min_injection" => injection_min,
        "max_injection" => injection_max,
        "nominal_injection" => injection_max,
        "min_pressure" => pressure_min, 
        "max_pressure" => pressure_max,
        "is_dispatchable" => is_dispatchable
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


function _get_manual_control_valve_entry(control_valve, density::Float64)
    fr_junction, to_junction = control_valve[:from], control_valve[:to]
    flow_min = density * _parse_gaslib_flow(control_valve["flowMin"])
    flow_max = density * _parse_gaslib_flow(control_valve["flowMax"])

    pressure_set = "pressureSet" in keys(control_valve) ? 
        _parse_gaslib_pressure(control_valve["pressureSet"]) : NaN 
    
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


    return Dict{String,Any}(
        "fr_node" => fr_junction,
        "to_node" => to_junction,
        "name" => control_valve[:id],
        "min_flow" => flow_min,
        "max_flow" => flow_max,
        "pressure_set" => pressure_set, 
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

function _get_automated_control_valve_entry(control_valve, density::Float64)
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
        "bypass_required" => bypass_required,
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

function _get_decision_group_entry(dg)
    components = Vector{Tuple{Any,Any}}()
    on_off_status = Dict{String,Any}()
    flow_direction = Dict{String,Any}()
    mode = Dict{String,Any}()
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


function _read_gaslib_deliveries(
    topology::XMLDict.XMLDictElement,
    nominations::XMLDict.XMLDictElement,
    density::Float64,
)
    node_xml = _get_component_dict(get(nominations["scenario"], "node", []))

    # Collect source nodes with negative injections.
    source_ids = [x[:id] for x in get(topology["nodes"], "source", [])]
    source_xml = filter(x -> x.second[:id] in source_ids, node_xml)
    source_data = Dict{String,Any}(i => _get_delivery_entry(x, density) for (i, x) in source_xml)
    source_data = filter(x -> x.second["max_withdrawal"] < 0.0, source_data)

    # Collect sink nodes with positive withdrawals.
    sink_ids = [x[:id] for x in get(topology["nodes"], "sink", [])]
    sink_xml = filter(x -> x.second[:id] in sink_ids, node_xml)
    sink_data = Dict{String,Any}(i => _get_delivery_entry(x, density) for (i, x) in sink_xml)
    sink_data = filter(x -> x.second["min_withdrawal"] > 0.0, sink_data)

    # For sink nodes with negative injections, negate the values.
    for (i, source) in source_data
        source["max_withdrawal"] *= -1.0
        source["max_withdrawal"] *= -1.0
        source["nominal_withdrawal"] *= -1.0
    end

    return merge(source_data, sink_data)
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


function _read_gaslib_receipts(
    topology::XMLDict.XMLDictElement,
    nominations::XMLDict.XMLDictElement,
    density::Float64,
)
    node_xml = _get_component_dict(get(nominations["scenario"], "node", []))

    # Collect source nodes with positive injections.
    source_ids = [x[:id] for x in get(topology["nodes"], "source", [])]
    source_xml = filter(x -> x.second[:id] in source_ids, node_xml)
    source_data = Dict{String,Any}(i => _get_receipt_entry(x, density) for (i, x) in source_xml)
    source_data = filter(x -> x.second["min_injection"] > 0.0, source_data)

    # Collect sink nodes with negative withdrawals.
    sink_ids = [x[:id] for x in get(topology["nodes"], "sink", [])]
    sink_xml = filter(x -> x.second[:id] in sink_ids, node_xml)
    sink_data = Dict{String,Any}(i => _get_receipt_entry(x, density) for (i, x) in sink_xml)
    sink_data = filter(x -> x.second["max_injection"] < 0.0, sink_data)

    # For sink nodes with negative withdrawals, negate the values.
    for (_, sink) in sink_data
        sink["min_injection"] *= -1.0
        sink["max_injection"] *= -1.0
        sink["nominal_injection"] *= -1.0
    end

    return merge(source_data, sink_data)
end

function _read_gaslib_manual_control_valves(topology::XMLDict.XMLDictElement, density::Float64)
    control_valve_xml = _get_component_dict(get(topology["connections"], "controlValve", []))
    control_valve_xml = filter(x -> "pressureSet" in keys(x.second), control_valve_xml)
    return Dict{String,Any}(
        i => _get_manual_control_valve_entry(x, density) for (i, x) in control_valve_xml
    )
end

function _read_gaslib_automated_control_valves(topology::XMLDict.XMLDictElement, density::Float64)
    control_valve_xml = _get_component_dict(get(topology["connections"], "controlValve", []))
    control_valve_xml = filter(x -> !("pressureSet" in keys(x.second)), control_valve_xml)
    return Dict{String,Any}(
        i => _get_automated_control_valve_entry(x, density) for (i, x) in control_valve_xml
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
    return Dict{String,Any}(i => _get_decision_group_entry(x) for (i, x) in cd_xml)
end 
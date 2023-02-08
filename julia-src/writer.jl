using JSON

get_params(; T::Float64 = 288.70599999999996) = 
    Dict{String,Any}(
    "Temperature (K):" => T,
    "Gas specific gravity (G):" => 0.6,
    "Specific heat capacity ratio" => 1.4,
    "units (SI = 0, standard = 1)" => 0.0
)

function write_params_data(data, outputfolder, zip_file) 
    output_folder = outputfolder * split(zip_file, ".")[1]
    (!isdir(output_folder)) && (mkdir(output_folder)) 

    open(output_folder * "/params.json", "w") do f 
        JSON.print(f, Dict("params" => get_params(T = data["temperature"])), 2)
    end
end 

function write_network_data(data, outputfolder, zip_file) 
    output_folder = outputfolder * split(zip_file, ".")[1]
    (!isdir(output_folder)) && (mkdir(output_folder)) 


    network_data = Dict{String,Any}(
        "nodes" => Dict{String,Any}(), 
        "pipes" => Dict{String,Any}(), 
        "compressors" => Dict{String,Any}(),
        "short_pipes" => Dict{String,Any}(), 
        "control_valves" => Dict{String,Any}(), 
        "valves" => Dict{String,Any}(), 
        "resistors" => Dict{String,Any}(), 
        "loss_resistors" => Dict{String,Any}(), 
        "entries" => Dict{String,Any}(), 
        "exits" => Dict{String,Any}()
    )

    slack_node_data = get_slack_nodes(data)

    for (i, node) in data["nodes"]
        network_data["nodes"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => node["name"],
            "x_coord" => node["x_coord"],
            "y_coord" => node["y_coord"], 
            "elevation" => node["elevation"],
            "min_pressure" => node["min_pressure"],
            "max_pressure" => node["max_pressure"]
        )
    end

    for (i, pipe) in data["pipes"]
        network_data["pipes"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => pipe["name"],
            "fr_node" => pipe["fr_node"],
            "to_node" => pipe["to_node"],
            "diameter" => pipe["diameter"],
            "length" => pipe["length"],
            "friction_factor" => pipe["friction_factor"],
            "min_pressure" => pipe["min_pressure"],
            "max_pressure" => pipe["max_pressure"],
            "min_flow" => round(pipe["min_flow"], digits=4),
            "max_flow" => round(pipe["max_flow"], digits=4)
        )
    end 

    for (i, compressor) in get(data, "compressors", []) 
        network_data["compressors"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => compressor["name"],
            "fr_node" => compressor["fr_node"], 
            "to_node" => compressor["to_node"], 
            "fuel_node" => compressor["fuel_node"],
            "internal_bypass_required" => compressor["internal_bypass_required"],
            "min_flow" => round(compressor["min_flow"], digits=4),
            "max_flow" => round(compressor["max_flow"], digits=4),
            "min_c_ratio" => compressor["c_ratio_min"],
            "max_c_ratio" => compressor["c_ratio_max"],
            "max_power" => compressor["power_max"]
        )
    end

    for (i, short_pipe) in get(data, "short_pipes", [])
        network_data["short_pipes"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => short_pipe["name"], 
            "fr_node" => short_pipe["fr_node"],
            "to_node" => short_pipe["to_node"],
            "min_flow" => round(short_pipe["min_flow"], digits=4), 
            "max_flow" => round(short_pipe["max_flow"], digits=4)
        )
    end 

    
    for (i, valve) in get(data, "valves", []) 
        network_data["valves"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => valve["name"], 
            "fr_node" => valve["fr_node"],
            "to_node" => valve["to_node"],
            "min_flow" => round(valve["min_flow"], digits=4), 
            "max_flow" => round(valve["max_flow"], digits=4), 
            "max_pressure_differential" => valve["max_pressure_differential"]
        )
    end 

    for (i, control_valve) in get(data, "control_valves", []) 
        network_data["control_valves"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => control_valve["name"], 
            "fr_node" => control_valve["fr_node"],
            "to_node" => control_valve["to_node"],
            "min_flow" => round(control_valve["min_flow"], digits=4), 
            "max_flow" => round(control_valve["max_flow"], digits=4), 
            "internal_bypass_required" => control_valve["internal_bypass_required"],
            "min_pressure_differential" => control_valve["min_pressure_differential"],
            "max_pressure_differential" => control_valve["max_pressure_differential"]
        )
    end 

    for (i, resistor) in get(data, "resistors", [])
        network_data["resistors"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => resistor["name"], 
            "fr_node" => resistor["fr_node"],
            "to_node" => resistor["to_node"],
            "min_flow" => round(resistor["min_flow"], digits=4),
            "max_flow" => round(resistor["max_flow"], digits=4),
            "drag" => resistor["drag"],
            "diameter" => resistor["diameter"]
        )
    end 

    for (i, resistor) in get(data, "loss_resistors", [])
        network_data["loss_resistors"][i] = Dict{String,Any}(
            "id" => parse(Int, i),
            "name" => resistor["name"], 
            "fr_node" => resistor["fr_node"],
            "to_node" => resistor["to_node"],
            "min_flow" => round(resistor["min_flow"], digits=4),
            "max_flow" => round(resistor["max_flow"], digits=4),
            "pressure_loss" => resistor["pressure_loss"]
        )
    end 

    all_entries = get(data, "entries", []) 
    if !isempty(all_entries)
        entries = all_entries |> first |> last
        for (i, entry) in entries
            network_data["entries"][i] = Dict{String,Any}(
                "id" => parse(Int, i), 
                "node_id" => entry["node_id"], 
                "name" => entry["name"]
            )
        end 
    end

    all_exits = get(data, "exits", []) 
    if !isempty(all_exits)
        exits = all_exits |> first |> last
        for (i, exit) in exits 
            network_data["exits"][i] = Dict{String,Any}(
                "id" => parse(Int, i), 
                "node_id" => exit["node_id"],
                "name" => exit["name"]
            )
        end 
    end 

    open(output_folder * "/network.json", "w") do f 
        JSON.print(f, network_data, 2)
    end


    open(output_folder * "/slack_nodes.json", "w") do f 
        JSON.print(f, slack_node_data, 2)
    end

end 

function get_slack_nodes(data::Dict)

    withdrawal = Dict{String,Any}( 
        k => Dict{String,Any}(
            i => 0.0 for (i, _) in data["nodes"]
        ) for k in data["entries"] |> keys)
    
    slack_node = Dict{String,Any}()
    
    for p in data["entries"] |> keys 
        k = split(p, ['/', '.'])[end-1]
        for (_, entry) in data["entries"][p]
            node = string(entry["node_id"])
            withdrawal[p][node] -= entry["max_injection"]
        end 
    
        for (_, exit) in data["exits"][p]
            node = string(exit["node_id"])
            withdrawal[p][node] += exit["max_withdrawal"] 
        end 

        withdrawal_values = [v for (_, v) in withdrawal[p]]
        max_injection = minimum(withdrawal_values)
        nodes_with_max_injection = sort([k for (k,v) in withdrawal[p] if v == max_injection])
        slack_node[k] = nodes_with_max_injection |> first
    end 

    return slack_node
end 

function write_decision_group_data(data, outputfolder, zip_file) 
    output_folder = outputfolder * split(zip_file, ".")[1]
    (!isdir(output_folder)) && (mkdir(output_folder)) 

    if !isempty(data["decision_groups"])
        open(output_folder * "/decision_groups.json", "w") do f 
            JSON.print(f, Dict("decision_groups" => data["decision_groups"]), 2)
        end
    end
end 

function write_nomination_data(data, outputfolder, zip_file) 
    output_folder = outputfolder * split(zip_file, ".")[1]
    (!isdir(output_folder)) && (mkdir(output_folder)) 

    nomination_data = Dict{String,Any}(
        split(p, ['/', '.'])[end-1] => 
        Dict{String,Any}() 
        for p in keys(data["entries"])
    )

    for (p, entries) in get(data, "entries", [])
        k = split(p, ['/', '.'])[end-1]
        nomination_data[k] = Dict{String,Any}(
            "entry_nominations" => Dict{String,Any}(), 
            "exit_nominations" => Dict{String,Any}()
        )
        for (i, entry) in entries
            nomination_data[k]["entry_nominations"][i] = Dict{String,Any}(
                "cost" => 1.0,
                "max_injection" => entry["max_injection"],
                "min_injection" => entry["min_injection"]
            )
        end 
    end 

    for (p, exits) in get(data, "exits", [])
        k = split(p, ['/', '.'])[end-1]
        for (i, exit) in exits 
            nomination_data[k]["exit_nominations"][i] = Dict{String,Any}(
                "cost" => 1.0, 
                "max_withdrawal" => exit["max_withdrawal"],
                "min_withdrawal" => exit["min_withdrawal"]
            )
        end
    end 

    open(output_folder * "/nominations.json", "w") do f 
        JSON.print(f, nomination_data, 2)
    end
end 
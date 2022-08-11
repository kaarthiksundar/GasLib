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

function write_network_data(data, outputfolder, zip_file; bc_flag::Bool=false) 
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
        "receipts" => Dict{String,Any}(), 
        "deliveries" => Dict{String,Any}()
    )

    slack_pressure_data = Dict{String,Any}(
        "slack_pressures" => Dict{String,Any}()
    )

    bc = Dict{String,Any}(
        "boundary_compressor" => Dict{String,Any}(),
        "boundary_nonslack_flow" => Dict{String,Any}(),
        "boundary_pslack" => Dict{String,Any}()

    )

    withdrawal = Dict{String,Any}()

    for (i, _) in data["nodes"]
        withdrawal[i] = 0.0
    end 

    for (_, receipt) in data["receipts"]
        node = string(receipt["node_id"])
        withdrawal[node] -= receipt["nominal_injection"]
    end 

    for (_, delivery) in data["deliveries"]
        node = string(delivery["node_id"])
        withdrawal[node] += delivery["nominal_withdrawal"] 
    end 

    withdrawal_values = [v for (_, v) in withdrawal]
    max_injection = minimum(withdrawal_values)
    nodes_with_max_injection = sort([k for (k,v) in withdrawal if v == max_injection])
    slack_node = nodes_with_max_injection[1]

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
        
        if i == string(slack_node)
            network_data["nodes"][i]["slack_bool"] = 1
            slack_pressure_data["slack_pressures"][i] = 5000000
            bc["boundary_pslack"][i] = 5000000
        else 
            network_data["nodes"][i]["slack_bool"] = 0
            if !(withdrawal[i] == 0)
                bc["boundary_nonslack_flow"][i] = withdrawal[i]
            end
        end
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
        bc["boundary_compressor"][i] = Dict{String,Any}(
            "control_type" => 0, "value" => 1.5
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

    bc["boundary_valve"] = Dict{String,Any}("on" => [], "off" => [])
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
        push!(bc["boundary_valve"]["on"], parse(Int, i))
    end 

    bc["boundary_control_valve"] = Dict{String,Any}("on" => [], "off" => [])
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
        push!(bc["boundary_control_valve"]["on"], parse(Int, i))
        bc["boundary_control_valve"][i] = Dict{String,Any}(
            "control_type" => 0, "value" => 1.0
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

    for (i, receipt) in get(data, "receipts", []) 
        network_data["receipts"][i] = Dict{String,Any}(
            "id" => parse(Int, i), 
            "node_id" => receipt["node_id"], 
            "name" => receipt["name"]
        )
    end 

    for (i, delivery) in get(data, "deliveries", []) 
        network_data["deliveries"][i] = Dict{String,Any}(
            "id" => parse(Int, i), 
            "node_id" => delivery["node_id"],
            "name" => delivery["name"]
        )
    end 

    open(output_folder * "/network.json", "w") do f 
        JSON.print(f, network_data, 2)
    end


    open(output_folder * "/slack_pressures.json", "w") do f 
        JSON.print(f, slack_pressure_data, 2)
    end

    if bc_flag
        open(output_folder * "/bc.json", "w") do f 
            JSON.print(f, bc, 2)
        end
    end
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
        "deliveries" => Dict{String,Any}(), 
        "receipts" => Dict{String,Any}()
    )

    for (i, receipt) in get(data, "receipts", []) 
        nomination_data["receipts"][i] = Dict{String,Any}(
            "cost" => 1.0,
            "max_injection" => round(receipt["max_injection"], digits=4),
            "min_injection" => round(receipt["min_injection"], digits=4)
        )
    end 

    for (i, delivery) in get(data, "deliveries", []) 
        nomination_data["deliveries"][i] = Dict{String,Any}(
            "cost" => 1.0, 
            "max_withdrawal" => round(delivery["max_withdrawal"], digits=4),
            "min_withdrawal" => round(delivery["min_withdrawal"], digits=4)
        )
    end 

    open(output_folder * "/nominations.json", "w") do f 
        JSON.print(f, nomination_data, 2)
    end
end 
using ArgParse
include("gaslib.jl")
include("writer.jl")

function parse_cli_args(args)
    s = ArgParseSettings(description = "GasLib")
    @add_arg_table! s begin
        "--datafolder"               
            arg_type = String 
            default = "./data/"
            help = "folder where GasLib data is available"
        "--file"
            help = "GasLib zipped directory" 
            arg_type = String
            default = "GasLib-11.zip"
        "--outputfolder"
            arg_type = String
            default = "./json/"
            help = "folder to save the output json"
    end
    return parse_args(s) # the result is a Dict{String,Any}
end

function main(ARGS)
    args = parse_cli_args(ARGS)
    data = parse_gaslib(args["datafolder"] * args["file"])
    write_network_data(data, args["outputfolder"], args["file"])
    write_params_data(data, args["outputfolder"], args["file"])
    write_decision_group_data(data, args["outputfolder"], args["file"])
    write_nomination_data(data, args["outputfolder"], args["file"])
    return data
end 

main(ARGS)
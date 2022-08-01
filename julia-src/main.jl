using ArgParse
include("gaslib.jl")

function parse_cli_args(args)
    s = ArgParseSettings(description = "GasLib")
    @add_arg_table! s begin
        "--datafolder"               
            arg_type = String 
            default = "./data/"
            help = "folder where GasLib data is available"
        "--file"
            help="GasLib zipped directory" 
            arg_type = String
            default="GasLib-4197.zip"
    end
    return parse_args(s) # the result is a Dict{String,Any}
end

function main(ARGS)
    args = parse_cli_args(ARGS)
    data = parse_gaslib(args["datafolder"] * args["file"])
    return data
end 
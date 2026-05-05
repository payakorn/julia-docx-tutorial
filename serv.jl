using LiveServer

# Parse port from command line arguments, default to 9000
port_arg = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 9000

# Directory to serve
dir_to_serve = joinpath(@__DIR__, "site")

println("Starting LiveServer on port $port_arg")
println("Serving directory: $dir_to_serve")
println("→  http://localhost:$port_arg/   (Ctrl-C to stop)")

# Start the server
serve(dir=dir_to_serve, port=port_arg, host="0.0.0.0")

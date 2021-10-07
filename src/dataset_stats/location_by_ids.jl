using Pkg
Pkg.activate(".")
using HGEpidemics
using DataFrames
using Dates
using CSV
using Query
using JSON
using JSON3
using JSONTables
using Serialization

"""
   Pre-compute the most crowded locations
"""

############################
# Loading simulation params
############################
# path = "src/dataset_stats/ble/configs/ble_params.json"
# path = "src/dataset_stats/foursquare/configs/fq_params.json"
path = "/home/lcrisci/workspace/Hypergraphs-Epidemics/deps/HGEpidemics/src/experiments/NPIs/SA/configs/no_npis/sa_params.json"
dataset_path = "/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/dataset.csv"
input_data = JSON.parse((open(path, "r")))

output_path = input_data["output_path"]
fdata_params = input_data["data_params"]
fparams = input_data["sim_params"]

jtable = jsontable(read(open(fparams, "r")))
paramsdf = DataFrame(jtable)
data_params = JSON3.read(read(open(fdata_params, "r")))

# The choice of the interval within which
# either an indirect (Δ) or direct (δ) contact
# may occur influences the data the
# simulation is run on.
# For this reason, it is necessary to store
# diffent information according to the
# values of both Δ and δ.
intervals_data = Dict{String, Dict{Symbol, Any}}()

header = [Symbol(col) for col in data_params.header]


# evaluating new dataset
# removing people
df, intervals, user2vertex, loc2he =
        generate_model_data(
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            datarow = data_params.datarow,
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

id_offsets = JSON.parsefile("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/info.json")

offset(start::String,last::String) = return id_offsets[start]:id_offsets[last]


transport_ids = filter(id -> haskey(loc2he,string(id)),id_offsets["transports_start_id"]:id_offsets["transports_start_id"]+500)
transports = [loc2he[string(id)] for id in transport_ids]
schools = [loc2he[string(id)] for id in offset("schools_start_id","households_start_id")]
workplaces = [loc2he[string(id)] for id in offset("workplaces_start_id","schools_start_id")]
schools_workplaces = schools ∪ workplaces ∪ transports
schools_100_workplaces_50 = schools ∪ rand(workplaces, floor(Int,length(workplaces) * 50/100)) ∪ transports

serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/schools.bin",schools)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/workplaces.bin",workplaces)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/schools_workplaces.bin",schools_workplaces)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/schools_100_workplaces_50.bin",schools_100_workplaces_50)
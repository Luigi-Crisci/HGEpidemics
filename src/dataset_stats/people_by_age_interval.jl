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


original_df = CSV.File(dataset_path) |> DataFrame

old_people_vect = @from row in original_df begin
    @where row.age >= 65
    @select row.person_id
    @collect
end

middle_age = @from row in original_df begin
	@where row.age >=20 && row.age < 65
	@select row.person_id
	@collect
end

young_people = @from row in original_df begin
	@where row.age < 20
	@select row.person_id
	@collect
end


middle_age_nodes = [user2vertex[string(id)] for id in middle_age]
young_people_nodes = [user2vertex[string(id)] for id in young_people]
old_people_vect_nodes = [user2vertex[string(id)] for id in old_people_vect]

# 50% young - 50% middle - 50% old
mixed_50_50_50 = rand(middle_age_nodes,floor(Int,length(middle_age_nodes)/2)) ∪ 
           rand(young_people_nodes,floor(Int,length(young_people_nodes)/2)) ∪ 
           rand(old_people_vect_nodes,floor(Int,length(old_people_vect_nodes)/2))

# 30% young - 50% middle - 80% old
mixed_30_50_80 = rand(middle_age_nodes,floor(Int,length(middle_age_nodes) * 30/100)) ∪ 
rand(young_people_nodes,floor(Int,length(young_people_nodes) * 50/100)) ∪ 
rand(old_people_vect_nodes,floor(Int,length(old_people_vect_nodes) * 80/100))

# 0% young - 100% middle - 100% old
mixed_0_100_100 = middle_age_nodes ∪ old_people_vect_nodes

serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/people_over_65.bin",old_people_vect_nodes)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/middle_age.bin",middle_age_nodes)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/people_under_20.bin",young_people_nodes)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/mixed_50_50_50.bin",mixed_50_50_50)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/mixed_30_50_80.bin",mixed_30_50_80)
serialize("/home/lcrisci/workspace/Hypergraphs-Epidemics/resources/datasets/Salerno/id_subsets/mixed_0_100_100.bin",mixed_0_100_100)
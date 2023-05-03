include("utils.jl")
include("functions.jl")


actT = collect((1:4:9)/10)
divT = collect((1:4:9)/10)
nDivs = collect(1:2:6)
ns = collect(1:2:6)
Rs = [2,4,7]
divTimes = 1
divFCs = 8
iter = collect(Iterators.product(actT, divT, nDivs, ns, Rs, divTimes, divFCs))
df = DataFrame(iter)
rename!(df, [:1, :2, :3, :4, :5, :6, :7] .=> [:actThresh, :divThresh, :nDiv, :nAct, :R, :divTime, :divFC])
CSV.write("parameterKey.csv", df)
# actDiv sims
Threads.@threads for x in 1:length(iter)
    # print(ad)
    ad = iter[x]
    mkpath(string("SimulationDataBoundaries/allParams/Case_", x))
    for i in 1:3
        fil = string("SimulationDataBoundaries/allParams/Case_", x,"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            actThresh = ad[1], divThresh = ad[2], nDiv = ad[3], n = ad[4],
            R = ad[5], divTime = 1, divTimeFold = 8)
        data, _ = run!(model, agent_step!, model_step!, 40; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end
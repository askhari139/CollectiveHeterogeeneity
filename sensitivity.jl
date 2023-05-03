 include("utils.jl")
 include("functions.jl")

#=
Parameters:
1. R
2. n
3. actThresh
4. divThresh
5. activation linearRate
=#


## Sensitivity all which where

#=
Default Parameters
1. divTime : 1.5
2. divTimeFold : 8
3. n : 5
4. nDiv : 1
5. R : 2
6. actThresh : 0.9
7. divThresh : 0.1

=#

actT = collect((1:9)/10)
divT = collect((1:9)/10)
nDivs = collect(1:6)
ns = collect(1:6)
Rs = collect(2:9)
divTimes = collect((1:6)/2)
divFCs = collect((1:5)*2)

# actDiv sims
Threads.@threads for ad in collect(Iterators.product(actT, divT))
    # print(ad)
    mkpath(string("SimulationDataBoundaries/actDiv/act_", ad[1], "_div_", ad[2]))
    for i in 1:3
        fil = string("SimulationDataBoundaries/actDiv/act_", ad[1], "_div_", ad[2],"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            actThresh = ad[1], divThresh = ad[2])
        data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end
# R 
Threads.@threads for x in 2:9
    mkpath(string("SimulationDataBoundaries/R/R_", x))
    for i in 1:3
        fil = string("SimulationDataBoundaries/R/R_",x,"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            R = x)
        data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end

# n activation
Threads.@threads for x in 1:6
    mkpath(string("SimulationDataBoundaries/nAct/n_", x))
    for i in 1:3
        fil = string("SimulationDataBoundaries/nAct/n_",x,"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            n = x)
        data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end

# # n division
Threads.@threads for x in 1:6
    mkpath(string("SimulationDataBoundaries/nDiv/nDiv_", x))
    for i in 1:3
        fil = string("SimulationDataBoundaries/nDiv/nDiv_",x,"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            nDiv = x)
        data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end

# min divtime
Threads.@threads for dFC in collect(Iterators.product(divTimes, divFCs))
    mkpath(string("SimulationDataBoundaries/divTime/divTime_", dFC[1], 
                "_divFC_", dFC[2]))
    for i in 1:3
        fil = string("SimulationDataBoundaries/divTime/divTime_", dFC[1], 
        "_divFC_", dFC[2],"/data_", i, ".csv")
        model = initialize(numagents = 2000, griddims = (100,100), 
            divTime = dFC[1], divTimeFold = dFC[2])
        data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
        data = @pipe data|>
            transform(_, :pos => ByRow(collect) => [:X, :Y])
        CSV.write(fil, data)
    end
end


# # random Parameters
# r = sample(Rs, 100)
# actTs = sample(actT, 100)
# divTs = sample(divT, 100)
# nAct = sample(ns, 100)
# nDivL = sample(nDivs, 100)
# divTm = sample(divTimes, 100)
# divFold = sample(divFCs, 100)
# df = DataFrame(R = r, actThresh = actTs, divThresh = divTs, 
#         n = nAct, nDiv = nDivL, divTime = divTm, divTimeFold = divFold)
# CSV.write(string("SimulationDataBoundaries/random_parameters.csv"), df)
# # df = CSV.read("SimulationDataBoundaries/random_parameters.csv", DataFrame)
# # r = df.R
# # actT = df.actThresh
# # divT = df.divThresh
# # nAct = df.n
# # nDivs = df.nDiv
# # divTm = df.divTime
# # divFold = df.divTimeFold
# # i = 93
# Threads.@threads for x in 9:100
#     mkpath(string("SimulationDataBoundaries/Randoms/random_", x))
#     for i in 1:3
#         fil = string("SimulationDataBoundaries/Randoms/random_",x,"/data_", i, ".csv")
#         model = initialize(numagents = 2000, griddims = (100,100), 
#             R = r[x], n = nAct[x], nDiv = nDivL[x], 
#             actThresh = actTs[x], divThresh = divTs[x], divTimeFold = divFold[x],
#             divTime = divTm[x])
#         data, _ = run!(model, agent_step!, model_step!, 80; adata = adata)
#         data = @pipe data|>
#             transform(_, :pos => ByRow(collect) => [:X, :Y])
#         CSV.write(fil, data)
#     end
# end





# parange = Dict(:actThresh => 0.0:0.1:5, :divThresh => 0.0:0.1:5, :n => 1:6, :R => 2:10)
# model = initialize(numagents = 200, griddims = (50,50))
# agentColor(b) = RGBf0(b.activationLevel, 0.2, 0.2)
# # 

# figure, adf, mdf = abm_data_exploration(
#     model, agent_step!, dummystep, parange;
#     ac=agentColor
#     # adata, alabels
# )


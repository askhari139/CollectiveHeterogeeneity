 #=
Model assumptions
1. We start with a seed population
2. All cells are inactive
3. Each cell has a proliferative capacity.
4. Proliferative capacity decreases with increasing local density
5. mitchondrial potential increases with increasing density, after a threshold density has been reached.
6. Cell division normally happens once in every 20 time steps. 
=#

function getTotalNeighbours(agent)
    pos = agent.pos

end

function agentAttrs!(agent, model)
    # print(agent.pos)
    dMax = length(collect(nearby_positions(agent.pos, model, model.R))) + 1
    # print(dMax)
    agent.density = length(collect(nearby_ids(agent.pos, model, model.R)))/dMax
    agent.activationLevel = agent.density^model.n/(agent.density^model.n + model.actThresh^model.n) + agent.lifeTime*model.linearRate
    agent.divTime = Int(round(model.divTime +  (model.divMaxTime - model.divTime)*(agent.density^model.nDiv/(agent.density^model.nDiv + model.divThresh^model.nDiv))))
end 

function pathFinder(p1, p2, model)
    x = p1[1]- p2[1]
    y = p1[2] - p2[2]
    path = []
    for i in 0:abs(x)
        for j in 0:abs(y)
            ps = (x>0 ? p1[1]-i : p1[1] + i, y>0 ? p1[2]-j : p1[2]+j)
            push!(path, ps)
        end
    end
    empty = maximum(findall(p -> isempty(p, model), path))
    path = path[empty:end]
    return path
end

function multiAgent(model)
    pos = allagents(model)
    multiPos = []
    for p in pos
        x = collect(agents_in_position(p, model))
        if length(x) > 1
            # for i in 2:length(x)
            #     kill_agent!(x[i], model)
            # end
            println(x.id)
        end
    end
    return multiPos
end


function initialize(; numagents = 2000, griddims = (100, 100),  
    divTime = 1.5, divTimeFold = 8, actThresh = 0.9, divThresh = 0.1, n=5, R = 2, nDiv = 1,
    linearRate = 0)
    # Random.seed!(seed)
    space = GridSpace(griddims)
    properties = Dict(:modelTime => 0, :divTime => divTime, :divMaxTime => divTime*divTimeFold,
    :divTimeFold => divTimeFold, :nDiv => nDiv,
    :actThresh => actThresh, :divThresh => divThresh,:n=>n, :R=>R, :linearRate=>linearRate)
    model = ABM(divisible, space;
                properties = properties, scheduler = Schedulers.fastest)
    # populate the model with agents, adding equal amount of the two types of agents
    # at random positions in the model
    for n in 1:numagents
        agent = divisible(n, (1, 1),
        false, 0.0, 0.01, rand(1:5, 1)[1], 0, 5, 0.0,0)
        # agent.activationLevel = agent.mode == false ? 0 : 1
        add_agent_single!(agent, model)
        agentAttrs!(agent, model)
    end
    return model
end

# Algorithm:
    #     1. get the empty positions around the agent in a radius of two cells by getting neighbour agent ids, 
    #         their positions and then comparing with all positions around the agent (empty or not)
    #     2. get the empty ones in the immediate neighbourhood
    #     3. if there are empty ones in the immediate neighbourhood, chose to fill one of those (emptyPos) with a probability of 2/3
    #     4. if not, get the empty ones in the second layer, chose one of them 
    #     5. get the positions neighbouring this newly chosen position that are also neighbours of the agent. Choose 1
    #     6. If this position is not empty, move the agent in this position to the 2nd layer. Mark this as emptyPos and fill it.

function divide_agent!(agent, model)
    pos = agent.pos
    nbs = collect(nearby_ids(agent.pos, model, model.R))
    nbPos = [model.agents[i].pos for i in nbs]
    allPos = collect(nearby_positions(pos, model, model.R))
    push!(allPos, pos)
    if length(allPos) == length(nbs) 
        return 
    end
    empty = allPos[findall(x->x ∉ nbPos, allPos)]
    emptyImmediate = empty[collect((abs(p1[1]-pos[1])<=1) && (abs(p1[2]-pos[2])<=1) for p1 in empty)]
    if length(emptyImmediate)!=0 && rand() > 0.5
        emptyPos = rand(emptyImmediate, 1)[1]
    else
        emptyExternal = empty[findall(x->x ∉ emptyImmediate, empty)]
        if isempty(emptyExternal) && isempty(emptyImmediate)
            return
        elseif isempty(emptyExternal)
            emptyPos = rand(emptyImmediate, 1)[1]
        else
            emptyPos = rand(emptyExternal,1)[1]
            path = pathFinder(emptyPos, pos, model)
            for i in 2:length(path)
                move_agent!(collect(agents_in_position(path[i], model))[1], path[i-1], model)
            end
        end
    end    
    n = nagents(model)
    agent.lifeTime = 0
    agent.time = rand(1:5,1)[1]
    daughter = divisible(n+1, emptyPos, agent.mode, 0.0,0.01, rand(1:5, 1)[1], 0, 5, 0.0, agent.id)
    add_agent!(daughter,emptyPos, model)
end



function agent_step!(agent, model)
    
    if agent.lifeTime >= agent.divTime && agent.density != 1
        divide_agent!(agent, model)
    end
    agentAttrs!(agent, model)
    agent.lifeTime += 1
    return
end

function model_step!(model)
    # println(model.modelTime)
    model.modelTime += 1
    # multiAgent(model)
end

function mean(x)
    sum(x)/length(x)
end

function intChane(x)
    x = replace(string(x), "\\" => "")
    parse(Int64, x)
end

adata = [:pos, :activationLevel, :lifeTime, :divTime, :density, :parentID]
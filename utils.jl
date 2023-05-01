using Agents,
# AgentsPlots,
DataFrames,
 CSV,
 Random,
#  CairoMakie,
 InteractiveDynamics,
#  GLMakie,
 Base.Threads,
 StatsBase,
 Pipe

function mean(x)
    sum(x)/length(x)
end

function sd(x)
    (sum((x.-mean(x)).^2)/length(x))^(1/2)
end

mutable struct metabolism <: AbstractAgent
    id::Int # The identifier number of the agent
    pos::Tuple{Int,Int} # The x, y location of the agent on a 2D grid
    mode::Bool # whether the agent is active. (true => active)
    activationLevel::Float64
    time::Int
    lifeTime::Int
    parentID::Int
end


@agent divisible GridAgent{2} begin
    mode::Bool # whether the agent is active. (true => active)
    activationLevel::Float64
    actInit::Float64
    time::Int
    lifeTime::Int
    divTime::Int
    density::Float64
    parentID::Int
end

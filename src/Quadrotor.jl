module Quadrotor

using Plots

# Write your package code here.
export QuadrotorParameters, quadrotor!, controller, animate
include("dynamics.jl")

end

module Quadrotor

using Plots

# Write your package code here.
export QuadrotorParameters, quadrotor_dynamics!, controller, animate, linear_quadrotor_dynamics!
include("dynamics.jl")

end

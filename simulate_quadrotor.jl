using Revise
using Quadrotor
using Plots
using DifferentialEquations
using Sundials

plotly()


# initialize the initial value problem
qp = Quadrotor.QuadrotorParameters(1.0, 1.0, 0.2, 0.15, 9.81)
dt = 0.05 # interpolate / snapshot simulation every 0.05 (s) in sim time
tspan = (0.0, 20.0) # simulate from 0 to 100 (s)
init_conds = zeros(12, 1) # set initial conditions of all state vars to 0
ic_prob = ODEProblem(Quadrotor.quadrotor_dynamics!, init_conds, tspan, qp)

# solve the initial value problem using by numerically integrating the ODE
sol = solve(ic_prob, AutoVern7(Rodas5()), saveat=dt)

# compute the control over time since it's not in the solution return struct
# since the control signals aren't state variables
u = zeros(length(sol.t), 4) # 4 because there are 4 control signals
for i in 1:length(sol.t)
    # for each saved state vector, compute what we the control, then add to a matrix
    u[i,:] = Quadrotor.controller(sol[:,1], qp)'
end

# plot state variables over time
plot(sol.t, sol[1,:], label="x"); plot!(sol.t, sol[2,:], label="y"); plot!(sol.t, sol[3,:], label="z", xlabel="time (s)", ylabel="position (m)", title="position over time", ylim=(-1.5, 1.5))
plot(sol.t, sol[7,:], label="x"); plot!(sol.t, sol[8,:], label="y"); plot!(sol.t, sol[9,:], label="z", xlabel="time (s)", ylabel="velocity (m/s)", title="velocity over time", ylim=(-1.5, 1.5))
plot(sol.t, sol[4,:], label="yaw"); plot!(sol.t, sol[5,:], label="pitch"); plot!(sol.t, sol[6,:], label="roll", xlabel="time (s)", ylabel="angles (rad)", title="attitude over time", ylim=(-pi, pi))
plot(sol.t, sol[10,:], label="yawrate"); plot!(sol.t, sol[11,:], label="pitchrate"); plot!(sol.t, sol[12,:], label="rollrate", xlabel="time (s)", ylabel="angle rates (rad/s)", title="attitude rate over time", ylim=(-pi, pi))

# plot control signals over time
plot(sol.t, u[:,1], label="rotor 1"); plot!(sol.t, u[:,2], label="rotor 2"); plot!(sol.t, u[:,3], label="rotor 3"); plot!(sol.t, u[:,4], label="rotor 4", xlabel="time (s)", ylabel="angular velocity (rad/s)", title="control over time")

anim = Quadrotor.animate(sol)
gif(anim, fps=1/dt)
"""
Stores the model parameters of the quadrotor dynamical system

M ~ mass of the fuselage (kg)
m ~ mass of the rotors (kg)
L ~ distance between center of mass and center of each rotor (m)
l ~ radius of each rotor (m)
g ~ acceleration due to gravity (m/s^2)
"""
struct QuadrotorParameters
    M::Real
    m::Real
    L::Real
    l::Real
    g::Real
end


"""
    quadrotor_dynamics!(dx, state, params, t)

implements the update equation for an ODE solver representing the nonlinear 
dynamics of a quadrotor

# Examples
```julia-repl
julia> # initialize the initial value problem

julia> qp = Quadrotor.QuadrotorParameters(1.0, 1.0, 0.2, 0.15, 9.81);

julia> dt = 0.05; # interpolate / snapshot simulation every 0.05 (s) in sim time

julia> tspan = (0.0, 100.0); # simulate from 0 to 100 (s)

julia> init_conds = zeros(12, 1); # set initial conditions of all state vars to 0

julia> ic_prob = ODEProblem(Quadrotor.quadrotor_dynamics!, init_conds, tspan, qp);

julia> # solve the initial value problem using by numerically integrating the ODE

julia> sol = solve(ic_prob, AutoVern7(Rodas5()), saveat=dt)
retcode: Success
Interpolation: 1st order linear
t: 2001-element Vector{Float64}:
   0.0
   0.05
   0.1
   0.15
   0.2
   0.25
   0.3
   0.35
   0.4
   ⋮
  99.6
  99.65
  99.7
  99.75
  99.8
  99.85
  99.9
  99.95
 100.0
u: 2001-element Vector{Matrix{Float64}}:
 [0.0; 0.0; … ; 0.0; 0.0;;]
 [0.0; 2.3648401e-316; … ; 2.3648401e-316; 2.3646365e-316;;]
 [0.0; 1.496431937e-315; … ; 1.496431937e-315; 1.496303134e-315;;]
 [0.0; -5.17187286e-316; … ; -5.17187286e-316; -5.1714277e-316;;]
 [0.0; -1.135799885e-315; … ; -1.135799885e-315; -1.135702124e-315;;]
 [0.0; -5.7802178e-316; … ; -5.7802178e-316; -5.77972024e-316;;]
 [0.0; 9.1689063e-316; … ; 9.1689063e-316; 9.16811715e-316;;]
 [0.0; 2.51756253e-315; … ; 2.51756253e-315; 2.517345843e-315;;]
 [0.0; 3.450161674e-315; … ; 3.450161674e-315; 3.44986471e-315;;]
 ⋮
 [0.0; -3.841254157e-315; … ; -3.841254157e-315; -3.84451148e-315;;]
 [0.0; -3.971835065e-315; … ; -3.971835065e-315; -3.975202617e-315;;]
 [0.0; -3.932100325e-315; … ; -3.932100325e-315; -3.9354345e-315;;]
 [0.0; -3.71997618e-315; … ; -3.71997618e-315; -3.723129855e-315;;]
 [0.0; -3.33337525e-315; … ; -3.33337525e-315; -3.336201747e-315;;]
 [0.0; -2.770197e-315; … ; -2.770197e-315; -2.772545675e-315;;]
 [0.0; -2.028328575e-315; … ; -2.028328575e-315; -2.03004836e-315;;]
 [0.0; -1.105643516e-315; … ; -1.105643516e-315; -1.106579385e-315;;]
 [0.0; 0.0; … ; 0.0; 0.0;;]
```
"""
function quadrotor_dynamics!(dx, state, params, t)

    # unpack model parameters into individual local variables
    M = params.M
    m = params.m
    L = params.L
    l = params.l
    g = params.g

    # unpack state vector into individual local variables
    x = state[1] # x-pos wrt some fixed point
    y = state[2] # y-pos wrt some fixed point
    z = state[3] # z-pos wrt some fixed point
    alpha = state[4] # yaw of the quadrotor in body-frame
    beta = state[5] # pitch of the quadrotor in body-frame
    gamma = state[6] # roll of the quadrotor in body-frame
    xdot = state[7] # vel along x-axis wrt fixed space frame
    ydot = state[8] # vel along y-axis wrt fixed space frame
    zdot = state[9] # vel along z-axis wrt fixed space frame
    alphadot = state[10] # yawspeed of quadrotor in body-frame
    betadot = state[11] # pitchspeed of quadrotor in body-frame
    gammadot = state[12] # rollspeed of quadrotor in body-frame

    # get control input from controller
    u = controller(state, params)

    # forward acceleration of the quadrotor (forward being upwards if level)
    uu = u' * u
    @assert typeof(uu) == Matrix{Float64}
    @assert size(uu) == (1, 1)
    @assert typeof(uu[1]) == Float64
    a = 1 / (M + 4 * m) * uu[1]

    # quadrotor dynamics
    xddot = a * (cos(alpha) * sin(beta) * cos(gamma) - sin(alpha) * sin(gamma))
    yddot = a * (sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma))
    zddot = a * (cos(beta) * cos(gamma)) - g

    alphaddot = m * l / ((M + 4 * m) * L) * sum(u)
    betaddot = 1 / m * (u[1]^2 - u[3]^2)
    gammaddot = 1 / m * (u[2]^2 - u[4]^2)

    # update equation
    dx[1:12] = [
        xdot; ydot; zdot; alphadot; betadot; gammadot; 
        xddot; yddot; zddot; alphaddot; betaddot; gammaddot
    ]
end

"""
    linear_quadrotor_dynamics!(dx, state, params, t)


"""
function linear_quadrotor_dynamics!(dx, state, params, t)

    # unpack model parameters into individual local variables
    M = params.M
    m = params.m
    L = params.L
    l = params.l
    g = params.g

    # unpack state vector into individual local variables
    x = state[1] # x-pos wrt some fixed point
    y = state[2] # y-pos wrt some fixed point
    z = state[3] # z-pos wrt some fixed point
    alpha = state[4] # yaw of the quadrotor in body-frame
    beta = state[5] # pitch of the quadrotor in body-frame
    gamma = state[6] # roll of the quadrotor in body-frame
    xdot = state[7] # vel along x-axis wrt fixed space frame
    ydot = state[8] # vel along y-axis wrt fixed space frame
    zdot = state[9] # vel along z-axis wrt fixed space frame
    alphadot = state[10] # yawspeed of quadrotor in body-frame
    betadot = state[11] # pitchspeed of quadrotor in body-frame
    gammadot = state[12] # rollspeed of quadrotor in body-frame

    # get control input from controller
    u = controller(state, params)
    u0 = sqrt(1 / 4 * (M + 4 * m) * g) # equilibrium control

    # quadrotor dynamics
    xddot = 1 / (M + 4 * m) * (4 * u0^2 * beta)
    yddot = 1 / (M + 4 * m) * (4 * u0^2 * gamma)
    zddot = 1 / (M + 4 * m) * (2 * u0 * u[1] - 2 * u0 * u[2] + 2 * u0 * u[3] - 2 * u0 * u[4] - 8 * u0^2)

    alphaddot = m * l / ((M + 4 * m) * L) * sum(u)
    betaddot = 1 / m * (2 * u0 * u[1] - 2 * u0 * u[3])
    gammaddot = 1 / m * (2 * u0 * u[2] - 2 * u0 * u[4])

    # update equation
    dx[1:12] = [
        xdot; ydot; zdot; alphadot; betadot; gammadot; 
        xddot; yddot; zddot; alphaddot; betaddot; gammaddot
    ]
end

"""
    controller(state, params)

Computes the control signals going to the four rotors given the system state.

Specifically, it returns a type "Vector{Float64}" and size (4,)

# Examples
```julia-repl
julia> controller(zeros(12,1), Quadrotor.QuadrotorParameters(1.0, 1.0, 0.2, 0.15, 9.81))
4x1 Matrix{Float64}:
 3.5017852589786256
 3.5017852589786256
 3.5017852589786256
 3.5017852589786256
```
"""
function controller(state, params)
    # unpack model parameters
    m = params.m
    M = params.M
    l = params.l
    L = params.L
    g = params.g

    # equilibrium control (manually computed so a bit of a magic expression)
    u0 = sqrt(1 / 4 * (M + 4 * m) * g)

    # Control Matrix (manually computer so a bit of a magic number)
    K = [
        0.7071     0.0000    0.5000   15.8114   24.0083   -0.0000    1.9902    0.0000    0.6546    7.2770    1.9819   -0.0000;
        0.0000     0.7071   -0.5000   15.8114   -0.0000   24.0083    0.0000    1.9902   -0.6546    7.2770    0.0000    1.9819;
        -0.7071    0.0000    0.5000   15.8114  -24.0083    0.0000   -1.9902    0.0000    0.6546    7.2770   -1.9819    0.0000;
        -0.0000   -0.7071   -0.5000   15.8114   -0.0000  -24.0083   -0.0000   -1.9902   -0.6546    7.2770   -0.0000   -1.9819
    ];

    # standard plain state-feedback control law
    delta_u = -K * state;

    # add the equilibrium since this is a linearization with an offset
    return delta_u .+ u0
end

"""
"""
function animate(sol)
    # definte animation limits to be +/-2.5 meters in each direction (x,y,z)
    axis_limits = 2.5
    lim = (-axis_limits, axis_limits)
    anim = Plots.Animation()

    # switch (if not already using it) to GR backend or Plots
    gr()

    for i in 1:length(sol.t)
        anim_title = string("t = ", round(sol.t[i], digits=1), " (s)")

        p1 = plot3d([sol[1,i]], [sol[2,i]], [sol[3,i]], size=(400, 400), xlim=lim, ylim=lim, zlim=lim, markershape=:circle, camera=(50,55), legend=:none)

        plt = plot(p1, legend=:none, title=anim_title)
        frame(anim, plt)
    end


    # switch back to plotly
    plotly()

    return anim
end
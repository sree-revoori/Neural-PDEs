using DelimitedFiles,Plots
using DiffEqSensitivity, OrdinaryDiffEq, Zygote, Flux, DiffEqFlux, Optim

# Problem setup parameters:
Lx = 10.0
x  = 0.0:0.01:Lx
dx = x[2] - x[1]
Nx = size(x)

u0 = exp.(-(x.-3.0).^2) # I.C

## Problem Parameters
p        = [1.0,1.0]    # True solution parameters
xtrs     = [dx,Nx]      # Extra parameters
dt       = 0.40*dx^2    # CFL condition
t0, tMax = 0.0 ,1000*dt
tspan    = (t0,tMax)
t        = t0:dt:tMax;

## Definition of Auxiliary functions
function ddx(u,dx)
    """
    2nd order Central difference for 1st degree derivative
    """
    return [[zero(eltype(u))] ; (u[3:end] - u[1:end-2]) ./ (2.0*dx) ; [zero(eltype(u))]]
end

function d2dx(u,dx)
    """
    2nd order Central difference for 2nd degree derivative
    """
    return [[zero(eltype(u))]; (u[3:end] - 2.0.*u[2:end-1] + u[1:end-2]) ./ (dx^2); [zero(eltype(u))]]
end

## ODE description of the Physics:
function burgers(u,p,t)


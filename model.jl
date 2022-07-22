using DataDrivenDiffEq, ModelingToolkit, OrdinaryDiffEq,  Plots, DelimitedFiles

time_ypos = readdlm("/Users/luth/Downloads/HeavyRedDrop1.txt", ',') # read in the data
time_ypos = time_ypos[setdiff(1:end, (1, 2)), :] #deletes headers
time_ypos = time_ypos[:, 1:end .!= 3] #deletes column 3
ypos = time_ypos[:, 1:end .!= 1]
ypos = vec(ypos)
time = time_ypos[:, 1:end .!= 2]
time = vec(time)

function absolute_ypos(ypos)
    for i in 1:length(ypos)
        ypos[i] = abs(ypos[i])
    end
end

absolute_ypos(ypos)

function velocitycalc(variable::Vector, ypos::Vector, time::Vector) #appends the velocity as calculated from given y-position and time
    variable = "1" ^ length(ypos)
    for i in 1:length(ypos)
        if i == 1
            variable[i] = 0
        else1
            variable[i] = (ypos[i] - ypos[i-1])/(time[i] - time[i-1])
        end
    end
end

velocitycalc(X, ypos, time)

X = Array(sol) # this is velocity, should be replaced by the measurement later on
T = time
ddprob = ContinuousDataDrivenProblem(X,T) # define a datadriven problem

@variables u(t) # u here is velocity, diferent 
order = 4
Ψ = Basis(monomial_basis([u], order), [u], independent_variable = t)
println(Ψ)
threshold = 0.01
res = solve(ddprob, Ψ, STLSQ(threshold))

println(result(res))
println(parameters(res))

res_sys = result(res)
@named reconstructed_sys = ODESystem(equations(res_sys), get_iv(res_sys), states(res_sys), parameters(res_sys));
x0 = [u => X[1]]
ps = parameter_map(res)
reconstructed_prob = ODEProblem(reconstructed_sys, x0, tspan, ps)
estimate = solve(reconstructed_prob, Tsit5(), saveat = T);

step=floor(Int, length(T)/20)
scatter(T[1:step:end], transpose(X)[1:step:end], legend=:bottom)
plot!(estimate)
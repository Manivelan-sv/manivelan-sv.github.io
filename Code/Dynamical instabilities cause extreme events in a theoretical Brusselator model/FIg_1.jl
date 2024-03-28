#=
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
=#

using Plots,DelimitedFiles

x_values = 0:0.005:1
y_values = 0:0.005:2
results = []
for x in x_values
    for y in y_values
        lambda_1 = 0.5 * ((-(x*x-y+1)) - sqrt(Complex((x*x-y+1)^2 - (4 *x*x) ))) # Calculated analytically
        lambda_2 = 0.5 * ((-(x*x-y+1)) + sqrt(Complex((x*x-y+1)^2 - (4 *x*x) )))
        
        # Determine stability
        if real(lambda_1) < 0 && real(lambda_2) < 0 && imag(lambda_1) == 0 && imag(lambda_2) == 0
            stability = "Stable Node"
        elseif real(lambda_1) > 0 && real(lambda_2) > 0 && imag(lambda_1) == 0 && imag(lambda_2) == 0
            stability = "Unstable Node"
        elseif real(lambda_1) > 0 && real(lambda_2) < 0
            stability = "Saddle Point"
        elseif real(lambda_1) < 0 && real(lambda_2) > 0
            stability = "Saddle Point"
        elseif real(lambda_1) > 0 && real(lambda_2) > 0 && imag(lambda_1) == -imag(lambda_2) != 0
            stability = "Unstable Focus"
        elseif real(lambda_1) < 0 && real(lambda_2) < 0 && imag(lambda_1) == -imag(lambda_2) != 0
            stability = "Stable Focus"
        elseif real(lambda_1) == 0 && real(lambda_2) == 0 && imag(lambda_1) == -imag(lambda_2) != 0
            stability = "Centre"
        else
            stability = "Other"
        end
        push!(results, (x,y, stability))
    end
end

stable_node = [(x,y) for (x,y, stability) in results if stability == "Stable Node"]
unstable_node = [(x,y) for (x,y, stability) in results if stability == "Unstable Node"]
saddle_point = [(x,y) for (x,y, stability) in results if stability == "Saddle Point"]
unstable_focus = [(x,y) for (x,y, stability) in results if stability == "Unstable Focus"]
stable_focus = [(x,y) for (x,y, stability) in results if stability == "Stable Focus"]
centre = [(x,y) for (x,y, stability) in results if stability == "Centre"]

writedlm("sn.dat",stable_node)
writedlm("usn.dat",unstable_node)
writedlm("sp.dat",saddle_point)
writedlm("usf.dat",unstable_focus)
writedlm("sf.dat",stable_focus)
writedlm("c.dat",centre)

scatter(stable_node, color=:green,markershape=:circle,markersize=:1,markerstrokecolor=:green, label="1")
scatter!(unstable_node, color=:orange,markersize=:1, markershape=:circle,markerstrokecolor=:orange, label="2")
scatter!(saddle_point, color=:red, markersize=:1,markershape=:circle,markerstrokecolor=:red, label="3")
scatter!(unstable_focus, color=:blue,markersize=:1, markershape=:circle,markerstrokecolor=:blue, label="4")
scatter!(stable_focus, color=:purple,markersize=:1, markershape=:circle,markerstrokecolor=:purple, label="5")
scatter!(centre, color=:yellow, markersize=:1,markershape=:circle,markerstrokecolor=:yellow, label="6")
xlabel!("x")
ylabel!("y")

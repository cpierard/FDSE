
# This script reads in output from KH.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Plots, Printf, Statistics

# Set the filename (without the extension)
filename = "KH"

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)
b_ic = FieldTimeSeries(filename * ".jld2", "b", iterations = 0)
ω_ic = FieldTimeSeries(filename * ".jld2", "ω", iterations = 0)
χ_ic = FieldTimeSeries(filename * ".jld2", "χ", iterations = 0)
ϵ_ic = FieldTimeSeries(filename * ".jld2", "ϵ", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xb, yb, zb = nodes(b_ic)
xω, yω, zω = nodes(ω_ic)
xχ, yχ, zχ = nodes(χ_ic)
xϵ, yϵ, zϵ = nodes(ϵ_ic)

## Now, open the file with our data
file_xz = jldopen(filename * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xz["timeseries/t"]))

@info "Making an animation from saved data..."


mean_KE_perturbation = zeros(length(iterations))
time = zeros(length(iterations))

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    u_xz = file_xz["timeseries/u/$iter"][:, 1, :];
    v_xz = file_xz["timeseries/v/$iter"][:, 1, :];
    w_xz = file_xz["timeseries/w/$iter"][:, 1, :];
    b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
    ω_xz = file_xz["timeseries/ω/$iter"][:, 1, :];
    χ_xz = file_xz["timeseries/χ/$iter"][:, 1, :];
    ϵ_xz = file_xz["timeseries/ϵ/$iter"][:, 1, :];

    u_avg = mean(u_xz, dims = 1);
    w_avg = mean(w_xz, dims = 1);

    u_pert = u_xz .- u_avg;
    w_pert = w_xz .- w_avg;

    KE = 0.5 .*(u_pert .^ 2) .+ (w_pert[:,:64] .^ 2)
    mean_KE_perturbation[i] = mean(KE)

    t = file_xz["timeseries/t/$iter"];
    time[i] = t


        b_xz_plot = heatmap(xb, zb, b_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ω_xz_plot = heatmap(xω, zω, ω_xz'; color = :viridis, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        χ_xz_plot = heatmap(xχ, zχ, χ_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ϵ_xz_plot = heatmap(xϵ, zϵ, KE'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 

    u_title = @sprintf("u, t = %s", round(t));
    v_title = @sprintf("v, t = %s", round(t));
    w_title = @sprintf("w, t = %s", round(t));
    b_title = @sprintf("b, t = %s", round(t));
    ω_title = @sprintf("vorticity (ω), t = %s", round(t));
    ϵ_title = @sprintf("KE perturbation, t = %s", round(t));
    χ_title = @sprintf("buoyancy variance dissipation (χ), t = %s", round(t));

# Combine the sub-plots into a single figure
    plot(b_xz_plot, ω_xz_plot, ϵ_xz_plot, χ_xz_plot, layout = (4, 1), size = (1200, 800),
    title = [b_title ω_title ϵ_title χ_title])

    iter == iterations[end] && close(file_xz)
end

# Save the animation to a file
mp4(anim, "KH.mp4", fps = 20) # hide

plot(time, mean_KE_perturbation, xlabel = "time", ylabel = "mean KE perturbation", legend = false)

### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 82c494e0-b75e-474e-b08c-cde3791a121e
using Pkg; Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/julne")

# ╔═╡ 8aaabf9e-08f7-11ed-1880-e798d62d8295
# NEXT Julia with Pluto

begin
using HDF5
using DataFrames
using Plots
using Statistics
using PlutoUI
#using PlotlyJS
end

# ╔═╡ 3d81d403-fbcc-4bc9-85f7-c1f257dc7695
function ingredients(path::String)
# this is from the Julia source code (evalfile in base/loading.jl)
# but with the modification that it returns the module instead of the last object
name = Symbol(basename(path))
m = Module(name)
Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
m
end

# ╔═╡ 03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
ng =ingredients("../src/julne.jl")

# ╔═╡ 853758b5-fecf-4a02-a299-4b5d7f0d99de
PlutoUI.TableOfContents(title="NEXT-100 Event Gallery", indent=true)

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

# NEXT-100 Events Gallery


This NB displays NEXT-100 MC events using Beersheba-reconstructed and MC hits of $^{214}$Bi and $\beta\beta0\nu$ tracks in the $Q_{\beta\beta}$ energy range


J.A. Hernado, JJ. Gomez-Cadenas; 
Donostia, July 2022

---
"""

# ╔═╡ 4f6f8003-b9f0-4aba-a42c-c79184bec77d
md"""
## Code
----

@todo: locate into a module
"""

# ╔═╡ 6df2cbcb-28fc-477f-986f-cfc9efea28b8
function _getdf(dataset)

    dda = read(dataset)
    
    ddic = Dict()
    for (i, key) in enumerate(keys(dda[1]))
        #println(i, ',', key)
        ddic[key] = [ddi[i] for ddi in dda]
    end

    df = DataFrame(ddic)
end

# ╔═╡ c095ec23-9567-4418-95c9-e011e632ddd0
"""
    get the DataFrames of Beershaba Voxels & Isaura

    Parameters
    ----------
    filename : str, filename

    Returns
    -------
    dfs      : Dict(DF), dictionary of DataFrames

"""
function get_dfs(filename::String)

	fid     = h5open(filename, "r")
	
	dfs = Dict()

	for key in ["BeershebaVoxels", "MCHits", "BinsInfo"]
		dataset  = fid["DATASET"][key]		
	    dfs[key] = _getdf(dataset)
	end

	close(fid)
	
    return dfs
end

# ╔═╡ 6407e43e-fe73-4c94-a839-4a42f112549b
function get_event(df::DataFrame, event::Int64)

    sel  = Vector(df.dataset_id .== event)
    println(sum(sel))
    xbin = df[!, :xbin]
    idf  = df[sel, :]
    
end;

# ╔═╡ f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# NEXT-100 214Bi
begin
datadir   = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
filenames = Dict(:Bi214 => "Bi/label_beersheba_554mm_214Bi_ICS.h5", 
             :bb0nu => "bb0nu/v1/label_beersheba_554mm_0nubb.h5"); 
#filename = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/Bi/label_beersheba_554mm_214Bi_ICS.h5"
end;

# ╔═╡ 783ebc81-712b-4f4b-a048-43e04c528652
""" 
    returns dimensions: frame origin and vixel size

    Parameters
    ----------
    dfs : Bins dataFrame

    Returns
    -------
    x0         : np.array(float), origin of the frame
    voxel_size : np.array(float), voxel size

"""
function _get_dimensions(dfs)

	df         = dfs["BinsInfo"]
    
    voxel_size = [df[!, var][1] for var in ("size_x", "size_y", "size_z")]
    x0         = [df[!, var][1] for var in ("min_x" , "min_y" , "min_z")]
	x1         = [df[!, var][1] for var in ("max_x" , "max_y" , "max_z")]
     
    return x0, x1, voxel_size
end

# ╔═╡ e12c8812-a220-4451-a2e7-602e29838010
md"""
## Load the data
"""

# ╔═╡ 8fe43246-a006-4ff9-bb64-3601d4e03938
begin
biso = @bind isotope Select([:Bi214, :bb0nu])

md"""

Select an isotope: $(biso)
"""
end

# ╔═╡ 6e7e284e-88a1-47f9-9da8-4a79797d8ebd
filename = string(datadir, filenames[isotope]);

# ╔═╡ c83851f5-5854-4218-8904-8418073c9296
md"""

Selected isotope  : $(isotope)

Datafile          : $(filename)
"""

# ╔═╡ 0398c69e-e594-4adb-8099-6524e2db5a31
begin
dfs           = ng.julne.get_dfs(filename)
df            = dfs["BeershebaVoxels"]
x0, x1, steps = _get_dimensions(dfs)
end

# ╔═╡ dcab561e-a163-4541-bccf-7a4065d01697
keys(dfs)

# ╔═╡ f57ac756-9ce6-4378-8d52-010c13795040
begin
df[!, "x"] = x0[1] .+ df[!, "xbin"] * steps[1]
df[!, "y"] = x0[2] .+ df[!, "ybin"] * steps[2]
df[!, "z"] = x0[3] .+ df[!, "zbin"] * steps[3]
end

# ╔═╡ 00c68cd9-f39b-4272-9991-95fcb487b119
df.z

# ╔═╡ 03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
begin
events = sort(collect(Set(df.dataset_id)))

md"""

Number of events : $(length(events))
"""
end

# ╔═╡ ddae0102-eba2-438f-8eb6-62380e905684
md"""

total number of entries in Beersheba's hit table: $(size(df)[1])

The dimension of the table are:

|coordinate | min (mm) | max (mm) | step (mm)| 
| :-- | :-- | :-- | :-- |
| x | $(x0[1])| $(x1[1]) | $(steps[1]) |
| y | $(x0[2])| $(x1[2]) | $(steps[2]) |
| z | $(x0[3])| $(x1[3]) | $(steps[3]) |
"""

# ╔═╡ 7913b885-4f25-40ae-ab84-3b5a4d50967c
begin

bevent = @bind event Select(events[1:20])

md"""
Select event number : $(bevent)
"""
end

# ╔═╡ 14160a98-07a0-4efc-9e71-9aab36ed01b6
begin
idf = get_event(df, event)

md"""

Event $(event) has $(size(idf)[1]) hits

"""
end

# ╔═╡ 0300688a-6899-46d2-92bd-4b018085e6e8
begin
bi_ntracks     = [34, 43, 112, 142, 157, 166, 191, 264, 296, 317, 329, 335, 348, 354, 364, 409, 421, 474, 494, 503]
bi_graph_noroi = [2, 10, 23, 25, 32, 34, 40, 43, 86, 93, 98, 100, 101, 107, 112, 121, 134, 142, 157, 166]
bi_graph_nob2  = [12, 95, 138, 188, 237, 244, 261, 288, 302, 393]
md"""
Selection of problematic $^{214}$Bi events
"""
end

# ╔═╡ e4ff7a0a-243d-4b51-9a49-7172edca1126
idf.z

# ╔═╡ 99a11ca4-b453-4e0b-a737-d9733b0e8c59
md"### Plots histograms"

# ╔═╡ 718fbf98-1725-429b-9a1c-609cb38430ac
begin
xwidth, ywidth, zwidth = steps
zrange = range(minimum(idf.z), maximum(idf.z) + zwidth; step = zwidth) 
xrange = range(minimum(idf.x), maximum(idf.x) + xwidth; step = xwidth) 
yrange = range(minimum(idf.y), maximum(idf.y) + ywidth; step = ywidth) 
dx = maximum(idf.x) - minimum(idf.x)
dy = maximum(idf.y) - minimum(idf.y)
dz = maximum(idf.z) - minimum(idf.z)
	
md"""

Event window:

| coordinate | min (mm) | max (mm)| width (mm)
| :-- | :-- | :-- | :-- |
| x | $(minimum(idf.x)) | $(maximum(idf.x))| $(dx) |
| y | $(minimum(idf.y)) | $(maximum(idf.y))| $(dy) |
| z | $(minimum(idf.z)) | $(maximum(idf.z))| $(dz) |

"""
end

# ╔═╡ 65597330-86b4-4d17-abae-5df638603cfc
plotly();
#gr()

# ╔═╡ ca22eb02-282a-412b-9297-15a1e379a7b2
md"""

## S2 signal, and scatter proyections of the hits
"""

# ╔═╡ 5f913a6e-1798-42da-8fbe-3dad64fc7bb7
begin
	
hez_ = histogram(idf.z, bins = zrange, weights = idf.energy, 
	xlabel = "z (mm)", ylabel = "energy (keV)")
hzx_ = histogram2d(idf.z, idf.x, bins = (zrange, xrange),  xlabel = "z (mm)", ylabel = "x (mm)")
hzy_ = histogram2d(idf.z, idf.y, bins = (zrange, yrange),  xlabel = "z (mm)", ylabel = "y (mm)")
hxy_ = histogram2d(idf.x, idf.y, bins = (xrange, yrange),  xlabel = "y (mm)", ylabel = "x (MeV)")
plot(hez_, hzx_, hzy_, hxy_; layout = grid(2, 2), size = (700, 600))
end

# ╔═╡ 0a6a2157-b02f-478b-947b-f579ca1c78c0
plotattr("size")

# ╔═╡ d13e0a1a-8eb4-49eb-a8b7-b732cf9ffe31
md"""
---
"""

# ╔═╡ f547f778-8a74-4923-9373-cb96b2254c9b
scatter(idf.xbin, idf.zbin,  markersize = 1e2*idf.energy)

# ╔═╡ 93d7da9e-320e-4e01-928e-40d29b2b2be7


# ╔═╡ 96581e7a-ee05-48a5-9120-88a1e0d37e7c
begin
theme(:dark)
scatter3d(idf.xbin, idf.ybin, idf.zbin, markersize = 1e2*idf.energy, markertype = "triangle", label = false,
xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")
end

# ╔═╡ 94310819-6bf4-427e-85d8-82d11ceae3ce
md"""
----
----
"""

# ╔═╡ 2ae832a9-f98e-42f8-a6d5-999c2bca44cb
mc = dfs["MCHits"]
imc = _get_event(mc, event)

# ╔═╡ 8316a0b8-4e9f-47bd-abfd-3bba96a5f3a3
plotattr("markercolor")

# ╔═╡ 795b757c-b5fd-4dac-9950-7b690848319a
md""" 
### statistics
"""

# ╔═╡ 746a87a1-0910-4bcd-93af-e4ca870c8847
"""
main statistics of a vector, length, mean, std
"""
function vstats(v::Vector{T}) where {T <: Number}
	(
	entries = length(v), 
	mean    = mean(v), 
	std     = std(v)
	)
end

# ╔═╡ 1f36611e-53f1-4371-b121-c8ea8de15d18
"""
Convert into a string a named-tuple
"""
function strnamedvalues(v::NamedTuple; digits = 3)
	vround = NamedTuple{keys(v)}(round.(values(v); digits = digits))
	ss = reduce(*, [string(key, " = ", vi, ", \n") for (key, vi) in pairs(vround)])
	#println(ss)
end

# ╔═╡ 042f1233-6cfb-4cee-b858-79c970dd56a7
begin
v = vstats(1e3*idf.energy)
ss = strnamedvalues(v, digits = 4)
end

# ╔═╡ c9929a2b-9d35-4f25-95b9-01e79da646c6
begin
h1 = histogram(idf.energy, bins = 100, xlabel = "x (mm)", normalize = true,
	label = false)
plot(h1, annotations = (0.05, 150, ss), annotationfontsize = 10)
end

# ╔═╡ 069fdfda-23b2-46be-bbc7-a612944dfad2
plotattr("legend")

# ╔═╡ 16a0d7e4-65df-48f2-9f8e-c83d9370cac6
plot(h1; annotations = (165, 0.08, ss), annotationfontsize = 10) 
#legend_position = :topbottom)

# ╔═╡ 1231b113-242d-4ae0-acf7-72fc0600eb6a
@bind xxx range()

# ╔═╡ 6d89975b-844c-496b-aecd-c941ca596a83
md"""
# This is markdown

## This is second level of MD
- This is 
- A list

And look at me = $xxx
"""

# ╔═╡ Cell order:
# ╠═82c494e0-b75e-474e-b08c-cde3791a121e
# ╠═8aaabf9e-08f7-11ed-1880-e798d62d8295
# ╠═3d81d403-fbcc-4bc9-85f7-c1f257dc7695
# ╠═03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
# ╠═853758b5-fecf-4a02-a299-4b5d7f0d99de
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╠═4f6f8003-b9f0-4aba-a42c-c79184bec77d
# ╠═c095ec23-9567-4418-95c9-e011e632ddd0
# ╠═6df2cbcb-28fc-477f-986f-cfc9efea28b8
# ╠═6407e43e-fe73-4c94-a839-4a42f112549b
# ╟─f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# ╠═783ebc81-712b-4f4b-a048-43e04c528652
# ╠═e12c8812-a220-4451-a2e7-602e29838010
# ╠═8fe43246-a006-4ff9-bb64-3601d4e03938
# ╟─6e7e284e-88a1-47f9-9da8-4a79797d8ebd
# ╠═c83851f5-5854-4218-8904-8418073c9296
# ╠═0398c69e-e594-4adb-8099-6524e2db5a31
# ╠═dcab561e-a163-4541-bccf-7a4065d01697
# ╠═f57ac756-9ce6-4378-8d52-010c13795040
# ╠═00c68cd9-f39b-4272-9991-95fcb487b119
# ╟─03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
# ╟─ddae0102-eba2-438f-8eb6-62380e905684
# ╠═7913b885-4f25-40ae-ab84-3b5a4d50967c
# ╠═14160a98-07a0-4efc-9e71-9aab36ed01b6
# ╠═0300688a-6899-46d2-92bd-4b018085e6e8
# ╠═e4ff7a0a-243d-4b51-9a49-7172edca1126
# ╟─99a11ca4-b453-4e0b-a737-d9733b0e8c59
# ╠═718fbf98-1725-429b-9a1c-609cb38430ac
# ╟─65597330-86b4-4d17-abae-5df638603cfc
# ╠═ca22eb02-282a-412b-9297-15a1e379a7b2
# ╠═5f913a6e-1798-42da-8fbe-3dad64fc7bb7
# ╠═0a6a2157-b02f-478b-947b-f579ca1c78c0
# ╠═d13e0a1a-8eb4-49eb-a8b7-b732cf9ffe31
# ╠═f547f778-8a74-4923-9373-cb96b2254c9b
# ╠═93d7da9e-320e-4e01-928e-40d29b2b2be7
# ╠═96581e7a-ee05-48a5-9120-88a1e0d37e7c
# ╠═94310819-6bf4-427e-85d8-82d11ceae3ce
# ╠═2ae832a9-f98e-42f8-a6d5-999c2bca44cb
# ╠═8316a0b8-4e9f-47bd-abfd-3bba96a5f3a3
# ╠═795b757c-b5fd-4dac-9950-7b690848319a
# ╠═746a87a1-0910-4bcd-93af-e4ca870c8847
# ╠═1f36611e-53f1-4371-b121-c8ea8de15d18
# ╠═042f1233-6cfb-4cee-b858-79c970dd56a7
# ╠═c9929a2b-9d35-4f25-95b9-01e79da646c6
# ╠═069fdfda-23b2-46be-bbc7-a612944dfad2
# ╠═16a0d7e4-65df-48f2-9f8e-c83d9370cac6
# ╠═1231b113-242d-4ae0-acf7-72fc0600eb6a
# ╟─6d89975b-844c-496b-aecd-c941ca596a83

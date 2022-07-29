using HDF5
using DataFrames

"""
    _getdf(dataset)
"""
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

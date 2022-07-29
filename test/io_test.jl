using julne
fn = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/Bi/label_beersheba_554mm_214Bi_ICS.h5"
@testset "io" begin
    dfs  = julne.get_dfs(fn)
@test length(keys(dfs)) == 3
end

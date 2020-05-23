#= Graph Utilities
Low level graph functions for processing OxygenMols

=#

"""
    dfs(mol)

Runs depth-first search on an OxygenMol.
Returns pre- and post- orders for each atom in the mol.

# Examples
```
julia> dfs(smilestomol("CC(=O)OCC"))
6-element Array{Tuple{Int64,Int64},1}:
 (1, 12)
 (2, 11)
 (3, 4)
 (5, 10)
 (6, 9)
 (7, 8)
```
"""
function dfs(m::OxygenMol)::Array{Tuple{Int, Int}}
    visited = [false for _ in 1:length(m.adj)]
    pres = [0 for _ in 1:length(m.adj)]
    posts = [0 for _ in 1:length(m.adj)]
    clock = 1

    function pre(v)
        pres[v] = clock
        clock += 1
    end

    function post(v)
        posts[v] = clock
        clock += 1
    end

    function explore(v)
        visited[v] = true
        pre(v)

        for (w, _) in m.adj[v]
            if !visited[w]
                explore(w)
            end
        end

        post(v)
    end

    for v in 1:length(m.adj)
        if !visited[v]
            explore(v)
        end
    end

    collect(zip(pres, posts))
end

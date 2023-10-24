#= =#
@views function mapstoN!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping back to mesh
    @threads for dim ∈ 1:meD.nD
        @simd for p ∈ 1:mpD.nmp
            # accumulation
            if dim == 1 
                meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p] 
            end
            meD.pn[  mpD.p2n[:,p],dim].+=(mpD.ϕ∂ϕ[:,p,1].*mpD.m[p].*(mpD.v[p,dim].+(mpD.∇v[dim,:,p]'*mpD.δx[:,p,:]'))')
            meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
            meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p])
        end
    end
    return nothing
end
@views function mapstoP!(mpD,meD,Δt)
    # mapping back to mp's
    @simd for dim ∈ 1:meD.nD
        # flip update
        @threads for p ∈ 1:mpD.nmp        
            mpD.v[p,dim] =    (mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
            mpD.x[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end          
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,whereto)
    if whereto == "p->n"
        mapstoN!(mpD,meD,g)
    elseif whereto == "p<-n"
        mapstoP!(mpD,meD,Δt)
    end
    return nothing
end


























#=
@views function mapstoN!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping back to mesh
    for dim ∈ 1:meD.nD
        lk = ReentrantLock()
        @threads for p ∈ 1:mpD.nmp
            # accumulation
            lock(lk) do 
                if dim == 1 
                    meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p] 
                end
                meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*mpD.v[p,dim])
                meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
                meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p]) 
            end
        end
    end
    return nothing
end
=#
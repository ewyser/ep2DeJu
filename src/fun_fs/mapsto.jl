#--------------------------------------------------------------------------------------------
## doer functions
#--------------------------------------------------------------------------------------------
@views function mUSLmapstoN!(mpD,meD,g)
    # initialize nodal quantities
    meD.Mn  .= 0.0
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping back to mesh
    @threads for dim ∈ 1:meD.nD
        @simd for p ∈ 1:mpD.nmp
            # accumulation
            if dim == 1 
                # lumped mass matrix
                meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p]
                # consistent mass matrix
                meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p] 
            end
            meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*mpD.v[p,dim])
            meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
            meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p])
        end
    end
    return nothing
end
@views function tpicUSLmapstoN!(mpD,meD,g)

    return nothing
end
@views function mUSLmapstoP!(mpD,meD,Δt)
    # mapping back to mp's
    @simd for dim ∈ 1:meD.nD
        # flip update
        @threads for p ∈ 1:mpD.nmp        
            mpD.v[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.an[mpD.p2n[:,p],dim])
            mpD.x[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end          
    end
    return nothing
end
@views function tpicUSLmapstoP!(mpD,meD,Δt)

    return nothing
end
@views function DM!(mpD,meD,Δt)
    # initialize for DM
    meD.pn .= 0.0
    meD.Δun.= 0.0
    # accumulate material point contributions
    @threads for dim ∈ 1:meD.nD
        # accumulation
        @simd for p ∈ 1:mpD.nmp
            meD.pn[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p].*mpD.v[p,dim])
        end
    end    
    # solve for nodal incremental displacement
    @simd for dim ∈ 1:meD.nD
        @threads for n ∈ 1:meD.nno[meD.nD+1]
            if meD.mn[n]>0.0
                meD.Δun[n,dim] = (Δt*meD.pn[n,dim]*(1.0/meD.mn[n])*meD.bc[n,dim])
            end    
        end
    end
    # update material point's displacement
    @simd for dim ∈ 1:meD.nD
        @threads for p ∈ 1:mpD.nmp
            mpD.u[p,dim]+= (mpD.ϕ∂ϕ[:,p,1]'*meD.Δun[mpD.p2n[:,p],dim])
        end        
    end
    return nothing
end
#--------------------------------------------------------------------------------------------
## dispatcher functions
#--------------------------------------------------------------------------------------------
@views function mapstoN!(mpD,meD,g,trsfrAp)
    if trsfrAp == :mUSL
        mUSLmapstoN!(mpD,meD,g)
    elseif trsfrAp == :tpicUSL

    end
    # lumped mass matrix
    #meD.mn .= sum(meD.Mn,dims=2)
    return nothing
end
@views function mapstoP!(mpD,meD,Δt,trsfrAp)
    if trsfrAp == :mUSL
        mUSLmapstoP!(mpD,meD,Δt)
        DM!(         mpD,meD,Δt)
    elseif trsfrAp == :tpicUSL

    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,trsfrAp,whereto)
    if whereto == "p->n"
        mapstoN!(mpD,meD,g,trsfrAp)
    elseif whereto == "p<-n"
        mapstoP!(mpD,meD,Δt,trsfrAp)
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
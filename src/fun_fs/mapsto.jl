@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0 
    # accumulate material point contributions
    iD  = zeros(Int64,meD.nn)
    @simd for dim ∈ 1:meD.nD
        @simd for p ∈ 1:mpD.nmp
            # index & buffer
            iD .= mpD.p2n[:,p]
            # accumulation
            if dim == 1 meD.mn[iD].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p] end
            meD.pn[  iD,dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*mpD.v[p,dim])
            meD.fext[iD,dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
            meD.fint[iD,dim].+= (mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p]))
        end
    end
    return nothing
end
@views function flip!(mpD,meD,Δt)
    # flip update
    @simd for dim ∈ 1:meD.nD
        # mapping back to mp's
        @threads for p ∈ 1:mpD.nmp        
            mpD.v[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.an[mpD.p2n[:,p],dim])
            mpD.x[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end          
    end
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
@views function mapsto!(mpD,meD,g,Δt,whereto)
    if whereto == "p->N"
        accum!(mpD,meD,g)
    elseif whereto == "p<-N"
        flip!(mpD,meD,Δt)
        DM!(  mpD,meD,Δt)
    end
    return nothing
end
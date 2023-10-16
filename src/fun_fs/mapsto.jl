@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    iD,ϕmp  = zeros(Int64,meD.nn),zeros(Float64,meD.nn)
    @simd for nD ∈ 1:meD.nD
        @simd for p ∈ 1:mpD.nmp
            # index & buffer
            iD .= mpD.p2n[p,:]
            ϕmp.= mpD.ϕ∂ϕ[p,:,1]*mpD.m[p]
            # accumulation
            if nD == 1 meD.mn[iD].+= ϕmp end
            meD.pn[  iD,nD].+= (ϕmp.*mpD.v[p,nD])
            meD.fext[iD,nD].+= (ϕmp.*g[nD])
            meD.fint[iD,nD].+= (mpD.V[p].*(mpD.B[:,nD:meD.nD:end,p]'*mpD.σ[:,p]))
        end
    end
    return nothing
end
@views function flipDM!(mpD,meD,Δt)
    # flip update
    @simd for nD ∈ 1:meD.nD
        # mapping back to mp's
        @threads for p ∈ 1:mpD.nmp        
            mpD.v[p,nD]+= Δt*(mpD.ϕ∂ϕ[p,:,1]'*meD.an[mpD.p2n[p,:],nD])
            mpD.x[p,nD]+= Δt*(mpD.ϕ∂ϕ[p,:,1]'*meD.vn[mpD.p2n[p,:],nD])
        end          
    end
    # initialize for DM + BCs procedure
    meD.pn .= 0.0
    meD.Δun.= 0.0
    # accumulate material point contributions
    @threads for nD ∈ 1:meD.nD
        # accumulation
        @simd for p ∈ 1:mpD.nmp
            meD.pn[mpD.p2n[p,:],nD].+= (mpD.ϕ∂ϕ[p,:,1].*mpD.m[p].*mpD.v[p,nD])
        end
    end    
    # solve for nodal incremental displacement
    @simd for nD ∈ 1:meD.nD
        @threads for n ∈ 1:meD.nno[meD.nD+1]
            if meD.mn[n]>0.0
                meD.Δun[n,nD]= (Δt*meD.pn[n,nD]*(1.0/meD.mn[n])*meD.bc[n,nD])
            end    
        end
    end
    # update material point's displacement
    @simd for nD ∈ 1:meD.nD
        @threads for p ∈ 1:mpD.nmp
            mpD.u[p,nD]+= (mpD.ϕ∂ϕ[p,:,1]'*meD.Δun[mpD.p2n[p,:],nD])
        end        
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,whereto)
    if whereto == "p->N"
        accum!(mpD,meD,g)
    elseif whereto == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end
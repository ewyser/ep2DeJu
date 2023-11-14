#--------------------------------------------------------------------------------------------
## doer functions
#--------------------------------------------------------------------------------------------
@views function flipMapping!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p->n"
        # initialize nodal quantities
        #meD.Mn  .= 0.0
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping back to mesh
        if meD.nD == 2
            @threads for dim ∈ 1:meD.nD
                @simd for p ∈ 1:mpD.nmp
                    # accumulation
                    for nn ∈ 1:meD.nn
                        if dim == 1 
                            # lumped mass matrix
                            meD.mn[mpD.p2n[nn,p]]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                        end
                        meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                        meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                        if dim == 1
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p].*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[3,p])
                        elseif dim == 2
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p].*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p])
                        end
                    end
                end
            end
        elseif meD.nD == 3
            @threads for dim ∈ 1:meD.nD
                @simd for p ∈ 1:mpD.nmp
                    # accumulation
                    for nn ∈ 1:meD.nn
                        if dim == 1 
                            # lumped mass matrix
                            meD.mn[mpD.p2n[nn,p]]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                        end
                        meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                        meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                        if dim == 1
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p].*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[5,p])
                        elseif dim == 2
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p].*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[4,p])
                        elseif dim == 3
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p].*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[3,p])
                        end
                    end
                end
            end
        end

    elseif mapsto == "p<-n"
        # mapping back to mp's
        @simd for dim ∈ 1:meD.nD
            @threads for p ∈ 1:mpD.nmp        
                # flip update
                mpD.v[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.an[mpD.p2n[:,p],dim])
                mpD.x[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
            end          
        end        
    end
    return nothing
end
@views function DM!(mpD,meD,Δt)
    # initialize for DM
    meD.pn.= 0.0
    meD.vn.= 0.0
    # accumulate material point contributions
    @threads for dim ∈ 1:meD.nD
        # accumulation
        @simd for p ∈ 1:mpD.nmp
            meD.pn[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p].*mpD.v[p,dim])
        end
    end    
    # solve for nodal incremental displacement
    @simd for dim ∈ 1:meD.nD
        @threads for n ∈ 1:meD.nno[end]
            if meD.mn[n]>0.0
                meD.vn[n,dim] = (meD.pn[n,dim]*(1.0/meD.mn[n])*meD.bc[n,dim])
            end    
        end
    end
    # update material point's displacement
    @simd for dim ∈ 1:meD.nD
        @threads for p ∈ 1:mpD.nmp
            mpD.u[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end        
    end
    return nothing
end
#--------------------------------------------------------------------------------------------
## dispatcher functions
#--------------------------------------------------------------------------------------------
@views function mapstoN!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
    elseif trsfrAp == :picflipUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance modee"
        throw(error(err_msg))
    elseif trsfrAp == :tpicUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance modee"
        throw(error(err_msg))
    end
    return nothing
end
@views function mapstoP!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
        DM!(       mpD,meD,Δt)
    elseif trsfrAp == :picflipUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance modee"
        throw(error(err_msg))
    elseif trsfrAp == :tpicUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance modee"
        throw(error(err_msg))
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,trsfrAp,whereto)
    if whereto == "p->n"
        mapstoN!(mpD,meD,g,Δt,trsfrAp,whereto)
    elseif whereto == "p<-n"
        mapstoP!(mpD,meD,g,Δt,trsfrAp,whereto)
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
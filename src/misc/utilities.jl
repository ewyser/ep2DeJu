function getKwargs(kwargs)
    if isempty(kwargs)
        # ϵp2De(40,"P","mohr")
        ϕ∂ϕType,fwrkDeform,isΔFbar = :bsmpm,:finite,true
    else
        #ϵp2De(40,"P","MC";shpfun="bsmpm",fwrk="finite",vollock=true)
        kwargs0 = (:shpfun => :bsmpm, :fwrk => :finite, :vollock => true)
        arg     = [kwargs0[1][2],kwargs0[2][2],kwargs0[3][2]]
        for REF in enumerate(kwargs0), argin in enumerate(kwargs)  
            NUM,FIELD = REF
            num,field = argin
            if field[1]==FIELD[1]
                arg[NUM] = field[2]
            end
        end
        ϕ∂ϕType,fwrkDeform,isΔFbar = arg
        if ϕ∂ϕType != :bsmpm && ϕ∂ϕType != :gimpm
            err_msg = "$(ϕ∂ϕType): shape function undefined"
            throw(error(err_msg))
        elseif fwrkDeform != :finite && fwrkDeform != :infinitesimal
            err_msg = "$(fwrkDeform): deformation framework undefined"
            throw(error(err_msg))
        elseif eltype(isΔFbar) != Bool
            err_msg = "$(isΔFbar): not a valid boolean"
            throw(error(err_msg))
        end
    end
    return ϕ∂ϕType,fwrkDeform,isΔFbar
end
function getVersion()
    return string(Pkg.project().version)
end
@views function getVals(meD,mpD,it,ηmax,ηtot,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("nel,np",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),
            ("iteration(s)",it),
            ("ηmax,ηtot",(ηmax,ηtot)),
            (symb*" t/T",cmpl)]
    return vals
end
function msg(message)
    message = "│\n└ "*message
    try
        return printstyled(message,color=:red,bold=true,blink=true)
    catch
        return printstyled(message,color=:blink)
    end
end

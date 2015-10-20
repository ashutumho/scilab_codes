function output = findct(h)
    [lhs,rhs]=argn(0)
    if typeof(h)<>'rational' then
        error(msprintf(gettext("input argument is not rational/ transfer function type")))
    end
    if h.dt=='c' then
        output = %T
    else
        output = %F
    end
endfunction

// Order is a function that returns the highest degree in the case of transfer
//function and in case of state space it returns the number of states.
function [output]= order(h)
    [lhs,rhs]=argn(0)
    //disp(rhs)
    if rhs ~= 1 then
        error(msprintf("Wrong size for input argument"))
    end
// find the data type
    datatype = typeof(h)
    if datatype == 'rational' then                  //if transfred data is transfer function type
        output = length(coeff(h.den))-1             //it will show the order of the transfer function
    elseif datatype == 'state-space' 
        dim = size(h("A"))
        output = dim(1)                            // it will tranfer the number of states of the state-space equation
    else
         error(msprintf("Wrong data type"))
    end
   // output = 1
endfunction

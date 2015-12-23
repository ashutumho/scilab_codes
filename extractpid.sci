// extarct the pid data from the given function
function [Kp, Ki, Kd, Tf, Ts] = extractpid(varargin)
//store the dimension of the data
    [lhs,rhs]=argn(0)
    total_dim = []
//find out the data type     
    for ii =1:rhs
        if and(typeof(varargin(ii))<>[ "rational" "state-space" "constant"  "hypermat"]) then
            error(msprintf(gettext("%s: Wrong type for input argument \n #%d: Linear dynamical system expected.\n"),"extractPID",1))
        end
    end
//store the data in varargin_data and find out the dimension of it
    varargin_data = varargin(1)                                                 //storing the array of the transfer function
    current_dim = size(varargin_data)                                           //size of the array
    current_length = length(current_dim)                                        //total number of dimension
    total_num_elements = 1
//find the total number of elements in the array
    for ii = 1:current_length
        total_num_elements = total_num_elements*current_dim(ii) 
    end
    //varargin_data = hypermat([1 1 current_dim])
    
//***************************************************************************************************************    
// find the considered array element
    if rhs>2 then
        for ii =2:rhs
            if varargin(ii)<=current_dim(ii-1) then
                total_dim = [total_dim, varargin(ii)]
            end//end of the current dimsion store
        end//end of the for loop
// compare the size of the array and trnasfred index
         if current_length ~= length(total_dim) then
             error(msprintf(gettext("Please enter a right dimensions")))
         end
// find the corresponding index
        temp_length = 1
        total_num_elements_count = 0
        temp_dim = ones(1,current_length)
        temp_index_counter = 1
        dim_index = 1
        while current_length ~= temp_length & total_num_elements_count < total_num_elements
             total_num_elements_count = total_num_elements_count+1
             ii=1
           //-------------------------------------------------------------------
            while ii <= current_length & temp_dim(1,ii)== total_dim(ii)
                temp_index_counter = temp_index_counter + 1
                ii = ii+1
            end
           //-------------------------------------------------------------------

           if temp_index_counter == current_length+1 then
               temp_length = current_length
            end
            for kk = 1:temp_index_counter
                if kk == current_length+1 then
                    kk = current_length
                end
                jj = kk
                if temp_dim(1,jj) < current_dim(jj) then
                    temp_dim(1,jj)=temp_dim(1,jj)+1
                elseif temp_dim(1,jj) == current_dim(jj) then
                    temp_dim(1,jj)=1
                end
            end
            
       
            if temp_index_counter <= current_length then
               dim_index = dim_index +1
               temp_index_counter = 1
            end
    end
//**************************************************************************************************************
    elseif rhs == 2 then
        if varargin(2)> total_num_elements then
            error(msprintf(gettext("Please enter a right index")))
        end
        dim_index = varargin(2)
    end //end of considerd array element
//**************************************************************************************************************
//temporary variable is define here
//    temp_Kp = hypermat([1 1 current_dim])
//    temp_Ki = hypermat([1 1 current_dim])
//    temp_Kd = hypermat([1 1 current_dim])
    //
//    temp_Ts = 0
//**************************************************************************************************************
    if current_length == 2 & current_dim == [1 1] then
        varargin_data_index = varargin_data
    elseif current_length > 2 then
        varargin_data_index = varargin_data.entries
    else//if current_length == 2 & current_dim ~= [1 1] then
        temp_varargin_data = hypermat([1 1 current_dim])
        temp_varargin_data(:,:,:,:) = varargin_data
        varargin_data_index = temp_varargin_data.entries
    end
//***************************************************************************************************************
    if varargin_data.dt == 'c' then
        sysTs = 0
        //temp_Tf = 0
        for ii = 1:total_num_elements
         //----------------------------------------------------
         // take out the coeff from the transfer function
            temp_tf = varargin_data_index(ii,1)
            num_coeff = coeff(temp_tf.num) 
            den_coeff = coeff(temp_tf.den)
            length_num_coeff = length(num_coeff)
            length_den_coeff = length(den_coeff)
            last_num_coeff = num_coeff(length_num_coeff)
            last_den_coeff = den_coeff(length_den_coeff)
            num_coeff = num_coeff/last_num_coeff
            den_coeff = den_coeff/last_den_coeff
         //----------------------------------------------------
            if length_num_coeff == 1 & length_den_coeff == 1 then
                temp_Kp_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                temp_Ki_data(ii,1) = 0
                temp_Kd_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
            //I type of controller
            elseif length_num_coeff == 1 & length_den_coeff == 2 & den_coeff(1) == 0 then
                temp_Kp_data(ii,1) = 0
                temp_Ki_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                temp_Kd_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
            // D type of controller
            elseif length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == 0 then     
                temp_Kp_data (ii,1) = 0
                temp_Ki_data (ii,1) = 0
                temp_Kd_data (ii,1) = last_num_coeff
                temp_Tf_data(ii,1) = 0
            // PI type of controller
            elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 then
                temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*num_coeff(1)
                temp_Kd_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
            // PD type of controller
            elseif length_num_coeff == 2 & length_den_coeff == 1 then 
                temp_Ki_data(ii,1) = 0
                temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                temp_Kp_data(ii,1) = temp_Kd_data(ii,1)*num_coeff(1)
                temp_Tf_data(ii,1) = 0
            // PID type of controller
            elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == 0 then
                temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                temp_Kp_data(ii,1) = num_coeff(2)*temp_Kd_data(ii,1)
                temp_Ki_data(ii,1) = num_coeff(1)*temp_Kd_data(ii,1)
                temp_Tf_data(ii,1) = 0
            elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(1) == 0 then
                temp_Tf_data(ii,1) = 1/den_coeff(2)
                //disp(temp_Tf_data(ii,1))
                A = [temp_Tf_data(ii,1) 0 1;
                    (1- temp_Tf_data(ii,1)*num_coeff(2)) temp_Tf_data(ii,1) -num_coeff(2);
                    -temp_Tf_data(ii,1)*num_coeff(1) 1 -num_coeff(1)]
                b = [temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff 0 0]
                temp_pid = A\b'
                temp_Kp_data(ii,1) = temp_pid(1)
                temp_Ki_data(ii,1) = temp_pid(2)
                temp_Kd_data(ii,1) = temp_pid(3)
            else
                disp(varargin_data_index(ii,1))
                error(msprintf(gettext("Above equation is not in parallel PID form")))
            end
        end
        //end of continuous domain
    elseif varargin_data.dt =='d' | varargin_data.dt =='constant' | varargin_data.dt ~= 0  then
        printf('\n in discrete domain')
        if varargin_data.dt == 'd' then
            varargin_data.dt = 1
        end
        sysTs = varargin_data.dt
        for ii = 1:total_num_elements
         //----------------------------------------------------
         // take out the coeff from the transfer function
            temp_tf = varargin_data_index(ii,1)
            num_coeff = coeff(temp_tf.num) 
            den_coeff = coeff(temp_tf.den)
            length_num_coeff = length(num_coeff)
            length_den_coeff = length(den_coeff)
            last_num_coeff = num_coeff(length_num_coeff)
            last_den_coeff = den_coeff(length_den_coeff)
            num_coeff = num_coeff/last_num_coeff
            den_coeff = den_coeff/last_den_coeff
            disp(num_coeff)
            disp(den_coeff)
         //----------------------------------------------------
         if length_num_coeff == 1 & length_den_coeff == 1 then
                temp_Kp_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                temp_Ki_data(ii,1) = 0
                temp_Kd_data(ii,1) = 0
          // I type controller
          elseif length_num_coeff == 1 & length_den_coeff == 2 then
                temp_Kp_data(ii,1) = 0
                temp_Ki_data(ii,1) = last_num_coeff/sysTs
                temp_Kd_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
          // D type of controller (when Tf ~= 0)
          elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                temp_Ki_data(ii,1) = 0
                temp_Kp_data(ii,1) = 0
           // D type of controller (when Tf = 0)
           elseif length_num_coeff == 2 & num_coeff(1) == 0 & length_den_coeff == 1 then
                temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                temp_Kp_data(ii,1) = 0
                temp_Ki_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
          //PI type of controller
          elseif length_num_coeff ==2 & length_den_coeff ==2 & den_coeff(1)==-1 then
                temp_Kp_data(ii,1) = last_num_coeff
                temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                temp_Kd_data(ii,1) = 0
                temp_Tf_data(ii,1) = 0
          // PID type of controller(when Tf ~= 0)
          elseif length_num_coeff == 3 &  length_den_coeff == 3 & (sysTs/(1-den_coeff(1))) == (sysTs/(2+den_coeff(2))) then
                temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                A = [temp_Tf_data(ii,1),0,1;
                     ((sysTs-2*temp_Tf_data(ii,1))-num_coeff(2)*temp_Tf_data(ii,1)),(sysTs*temp_Tf_data(ii,1)),-(2+num_coeff(2));
                     ((temp_Tf_data(ii,1)-sysTs)-num_coeff(1)*temp_Tf_data(ii,1)),(sysTs*(sysTs-temp_Tf_data(ii,1))),(1-num_coeff(1))]
                b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0 0]
                temp_pid= A\b'
                //disp(temp_pid)
                temp_Kp_data(ii,1) = temp_pid(1)
                temp_Ki_data(ii,1) = temp_pid(2)
                temp_Kd_data(ii,1) = temp_pid(3)
            // PID type controller(when Tf = 0)
            elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1)== 0 then
                temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                temp_Kp_data(ii,1) = temp_Kd_data(ii,1)*num_coeff(2)
                temp_Ki_data(ii,1) = temp_Kd_data(ii,1)*num_coeff(1)
                temp_Tf_data(ii,1) = 0
          //PD type controller
            elseif length_num_coeff == 2 & length_den_coeff == 2 then
                //printf('\n bingo')
                temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                A = [temp_Tf_data(ii,1) 1;
                     (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0]
                temp_pid = A\b'
                temp_Kp_data(ii,1) = temp_pid(1)
                temp_Kd_data(ii,1) = temp_pid(2)
                temp_Ki_data(ii,1) = 0 
          else
                disp(varargin_data_index(ii,1))
                error(msprintf(gettext("Above equation is not in parallel PID form")))
         end
         //disp(temp_Kp_data(ii,1))
        end
        //end of discrete domain
    end // end of the data extractor
//**************************************************************************************************************
    
    //disp(varargin_data)
        if rhs >=2 & current_dim ~= [1 1] then
            Kp = temp_Kp_data(dim_index,1)
            Ki = temp_Ki_data(dim_index,1)
            Kd = temp_Kd_data(dim_index,1)
            Ts = sysTs
            Tf = temp_Tf_data(dim_index,1)
            
        elseif rhs == 1 then
            temp_Kp = hypermat([current_dim],temp_Kp_data)
            temp_Ki = hypermat([current_dim],temp_Ki_data)
            temp_Kd = hypermat([current_dim],temp_Kd_data)
            temp_Tf = hypermat([current_dim],temp_Tf_data)
            Kp = temp_Kp
            Ki = temp_Ki
            Kd = temp_Kd
            Tf = temp_Tf
            Ts = sysTs
        end
endfunction// end of the function

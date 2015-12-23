function [output, impdata ] =  pid(varargin)
        [lhs,rhs] = argn(0)
//        disp(rhs)
//        disp(lhs)
        s = poly(0,'s')
        z = poly(0,'z')
        pid1 = [];
        integral_data = [];
        diff_data = [];
        pid_function = [];
        store_data = [];
//______________________________________________________________________________
        // find the type of data that are transferred
        if typeof(varargin(1))=='constant' | typeof(varargin(1))=='hypermat'  then
            datatype = 1
        elseif typeof(varargin(1))=='rational' then 
            datatype =2
        end
        impdata = tlist(["listtype","Kp","Ki","Kd","Tf","Iformula","Dformula","InputDelay","OutputDelay","Ts","Name","Notes","UserData"]) //
//______________________________________________________________________________
        // find the count
        count_numb = 0
        for jj = 1:rhs
            if typeof(varargin(jj)) == 'constant' | typeof(varargin(jj)) == 'hypermat' | typeof(varargin(jj)) == 'rational' then
                count_numb = count_numb + 1
            end
            
        end
//------------------------------------------------------------------------------
        if count_numb <= 4 & datatype == 1 then                                                             //continuous time domain
            printf('\n continuous time domain\n')
            sysTs = 'c'
            kk = 0
            qq = 0
                for ii=1:rhs
                    // check dimension
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    //disp(size(varargin(ii)))
                    if size(varargin(ii))==[1 1] then
                        qq=qq+1
                        if qq == rhs then
                            current_length =1
                            current_size = [1 1]
                        end
                        
                    else
                        if kk==0 then 
                            temp_size = size(varargin(ii))
                            temp_length = length(temp_size)
                            kk=kk+1
                        end
                            current_size = size(varargin(ii))
                            current_length = length(current_size)
                        if temp_length ~= current_length then
                                error(msprintf(gettext("Please enter a right dimensions")))
                        else
                            max_size = max(current_length,temp_length)
                            for jj =1:max_size
                                 if current_size(jj)~=temp_size(jj) then
                                     error(msprintf(gettext("Please enter a right dimensions")))
                                 end
                            end
                        end
                    end
                   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            //******************************************************************
                //find the total number of array elements
                
                total_num_element = 1
                for ii=1:current_length
                    total_num_element = total_num_element*current_size(ii)
                end
                //disp(total_num_element)
                // data transfer
                varargin_data=varargin
                
                for ii =1:total_num_element
                    //##########################################################
                    for jj =1:4
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if jj<=rhs then
                            if typeof(varargin(jj))=='hypermat' then
                                if jj<=rhs then
                                    pid_data(jj)= varargin_data(jj).entries(ii)
                                elseif jj>rhs then
                                    pid_data(jj)=0
                                end
                            elseif typeof(varargin_data(jj))=='constant' then
                                if size(varargin_data(jj))== [1 1] then
                                    pid_data(jj)=varargin_data(jj)
                                else
                                    yup = varargin_data(jj)
                                    pid_data(jj)= yup(ii)
                                end
                            end
                        elseif jj>rhs then
                            pid_data(jj)= 0
                        end
                    
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    //##########################################################
                    //disp(pid_data)
                    dd(ii,1) = syslin('c', pid_data(1)+pid_data(2)/s+pid_data(3)*s/(pid_data(4)*s+1))
                    temp_Kp_data(ii,1) = pid_data(1)
                    temp_Ki_data(ii,1) = pid_data(2)
                    temp_Kd_data(ii,1) = pid_data(3)
                    temp_Tf_data(ii,1) = pid_data(4)
                end
                if current_length<=2 then
                    output_data = hypermat([current_size],dd(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)=dd
                end
                output = output_data
            //end
//*************************************************************************************************************************************************************************************
            elseif count_numb == 5 & datatype ==1 then                                                            //discrete time domain
                printf('\n Discrete time system\n')
                if size(varargin(5))~=[1 1] then
                 error(msprintf(gettext("Please enter valid sampling timing. " )))
                end
                if varargin(5)<=0 then
                    error(msprintf(gettext("0 or neagtive sampling is not valid sampling timing. " )))
                end
            kk = 0
            qq = 0
                for ii=1:5
                    // check dimension
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    //disp(size(varargin(ii)))
                    if size(varargin(ii))==[1 1] then
                        qq=qq+1
                        if qq == 5 then
                            current_length =1
                            current_size = [1 1]
                        end
                    else
                        if kk==0 then 
                            temp_size = size(varargin(ii))
                            temp_length = length(temp_size)
                            kk=kk+1
                        end
                            current_size = size(varargin(ii))
                            current_length = length(current_size)
                        if temp_length ~= current_length then
                                error(msprintf(gettext("Please enter a right dimensions")))
                        else
                            max_size = max(current_length,temp_length)
                            for jj =1:max_size
                                 if current_size(jj)~=temp_size(jj) then
                                     error(msprintf(gettext("Please enter a right dimensions")))
                                 end
                            end
                        end
                    end
                   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            //******************************************************************
                //find the total number of array elements
                
                total_num_element = 1
                for ii=1:current_length
                    total_num_element = total_num_element*current_size(ii)
                end
                //disp(total_num_element)
                
                varargin_data=varargin
                
                for ii =1:total_num_element
                    //##########################################################
                    for jj =1:5
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if jj<=rhs then
                            if typeof(varargin(jj))=='hypermat' then
                                if jj<=rhs then
                                    pid_data(jj)= varargin_data(jj).entries(ii)
                                elseif jj>rhs then
                                    pid_data(jj)=0
                                end
                            elseif typeof(varargin_data(jj))=='constant' then
                                if size(varargin_data(jj))== [1 1] then
                                    pid_data(jj)=varargin_data(jj)
                                else
                                    yup = varargin_data(jj)
                                    pid_data(jj)= yup(ii)
                                end
                            end
                        elseif jj>rhs then
                            pid_data(jj)= 0
                        end
                    
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    //##########################################################
                    pid1 = [pid1 pid_data]
                end
                temp_Kp_data = pid1(1,:)'
                temp_Ki_data = pid1(2,:)'
                temp_Kd_data = pid1(3,:)'
                temp_Tf_data = pid1(4,:)'
//********************************************************************************************************************
                //now find the type of discrete time integral and derivative formula
                Iformula_index = 0
                Dformula_index = 0
                if rhs >5 then
                    for ii = 6:rhs
                        if varargin_data(ii)=='Iformula' then
                            Iformula_index = ii
                        elseif varargin_data(ii)=='Dformula' then
                            Dformula_index = ii
                        end
                    end
                end
                //validity of sampling formula
                if Iformula_index~= 0 then
                    if Iformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Iformula_index+1)~='F'& varargin_data(Iformula_index+1)~= 'B'& varargin_data(Iformula_index+1)~= 'T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
                if Dformula_index~= 0 then
                    if Dformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Dformula_index+1)~='F' & varargin_data(Dformula_index+1)~='B' & varargin_data(Dformula_index+1)~='T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
//********************************************************************************************************************
                 //section of integral formula
                if Iformula_index==0 | varargin_data(Iformula_index+1)=='F' then
                    Iformula = varargin_data(5)/(z-1)
                    Iformula_method = 'F'
                elseif varargin_data(Iformula_index+1)=='B' then
                    Iformula = varargin_data(5)*z/(z-1)
                    Iformula_method = 'B'
                elseif varargin_data(Iformula_index+1)=='T' then
                    Iformula = varargin_data(5)*(z+1)/(2*(z-1))
                    Iformula_method = 'T'
                else//if Iformula_index+1>rhs|varargin_data(Iformula_index+1)~='T'|varargin_data(Iformula_index+1)~='B'|varargin_data(Iformula_index+1)~='F' then
                    Iformula = varargin_data(5)/(z-1)
                    Iformula_method = 'F'
                end
                 //selection of derivative formula
                if Dformula_index==0 |varargin_data(Dformula_index+1)=='F' then
                    temp_Tf = varargin_data(4)
                    temp_Kd = varargin_data(3)
                    size_Tf = size(temp_Tf)
                    if temp_Kd ~= 0 & temp_Tf ~= 0 then
                            if size_Tf == [1 1] then
                                temp_index = 1
                            else
                                temp_index = total_num_element 
                            end
                            for ii =1:temp_index
                                if varargin_data(5)>=2*temp_Tf(ii) then
                                    error(msprintf(gettext("Sample timing must be greater then two times of filter timing Ts>2*Tf")))
                                end
                            end
                     end
                    Dformula = varargin_data(5)/(z-1)
                    Dformula_method = 'F'
                elseif varargin_data(Dformula_index+1)=='B' then
                    Dformula = varargin_data(5)*z/(z-1)
                    Dformula_method = 'B'
                elseif varargin_data(Dformula_index+1)=='T' then
                    Dformula = varargin_data(5)*(z+1)/(2*(z-1))
                    Dformula_method = 'T'
                else//if Dformula_index+1>rhs|varargin_data(Dformula_index+1)~='T'|varargin_data(Dformula_index+1)~='B'|varargin_data(Dformula_index+1)~='F' then
                    temp_Tf = varargin_data(4)
                    temp_Kd = varargin_data(3)
                    size_Tf = size(temp_Tf)
                    if temp_Kd ~= 0 & temp_Tf ~= 0 then
                        if size_Tf == [1 1] then
                            temp_index = 1
                        else
                            temp_index = total_num_element 
                        end
                        for ii =1:temp_index
                            if varargin_data(5)>=2*temp_Tf(ii) then
                                error(msprintf(gettext("Sample timing must be greater then two times of filter timing Ts>2*Tf")))
                            end
                        end
                     end
                    Dformula = varargin_data(5)/(z-1)
                    Dformula_method = 'F'
                end
////********************************************************************************************************************
                  for ii = 1:total_num_element
                    pid_function(ii,1) = syslin('d',pid1(1,ii)+pid1(2,ii)*Iformula+pid1(3,ii)/(pid1(4,ii)+Dformula))
                  end
                if current_length<=2 then
                    output_data = hypermat([current_size], pid_function(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)= pid_function
                end
                output_data.dt = varargin_data(5)
                sysTs = output_data.dt
                output = output_data
//**************************************************************************************************************************************************************************************
            elseif datatype == 2 then
                printf('\n rational data type')
                varargin_data = varargin
// storing the dimension of the transfer matrix in current_size
                    current_size = size(varargin_data(1))
                    current_length = length(current_size)
// total number of elements
                    total_num_element = 1
                    for ii=1:current_length
                        total_num_element = total_num_element*current_size(ii)
                    end
                    systf_num = varargin_data(1).num
                    systf_den = varargin_data(1).den
                    sysTs = varargin_data(1).dt
                    disp(sysTs)
                if varargin_data(1).dt == 'c' then                                   // continuous domain
                    for ii = 1:total_num_element
//------------------------------------------------------------------------------
// finding out the transfer function 
                    //trans_func = systf(ii)
                    num_coeff = coeff(systf_num(ii))
                    den_coeff = coeff(systf_den(ii))
                    
                    length_num_coeff = length(num_coeff)
                    length_den_coeff = length(den_coeff)
                    last_num_coeff = num_coeff(length_num_coeff)
                    last_den_coeff = den_coeff(length_den_coeff)
                    num_coeff = num_coeff/last_num_coeff
                    den_coeff = den_coeff/last_den_coeff
                    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
                    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    end
                    // end of continuous domain
                elseif varargin_data(1).dt == 'd' | typeof(varargin_data(1).dt)=='constant' then                //discrete time domain
                    if varargin_data(1).dt == 'd' then
                        varargin_data(1).dt = 1
                    end
                    
//------------------------------------------------------------------------------
        // find the Iformula and Dformula index
                    Iformula_index = 0
                    Dformula_index = 0
                    for ii = 1:rhs
                        temp_data = varargin_data(ii)
                        I_index = isequal(temp_data,'Iformula')
                        D_index= isequal(temp_data,'Dformula')
                        if I_index== %T then
                            Iformula_index = ii
                        elseif D_index == %T then
                            Dformula_index = ii
                        end
                    end
            // end of finding of Iformula and Dformula
//------------------------------------------------------------------------------
    
//------------------------------------------------------------------------------
                    //selection of the type of Iformula and Dformula
                    //Iformula
                    if Iformula_index == 0|Iformula_index+1>rhs then
                        Iformula = 'F'
                    elseif varargin_data(Iformula_index+1)=='F' then
                        Iformula = 'F'
                    elseif varargin_data(Iformula_index+1)=='B' then
                        Iformula = 'B'
                    elseif varargin_data(Iformula_index+1)=='T' then
                        Iformula = 'T'
                    elseif varargin_data(Iformula_index+1)~='F' |varargin_data(Iformula_index+1)~='B'|varargin_data(Iformula_index+1)~='T' then
                        Iformula ='F'
                    end
                    //Dformula
                    if Dformula_index == 0|Dformula_index+1>rhs then
                        Dformula = 'F'
                    elseif varargin_data(Dformula_index+1)=='F' then
                        Dformula = 'F'
                    elseif varargin_data(Dformula_index+1)=='B' then
                        Dformula = 'B'
                    elseif varargin_data(Dformula_index+1)=='T' then
                        Dformula = 'T'
                    elseif varargin_data(Dformula_index+1)~='F' |varargin_data(Dformula_index+1)~='B'|varargin_data(Dformula_index+1)~='T' then
                        Dformula ='F'
                    end
                    //end of type of Iformula and Dformula
                    Iformula_method = Iformula
                    Dformula_method = Dformula
//extract the element from the given trnsfer function
//____________________________________________________________________________________________________________________________________________________________________________
                    for ii = 1:total_num_element
//------------------------------------------------------------------------------
// finding out the transfer function 
                    //trans_func = systf(ii)
                    num_coeff = coeff(systf_num(ii))
                    den_coeff = coeff(systf_den(ii))
                    length_num_coeff = length(num_coeff)
                    length_den_coeff = length(den_coeff)
                    last_num_coeff = num_coeff(length_num_coeff)
                    last_den_coeff = den_coeff(length_den_coeff)
                    num_coeff = num_coeff/last_num_coeff
                    den_coeff = den_coeff/last_den_coeff
//------------------------------------------------------------------------------
                            //type of PID equation 
                            //(Iformula Dformula)
                            //(F F),(F B) (F T)
                            //(B F),(B B) (B T)
                            //(T F),(T B) (T T)
                            if Iformula =='F' & Dformula == 'F' then
                    //******************************************************************************
                    //printf('\n do ***')
                    // when Ifomula = Forward Euler and Dformula = Forward Euler
                            // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                        temp_Kp_data(ii,1) = 0
                                        temp_Ki_data(ii,1) = last_num_coeff/sysTs
                                        temp_Kd_data(ii,1) = 0
                                        temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                             // for PI type Controller
                             //           Ts
                             //  Kp+Ki*---------
                             //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1              Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 &  length_den_coeff == 3  then //& () == (sysTs/(1-den_coeff(1)))
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then    
                                        A = [temp_Tf_data(ii,1),0,1;
                                             ((sysTs-2*temp_Tf_data(ii,1))-num_coeff(2)*temp_Tf_data(ii,1)),(sysTs*temp_Tf_data(ii,1)),-(2+num_coeff(2));
                                             ((temp_Tf_data(ii,1)-sysTs)-num_coeff(1)*temp_Tf_data(ii,1)),(sysTs*(sysTs-temp_Tf_data(ii,1))),(1-num_coeff(1))]
                                        b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0 0]
                                        temp_pid= A\b'
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                      else
                                          disp(systf_num(ii)/(systf_den(ii)))
                                          error(msprintf(gettext("Above equation is not in parallel PID form")))
                                      end
                                      
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form"))) 
                            end
                                
                    // when Ifomula = Forward Euler and Dformula = Backward Euler        
                    //******************************************************************************
                            elseif Iformula == 'F' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //           Ts
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff ==2 & length_den_coeff ==2 & den_coeff(1)==-1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           z*Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = -sysTs*(1+den_coeff(2))/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then 
                                        A = [(temp_Tf_data(ii,1)+sysTs) 0 1;
                                              -((2*temp_Tf_data(ii,1)+sysTs)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(temp_Tf_data(ii,1)+sysTs) -(2+num_coeff(2));
                                               temp_Tf_data(ii,1)-(temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -sysTs*temp_Tf_data(ii,1) 1-num_coeff(1)]
                                        b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                            elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))

                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                            else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                            end
                    // when Ifomula = Forward Euler and Dformula = Trapezoidal
                    //******************************************************************************
                            elseif Iformula == 'F' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = last_num_coeff/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //           Ts
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff ==2 & length_den_coeff ==2 & den_coeff(1)==-1 then
                                    //
                                    temp_Kp_data(ii,1) = last_num_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1          Ts*(z+1)
                            //                      Tf+---------
                            //                           2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = -sysTs*(1+den_coeff(1))/(2*(den_coeff(1)-1))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A =[(2*temp_Tf_data(ii,1)+sysTs) 0 2;
                                            -(4*temp_Tf_data(ii,1)+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(2*temp_Tf_data(ii,1)+sysTs) -(4+2*num_coeff(2));
                                            -((sysTs-2*temp_Tf_data(ii,1))+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) sysTs*(sysTs-2*temp_Tf_data(ii,1)) (2-2*num_coeff(1))]
                                        b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Forward Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'F' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2&length_den_coeff==2&num_coeff(1)==0&den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [temp_Tf_data(ii,1) temp_Tf_data(ii,1)*sysTs 1;
                                             ((sysTs-2*temp_Tf_data(ii,1))-(temp_Tf_data(ii,1)*num_coeff(2))) (sysTs*(sysTs-temp_Tf_data(ii,1))-temp_Tf_data(ii,1)*sysTs*num_coeff(2)) -(2+num_coeff(2));
                                             ((temp_Tf_data(ii,1)-sysTs)-(temp_Tf_data(ii,1)*num_coeff(1))) -temp_Tf_data(ii,1)*sysTs*num_coeff(1) 1-num_coeff(1)]
                                        b = [temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                     else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                    
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1) 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0 
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Backward Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2 & length_den_coeff==2&num_coeff(1)==0 & den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           z*Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = sysTs*(den_coeff(2)-1)/(2-den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A =[(temp_Tf_data(ii,1)+sysTs) sysTs*(temp_Tf_data(ii,1)+sysTs) 1;
                                            ((2*temp_Tf_data(ii,1)+sysTs)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) (temp_Tf_data(ii,1)*sysTs+sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(2)) (2+num_coeff(2));
                                            (temp_Tf_data(ii,1)-(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -sysTs*(temp_Tf_data(ii,1)+sysTs)*num_coeff(1) (1-num_coeff(1))]
                                        b = [(temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                                elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf)*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Trapezoidal Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2 & length_den_coeff==2 & num_coeff(1)==0 & den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1          Ts*(z+1)
                            //                      Tf+---------
                            //                           2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*(1+den_coeff(1))/(2*(1-den_coeff(1)))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [(2*temp_Tf_data(ii,1)+sysTs) sysTs*(sysTs+2*temp_Tf_data(ii,1)) 2;
                                             -(4*temp_Tf_data(ii,1)+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(sysTs-2*temp_Tf_data(ii,1))-sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(2) -(4+2*num_coeff(2));
                                             ((2*temp_Tf_data(ii,1)-sysTs)-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(1) (2-2*num_coeff(1))]
                                        b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Forward Euler
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'F' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1
                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)             Ts
                            //                           Tf+---------
                            //                                 z-1
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [2*temp_Tf_data(ii,1) temp_Tf_data(ii,1)*sysTs 2;
                                            ((2*sysTs-4*temp_Tf_data(ii,1))-2*temp_Tf_data(ii,1)*num_coeff(2)) (sysTs^2-temp_Tf_data(ii,1)*sysTs*num_coeff(2)) -(4+2*num_coeff(2));
                                            (2*(temp_Tf_data(ii,1)-sysTs)-2*temp_Tf_data(ii,1)*num_coeff(1)) (sysTs*(sysTs-temp_Tf_data(ii,1))-temp_Tf_data(ii,1)*sysTs*num_coeff(1)) (2-2*num_coeff(1))]
                                        b = [2*temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1) 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Backward Euler
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1
                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)             Ts
                            //                           Tf+---------
                            //                                 z-1
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = -sysTs*(den_coeff(2)+1)/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [2*(temp_Tf_data(ii,1)+sysTs) sysTs*(temp_Tf_data(ii,1)+sysTs) 2;
                                             -((2*sysTs+4*temp_Tf_data(ii,1))+2*(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) (sysTs^2 -sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(2)) -(4+2*num_coeff(2));
                                             (2*temp_Tf_data(ii,1)-2*(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(temp_Tf_data(ii,1)*sysTs+sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(1)) 2*(1-num_coeff(1))]
                                        b = [2*(temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                                elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)                
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Trapezoidal
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                             // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1
                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)          Ts*(z+1)
                            //                           Tf+---------
                            //                                2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*(1+den_coeff(1))/(2*(1-den_coeff(1)))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [(4*temp_Tf_data(ii,1)+2*sysTs) sysTs*(2*temp_Tf_data(ii,1)+sysTs) 4;
                                             -(8*temp_Tf_data(ii,1)+(4*temp_Tf_data(ii,1)+2*sysTs)*num_coeff(2)) (2*sysTs^2-sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(2)) -(8+4*num_coeff(2));
                                             ((4*temp_Tf_data(ii,1)-2*sysTs)-(4*temp_Tf_data(ii,1)+2*sysTs)*num_coeff(1)) (sysTs*(sysTs-2*temp_Tf_data(ii,1))-sysTs*(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) 4-4*num_coeff(1)]
                                        b = [(4*temp_Tf_data(ii,1)+2*sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                    
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                                
                    //******************************************************************************
                            end
                    //------------------------------------------------------------------------------
                    // when Iformula is Backward euler
                     end// extract the PID data from given transfer function
                    
//______________________________________________________________________________                    
                end// identifty that the given function is "c" type or "d" type
                output = varargin_data(1)
        end
//______________________________________________________________________________

                if current_length<=2 then
                    //output_data = hypermat([current_size], pid_function(:,1))
                    Kp = hypermat([current_size],temp_Kp_data(:,1))
                    Ki = hypermat([current_size],temp_Ki_data(:,1))
                    Kd = hypermat([current_size],temp_Kd_data(:,1))
                    Tf = hypermat([current_size],temp_Tf_data(:,1))
                elseif current_length > 2 then
                    //output_data = hypermat([current_size])
                    //output_data(:,:,:,:)= pid_function
                    Kp = hypermat([current_size])
                    Kp(:,:,:,:) = temp_Kp_data
                    Ki = hypermat([current_size])
                    Ki(:,:,:,:) = temp_Ki_data
                    Kd = hypermat([current_size])
                    Kd(:,:,:,:) = temp_Kd_data
                    Tf = hypermat([current_size])
                    Tf(:,:,:,:) = temp_Tf_data
                end
                impdata.Kp = Kp
                impdata.Ki = Ki
                impdata.Kd = Kd
                impdata.Tf = Tf
                impdata.Ts = sysTs
//------------------------------------------------------------------------------
                if sysTs == 'c' then
                    impdata.Iformula = '[]'
                else
                    if Iformula_method == 'F' then
                        impdata.Iformula = 'ForwardEuler'
                    elseif Iformula_method == 'B' then
                        impdata.Iformula = 'BackwardEuler'
                    elseif Iformula_method == 'T' then
                        impdata.Iformula = 'Trapezoidal'
                    end
                end
//------------------------------------------------------------------------------
                if sysTs == 'c' then
                    impdata.Dformula = '[]'
                else
                    if Dformula_method == 'F' then
                        impdata.Dformula = 'ForwardEuler'
                    elseif Dformula_method == 'B' then
                        impdata.Dformula = 'BackwardEuler'
                    elseif Dformula_method == 'T' then
                        impdata.Dformula = 'Trapezoidal'
                    end
                end
//------------------------------------------------------------------------------
                impdata.InputDelay = 0
                impdata.OutputDelay = 0
//------------------------------------------------------------------------------
                Namedata = 0
                Notesdata = 0
                Datauser = 0
                //Namedata = find( varargin == "Name" )
                //disp(Namedata)
                if count_numb < rhs then
                    for ii = count_numb:rhs
                        if typeof(varargin(ii)) == 'hypermat' | typeof(varargin(ii)) == 'rational' then
                        else
                            if varargin(ii) == 'Name' then
                                Namedata = ii
                            elseif varargin(ii) == 'Notes' then
                                Notesdata = ii
                            elseif varargin(ii) == 'UserData' then
                                Datauser = ii
                            end
                        end
                        
                    end
                end
//------------------------------------------------------------------------------                
                if Namedata == 0 then
                    impdata.Name = []
                elseif Namedata + 1 > rhs then
                    impdata.Name = []
                else
                    impdata.Name = varargin(Namedata+1)
                end
//------------------------------------------------------------------------------                
                if Notesdata == 0 then
                    impdata.Notes = []
                elseif Notesdata + 1 > rhs then
                    impdata.Notes = []
                else
                    impdata.Notes = varargin(Notesdata+1)
                end
//------------------------------------------------------------------------------                
                if Datauser == 0 then
                    impdata.UserData = []
                elseif Datauser + 1 > rhs then
                    impdata.UserData = []
                else
                    impdata.UserData = varargin(datauser)
                end
//------------------------------------------------------------------------------                
disp(Namedata)
disp(impdata.Name)
disp(Notesdata)
disp(impdata.Notes)
endfunction


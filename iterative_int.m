function itint=iterative_int(u,t,Ntrunc)
    
    u0=ones(1,length(t));
    u=[u0;u];

    dt=t(2)-t(1);
    
    length_t=length(t);
    
    num_input=size(u,1);
    
    total_iterint=num_input*(1-num_input^Ntrunc)/(1-num_input);
    
    
    Etemp=zeros(total_iterint,length_t);
    ctrEtemp=zeros(Ntrunc,1);
    
    
    for i=1:Ntrunc
        ctrEtemp(i+1)=ctrEtemp(i)+num_input^i;
    end
    
    
    
    sum_acc=cumsum(u,2)*dt;
    Etemp(1:num_input,:)=[zeros(num_input,1) sum_acc(:,1:end-1)];
    
    
    
    for i=2:Ntrunc
        start_prev_block=ctrEtemp(i-1)+1;
        end_prev_block=ctrEtemp(i);
        end_current_block=ctrEtemp(i+1);
        num_prev_block = end_prev_block - start_prev_block+1;
        num_current_block = end_current_block - end_prev_block;
        
        U_block=u(repelem(1:size(u, 1), num_prev_block), :);%reshape(repmat(u',num_prev_block,1),length_t,[])'
        prev_int_block=repmat(Etemp(start_prev_block:end_prev_block,:),num_input,1);
        current_int_block=cumsum(U_block.*prev_int_block,2)*dt;
        Etemp(end_prev_block+1:end_current_block,:)=[zeros(num_current_block,1) current_int_block(:,1:end-1)];
    end
    itint=Etemp;


end




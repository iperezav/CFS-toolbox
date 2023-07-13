function pws=power_series(h,vector_field,z,Ntrunc)

    num_vfield=size(vector_field,2);
    
    total_lderiv=num_vfield*(1-num_vfield^Ntrunc)/(1-num_vfield);
    
    
    Ltemp=zeros(total_lderiv,1,'sym');
    ctrLtemp=zeros(Ntrunc,1);
    
    
    for i=1:Ntrunc
        ctrLtemp(i+1)=ctrLtemp(i)+num_vfield^i;
    end
    
    
    Ltemp(1:num_vfield)=transpose(jacobian(h,z)*vector_field);
    
    
    
    for i=2:Ntrunc
        start_prev_block=ctrLtemp(i-1)+1;
        end_prev_block=ctrLtemp(i);
        end_current_block=ctrLtemp(i+1);
        num_prev_block = end_prev_block - start_prev_block+1;
        num_current_block = end_current_block - end_prev_block;
    
        LT=Ltemp(start_prev_block:end_prev_block);
        LT=jacobian(LT,z)*vector_field;
        LT=transpose(LT);
        LT=reshape(LT,[],1);
        Ltemp(end_prev_block+1:end_current_block,:)=LT;
    end
    pws=Ltemp;



end
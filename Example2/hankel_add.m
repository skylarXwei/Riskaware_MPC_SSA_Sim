function hanekl_out = hankel_add(hankel,newmeasurement,L,index_measurement)
    if index_measurement <= L
        hanekl_out = [hankel;newmeasurement];
    else
        hanekl_out = [hankel,[hankel(2:end,end);newmeasurement]];
    end
end
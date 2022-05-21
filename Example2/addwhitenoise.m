function noise_xtraj = addwhitenoise(xtraj, noiselevel,ttraj)
    [ncoln,nrow] = size(xtraj);
    noise_xtraj = [];
    for j = 1:nrow
        stds =1;
        noise_xtraj(:,j) = xtraj(:,j) + stds*noiselevel*(rand(ncoln,1)-0.5);
    end
    
%     figure();
%     subplot(3,1,1);
%     plot(ttraj,xtraj(:,1));
%     hold on
%     plot(ttraj,noise_xtraj(:,1));
%     subplot(3,1,2);
%     plot(ttraj,xtraj(:,2));
%         hold on
%     plot(ttraj,noise_xtraj(:,2));
%     subplot(3,1,3);
%     plot(ttraj,xtraj(:,3));
%         hold on
%     plot(ttraj,noise_xtraj(:,3));
end
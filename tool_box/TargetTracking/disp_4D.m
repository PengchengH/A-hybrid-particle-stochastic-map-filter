function [] = disp_4D(xf_2, Yt, xp, model, mean_t, n, Cycle_count, Acycles)
    figure(1);
    subplot(4,4,[1 5]);
    plot(xf_2(:,4),xf_2(:,1),'og','LineWidth',3,'MarkerSize',6);
%   hold on;
%   plot(Yt(1),Yt(2),'^b','LineWidth',3,'MarkerSize',10);
%   hold off;
    title('prior');
    
    subplot(4,4,[2 6]);
    plot(xf_2(:,4),xf_2(:,2),'og','LineWidth',3,'MarkerSize',6);
    title('prior');
    
    subplot(4,4,[3 7]);
    plot(xf_2(:,4),xf_2(:,3),'og','LineWidth',3,'MarkerSize',6);
    title('prior');
    
    
    subplot(4,4,[9 13]);
    plot(xp(:,4),xp(:,1),'og','LineWidth',3,'MarkerSize',6);
%     hold on;
%     plot(Yt(1),Yt(2),'^b','LineWidth',3,'MarkerSize',10);
%     plot(model.xt(1,Acycles(n),Cycle_count),model.xt(2,Acycles(n),Cycle_count),'pr','LineWidth',3,'MarkerSize',10);    
%     plot(mean_t(1,n,Cycle_count),mean_t(2,n,Cycle_count),'Or','LineWidth',3,'MarkerSize',10);
%     hold off;
    title('posterior');
    
    subplot(4,4,[10 14]);
    plot(xp(:,4),xp(:,2),'og','LineWidth',3,'MarkerSize',6);
    title('posterior');
    
    subplot(4,4,[11 15]);
    plot(xp(:,4),xp(:,3),'og','LineWidth',3,'MarkerSize',6);
    title('posterior');
    
    subplot(4,4,[4 8]);
    plot(model.xt(1,Acycles(1):Acycles(n),Cycle_count),model.xt(2,Acycles(1):Acycles(n),Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(1,1:n,Cycle_count),mean_t(2,1:n,Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 1,2');
    
    
    subplot(4,4,12);
    plot(model.xt(3,Acycles(1):Acycles(n),Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(3,1:n,Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 3');
    
    subplot(4,4,16);
    plot(model.xt(4,Acycles(1):Acycles(n),Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(4,1:n,Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 4');
end
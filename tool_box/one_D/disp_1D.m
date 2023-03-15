function [] = disp_1D(xf_2, Yt, xp, model, mean_t, n, Cycle_count, Acycles)

    K= size(xp,1);
    
    figure(1);
    subplot(2,2,1);
    plot(xf_2(1:K/2),xf_2(K/2+1:K),'og','LineWidth',3,'MarkerSize',6);
    hold on;
%     plot(Yt(1),Yt(1),'^b','LineWidth',3,'MarkerSize',10);
    hold off;
    title('prior');
    subplot(2,2,2);
    plot(xp(1:K/2),xp(K/2+1:K),'og','LineWidth',3,'MarkerSize',6);
    hold on;
%     plot(Yt(1),Yt(1),'^b','LineWidth',3,'MarkerSize',10);
    plot(model.xt(1,Acycles(n), Cycle_count),model.xt(1,Acycles(n), Cycle_count),'pr','LineWidth',3,'MarkerSize',10);    
    plot(mean_t(1,n, Cycle_count),mean_t(1,n, Cycle_count),'Or','LineWidth',3,'MarkerSize',10);
    hold off;
    title('posterior');
    
    subplot(2,2,[3 4]);
    plot(model.xt(1,Acycles(1):Acycles(n), Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(1,1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 1');
end
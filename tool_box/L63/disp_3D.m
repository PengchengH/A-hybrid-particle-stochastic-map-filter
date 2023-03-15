function [] = disp_3D(xf_2, Yt, xp, model, mean_t, n, Cycle_count, Acycles)
    figure(1);
    subplot(4,2,[1 3]);
    plot3(xf_2(:,1), xf_2(:,2), xf_2(:,3),'og','LineWidth',3,'MarkerSize',6);
    hold on;
    plot3(Yt(1), Yt(2), Yt(3),'^b','LineWidth',3,'MarkerSize',10);
    hold off;
    title('prior');
    subplot(4,2,[5 7]);
    plot3(xp(:,1), xp(:,2), xp(:,3),'og','LineWidth',3,'MarkerSize',6);
    hold on;
    plot3(Yt(1), Yt(2), Yt(3),'^b','LineWidth',3,'MarkerSize',10);
    plot3(model.xt(1, Acycles(n), Cycle_count), model.xt(2,Acycles(n)),model.xt(3,Acycles(n), Cycle_count),'pr','LineWidth',3,'MarkerSize',10);    
    plot3(mean_t(1,n, Cycle_count), mean_t(2,n, Cycle_count), mean_t(3,n, Cycle_count),'Or','LineWidth',3,'MarkerSize',10);
    hold off;
    title('posterior');
    
    subplot(4,2,2);
    plot(model.xt(1, Acycles(1):Acycles(n), Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(1, 1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 1');
    
    subplot(4,2,4);
    plot(model.xt(2, Acycles(1):Acycles(n), Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(2, 1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 2');
    
    subplot(4,2,6);
    plot(model.xt(3, Acycles(1):Acycles(n), Cycle_count),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(3, 1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 3');
end
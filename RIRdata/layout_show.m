% 房间尺寸
function layout_show(array)
    room_dim = array.roomSize;   % [x y z] in meters
    
    % 麦克风与扬声器位置
    mic_pos = array.bCtrPtsPositions;  % 1 x 3 明区控制
    mic_pos2 = array.dCtrPtsPositions;  % 暗区控制
    spk_pos = array.s;                % N x 3
    
    figure;
    hold on;
    
    % ==== 1. 绘制房间立方体 ====
    % 顶点顺序（共8个）
    corner = [0 0 0;
              room_dim(1) 0 0;
              room_dim(1) room_dim(2) 0;
              0 room_dim(2) 0;
              0 0 room_dim(3);
              room_dim(1) 0 room_dim(3);
              room_dim(1) room_dim(2) room_dim(3);
              0 room_dim(2) room_dim(3)];
    
    % 边连接关系（从 corner 的行索引）
    edges = [1 2; 2 3; 3 4; 4 1;    % 底面
             5 6; 6 7; 7 8; 8 5;    % 顶面
             1 5; 2 6; 3 7; 4 8];   % 竖边
    
    edges([5 6 10],:) = [];  % Remove lines that obstruct the view
    
    % 画框架
    for i = 1:size(edges,1)
        p1 = corner(edges(i,1),:);
        p2 = corner(edges(i,2),:);
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-');
    end
    
    % ==== 2. 麦克风 ====
    h1 = scatter3(mic_pos(:,1), mic_pos(:,2), mic_pos(:,3), 10, ...
        'MarkerEdgeColor','#5892E8','MarkerFaceColor','#5892E8');
        % for i = 1:size(mic_pos,1)
        %     text(mic_pos(i,1), mic_pos(i,2), mic_pos(i,3)+0.1,sprintf('M%d', i), 'Color', 'b');
        % end
    h2 = scatter3(mic_pos2(:,1), mic_pos2(:,2), mic_pos2(:,3), 10, ...
        'MarkerEdgeColor','#B2B9CB','MarkerFaceColor','#B2B9CB');
    
    % ==== 3. 扬声器 ====
    h3 = scatter3(spk_pos(:,1), spk_pos(:,2), spk_pos(:,3), 80, ...
        'MarkerEdgeColor','#000000','MarkerFaceColor','#CC4230');
    % for i = 1:size(spk_pos,1)
    %     text(spk_pos(i,1), spk_pos(i,2), spk_pos(i,3)+0.1, ...s
    %         sprintf('L%d', i), 'Color','r');
    % end
    
    % ==== 4. 设置坐标范围与标签 ====
    xlabel('X (m)','FontName','Times New Roman','FontWeight','bold'); 
    ylabel('Y (m)','FontName','Times New Roman','FontWeight','bold'); 
    zlabel('Z (m)','FontName','Times New Roman','FontWeight','bold');
    xlim([0 room_dim(1)]);
    ylim([0 room_dim(2)]);
    zlim([0 room_dim(3)]);
    
    % 只想在图例里显示 h1 和 h3
    lgd = legend([h1 h2 h3], {'BZ Control Points','DZ Control Points' ...
        ,'Loudspeakers'});
    lgd.NumColumns = 3;
    lgd.FontName   = 'Times New Roman'; % 字体
    lgd.FontSize   = 8;                % 字号
    lgd.FontWeight = 'bold';            % 加粗
    lgd.Position = [0.383,0.72681,0.2742,0.0312];
    grid on; axis equal;
    view(47, 20);  % 调整观察视角
    title('Room Layout with Microphone and Speakers','FontSize',12, ...
        'FontName','Times New Roman','FontWeight','bold');
end

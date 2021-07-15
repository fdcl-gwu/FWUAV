function [theta, W] = R2axang(R)
% Conversion from a Rotation matrix to axis-angle representation
    theta = acos((trace(R)-1)/2);
    if theta ~= 0
		W = (1/ (2*sin(theta))) * [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
    else
        W = zeros(3, 1);
    end
%     out = vrrotmat2vec(R);
%     W = out(1:3)';
%     theta = out(4);
%     % To keep the second axis positive
%     if W(2) < 0
%         W = -W;
%         theta = -theta;
%     end
end

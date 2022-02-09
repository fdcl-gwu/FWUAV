function [JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_6_9(JJ)
coder.inline('always');
% Decomposes the inertia matrix into sub-blocks of 6 and 9
JJ_11 = JJ(1:6,1:6);
JJ_12 = JJ(1:6,7:15);
JJ_21 = JJ(7:15,1:6);
JJ_22 = JJ(7:15,7:15);
end

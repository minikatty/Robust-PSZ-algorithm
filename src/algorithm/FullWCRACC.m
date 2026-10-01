function res = FullWCRACC(HB_ctrl, HD_ctrl, para, fre_indices)
%FULLWCRACC Worst-case robust ACC with the complete ATF expansion.
%
% This is the project-level algorithm entry for the Full-WCRACC ablation.
% The SDP implementation is kept in FullCrossTermRACC so existing response
% experiments remain backward compatible.
%
% The recommended common uncertainty setting is
%
%   para.full_racc.relative_eta = nu;
%
% which gives eta_Z(f)=nu*||H_Z(f)||_F. Both Full-WCRACC and its NoCT
% counterpart use the same eta_Z. The optional common loading is
%
%   para.full_racc.diagonal_loading_ratio = 1e-6;

    if nargin < 3 || isempty(para)
        para = struct();
    end
    if nargin < 4
        res = FullCrossTermRACC(HB_ctrl, HD_ctrl, para);
    else
        res = FullCrossTermRACC(HB_ctrl, HD_ctrl, para, fre_indices);
    end
    res.algorithm = 'Full-WCRACC';
end

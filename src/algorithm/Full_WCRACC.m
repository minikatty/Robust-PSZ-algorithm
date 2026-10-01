function res = Full_WCRACC(HB_ctrl, HD_ctrl, para, fre_indices)
%FULL_WCRACC Project-facing Full-WCRACC algorithm entry.
    if nargin < 4
        res = FullWCRACC(HB_ctrl, HD_ctrl, para);
    else
        res = FullWCRACC(HB_ctrl, HD_ctrl, para, fre_indices);
    end
end

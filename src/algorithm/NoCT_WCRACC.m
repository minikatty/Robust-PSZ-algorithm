function res = NoCT_WCRACC(HB_ctrl, HD_ctrl, para, fre_indices)
%NOCT_WCRACC Project-facing no-cross-term WCRACC algorithm entry.
    if nargin < 4
        res = NoCrossTermWCRACC(HB_ctrl, HD_ctrl, para);
    else
        res = NoCrossTermWCRACC(HB_ctrl, HD_ctrl, para, fre_indices);
    end
end

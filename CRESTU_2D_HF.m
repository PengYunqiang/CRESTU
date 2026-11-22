function r=CRESTU_2D_HF(width_m,draft_m,options)
%CRESTU_2D_HF Research-only rectangular high-frequency mass, separate from RC1.
    if nargin<3,options=struct;end
    root=fileparts(mfilename('fullpath'));CRESTU_SETUP('2d');addpath(fullfile(root,'research','p2_hf'));
    defaults=struct('rho',1025,'depth',420,'reference',[0 0],'quadrature_order',12,'panels_per_edge',128);
    fields=fieldnames(defaults);for j=1:numel(fields),f=fields{j};if ~isfield(options,f),options.(f)=defaults.(f);end,end
    assert(width_m>0 && draft_m>0 && draft_m<options.depth);
    nodes=hf_rectangle(width_m,draft_m,options.panels_per_edge);
    r=crestu2d.hf_solve(nodes,options);
end

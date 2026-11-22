function run_hf_validation(stage,outputRoot)
%RUN_HF_VALIDATION Bounded, reproducible stages. Old outputs are never reused silently.
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(root);CRESTU_SETUP('2d');addpath(fileparts(mfilename('fullpath')));
    if nargin<2,outputRoot=fullfile(root,'hf_results');end
    out=fullfile(outputRoot,stage);assert(~isfolder(out),'Choose a fresh stage: historical evidence is immutable.');mkdir(out);
    opt=struct('rho',1025,'depth',420,'reference',[0 0],'quadrature_order',12);
    if strcmp(stage,'smoke')
        opt.depth=Inf;theta=linspace(-pi,0,129)';nodes=[cos(theta),sin(theta)];nodes([1 end],2)=0;
        r=crestu2d.hf_solve(nodes,opt);save(fullfile(out,'circle.mat'),'r');
        fprintf('CIRCLE Aheave=%.12g exact=%.12g rel=%.6g recip=%.6g\n',r.A(2,2),1025*pi/2,abs(r.A(2,2)/(1025*pi/2)-1),r.reciprocity);
        opt.depth=420;r=CRESTU_2D_HF(7,7,opt);save(fullfile(out,'rectangle.mat'),'r');
        fprintf('RECT Aheave=%.12g recip=%.6g rcond=%.6g\n',r.A(2,2),r.reciprocity,r.matrix_rcond);return
    end
    if strcmp(stage,'stage_a'),widths=[7 10];drafts=[3 7 10];counts=[16 32 64 128];
    elseif strcmp(stage,'table'),widths=[7 10];drafts=[.2 .5 1 2 3 4 5 6 7 7.5 8 9 10 11 12 13 14 14.8 15];counts=[64 128];
    else,error('Unknown stage');end
    rows={};
    for b=widths,for d=drafts,for n=counts
        opt.panels_per_edge=n;r=CRESTU_2D_HF(b,d,opt);
        record=struct('width_m',b,'draft_m',d,'panels_per_edge',n,'body_panels',r.boundary.n,'water_depth_m',420, ...
          'A_heave_kg_per_m',r.A(2,2),'CA_inf',r.A(2,2)/(1025*b*d),'reciprocity',r.reciprocity, ...
          'minimum_scaled_eigenvalue',r.minimum_scaled_eigenvalue,'linear_residual',r.linear_residual,'matrix_rcond',r.matrix_rcond,'elapsed_seconds',r.elapsed_seconds);
        rows{end+1}=record;writetable(struct2table([rows{:}]),fullfile(out,'results.csv'));
        save(fullfile(out,sprintf('B%g_D%g_N%d.mat',b,d,n)),'r');
        fprintf('HF %s b=%g d=%g N=%d A=%.10f recip=%.3g time=%.3g\n',stage,b,d,n,r.A(2,2),r.reciprocity,r.elapsed_seconds);
    end,end,end
end

function results = HydroMain(config_file)
%HYDROMAIN Entry point for the CRESTU frequency-domain Rankine solver.
    here=fileparts(mfilename('fullpath')); root=fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'), ...
        fullfile(root,'6.MeanDriftLoads'),here);
    if nargin<1||isempty(config_file), config_file=fullfile(here,'CRESTU.cfg'); end
    results=run_frequency_domain_case(config_file);
end

function result = CRESTU_2D(input)
%CRESTU_2D Zero-speed, linear 2-D transverse-section hydrodynamics.
% result=CRESTU_2D('examples/2d/semicircle_rigid_lid.c2d.json')
% Coordinates: (y,z), z up. DOFs: sway, heave, positive roll about +x.
% Phasor exp(+i*omega*t); potentials are normalized to unit velocity.
% Independent of the CRESTU_V2 3-D input and numerical implementation.
    root = CRESTU_SETUP('2d');
    cfg = crestu2d.read_case(input);
    if cfg.output.save
        retained={'result.mat','radiation.csv','excitation.csv','coefficients.json'};
        for j=1:numel(retained)
            assert(~isfile(fullfile(cfg.output.directory,retained{j})), ...
                'CRESTU2D:OutputExists','Retained output exists. Supply a new output.directory: %s',cfg.output.directory);
        end
    end
    timer=tic;
    metadata=crestu2d.run_metadata(root,cfg);
    result = crestu2d.solve(cfg);
    metadata.elapsed_seconds=toc(timer);
    metadata.completed_utc=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
    metadata.radiation_normalization=result.radiation_normalization;
    ids=cfg.dof_indices;
    metadata.units=struct('A',{result.units.A(ids,ids)},'B',{result.units.B(ids,ids)}, ...
        'F_exc',{result.units.F_exc(ids)},'F_exc_per_amplitude',{result.units.F_exc_per_amplitude(ids)});
    metadata.geometry_diagnostics=result.geometry.diagnostics;
    result.metadata=metadata;
    result.module_version='0.4.0-rc1';
    if cfg.output.save
        crestu2d.save_result(result);
    end
end

function save_result(result)
%SAVE_RESULT Reproducible resolved input, MAT arrays, tidy CSV and static plots.
    out=result.input.output.directory;
    if ~isfolder(out),mkdir(out);end
    save(fullfile(out,'result.mat'),'result');
    fid=fopen(fullfile(out,'resolved_input.json'),'w'); assert(fid>=0);
    clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(result.input,PrettyPrint=true)); clear clean;
    if isfield(result,'metadata')
        fid=fopen(fullfile(out,'metadata.json'),'w','n','UTF-8'); assert(fid>=0);
        clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(result.metadata,PrettyPrint=true)); clear clean;
    end
    crestu2d.export_coefficients(result,fullfile(out,'coefficients.json'));
    nd=numel(result.dofs); nf=numel(result.omega_rad_s);
    rows=cell(nd*nd*nf,7); r=0;
    for f=1:nf
        for i=1:nd
            for j=1:nd
                r=r+1;rows(r,:)={result.omega_rad_s(f),result.dofs{i},result.dofs{j}, ...
                    result.A(i,j,f),result.B(i,j,f),result.A_hat(i,j,f),result.B_hat(i,j,f)};
            end
        end
    end
    writetable(cell2table(rows,'VariableNames',{'omega_rad_s','force_dof','motion_dof','A_SI','B_SI','A_hat','B_hat'}),fullfile(out,'radiation.csv'));
    if ~isempty(result.F_exc)
        rows=cell(nd*nf,8);r=0;
        for f=1:nf
            for i=1:nd
                r=r+1; F=result.F_exc(i,f);C=result.F_exc_per_amplitude(i,f); K=result.F_FK_per_amplitude(i,f);
                rows(r,:)={result.omega_rad_s(f),result.dofs{i},real(F),imag(F),real(C),imag(C),real(K),imag(K)};
            end
        end
        writetable(cell2table(rows,'VariableNames',{'omega_rad_s','dof','F_real_SI','F_imag_SI', ...
            'F_per_a_real_SI','F_per_a_imag_SI','FK_per_a_real_SI','FK_per_a_imag_SI'}),fullfile(out,'excitation.csv'));
    end
    writetable(struct2table(result.diagnostics),fullfile(out,'diagnostics.csv'));
    fig=figure('Visible','off','Color','w','Position',[100 100 1100 720]); cleaner=onCleanup(@()close(fig));
    tiledlayout(2,2);
    nexttile;plot(result.geometry.nodes(:,1),result.geometry.nodes(:,2),'k.-');axis equal;grid on;
    xlabel('y [m]');ylabel('z [m]');title([result.name ' wetted section'],'Interpreter','none');
    nexttile;hold on;
    for i=1:nd,plot(result.omega_rad_s,squeeze(result.A_hat(i,i,:)),'o-','DisplayName',result.dofs{i});end
    grid on;xlabel('\omega [rad/s]');ylabel('A_{ii} / (\rho L^2 D_i^2)');legend('Location','best');
    nexttile;hold on;
    for i=1:nd,plot(result.omega_rad_s,squeeze(result.B_hat(i,i,:)),'o-','DisplayName',result.dofs{i});end
    grid on;xlabel('\omega [rad/s]');ylabel('B_{ii} / (\rho \omega L^2 D_i^2)');legend('Location','best');
    nexttile;hold on;
    if ~isempty(result.F_exc)
        scale=[1 1 result.input.normalization_length_m];scale=scale(result.input.dof_indices);
        for i=1:nd
            plot(result.omega_rad_s,abs(result.F_exc_per_amplitude(i,:))/(result.input.fluid.rho_kg_m3* ...
                result.input.fluid.g_m_s2*result.input.normalization_length_m*scale(i)),'o-','DisplayName',result.dofs{i});
        end
        ylabel('|F/a| / (\rho g L D_i)');legend('Location','best');
    end
    grid on;xlabel('\omega [rad/s]');
    exportgraphics(fig,fullfile(out,'coefficients.png'),'Resolution',150);
end

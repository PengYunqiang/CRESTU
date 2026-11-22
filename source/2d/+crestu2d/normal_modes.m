function modes = normal_modes(yz,normal,ref)
%NORMAL_MODES Unit sway/heave velocity and roll angular velocity about +x.
    modes=[normal(:,1),normal(:,2), ...
        (yz(:,1)-ref(1)).*normal(:,2)-(yz(:,2)-ref(2)).*normal(:,1)];
end

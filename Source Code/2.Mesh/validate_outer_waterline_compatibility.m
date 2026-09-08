function maximumDistanceM = validate_outer_waterline_compatibility(bodyLine, outerLine, toleranceM)
% VALIDATE_OUTER_WATERLINE_COMPATIBILITY Prevent a gap at body/FS intersection.
% Different node counts are allowed on the same piecewise-linear boundary.
% Vertices and edge midpoints are tested in both directions against segments.
    a=bodyLine.nodes; b=outerLine.nodes;
    maximumDistanceM=max(one_direction(a,bodyLine.is_closed,b,outerLine.is_closed), ...
        one_direction(b,outerLine.is_closed,a,bodyLine.is_closed));
    assert(maximumDistanceM<=toleranceM,'CRESTU:OuterWaterlineMismatch', ...
        ['Body and fixed outer reference waterlines differ by %.6g m (tolerance %.6g m). ', ...
         'Using that reference leaves a body/FS seam gap. Use NONE or a compatible boundary.'], ...
         maximumDistanceM,toleranceM);
end

function d=one_direction(a,closedA,b,closedB)
    if closedA
        queries=[a;0.5*(a+a([2:end,1],:))];
    else
        queries=[a;0.5*(a(1:end-1,:)+a(2:end,:))];
    end
    if closedB, b1=b; b2=b([2:end,1],:); else, b1=b(1:end-1,:); b2=b(2:end,:); end
    delta=b2-b1; length2=sum(delta.^2,2);d=0;
    for j=1:size(queries,1)
        t=sum((queries(j,:)-b1).*delta,2)./max(length2,eps);
        t=max(0,min(1,t)); projection=b1+t.*delta;
        d=max(d,min(sqrt(sum((projection-queries(j,:)).^2,2))));
    end
end

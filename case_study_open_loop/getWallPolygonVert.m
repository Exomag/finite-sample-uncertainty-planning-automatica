function wallPolygonVert = getWallPolygonVert(wallCoefsMean, yMin, xMax)

    % Wall polygon A/B matrices
    nWalls = length(wallCoefsMean);
    wallPolygonA = zeros(nWalls, 2);
    wallPolygonB = zeros(nWalls, 1);
    for i = 1 : nWalls
        wallPolygonA(i,:) = wallCoefsMean{i}(1:2)';
        wallPolygonB(i) = - wallCoefsMean{i}(3);
    end
    wallPolygonA = [0 1; wallPolygonA; 1 0; 0 1];
    wallPolygonB = [-yMin; wallPolygonB; xMax; -yMin];
    
    % Wall polygon vertices
    wallPolygonVert = zeros(length(wallPolygonB)-1, 2);
    for i = 1 : size(wallPolygonVert, 1)
        wallPolygonVert(i,:) = linsolve([wallPolygonA(i,:); wallPolygonA(i+1,:)], ...
            [wallPolygonB(i); wallPolygonB(i+1)]);
    end

end


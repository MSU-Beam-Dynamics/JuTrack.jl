function translate_AT_to_Julia1(THERING, filename)
    % The definition of multipole strength is the same as ELEGANT which is different from AT.
    % *2 for sextupole and *6 for octupole

    fid = fopen(filename, 'w');

    juliaElements = {};

    for idx = 1:length(THERING)
        element = THERING{idx};
        className = element.Class;  

        switch className
            case 'RFCavity'
                eletype = 'RFCA';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, volt=%.10f, h=%.10f, freq=%.10e, energy=%.10e)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.Voltage, element.HarmNumber, element.Frequency, element.Energy);
            case 'Bend'
                eletype = 'SBEND';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, angle=%.10f, e1=%.10f, e2=%.10e, PolynomB=[0.0, %.10e, 0.0, 0.0])\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.BendingAngle, element.EntranceAngle, element.ExitAngle, element.K);
            case 'Quadrupole'
                eletype = 'KQUAD';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k1=%.10f)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.K);
            case 'Sextupole'
                eletype = 'KSEXT';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k2=%.10f)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.PolynomB(3)*2);
            case 'Octupole'
                eletype = 'KOCT';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k3=%.10f)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.PolynomB(4)*6);
            case 'Drift'
                eletype = 'DRIFT';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f)\n', ...
                element.FamName, eletype, element.FamName, element.Length);
            case 'BPM'
                eletype = 'MARKER';
                juliaStr = sprintf('%s = %s(name="%s")\n', ...
                element.FamName, eletype, element.FamName);
            case 'Corrector'
                eletype = 'CORRECTOR';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, xkick=%.10f, ykick=%.10f)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.KickAngle(1), element.KickAngle(2));
            case 'Aperture'
                disp('Aperture is skipped');
                continue;
            case 'Injection'
                disp('Aperture is skipped');
                continue;
        % Add more cases as needed
        otherwise
            error(['unknown element type: ', className]);
        end
        
        fprintf(fid, juliaStr);
        % Append the current element string to the cell array
        juliaElements{end+1} = element.FamName;
    end

    % fprintf(fid, 'elements = [%s]\n', strjoin(juliaElements, ','))
    fclose(fid);  

    remove_repeated_lines(filename, filename)
    
    fid = fopen(filename, 'a');
    fprintf(fid, ['line = Lattice(nelems=',num2str(length(juliaElements)),')\n']);
    for i = 1:length(juliaElements)
        fprintf(fid, ['add!(line, ',juliaElements{i},')\n']);
    end
    fclose(fid);
    disp(['Julia lattice file created: julia_lattice.jl']);
end

function remove_repeated_lines(inputFile, outputFile)
    % Check if the input file exists
    if ~exist(inputFile, 'file')
        error('Input file does not exist');
    end

    % Read lines from the file
    fid = fopen(inputFile, 'r');
    fileContents = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = fileContents{1};

    % Find unique lines and preserve their original order
    [uniqueLines, idx] = unique(lines, 'stable');

    % Write unique lines to the output file
    fid = fopen(outputFile, 'w');
    fprintf(fid, '%s\n', uniqueLines{:});
    fclose(fid);

    disp(['Output file created: ', outputFile]);
end


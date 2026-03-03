function translate_AT_to_Julia(THERING, filename)
    % The definition of multipole strength is the same as ELEGANT which is different from AT.
    % *2 for sextupole and *6 for octupole

    fid = fopen(filename, 'w');

    juliaElements = {};

    for idx = 1:length(THERING)
        element = THERING{idx};
        className = element.Class;  

        switch className
            case 'Marker'
                eletype = 'MARKER';
                juliaStr = sprintf('%s = %s(name="%s")\n', ...
                element.FamName, eletype, element.FamName);
            case 'RFCavity'
                eletype = 'RFCA';
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, volt=%.10f, h=%.10f, freq=%.10e, energy=%.10e)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.Voltage, element.HarmNumber, element.Frequency, element.Energy);
            case 'Bend'
                eletype = 'SBEND';
                T1 = get_element_field(element, 'T1', zeros(6,1));
                R1 = get_element_field(element, 'R1', eye(6));
                T2 = get_element_field(element, 'T2', zeros(6,1));
                R2 = get_element_field(element, 'R2', eye(6));
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, angle=%.10f, e1=%.10f, e2=%.10e, PolynomB=[0.0, %.10e, 0.0, 0.0], T1=%s, R1=%s, T2=%s, R2=%s)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.BendingAngle, element.EntranceAngle, element.ExitAngle, element.K, ...
                matrix_to_julia_str(T1), matrix_to_julia_str(R1), matrix_to_julia_str(T2), matrix_to_julia_str(R2));
            case 'Quadrupole'
                eletype = 'KQUAD';
                T1 = get_element_field(element, 'T1', zeros(6,1));
                R1 = get_element_field(element, 'R1', eye(6));
                T2 = get_element_field(element, 'T2', zeros(6,1));
                R2 = get_element_field(element, 'R2', eye(6));
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k1=%.10f, T1=%s, R1=%s, T2=%s, R2=%s)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.K, matrix_to_julia_str(T1), matrix_to_julia_str(R1), matrix_to_julia_str(T2), matrix_to_julia_str(R2));
            case 'Sextupole'
                eletype = 'KSEXT';
                T1 = get_element_field(element, 'T1', zeros(6,1));
                R1 = get_element_field(element, 'R1', eye(6));
                T2 = get_element_field(element, 'T2', zeros(6,1));
                R2 = get_element_field(element, 'R2', eye(6));
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k2=%.10f, T1=%s, R1=%s, T2=%s, R2=%s)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.PolynomB(3)*2, matrix_to_julia_str(T1), matrix_to_julia_str(R1), matrix_to_julia_str(T2), matrix_to_julia_str(R2));
            case 'Octupole'
                eletype = 'KOCT';
                T1 = get_element_field(element, 'T1', zeros(6,1));
                R1 = get_element_field(element, 'R1', eye(6));
                T2 = get_element_field(element, 'T2', zeros(6,1));
                R2 = get_element_field(element, 'R2', eye(6));
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, k3=%.10f, T1=%s, R1=%s, T2=%s, R2=%s)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                element.PolynomB(4)*6, matrix_to_julia_str(T1), matrix_to_julia_str(R1), matrix_to_julia_str(T2), matrix_to_julia_str(R2));
            case 'Drift'
                eletype = 'DRIFT';
                T1 = get_element_field(element, 'T1', zeros(6,1));
                R1 = get_element_field(element, 'R1', eye(6));
                T2 = get_element_field(element, 'T2', zeros(6,1));
                R2 = get_element_field(element, 'R2', eye(6));
                juliaStr = sprintf('%s = %s(name="%s", len=%.10f, T1=%s, R1=%s, T2=%s, R2=%s)\n', ...
                element.FamName, eletype, element.FamName, element.Length, ...
                matrix_to_julia_str(T1), matrix_to_julia_str(R1), matrix_to_julia_str(T2), matrix_to_julia_str(R2));
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

    fprintf(fid, 'elements = [%s]\n', strjoin(juliaElements, ','));

    fclose(fid);  

    remove_repeated_lines(filename, filename)
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

function value = get_element_field(element, fieldname, default_value)
    % Helper function to safely get a field from element
    % Returns default_value if field doesn't exist
    if isfield(element, fieldname)
        value = element.(fieldname);
    else
        value = default_value;
    end
end

function str = matrix_to_julia_str(mat)
    % Helper function to convert a matrix or vector to Julia format string
    if isvector(mat)
        % For vectors, format as [val1, val2, val3, ...]
        str = sprintf('[%.10e', mat(1));
        for i = 2:length(mat)
            str = sprintf('%s, %.10e', str, mat(i));
        end
        str = sprintf('%s]', str);
    else
        % For matrices, format as [row1; row2; ...]
        [m, n] = size(mat);
        str = sprintf('[');
        for i = 1:m
            if i > 1
                str = sprintf('%s; ', str);
            end
            for j = 1:n
                if j > 1
                    str = sprintf('%s ', str);
                end
                str = sprintf('%s%.10e', str, mat(i,j));
            end
        end
        str = sprintf('%s]', str);
    end
end

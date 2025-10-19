function c = temp2speed(T_celsius, method)
%SOUND_SPEED Compute speed of sound in air based on temperature.
%   c = sound_speed(T_celsius)
%       Computes speed of sound (in m/s) using the empirical formula.
%
%   c = sound_speed(T_celsius, method)
%       Allows selection of computation method:
%       'empirical'     - uses c ≈ 331 + 0.6*T (default)
%       'thermodynamic' - uses c = sqrt(gamma * R * T)

    % Set default method
    if nargin < 2
        method = 'empirical';
    end

    % Convert Celsius to Kelvin
    T_kelvin = T_celsius + 273.15;

    switch lower(method)
        case 'empirical'
            % Empirical formula (approximation)
            % Valid near room temperature in dry air
            c = 331 + 0.6 * T_celsius;

        case 'thermodynamic'
            % Thermodynamic formula: c = sqrt(gamma * R * T)
            gamma = 1.4;        % Adiabatic index for dry air
            R = 287.05;         % Specific gas constant for dry air (J/kg·K)
            c = sqrt(gamma * R * T_kelvin);

        otherwise
            error('Unknown method. Choose ''empirical'' or ''thermodynamic''.');
    end
end

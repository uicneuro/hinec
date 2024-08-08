function colors = vector_to_color(Vx, Vy, Vz)
% Normalize the direction vectors
mag = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
Vx = Vx ./ mag;
Vy = Vy ./ mag;
Vz = Vz ./ mag;

% Handle NaNs resulting from zero magnitudes
Vx(isnan(Vx)) = 0;
Vy(isnan(Vy)) = 0;
Vz(isnan(Vz)) = 0;

% Convert the direction to a colormap value
% Red for left-right, Green for anterior-posterior, Blue for superior-inferior
R = abs(Vx); % Map x-component to Red
G = abs(Vy); % Map y-component to Green
B = abs(Vz); % Map z-component to Blue

colors = [R(:), G(:), B(:)]; % Combine into an RGB matrix
end

function [xb, yb] = naca4_boundary(naca, N_per_side, c, spacing)
% NACA4_BOUNDARY  NACA 4-digit airfoil boundary points (clockwise, TE→TE).
%
%   [xb, yb] = naca4_boundary(naca, N_per_side, c, spacing)
%
% Inputs
%   naca         : char/string like '2412'  (m-p-tt → m/100, p/10, t/100)
%   N_per_side   : number of points on EACH surface (>= 2). Includes LE and TE.
%   c            : chord length (default 1.0)
%   spacing      : 'cosine' (default) or 'linear' distribution along chord
%
% Outputs
%   xb, yb : column vectors of boundary coordinates ordered clockwise,
%            starting and ending at the trailing edge (TE upper → LE → TE lower).
%
% Notes
%   Thickness distribution (standard NACA 4-digit; closed TE):
%     yt = (t*c/0.2)*( 0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c)^2 ...
%                      + 0.2843*(x/c)^3 - 0.1036*(x/c)^4 )
%   Mean camber line (0 ≤ x ≤ c), with ξ = x/c:
%     For ξ < p:  yc = c*(m/p^2)*(2pξ - ξ^2),      dyc/dx = (m/p^2)*(2p - 2ξ)
%     For ξ ≥ p:  yc = c*(m/(1-p)^2)*[(1-2p)+2pξ-ξ^2],  dyc/dx = (m/(1-p)^2)*(2p - 2ξ)
%   Upper/lower surfaces applied normal to camberline with θ = atan(dyc/dx):
%     xU = x - yt*sinθ,   yU = yc + yt*cosθ
%     xL = x + yt*sinθ,   yL = yc - yt*cosθ
%
%   Output order is:  TE (upper) → ... → LE → ... → TE (lower) → TE (repeat)
%   so the first and last points are both the trailing edge.

    if nargin < 4 || isempty(spacing), spacing = 'cosine'; end
    if nargin < 3 || isempty(c),       c = 1.0;            end
    if nargin < 2 || isempty(N_per_side), N_per_side = 101; end

    code = char(naca);
    assert(numel(code) == 4, 'naca must be a 4-character code like ''2412''.');
    m = str2double(code(1))   / 100;   % max camber
    p = str2double(code(2))   / 10;    % location of max camber
    t = str2double(code(3:4)) / 100;   % thickness ratio

    % ---- x distribution on [0,c] (LE at 0, TE at c) ----
    switch lower(spacing)
        case 'cosine'
            beta = linspace(0, pi, N_per_side).';         % 0..pi
            x    = 0.5*c*(1 - cos(beta));                  % LE->TE
        case 'linear'
            x    = linspace(0, c, N_per_side).';
        otherwise
            error('spacing must be ''cosine'' or ''linear''.');
    end
    xi = x / c;  % nondimensional chord coordinate

    % ---- thickness distribution (closed TE: -0.1036) ----
    yt = (t*c/0.2) .* ( 0.2969*sqrt(xi) - 0.1260*xi - 0.3516*xi.^2 ...
                        + 0.2843*xi.^3 - 0.1036*xi.^4 );

    % ---- mean camber line and slope ----
    yc    = zeros(size(x));
    dycdx = zeros(size(x));
    if m > 0 && p > 0 && p < 1
        i1 = xi <  p;
        i2 = xi >= p;

        yc(i1)    = c*(m/p^2)      .* ( 2*p*xi(i1) - xi(i1).^2 );
        yc(i2)    = c*(m/(1-p)^2)  .* ( (1 - 2*p) + 2*p*xi(i2) - xi(i2).^2 );

        dycdx(i1) = (m/p^2)        .* ( 2*p - 2*xi(i1) );
        dycdx(i2) = (m/(1-p)^2)    .* ( 2*p - 2*xi(i2) );
    end
    theta = atan(dycdx);

    % ---- upper & lower surfaces (applied normal to camberline) ----
    xU = x - yt .* sin(theta);
    yU = yc + yt .* cos(theta);
    xL = x + yt .* sin(theta);
    yL = yc - yt .* cos(theta);

    % ---- build clockwise boundary starting and ending at TE ----
    % Upper: go TE->LE (reverse), Lower: LE->TE (forward)
    % Include TE as first AND last point.
    TE_xU = xU(end); TE_yU = yU(end);      % trailing edge (upper)
    TE_xL = xL(end); TE_yL = yL(end);      % trailing edge (lower), should equal TE_xU,TE_yU

    % Assemble (avoid duplicating LE in the middle more than once)
    xb = [ TE_xU;
           xU(end-1:-1:1);     % TE (excluded) back to LE
           xL(2:end-1);          % LE (excluded on lower) forward to TE
           TE_xL ];            % repeat TE as final point
    yb = [ TE_yU;
           yU(end-1:-1:1);
           yL(2:end-1);
           TE_yL ];

    % Ensure column vectors
    xb = xb(:);
    yb = yb(:);
end
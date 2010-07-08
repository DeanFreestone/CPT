
function phi_unwrapped = unwrap_phi(phi)

NSamples = length(phi);
phi_unwrapped = phi+pi;     % initialize, shift up by phi so there is no negative values

PreviousGoodPhase = 0;

for n=2:NSamples
    PreviousPhase = phi_unwrapped(n-1);
    if ~isnan(PreviousPhase)
        PreviousGoodPhase = PreviousPhase;  % if it is not NaN than we can use it
    end
    
    CurrentPhase = phi_unwrapped(n);
    if ~isnan(CurrentPhase)
        if CurrentPhase <= PreviousGoodPhase    % then we must have made a 2 \pi phase jump
            phi_unwrapped(n:end) = phi_unwrapped(n:end) + 2*pi;       % this holds as we assume the phase is monotonically increasing
        end
    end
end
phi_unwrapped=phi_unwrapped-pi;
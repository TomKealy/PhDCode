function impulse_noise_elimination(file_1, file_2, model_order, n, lambda, p, save_ar)

%
%function impulse_noise_elimination(file_1, file_2, model_order, n, lambda,
%p, save_ar)
%
%Function Impulse Noise Elimination input parametres:
% * file_1 is a name of input track, ex. 'music.wav';
% * file_2 is a name of output track, ex. 'music_mod.wav';
% * model_order is a choosen number of estimated model parameters; can be
%   choosen from range: 1 <= model_order <= 25; mostly the best choice is
%   model_order = 4;
% * n is a scaling coefficient of local variance in detecting algorithm;
%   its value decides of the level of impulse noise samples; mostly n = 3
%   is the best choice, but n is from range 2 <= n <= 20;
% * lambda is a memory coefficient in EW - LS algorithm; its value
%   dependends on track dynamic; lambda is from the range: 
%   0.001 (low dynamic) <= lambda <= 0.999 (high dynamic);
% * p is a scaling coefficient of P matrix (eye matrix rank model_order)
%   starting value; p is from range 10^3 <= p <= 10^6;
% * save_ar = 0 or 1; 0 if model parametres should not be saved, 1 if model
%   parameters should be saved.
%
%Function Impulse Noise Elimination is using an EW - LS algorithm
%(exepotencialy weigthed least squares) to estimate signal parametres of
%a music track. Created that way model is dynamic - changes every process
%iteration. Using created model, algorithm is gives an one - step
%prediction error. Impulse noise detecting algorithm is checking whether
%prediction error is not too large. If yes, it means that this sample is
%noise.
%
%--------------------------------------------------------------------------
%created by:
%Jakub Marek Szymanowski (jakub.szymanowski@gmail.com)
%--------------------------------------------------------------------------


%Reading and checking input variebles:
[track fs tmp] = wavread(file_1);

round(model_order);
round(n);

if model_order > 25
    model_order = 25;
elseif model_order < 1
    model_order = 1;
end

if n > 20
    n = 20;
elseif n < 2
    n = 1;
end

if lambda > 1
    lambda = 0.999;
elseif lambda < 0
    lambda = 0.001;
end

if p < 1000
    p = 1000;
elseif p > 1000000
    p = 1000000;
end

if save_ar ~= 1
    save_ar = 0;
elseif save_ar == 1
    model_vector = zeros(max(size(track)), model_order);
end

%Starting parametres of EW - LS algorithm:
org_track = track;
model = zeros(model_order, 1);
K = 0;
P = p * eye(model_order);
variance = 0;


%EW - LS algorithm & detecting algorithm:
for t = (model_order + 1) : (max(size(track)) - 4)
    
    signal = track(t-1 : -1 : t-model_order);
    prediction_error = track(t) - model' * signal;
        
    if save_ar == 1
        model_vector(t,:) = model;
    end
    
    if t == model_order + 1
        variance = prediction_error^2;
    end
       
    if abs(prediction_error) > n * sqrt(variance) %detecting test
        
        noise_samples = 1;
        signal_1 = track(t : -1 : t-model_order+1);
        prediction_error_1 = track(t+1) - model' * signal_1;

        if abs(prediction_error_1) > n * sqrt(variance)
            
            noise_samples = 2;
            signal_2 = track(t+1 : -1 : t-model_order+2);
            prediction_error_2 = track(t+2) - model' * signal_2;
            
            if abs(prediction_error_2) > n * sqrt(variance)
                            
                noise_samples = 3;
                signal_3 = track(t+2 : -1 : t-model_order+3);
                prediction_error_3 = track(t+3) - model' * signal_3;
                
                if abs(prediction_error_3) > n * sqrt(variance)
                    
                    noise_samples = 4;
                    
                end
            
            end
            
        end
  
    else
        noise_samples = 0;
        variance = lambda*variance + (1 - lambda)*(prediction_error)^2;
    end
    
    if noise_samples ~= 0 %interpolation
        for i = 0 : noise_samples - 1
            track(t + i) = track(t-1) + (-1 / (noise_samples + 1))*(track(t-1) - track(t + noise_samples))*(i+1);
        end
        noise_samples = 0;
    end
    
    K = (P * signal) / (lambda + signal'*P*signal);
    P = (1 / lambda)*(P - (P*signal*signal'*P)/(lambda + signal'*P*signal));
    model = model + K * prediction_error;
    
end

wavwrite(track, fs, file_2);

if save_ar == 1
        save model_vector.mat model_vector -MAT
end

%graphs:
clf;
wykres = figure(1);
subplot(221); plot(org_track,'g'); grid; title('Original track'); xlabel('sample [n]'); ylabel('y_o(n)'); axis auto;
subplot(223); plot(track,'r'); grid; title('Modified track'); xlabel('sample [n]'); ylabel('y_m(n)'); axis auto;
subplot(222); plot(abs(org_track - track), 'b'); title('Differences in tracks'); xlabel('sample [n]'); ylabel('|y_o(n) - y_m(n)|'); axis auto;

if save_ar == 1
    subplot(224); plot(model_vector); title('Estimated model parameters'); xlabel('sample [n]'); ylabel('vector of parameters'); axis auto;
end

end
